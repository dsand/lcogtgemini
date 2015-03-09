#!/usr/bin/env python
'''
Created on Nov 7, 2014

@author: cmccully
'''
import os, shutil
from glob import glob
import pyfits
import numpy as np
import lacosmicx
from pyraf import iraf
from scipy import interpolate, ndimage, signal, optimize
import pf_model as pfm
import statsmodels as sm

########## in iraf:
iraf.cd(os.getcwd())
iraf.gemini()
iraf.gmos()
iraf.onedspec()


iraf.gmos.logfile = "log.txt"
iraf.gmos.mode = 'h'
iraf.set(clobber='yes')

iraf.set(stdimage='imtgmos')

dooverscan = True
def sanitizeheader(hdr):
    # Remove the mandatory keywords from a header so it can be copied to a new
    # image.
    hdr = hdr.copy()

    # Let the new data decide what these values should be
    for i in ['SIMPLE', 'BITPIX', 'BSCALE', 'BZERO']:
        if i in hdr.keys():
            hdr.pop(i)

#    if hdr.has_key('NAXIS'):
    if 'NAXIS' in hdr.keys():
        naxis = hdr.pop('NAXIS')
        for i in range(naxis):
            hdr.pop('NAXIS%i' % (i + 1))

    return hdr

def tofits(filename, data, hdr=None, clobber=False):
    """simple pyfits wrapper to make saving fits files easier."""
    from pyfits import PrimaryHDU, HDUList
    hdu = PrimaryHDU(data)
    if not (hdr is None):
        hdu.header += hdr
    hdulist = HDUList([hdu])
    hdulist.writeto(filename, clobber=clobber, output_verify='ignore')

def mad(d):
    return np.median(np.abs(np.median(d) - d))

def magtoflux(wave, mag, zp):
    # convert from ab mag to flambda
    # 3e-19 is lambda^2 / c in units of angstrom / Hz
    return zp * 10 ** (-0.4 * mag) / 3.33564095e-19 / wave / wave

def fluxtomag(flux):
    return -2.5 * np.log10(flux)

def spectoascii(infilename, outfilename):
    hdu = pyfits.open(infilename)
            
    lam = fitshdr_to_wave(hdu['SCI'].header.copy())
    flux = hdu['SCI'].data.copy()
    hdu.close()
    d = np.zeros((2, len(lam)))
    d[0] = lam
    d[1] = flux
    np.savetxt(outfilename, d.transpose())
    
def specsens(specfile, outfile, stdfile, extfile, airmass=None, exptime=None,
             stdzp=3.68e-20, thresh=8, clobber=True):

    # read in the specfile and create a spectrum object
    obs_hdu = pyfits.open(specfile)
    obs_flux = obs_hdu[2].data.copy()[0]
    obs_hdr = obs_hdu[2].header.copy()
    obs_hdu.close()
    obs_wave = fitshdr_to_wave(obs_hdr)
    # Mask out everything below 3350 where there is no signal
    obs_flux = obs_flux[obs_wave >= 3350.0]
    obs_wave = obs_wave[obs_wave >= 3350.0]
    # smooth the observed spectrum
    # read in the std file and convert from magnitudes to fnu
    # then convert it to fwave (ergs/s/cm2/A)
    std_wave, std_mag, _stdbnd = np.genfromtxt(stdfile).transpose()
    std_flux = magtoflux(std_wave, std_mag, stdzp)
 
    # Get the typical bandpass of the standard star,
    std_bandpass = np.diff(std_wave).mean()
    # Smooth the observed spectrum to that bandpass
    obs_flux = boxcar_smooth(obs_wave, obs_flux, std_bandpass)
    # read in the extinction file (leave in magnitudes)
    ext_wave, ext_mag = np.genfromtxt(extfile).transpose()

    # calculate the calibrated spectra
    cal_flux = cal_std(obs_wave, obs_flux, std_wave, std_flux, ext_wave,
                             ext_mag, airmass, exptime)

    
    # Normailize the fit variables so the fit is well behaved
    fitme_x = (obs_wave - obs_wave.min()) / (obs_wave.max() - obs_wave.min())
    fitme_y = cal_flux / np.median(cal_flux)
    coeffs = pfm.pffit(fitme_x, fitme_y, 15, 7, robust=True,
                       M=sm.robust.norms.AndrewWave())
 
    cal_flux = pfm.pfcalc(coeffs, fitme_x) * np.median(cal_flux)

    cal_mag = -1.0 * fluxtomag(cal_flux)
    # write the spectra out
    cal_hdr = sanitizeheader(obs_hdr.copy())
    cal_hdr['OBJECT'] = 'Sensitivity function for all apertures'
    cal_hdr['CRVAL1'] = obs_wave.min()
    tofits(outfile, cal_mag, hdr=cal_hdr, clobber=True)


def cal_std(obs_wave, obs_flux, std_wave, std_flux, ext_wave, ext_mag, airmass, exptime):
    """Given an observe spectra, calculate the calibration curve for the
       spectra.  All data is interpolated to the binning of the obs_spectra.
       The calibrated spectra is then calculated from
       C =  F_obs/ F_std / 10**(-0.4*A*E)/T/dW
       where F_obs is the observed flux from the source,  F_std  is the
       standard spectra, A is the airmass, E is the
       extinction in mags, T is the exposure time and dW is the bandpass
    Parameters
    -----------
    obs_spectra--spectrum of the observed star (counts/A)
    std_spectra--know spectrum of the standard star (ergs/s/cm2/A)
    ext_spectra--spectrum of the extinction curve (in mags)
    airmass--airmass of the observations
    exptime--exposure time of the observations
    function
    """
    
    # re-interpt the std_spectra over the same wavelength
    std_flux = np.interp(obs_wave, std_wave, std_flux)

    # re-interp the ext_spetra over the same wavelength
    ext_mag = np.interp(obs_wave, ext_wave, ext_mag)
    
    # create the calibration spectra
    # set up the bandpass
    bandpass = np.diff(obs_wave).mean()

    # correct for extinction
    cal_flux = obs_flux / 10 ** (-0.4 * airmass * ext_mag)

    # correct for the exposure time and calculation the sensitivity curve
    cal_flux = cal_flux / exptime / bandpass / std_flux

    return cal_flux

def boxcar_smooth(spec_wave, spec_flux, smoothwidth):
    # get the average wavelength separation for the observed spectrum
    # This will work best if the spectrum has equal linear wavelength spacings
    wavespace = np.diff(spec_wave).mean()
    # kw
    kw = int(smoothwidth / wavespace)
    # make sure the kernel width is odd
    if kw % 2 == 0:
        kw += 1
    kernel = np.ones(kw)
    # Conserve flux
    kernel /= kernel.sum()
    smoothed = spec_flux.copy()
    smoothed[(kw / 2):-(kw / 2)] = np.convolve(spec_flux, kernel, mode='valid')
    return smoothed

telluricWaves = [(2000., 3190.), (3216., 3420.), (5500., 6050.), (6250., 6360.),
                 (6450., 6530.), (6840., 7410.), (7550., 8410.), (8800., 9900.)]
def combine_spec_chi2(p, lam, specs, specerrs):
    # specs should be an array with shape (nspec, nlam)
    nspec = specs.shape[0]
    # scale each spectrum by the given value

    scaledspec = (specs.transpose() * p).transpose()
    scaled_spec_err = (specerrs.transpose() * p).transpose()

    chi = 0.0
    # loop over each pair of spectra
    for i in range(nspec):
        for j in range(i + 1, nspec):
            # Calculate the chi^2 for places that overlap
            # (i.e. spec > 0 in both)
            w = np.logical_and(scaledspec[i] != 0.0, scaledspec[j] != 0)
            if w.sum() > 0:
                residuals = scaledspec[i][w] - scaledspec[j][w]
                errs2 = scaled_spec_err[i][w] ** 2.0
                errs2 += scaled_spec_err[j][w] ** 2.0
                chi += (residuals ** 2.0 / errs2).sum()
    return chi

def speccombine(fs, outfile):
    nsteps = 8001
    lamgrid = np.linspace(3000.0, 11000.0, nsteps)

    nfs = len(fs)
    # for each aperture
    # get all of the science images
    specs = np.zeros((nfs, nsteps))
    specerrs = np.zeros((nfs, nsteps))

    for i, f in enumerate(fs):
        hdu = pyfits.open(f)
        lam = fitshdr_to_wave(hdu['SCI'].header.copy())
        # interpolate each spectrum onto a comman wavelength scale

        specs[i] = np.interp(lamgrid, lam, hdu['SCI'].data[0],
                          left=0.0, right=0.0)
        # Also calculate the errors. Right now we assume that the variances
        # interpolate linearly. This is not stricly correct but it should be
        # close. Also we don't include terms in the variance for the
        # uncertainty in the wavelength solution.
        specerrs[i] = 0.1 * specs[i]

    # minimize the chi^2 given free parameters are multiplicative factors
    # We could use linear or quadratic, but for now assume constant
    p0 = np.ones(nfs)

    results = optimize.minimize(combine_spec_chi2, p0,
                                args=(lamgrid, specs, specerrs),
                                method='Nelder-Mead',
                                options={'maxfev': 1e5, 'maxiter': 1e5})

    # write the best fit parameters into the headers of the files
    # Dump the list of spectra into a string that iraf can handle
    iraf_filelist = str(fs).replace('[', '').replace(']', '').replace("'", '').replace(',', '[SCI],')
    iraf_filelist += '[SCI]'
    print(iraf_filelist)
    # write the best fit results into a file
    lines = []
    for p in results['x']:
        lines.append('%f\n' % p)
    f = open('scales.dat', 'w')
    f.writelines(lines)
    f.close()
    # run scombine after multiplying the spectra by the best fit parameters
    if os.path.exists(outfile):
        os.remove(outfile)
    iraf.unlearn(iraf.scombine)
    iraf.scombine(iraf_filelist, outfile, scale='@scales.dat',
                  reject='avsigclip', lthreshold=1e-19, w1=3350)

def fitshdr_to_wave(hdr):
    crval = float(hdr['CRVAL1'])
    if 'CDELT1' in hdr.keys():
        cdelt = float(hdr['CDELT1'])
    else:
        cdelt = float(hdr['CD1_1'])
    nlam = float(hdr['NAXIS1'])
    lam = np.arange(crval, crval + cdelt * nlam - 1e-4, cdelt)
    return lam

def telluric_mask(waves):
    # True where not telluric contaminated
    not_telluric = np.ones(waves.shape, dtype=np.bool)
    for wavereg in telluricWaves:
        in_telluric_region = np.logical_and(waves >= wavereg[0],
                                            waves <= wavereg[1])
        not_telluric = np.logical_and(not_telluric,
                                         np.logical_not(in_telluric_region))
    return not_telluric

def mktelluric(filename):
    # if it is a standard star combined file
    # read in the spectrum and calculate the wavelengths of the pixels
    hdu = pyfits.open(filename)
    spec = hdu[0].data.copy()
    hdr = hdu[0].header.copy()
    hdu.close()
    waves = fitshdr_to_wave(hdr)
    
    template_spectrum = signal.savgol_filter(spec, 21, 3)
    noise = np.abs(spec - template_spectrum)
    noise = ndimage.filters.gaussian_filter1d(noise, 100.0)
    not_telluric = telluric_mask(waves)
        
    # Smooth the spectrum so that the spline doesn't go as crazy
    # Use the Savitzky-Golay filter to presevere the edges of the
    # absorption features (both atomospheric and intrinsic to the star)
    sgspec = signal.savgol_filter(spec, 31, 3)
    # Get the number of data points to set the smoothing criteria for the 
    # spline
    m = not_telluric.sum()
    intpr = interpolate.splrep(waves[not_telluric], sgspec[not_telluric],
                               w=1 / noise[not_telluric], k=2, s=20 * m)

    # Replace the telluric with the smoothed function
    smoothedspec = interpolate.splev(waves, intpr)
    
 
    # Extrapolate the ends linearly
    # Blue side
    w = np.logical_and(waves > 3420, waves < 3600)
    bluefit = np.poly1d(np.polyfit(waves[w], spec[w], 1))
    bluewaves = waves < 3420
    smoothedspec[bluewaves] = bluefit(waves[bluewaves])
     
    # Red side
    w = np.logical_and(waves > 8410, waves < 8800)
    redfit = np.poly1d(np.polyfit(waves[w], spec[w], 1))
    redwaves = waves > 8800
    smoothedspec[redwaves] = redfit(waves[redwaves])
    smoothedspec[not_telluric] = spec[not_telluric]
    # Divide the original and the telluric corrected spectra to
    # get the correction factor
    correction = spec / smoothedspec

    airmass = float(hdr['AIRMASS'])
    correction = correction ** (1.0 / airmass ** 0.55)
    # Save the correction
    dout = np.ones((2, len(waves)))
    dout[0] = waves
    dout[1] = correction
    np.savetxt('telcor.dat', dout.transpose())


def telluric(filename, outfile):

    # Get the standard to use for telluric correction
    stdfile = 'telcor.dat'
    
    hdu = pyfits.open(filename)
    spec = hdu[0].data.copy()
    hdr = hdu[0].header.copy()
    hdu.close()
    waves = fitshdr_to_wave(hdr)
    
    telwave, telspec = np.genfromtxt(stdfile).transpose()
    # Cross-correlate the standard star and the sci spectra
    # to find wavelength shift of standard star.
    w = np.logical_and(waves > 7550., waves < 8410.)
    tw = np.logical_and(telwave > 7550., telwave < 8410.)
    p = fitxcor(waves[w], spec[w], telwave[tw], telspec[tw])
    # shift and stretch standard star spectrum to match science
    # spectrum.
    telcorr = np.interp(waves, p[0] * telwave + p[1], telspec)

    # Correct for airmass
    airmass = float(hdr['AIRMASS'])
    telcorr = telcorr ** (airmass ** 0.55)

    # Divide science spectrum by transformed standard star sub-spectrum
    correct_spec = spec / telcorr

    # Copy telluric-corrected data to new file.
    tofits(outfile, correct_spec, hdr=hdr)

def ncor(x, y):
    """Calculate the normalized correlation of two arrays"""
    d = np.correlate(x, x) * np.correlate(y, y)

    return np.correlate(x, y) / d ** 0.5

def xcorfun(p, warr, farr, telwarr, telfarr):
    # Telluric wavelengths and flux
    # Observed wavelengths and flux
    # resample the telluric spectrum at the same wavelengths as the observed
    # spectrum
    # Make the artifical spectrum to cross correlate
    asfarr = np.interp(warr, p[0] * telwarr + p[1], telfarr, left=1.0, right=1.0)
    return abs(1.0 / ncor(farr, asfarr))

def fitxcor(warr, farr, telwarr, telfarr):
    """Maximize the normalized cross correlation coefficient for the telluric
    correction
    """
    res = optimize.minimize(xcorfun, [1.0, 0.0], method='Nelder-Mead',
                   args=(warr, farr, telwarr, telfarr))

    return res['x']

def sort():
    if not os.path.exists('raw'):
        os.mkdir('raw')
    fs = glob('*.fits')
    for f in fs:
        shutil.move(f, 'raw/')
    
    sensfs = glob('raw/sens*.fits')
    if len(sensfs) != 0:
        for f in sensfs:
            shutil.move(f, './')
    # Make a reduction directory
    if not os.path.exists('work'):
        os.mkdir('work')

    sensfs = glob('sens*.fits')
    if len(sensfs) != 0:
        for f in sensfs:
            shutil.copy(f, 'work/')
    
    if os.path.exists('telcor.dat'):
        shutil.copy('telcor.dat', 'work/')
    
    if os.path.exists('raw/bias.fits'):
        shutil.copy('raw/bias.fits', 'work/')
    
    # make a list of the raw files
    fs = glob('raw/*.fits')
    # Add a ../ in front of all of the file names
    for i in range(len(fs)):
        fs[i] = '../' + fs[i]
    return np.array(fs)

def init_northsouth(fs, topdir, rawpath):
    # Get the correct directory for the standard star
    base_stddir = 'spec50cal/'
    extfile = iraf.osfn('gmisc$lib/onedstds/kpnoextinct.dat') 
    observatory = 'Gemini-North'
    
    is_GS = pyfits.getval(fs[0], 'DETECTOR') == 'GMOS + Hamamatsu'
    if is_GS:
        if not os.path.exists(topdir + '/raw_fixhdr'):
            os.mkdir(topdir + '/raw_fixhdr')
        rawpath = '%s/raw_fixhdr/' % topdir
        os.system('gmoss_fix_headers.py --files="%s/raw/*.fits" --destination=%s' % (topdir, rawpath))
        base_stddir = 'ctionewcal/'
        observatory = 'Gemini-South'
        extfile = iraf.osfn('gmisc$lib/onedstds/ctioextinct.dat') 
    return extfile, observatory, base_stddir, rawpath

def getobstypes(fs):
    # get the type of observation for each file
    obstypes = []
    obsclasses = []
    for f in fs: 
        obstypes.append(pyfits.getval(f, 'OBSTYPE', ext=0))
        obsclasses.append(pyfits.getval(f, 'OBSCLASS', ext=0))
        
    obstypes = np.array(obstypes)
    obsclasses = np.array(obsclasses)
    return obstypes, obsclasses

def makebias(fs, obstypes, rawpath):
    if not os.path.exists('bias.fits'):
        # Make the master bias
        biastxtfile = open('bias.txt', 'w')
        biasfs = fs[obstypes == 'BIAS']
        for f in biasfs:
            biastxtfile.writelines(f.split('/')[-1] + '\n')
        biastxtfile.close()
        iraf.gbias('@%s/bias.txt' % os.getcwd(), 'bias', rawpath=rawpath, fl_over=dooverscan)

def getobjname(fs, obstypes):
    objname = pyfits.getval(fs[obstypes == 'OBJECT'][0], 'OBJECT', ext=0).lower()
    
    # get rid of nonsense in the name (like the plus and whitespace
    objname = objname.replace('+', '')
    objname = ''.join(objname.split())
    return objname

def maketxtfiles(fs, obstypes, obsclasses, objname):
    # go through each of the files (Ignore bias and aquisition files)
    goodfiles = np.logical_and(obsclasses != 'acqCal', obsclasses != 'acq')
    goodfiles = np.logical_and(goodfiles, obstypes != 'BIAS')
    
    for f in fs[goodfiles]:
        # put the filename in the correct text file.
        obsstr = ''
        obstype = pyfits.getval(f, 'OBSTYPE', ext=0)
        if obstype != 'OBJECT':
            obsstr = '.' + obstype.lower()
            expnum = ''
        else:
            expnum = 1
        
        # Drop the raw/
        fname = f.split('/')[-1]
        # red or blue setting
        redblue = pyfits.getval(f, 'GRATING')[0].lower()
        # central wavelength
        lamcentral = pyfits.getval(f, 'CENTWAVE')

        txtname = '%s.%s%s%i%s.txt' % (objname, str(expnum), redblue, lamcentral, obsstr)
        # If more than one arc or flat, append to the text file 
        if os.path.exists(txtname):
            if obsstr == '.flat' or obsstr == '.arc':
                # write to a text file
                txtfile = open(txtname, 'a')
            else:
                # We need to increment the exposure number
                moreimages = True
                expnum += 1
                while moreimages:
                    txtname = '%s.%s%s%i%s.txt' % (objname, str(expnum), redblue, lamcentral, obsstr)
                    if not os.path.exists(txtname):
                        txtfile = open(txtname, 'w')
                        moreimages = False
                    else:
                        expnum += 1
        else:
            txtfile = open(txtname, 'w')
        
        txtfile.write(fname + '\n')
        txtfile.close()

def gettxtfiles(fs, objname):
       
    flatfiles = np.array(glob('*.flat.txt'))
    
    # reduce the CuAr arcfiles.  Not flat fielded, gaps are not fixpixed
    arcfiles = np.array(glob('*.arc.txt'))

    # reduce the science files
    scifiles = glob(objname + '*.txt')
    
    nonscifiles = []
    # remove the arcs and flats
    for f in scifiles:
        if 'arc' in f or 'flat' in f: nonscifiles.append(f)
        
    for f in nonscifiles: 
        scifiles.remove(f)
    scifiles = np.array(scifiles)

    return flatfiles, arcfiles, scifiles

def makemasterflat(flatfiles, rawpath):
    # normalize the flat fields
    for f in flatfiles:
        iraf.unlearn(iraf.gsflat)
        # Use IRAF to get put the data in the right format and subtract the
        # bias
        iraf.gsflat('@' + f, f[:-4], order=51, rawpath=rawpath, fl_fixpix='yes',
                    bias="bias", fl_inter='no', niterate=3, low_reject=5.0, high_reject=5.0, fl_over=dooverscan, function='legendre')


def wavesol(arcfiles, rawpath):
    for f in arcfiles:
        iraf.unlearn(iraf.gsreduce)
        iraf.gsreduce('@' + f, outimages=f[:-4], rawpath=rawpath,
                      fl_flat=False, bias="bias", fl_fixpix=False,
                      fl_over=dooverscan)
        
        # determine wavelength calibration -- 1d and 2d
        iraf.unlearn(iraf.gswavelength)
        iraf.gswavelength(f[:-4], fl_inter='yes', low_reject=2.0,
                          high_reject=2.0, step=10, nsum=10, gsigma=3.0,
                          cradius=25.0, match=-12, order=7, fitcxord=7,
                          fitcyord=7)
    
        # transform the CuAr spectrum, for checking that the transformation is OK
        # output spectrum has prefix t
        iraf.unlearn(iraf.gstransform)
        iraf.gstransform(f[:-4], wavtran=f[:-4])
    
def getsetupname(f):
    # Get the setup base name by removing the exposure number
    return f.split('.')[0] + '.' + f.split('.')[1][1:]

def getredorblue(f):
    return f.split('.')[1][1]

def scireduce(scifiles, rawpath):
    for f in scifiles:
        setupname = getsetupname(f)
        # gsreduce subtracts bias, mosaics detectors, flat fields
        iraf.unlearn(iraf.gsreduce)
        iraf.gsreduce('@' + f, outimages=f[:-4], rawpath=rawpath, bias="bias",
                      flat=setupname + '.flat', fl_over=dooverscan)

        # Transform the data based on the arc  wavelength solution 
        iraf.unlearn(iraf.gstransform)
        iraf.gstransform(f[:-4], wavtran=setupname + '.arc')

def skysub(scifiles, rawpath):
    for f in scifiles:
        # sky subtraction
        # output has an s prefixed on the front
        # This step is currently quite slow for Gemini-South data
        iraf.unlearn(iraf.gsskysub)
        iraf.gsskysub('t' + f[:-4], long_sample='*', fl_inter='no',
                      naverage=-10, order=1, low_reject=2.0, high_reject=2.0,
                      niterate=10, mode='h')
    

def crreject(scifiles):
    for f in scifiles:
        # run lacosmicx
        hdu = pyfits.open('st' + f.replace('.txt', '.fits'))

        readnoise = 3.5
        # figure out what pssl should be approximately
        d = hdu[2].data.copy()
        dsort = np.sort(d.ravel())
        nd = dsort.shape[0]
        # Calculate the difference between the 16th and 84th percentiles to be 
        # robust against outliers
        dsig = (dsort[0.84 * nd] - dsort[0.16 * nd]) / 2.0
        pssl = (dsig * dsig - readnoise * readnoise)
        
        crmask, _cleanarr = lacosmicx.lacosmicx(d, sigclip=4.0, objlim=1.0, sigfrac=0.05,
                                     gain=1.0, readnoise=readnoise, pssl=pssl)
        
        tofits(f[:-4] + '.lamask.fits', np.array(crmask, dtype=np.uint8), hdr=hdu['SCI'].header.copy())

def fixpix(scifiles):
    # Run fixpix to interpolate over cosmic rays
    for f in scifiles:
        # run fixpix
        iraf.unlearn(iraf.fixpix)
        iraf.fixpix('t' + f[:-4] + '.fits[2]', f[:-4] + '.lamask.fits', mode='h')
        
def extract(scifiles):
    for f in scifiles:    
        iraf.unlearn(iraf.gsextract)
        # Extract the specctrum
        iraf.gsextract('t' + f[:-4], fl_inter='yes', bfunction='legendre',
                       border=2, bnaverage=-3, bniterate=2, blow_reject=2.0,
                       bhigh_reject=2.0, long_bsample='-100:-40,40:100',
                       background='fit', weights='variance',
                       lsigma=3.0, usigma=3.0, mode='h')

def makesensfunc(scifiles, objname, base_stddir, extfile):   
    for f in scifiles:
        redorblue = getredorblue(f)
        # If this is a standard star, run standard
        # Standards will have an observation class of either progCal or partnerCal
        obsclass = pyfits.getval(f[:-4] + '.fits', 'OBSCLASS')
        if obsclass == 'progCal' or obsclass == 'partnerCal':
            # Figure out which directory the stardard star is in
            stddir = iraf.osfn('gmisc$lib/onedstds/') + base_stddir
            
            # iraf.gsstandard('est' + f[:-4], 'std' + redorblue,
            #                'sens' + redorblue, starname=objname.lower(),
            #                caldir='gmisc$lib/onedstds/'+stddir, fl_inter=True)

            specsens('et' + f[:-4] + '.fits', 'sens' + redorblue + '.fits',
                     stddir + objname + '.dat' , extfile,
                     float(pyfits.getval(f[:-4] + '.fits', 'AIRMASS')),
                     float(pyfits.getval(f[:-4] + '.fits', 'EXPTIME')))

def calibrate(scifiles, extfile, observatory):
    for f in scifiles:
        redorblue = getredorblue(f)
        iraf.unlearn(iraf.gscalibrate)
        iraf.gscalibrate('et' + f[:-4] + '.fits',
                         sfunc='sens' + redorblue + '.fits', fl_ext=True,
                         extinction=extfile, observatory=observatory)
        
        if os.path.exists('cet' + f[:-4] + '.fits'):
            iraf.unlearn(iraf.splot)
            iraf.splot('cet' + f.replace('.txt', '.fits') + '[sci]')  # just to check
            
         
def updatecomheader(extractedfiles, objname):
    airmasses = []
    exptimes = []
    for f in extractedfiles:
        airmasses.append(float(pyfits.getval(f, 'AIRMASS')))
        exptimes.append(float(pyfits.getval(f, 'EXPTIME')))
    
    pyfits.setval(objname + '_com.fits', 'AIRMASS', value=np.mean(airmasses))
    pyfits.setval(objname + '_com.fits', 'SLIT', value=pyfits.getval(extractedfiles[0], 'MASKNAME').replace('arcsec', ''))
    
    comhdu = pyfits.open(objname + '_com.fits', mode='update')
    
    extractedhdu = pyfits.open(extractedfiles[0])
    for k in extractedhdu[0].header.keys():
        if not k in comhdu[0].header.keys():
            extractedhdu[0].header.cards[k].verify('fix')
            comhdu[0].header.append(extractedhdu[0].header.cards[k])
    
    comhdu.flush(output_verify='fix')
    comhdu.close()
    extractedhdu.close()
    dateobs = pyfits.getval(objname + '_com.fits', 'DATE-OBS')
    dateobs += 'T' + pyfits.getval(objname + '_com.fits', 'TIME-OBS')
    pyfits.setval(objname + '_com.fits', 'DATE-OBS', value=dateobs)

def cleanfinal(filename):
    # Clean the data of infs and nans
    hdu = pyfits.open(filename, mode='update')
    hdu[0].data[np.isnan(hdu[0].data)] = 0.0
    hdu[0].data[np.isinf(hdu[0].data)] = 0.0
    hdu.flush()
    hdu.close()
    
    
def rescale1e15(filename):
    hdu = pyfits.open(filename, mode='update')
    hdu[0].data *= 1e-15
    hdu.flush()
    hdu.close()
    
    
if __name__ == "__main__":
    # copy over sensr.fits, sensb.fits files
    # before running this script
    
    # launch the image viewer
    os.system('ds9 &')
    
    topdir = os.getcwd()
    # Get the raw directory
    rawpath = '%s/raw/' % topdir
    
    # Sort the files into the correct directories
    fs = sort()
    
    # Change into the reduction directory
    iraf.chdir('work')

    # Initialize variables that depend on which site was used
    extfile, observatory, base_stddir, rawpath = init_northsouth(fs, topdir, rawpath)
       
    # Get the observation type
    obstypes, obsclasses = getobstypes(fs)
    
    # Make the bias frame
    makebias(fs, obstypes, rawpath)
    
    # get the object name
    objname = getobjname(fs, obstypes)
    
    # Make the text files for the IRAF tasks
    maketxtfiles(fs, obstypes, obsclasses, objname)
                
    # remember not to put ".fits" on the end of filenames!
    flatfiles, arcfiles, scifiles = gettxtfiles(fs, objname)
    
    # Make the master flat field image
    makemasterflat(flatfiles, rawpath)
    
    # Get the wavelength solution
    wavesol(arcfiles, rawpath)
    
    # Flat field and rectify the scienc images
    scireduce(scifiles, rawpath)
    
    # Run sky subtraction
    skysub(scifiles, rawpath)
    
    # Run LA Cosmic
    crreject(scifiles)
    
    # Fix the cosmic ray pixels
    fixpix(scifiles)
    
    # Extract the 1D spectrum
    extract(scifiles)
    
    # If standard star, make the sensitivity function
    makesensfunc(scifiles, objname, base_stddir, extfile)
    
    # Flux calibrate the spectrum
    calibrate(scifiles, extfile, observatory)
    
    # Write the spectra to ascii
    for f in scifiles:
        if os.path.exists('cet' + f[:-4] + '.fits'):
            # Make the ascii file            
            spectoascii('cet' + f[:-4] + '.fits',f[:-4] + '.dat')
            
    # Get all of the extracted files
    extractedfiles = glob('cet*.fits')
    
    # Combine the spectra
    speccombine(extractedfiles, objname + '_com.fits')
    
    # write out the ascii file
    spectoascii(objname + '_com.fits', objname + '_com.dat')
    
    # Update the combined file with the necessary header keywords
    updatecomheader(extractedfiles, objname)
    
    # If a standard star, make a telluric correction
    obsclass = pyfits.getval(objname + '_com.fits', 'OBSCLASS')
    if obsclass == 'progCal' or obsclass == 'partnerCal':
        # Make telluric
        mktelluric(objname + '_com.fits')
        
    # Telluric Correct
    telluric(objname + '_com.fits', objname + '.fits')
    
    #Clean the data of nans and infs
    cleanfinal(objname + '.fits')
    
    # Write out the ascii file
    spectoascii(objname + '.fits', objname + '.dat')
    
    # Multiply by 1e-15 so the units are correct in SNEx:
    rescale1e15(objname + '.fits')
    
    # Change out of the reduction directory
    iraf.chdir('..')
