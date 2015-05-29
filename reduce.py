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
from astropy.modeling import models, fitting

iraf.cd(os.getcwd())
iraf.gemini()
iraf.gmos()
iraf.onedspec()


iraf.gmos.logfile = "log.txt"
iraf.gmos.mode = 'h'
iraf.set(clobber='yes')

iraf.set(stdimage='imtgmos')

dooverscan = False
is_GS = False

def normalize_fitting_coordinate(x):
    xrange = x.max() - x.min()
    return (x - x.min()) / xrange


# Iterative reweighting linear least squares
def irls(x, data, errors, model, tol=1e-6, M=sm.robust.norms.AndrewWave(), maxiter=10):
    fitter = fitting.LinearLSQFitter()

    if x is None:
        # Make x and y arrays out of the indicies
        x = np.indices(data.shape, dtype=np.float)

        if len(data.shape) == 2:
            y, x = x
        else:
            x = x[0]

        #Normalize to make fitting easier
        x = normalize_fitting_coordinate(x)
        if len(data.shape) == 2:
            y = normalize_fitting_coordinate(y)

    scatter = errors
    # Do an initial fit of the model
    # Use 1 / sigma^2 as weights
    weights = (errors ** -2.0).flatten()

    if len(data.shape) == 2:
        fitted_model = fitter(model, x, y, data, weights=weights)
    else:
        fitted_model = fitter(model, x, data, weights=weights)

    notconverged=True
    last_chi = np.inf
    iter = 0
    # Until converged
    while notconverged:
        # Update the weights
        if len(data.shape) == 2:
            residuals = data - fitted_model(x, y)
        else:
            residuals = data - fitted_model(x)
        # Save the chi^2 to check for convergence
        chi = ((residuals / scatter) ** 2.0).sum()

        # update the scaling (the MAD of the residuals)
        scatter = mad(residuals)  * 1.4826 # To convert to standard deviation
        weights = M.weights(residuals / scatter).flatten()

        # refit
        if len(data.shape) == 2:
            fitted_model = fitter(model, x, y, data, weights=weights)
        else:
            fitted_model = fitter(model, x, data, weights=weights)
        # converged when the change in the chi^2 (or l2 norm or whatever) is
        # less than the tolerance. Hopefully this should converge quickly.
        if iter >= maxiter or np.abs(chi - last_chi) < tol:
            notconverged = False
        else:
            last_chi = chi
            iter += 1

    return fitted_model

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
    try:
        lam = fitshdr_to_wave(hdu['SCI'].header.copy())
        flux = hdu['SCI'].data.copy()
    except:
        lam = fitshdr_to_wave(hdu[0].header.copy())
        flux = hdu[0].data.copy()
    hdu.close()
    d = np.zeros((2, len(lam)))
    d[0] = lam
    d[1] = flux
    np.savetxt(outfilename, d.transpose())
    
def specsens(specfile, outfile, stdfile, extfile, airmass=None, exptime=None,
             stdzp=3.68e-20, thresh=8, clobber=True):

    # read in the specfile and create a spectrum object
    obs_hdu = pyfits.open(specfile)
    try:
        obs_flux = obs_hdu[2].data.copy()[0]
        obs_hdr = obs_hdu[2].header.copy()
    except:
        obs_flux = obs_hdu[0].data.copy()
        obs_hdr = obs_hdu[0].header.copy()
    obs_hdu.close()
    obs_wave = fitshdr_to_wave(obs_hdr)

    # Mask out everything below 3350 where there is no signal
    obs_flux = obs_flux[obs_wave >= 3350.0]
    obs_wave = obs_wave[obs_wave >= 3350.0]

    # Figure out where the chip gaps are
    chip_edges = get_chipedges(obs_flux)


    try:
        chip_gaps = np.ones(obs_flux.size, dtype=np.bool)
        for edge in chip_edges:
            chip_gaps[edge[0]: edge[1]] = False
    except:
        chip_gaps = np.zeros(obs_flux.size, dtype=np.bool)

    template_spectrum = signal.savgol_filter(obs_flux, 21, 3)
    noise = np.abs(obs_flux - template_spectrum)
    noise = ndimage.filters.gaussian_filter1d(noise, 100.0)

    if chip_gaps.sum() != len(chip_gaps):
        # Smooth the chip gaps
        intpr = interpolate.splrep(obs_wave[np.logical_not(chip_gaps)],
                                   obs_flux[np.logical_not(chip_gaps)],
                                   w=1 / noise[np.logical_not(chip_gaps)], k=2,
                                   s=20 * np.logical_not(chip_gaps).sum())
        obs_flux[chip_gaps] = interpolate.splev(obs_wave[chip_gaps], intpr)
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

    # Normalize the fit variables so the fit is well behaved
    fitme_x = (obs_wave - obs_wave.min()) / (obs_wave.max() - obs_wave.min())
    fitme_y = cal_flux / np.median(cal_flux)
    coeffs = pfm.pffit(fitme_x, fitme_y, 5 , 7, robust=True,
                    M=sm.robust.norms.AndrewWave())

    fitted_flux = pfm.pfcalc(coeffs, fitme_x) * np.median(cal_flux)

    cal_mag = -1.0 * fluxtomag(fitted_flux)
    # write the spectra out
    cal_hdr = sanitizeheader(obs_hdr.copy())
    cal_hdr['OBJECT'] = 'Sensitivity function for all apertures'
    cal_hdr['CRVAL1'] = obs_wave.min()
    cal_hdr['CRPIX1'] = 1
    tofits(outfile, cal_mag, hdr=cal_hdr, clobber=True)


def hdr_pixel_range(x0, x1, y0, y1):
    return '[{0:d}:{1:d},{2:d}:{3:d}]'.format(x0, x1, y0, y1)

def cut_gs_image(filename, output_filename, pixel_range):
    """

    :param filename:
    :param output_filename:
    :param pixel_range: array-like, The range of pixels to keep, python indexed,
                        given in binned pixels
    :return:
    """
    hdu = pyfits.open(filename, unit16=True)
    for i in range(1, 13):
        ccdsum = hdu[i].header['CCDSUM']
        ccdsum = np.array(ccdsum.split(), dtype=np.int)

        y_ccdsec = [(pixel_range[0]  * ccdsum[1]) + 1,
                    (pixel_range[1]) * ccdsum[1]]

        detsec = hdr_pixel_range(512 * (i - 1) + 1, 512 * i,
                                 y_ccdsec[0], y_ccdsec[1])
        hdu[i].header['DETSEC'] = detsec

        ccdsec = hdr_pixel_range(512 * ((i - 1) % 4) + 1, 512 * ((i - 1) % 4 + 1),
                                 y_ccdsec[0], y_ccdsec[1])

        hdu[i].header['CCDSEC'] = ccdsec

        numpix = pixel_range[1] - pixel_range[0]

        xsize = 512 / ccdsum[0]
        # Add a 4 pixel buffer to the overscan region because of bleed over
        if i % 2 == 1:
            hdu[i].header['BIASSEC'] = hdr_pixel_range(xsize + 5, xsize + 32, 1, numpix)
            hdu[i].header['DATASEC'] = hdr_pixel_range(1, xsize, 1, numpix)
        if i % 2 == 0:
            hdu[i].header['BIASSEC'] = hdr_pixel_range(1, 28, 1, numpix)
            hdu[i].header['DATASEC'] = hdr_pixel_range(33, xsize + 32, 1, numpix)

        hdu[i].data = hdu[i].data[pixel_range[0]:pixel_range[1], :]

    hdu.writeto(output_filename)
    hdu.close()

def get_chipedges(data):
        # Get the x coordinages of all of the chip gap pixels
        # recall that pyfits opens images with coordinates y, x
        if len(data.shape) > 1:
            data = data[0]

        try:
            w = np.where(data == 0.0)[0]

            # Note we also grow the chip gap by 8 pixels on each side
            # Chip 1
            chip_edges = []
            left_chipedge = 10
            morechips = True
            while morechips:
                try:
                    right_chipedge = np.min(w[w > left_chipedge]) - 10
                except:
                    right_chipedge = data.size - 10
                    morechips = False
                chip_edges.append((left_chipedge, right_chipedge))

                left_chipedge = np.max(w[w < right_chipedge + 100]) + 10
        except:
            chip_edges = []
        return chip_edges


def split1d(filename):

    hdu = pyfits.open(filename)
    chipedges = get_chipedges(hdu['SCI'].data[0])
    lam = fitshdr_to_wave(hdu['SCI'].header)
    # Copy each of the chips out seperately. Note that iraf is 1 indexed
    for i in range(3):
        # get the wavelengths that correspond to each chip
        w1 = lam[chipedges[i][0]]
        w2 = lam[chipedges[i][1]]
        iraf.scopy(filename+ '[SCI]', output=filename[:-5] + 'c%i.fits' % (i + 1), w1=w1,
                   w2=w2, rebin='no')
        hdu.close()


def mask_chipedges(filename):
    """
    Mask the edges of the chips with zeros to minimize artifacts.
    :param filename: Name of file that contains the spectrum
    :return:
    """
    hdu = pyfits.open(filename, mode='update')
    chip_edges = get_chipedges(hdu['SCI'].data[0])
    print(chip_edges)
    # Set the data = 0 in the chip gaps
    # Assume 3 chips for now.
    for i in range(2):
        hdu['SCI'].data[0, chip_edges[i][1]:chip_edges[i+1][0]] = 0.0

    hdu.flush()
    hdu.close()


def cal_std(obs_wave, obs_flux, std_wave, std_flux, ext_wave, ext_mag, airmass, exptime):
    """Given an observe spectra, calculate the calibration curve for the
       spectra.  All data is interpolated to the binning of the obs_spectra.
       The calibrated spectra is then calculated from
       C =  F_obs/ F_std / 10**(-0.4*A*E)/T/dW
       where F_obs is the observed flux from the source,  F_std  is the
       standard spectra, A is the airmass, E is the
       extinction in mags, T is the exposure time and dW is the bandpass
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
    # Assume 3 chips here
    scales = np.repeat(p, 3)

    scaledspec = (specs.transpose() * scales).transpose()
    scaled_spec_err = (specerrs.transpose() * scales).transpose()

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
        lam = fitshdr_to_wave(hdu[0].header.copy())

        # interpolate each spectrum onto a common wavelength scale

        specs[i] = np.interp(lamgrid, lam, hdu[0].data,
                             left=0.0, right=0.0)
        # Also calculate the errors. Right now we assume that the variances
        # interpolate linearly. This is not strictly correct but it should be
        # close. Also we don't include terms in the variance for the
        # uncertainty in the wavelength solution.
        specerrs[i] = 0.1 * specs[i]

    # minimize the chi^2 given free parameters are multiplicative factors
    # We could use linear or quadratic, but for now assume constant
    # Assume 3 chips for now
    p0 = np.ones(nfs / 3)

    results = optimize.minimize(combine_spec_chi2, p0,
                                args=(lamgrid, specs, specerrs),
                                method='Nelder-Mead',
                                options={'maxfev': 1e5, 'maxiter': 1e5, 'ftol':1e-5})

    # write the best fit parameters into the headers of the files
    # Dump the list of spectra into a string that iraf can handle
    iraf_filelist = str(fs).replace('[', '').replace(']', '').replace("'", '') #.replace(',', '[SCI],')
    #iraf_filelist += '[SCI]'

    # write the best fit results into a file
    lines = []

    for p in np.repeat(results['x'], 3):
        lines.append('%f\n' % (1.0 / p))
    f = open('scales.dat', 'w')
    f.writelines(lines)
    f.close()
    # run scombine after multiplying the spectra by the best fit parameters
    if os.path.exists(outfile):
        os.remove(outfile)
    iraf.unlearn(iraf.scombine)
    iraf.scombine(iraf_filelist, outfile, scale='@scales.dat',
                  reject='avsigclip', lthreshold=1e-4, w1=3350)

def fitshdr_to_wave(hdr):
    crval = float(hdr['CRVAL1'])
    crpix = float(hdr['CRPIX1'])
    # Convert crpix to be zero indexed
    crpix -= 1
    if 'CDELT1' in hdr.keys():
        cdelt = float(hdr['CDELT1'])
    else:
        cdelt = float(hdr['CD1_1'])
    npix = float(hdr['NAXIS1'])
    lam = np.arange(crval - cdelt * crpix ,
                    crval + cdelt * (npix - crpix) - 1e-4,
                    cdelt)
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

    # Start by interpolating over the chip gaps
    chip_edges = get_chipedges(spec)
    chip_gaps = np.ones(spec.size, dtype=np.bool)
    for edge in chip_edges:
        chip_gaps[edge[0]: edge[1]] = False

    template_spectrum = signal.savgol_filter(spec, 21, 3)
    noise = np.abs(spec - template_spectrum)
    noise = ndimage.filters.gaussian_filter1d(noise, 100.0)

    # Smooth the chip gaps
    intpr = interpolate.splrep(waves[np.logical_not(chip_gaps)],
                               spec[np.logical_not(chip_gaps)],
                               w=1 / noise[np.logical_not(chip_gaps)], k=2,
                               s=20 * np.logical_not(chip_gaps).sum())
    spec[chip_gaps] = interpolate.splev(waves[chip_gaps], intpr)

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
    return np.abs(1.0 / ncor(farr, asfarr))

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

    global is_GS
    is_GS = pyfits.getval(fs[0], 'DETECTOR') == 'GMOS + Hamamatsu'
    if is_GS:
        global dooverscan
        dooverscan = True
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
    for f in fs:
        if f[-10:] == '_bias.fits':
            shutil.copy(f, 'bias.fits')

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
        # This will currently break if multiple flats are used for a single setting
        iraf.unlearn(iraf.gsreduce)
        iraf.gsreduce('@' + f, outimages = f[:-4]+'.mef.fits',rawpath=rawpath,
                      bias="bias", fl_over=dooverscan, fl_flat=False, fl_gmosaic=False,
                      fl_fixpix=False, fl_gsappwave=False, fl_cut=False, fl_title=False,
                      fl_oversize=False)

        # Divide the flat by the median over the spatial direction
        hdunorm = pyfits.open(f[:-4] + '.mef.fits')

        for i in range(1, 13):
            meddata = np.median(hdunorm[i].data, axis=0)

            # Fit out the ramp up on the edge of the chip
            x = np.arange(len(meddata), dtype=np.float)
            x /= x.max()
            if i in [1, 5, 9]:
                fit = pfm.pffit(x[50:100], meddata[50:100], 0, 1,
                                robust=True, M=sm.robust.norms.AndrewWave())
                meddata[:50] = pfm.pfcalc(fit, x[:50])
            elif i in [4, 8, 12]:
                fit = pfm.pffit(x[-100:-50], meddata[-100:-50], 0, 1,
                                robust=True, M=sm.robust.norms.AndrewWave())
                meddata[-50:] = pfm.pfcalc(fit, x[-50:])

            hdunorm[i].data /= meddata


        # Now we want to fit the jumps between chips
        hdu = pyfits.open(f[:-4] + '.mef.fits')
        # Divide out the pixel to pixel sensitivity
        for i in range(1, 13):
            hdu[i].data /= hdunorm[i].data
        hdu.writeto(f[:-4]+'.flt.fits', clobber=True)
        hdu.close()

        # Mosaic the combined flat
        iraf.unlearn(iraf.gmosaic)
        iraf.gmosaic(f[:-4] + '.flt.fits', outimages=f[:-4] + '.mos.fits', fl_clean=False)

        # Renormalize the chips to remove the discrete jump in the
        # sensitivity due to differences in the QE for different chips

        flat_hdu = pyfits.open(f[:-4] + '.mos.fits')

        data = np.median(flat_hdu['SCI'].data, axis=0)
        chip_edges = get_chipedges(data)
        fitme_x = np.arange(len(data), dtype=np.float)
        fitme_x /= fitme_x.max()
        fitme_y = data / np.median(data)

        leftchip = slice(chip_edges[0][1] - 200, chip_edges[0][1])
        leftchip_fit = pfm.pffit(fitme_x[leftchip], fitme_y[leftchip], 0 , 1,
                                 robust=True, M=sm.robust.norms.AndrewWave())

        # Normalize to the middle chip for now
        left_middlechip = slice(chip_edges[1][0], chip_edges[1][0] + 150)

        left_middlechip_fit = pfm.pffit(fitme_x[left_middlechip], fitme_y[left_middlechip], 0 , 1, robust=True,
                                   M=sm.robust.norms.AndrewWave())


        # Normalize to the middle chip for now
        rightchip = slice(chip_edges[2][0], chip_edges[2][0] + 200)
        rightchip_fit = pfm.pffit(fitme_x[rightchip], fitme_y[rightchip], 0 , 1, robust=True,
                                  M=sm.robust.norms.AndrewWave())

        # Then normalize the right chip to the middle chip
        right_middlechip = slice(chip_edges[1][1] - 200, chip_edges[1][1])
        right_middlechip_fit = pfm.pffit(fitme_x[right_middlechip], fitme_y[right_middlechip], 0 , 1,
                                   robust=True, M=sm.robust.norms.AndrewWave())

        # Save the normalized flat as a MEF file
        hdunorm.writeto(f[:-4]+'.norm.mef.fits', clobber=True)

        hdunorm.close()
        # mosaic the normalized file
        iraf.unlearn(iraf.gmosaic)
        iraf.gmosaic(f[:-4] + '.norm.mef.fits', outimages=f[:-4] + '.fits', fl_clean=False)
        #print(leftchip_fit, left_middlechip_fit)
        #print(rightchip_fit, right_middlechip_fit)
        #from matplotlib import pyplot
        #pyplot.plot(fitme_x, fitme_y)
        #pyplot.plot(fitme_x, pfm.pfcalc(leftchip_fit, fitme_x))
        #pyplot.plot(fitme_x, pfm.pfcalc(left_middlechip_fit, fitme_x))
        #pyplot.plot(fitme_x, pfm.pfcalc(rightchip_fit, fitme_x))
        #pyplot.plot(fitme_x, pfm.pfcalc(right_middlechip_fit, fitme_x))
        #pyplot.plot(fitme_x, pfm.pfcalc(rightchip_fit, fitme_x) / pfm.pfcalc(right_middlechip_fit, fitme_x))
        #pyplot.show()
        hdunorm = pyfits.open(f[:-4]+'.fits', mode='update')
        # Renormalize the left and right chips from the fits above.
        # Get the value of the fits at the left edge of the middle chip
        leftedge = chip_edges[1][0]
        #if pfm.pfcalc(leftchip_fit, fitme_x[leftedge]) <  pfm.pfcalc(left_middlechip_fit, fitme_x[leftedge]):
        #     # Raise the left chip
        #     leftx = fitme_x[:chip_edges[0][1]]
        #     #leftcorrection = pfm.pfcalc(leftchip_fit, leftx) / pfm.pfcalc(left_middlechip_fit, leftx)
        #     leftcorrection = pfm.pfcalc(leftchip_fit, fitme_x[leftedge]) / pfm.pfcalc(left_middlechip_fit, fitme_x[leftedge])
        #     hdunorm['SCI'].data[:,:chip_edges[0][1]] *= leftcorrection
        # else:
        #     middlex = fitme_x[chip_edges[1][0]:chip_edges[1][1]+1]
        #     #middlecorrection = pfm.pfcalc(left_middlechip_fit, middlex) / pfm.pfcalc(leftchip_fit, middlex)
        #     middlecorrection = pfm.pfcalc(left_middlechip_fit, fitme_x[leftedge]) / pfm.pfcalc(leftchip_fit, fitme_x[leftedge])
        #
        #     hdunorm['SCI'].data[:,chip_edges[1][0]:chip_edges[1][1]+1] /= middlecorrection

        leftcorrection = pfm.pfcalc(leftchip_fit, fitme_x[leftedge]) / pfm.pfcalc(left_middlechip_fit, fitme_x[leftedge])
        hdunorm['SCI'].data[:,:chip_edges[0][1]] *= leftcorrection
        rightedge = chip_edges[1][1]
        # if pfm.pfcalc(rightchip_fit, fitme_x[rightedge]) <  pfm.pfcalc(right_middlechip_fit, fitme_x[rightedge]):
        #     rightx = fitme_x[chip_edges[2][0]:]
        #     #rightcorrection = pfm.pfcalc(rightchip_fit, rightx) / pfm.pfcalc(right_middlechip_fit, rightx)
        #     rightcorrection = pfm.pfcalc(rightchip_fit, fitme_x[rightedge]) / pfm.pfcalc(right_middlechip_fit, fitme_x[rightedge])
        #     hdunorm['SCI'].data[:,chip_edges[2][0]:] *= rightcorrection
        # else:
        #     # Raise both left and middle
        #     leftcenterx = fitme_x[0: chip_edges[1][1] + 1]
        #     #leftcorrection = pfm.pfcalc(right_middlechip_fit, leftcenterx) / pfm.pfcalc(rightchip_fit, leftcenterx)
        #     leftcorrection = pfm.pfcalc(right_middlechip_fit, fitme_x[rightedge]) / pfm.pfcalc(rightchip_fit, fitme_x[rightedge])
        #     hdunorm['SCI'].data[:,:chip_edges[1][1]+1] *= leftcorrection
        #
        rightcorrection = pfm.pfcalc(rightchip_fit, fitme_x[rightedge]) / pfm.pfcalc(right_middlechip_fit, fitme_x[rightedge])
        hdunorm['SCI'].data[:,chip_edges[2][0]:] *= rightcorrection

        # Replace the chip gaps with ones
        hdunorm['SCI'].data[:, chip_edges[0][1]+1:chip_edges[1][0]] = 1
        hdunorm['SCI'].data[:, chip_edges[1][1]+1:chip_edges[2][0]] = 1
        # Save the final normalized flat
        hdunorm.flush()
        hdunorm.close()

def wavesol(arcfiles, rawpath):
    for f in arcfiles:
        iraf.unlearn(iraf.gsreduce)
        iraf.gsreduce('@' + f, outimages=f[:-4], rawpath=rawpath,
                      fl_flat=False, bias="bias",
                      fl_fixpix=False, fl_over=dooverscan, fl_cut=False, fl_gmosaic=True,
                      fl_gsappwave=True, fl_oversize=False)

        # Divide out the flat file
        flat_hdu = pyfits.open(f[:-4].replace('arc', 'flat')+'.fits')
        arc_hdu = pyfits.open(f[:-4]+'.fits', mode='update')
        arc_hdu['SCI'].data /= flat_hdu['SCI'].data
        arc_hdu.flush()
        arc_hdu.close()
        flat_hdu.close()

        # determine wavelength calibration -- 1d and 2d
        iraf.unlearn(iraf.gswavelength)
        iraf.gswavelength(f[:-4], fl_inter='yes', fwidth=15.0, low_reject=2.0,
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
        # gsreduce subtracts bias and mosaics detectors
        iraf.unlearn(iraf.gsreduce)
        iraf.gsreduce('@' + f, outimages=f[:-4], rawpath=rawpath, bias="bias",
                      fl_over=dooverscan, fl_fixpix='no', fl_flat=False,
                      fl_gmosaic=True, fl_cut=False, fl_gsappwave=False, fl_oversize=False)

        # Divide out the flat file
        flat_hdu = pyfits.open(setupname + '.flat.fits')
        sci_hdu = pyfits.open(f[:-4]+'.fits', mode='update')
        sci_hdu['SCI'].data /= flat_hdu['SCI'].data
        sci_hdu.flush()
        sci_hdu.close()
        flat_hdu.close()

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

        mask = d == 0.0
        crmask, _cleanarr = lacosmicx.lacosmicx(d, inmask=mask, sigclip=4.0,
                                                objlim=1.0, sigfrac=0.05, gain=1.0,
                                                readnoise=readnoise, pssl=pssl)
        
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
    # os.system('ds9 &')
    
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

    extractedfiles = glob('cet*.fits')

    # Write the spectra to ascii
    for f in scifiles:
        if os.path.exists('cet' + f[:-4] + '.fits'):
            split1d('cet' + f[:-4] + '.fits')
            # Make the ascii file
            spectoascii('cet' + f[:-4] + '.fits',f[:-4] + '.dat')
            
    # Get all of the extracted files
    splitfiles = glob('cet*c[1-9].fits')
    
    # Combine the spectra
    speccombine(splitfiles, objname + '_com.fits')
    
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
