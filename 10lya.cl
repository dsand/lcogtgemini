########## in a terminal 
copy over sensr.fits, sensb.fits, and bias.fits files

ds9 & #launch image viewer

########## in iraf:

gemini
gmos

# cd to the data reduction directory (not the raw)

gmos.logfile="log.txt"
set stdimage=imt4096
set raw=/Users/ahowell/data/ptf/gemini/GS10A/10lya/raw/

# I can't get gdisplay to work, so I don't set stdimage=imtgmos

#determine what's what
hselect raw/*.fits[0] $I,OBJECT,DATE_OBS,MASKNAME,FILTER1,GRATING,EXPTIME,CENTWAVE,OBSTYPE,OBSCLASS

# edit this to give the right files and copy it into a terminal 
#echo S20100628S0137.fits >  ptf10lya.b450.txt
#echo S20100628S0136.fits >  ptf10lya.b450.flat.txt
echo S20100629S0091.fits >  ptf10lya.1b450.txt
echo S20100629S0092.fits >  ptf10lya.2b450.txt
echo S20100629S0093.fits >  ptf10lya.1b455.txt
echo S20100629S0094.fits >  ptf10lya.2b455.txt
echo S20100629S0083.fits >  ptf10lya.1r750.txt
echo S20100629S0084.fits >  ptf10lya.2r750.txt
echo S20100629S0087.fits >  ptf10lya.1r755.txt
echo S20100629S0088.fits >  ptf10lya.2r755.txt

echo S20100629S0090.fits >  ptf10lya.1b450.flat.txt
echo S20100629S0109.fits >  ptf10lya.1b450.arc.txt
echo S20100629S0095.fits >  ptf10lya.1b455.flat.txt
echo S20100629S0110.fits >  ptf10lya.1b455.arc.txt
echo S20100629S0086.fits >  ptf10lya.r750.flat.txt
echo S20100629S0107.fits >  ptf10lya.r750.arc.txt
echo S20100629S0089.fits >  ptf10lya.r755.flat.txt
echo S20100629S0108.fits >  ptf10lya.r755.arc.txt

#note arcs used in this reduction taken on different day

#remember not to put ".fits" on the end of filenames below!

##### do in iraf, in reduction (not raw) directory
######## the object spectrum B450 ###############

#normalize the flat field
gsflat @ptf10lya.1b450.flat.txt ptf10lya.1b450.flat order=23 rawpath="raw$" bias="bias" fl_inter-

#reduce CuAr.  Not flat fielded, gaps are not fixpixed
gsreduce @ptf10lya.1b450.arc.txt outimages="ptf10lya.1b450.arc" rawpath="raw$" fl_flat- bias="bias" fl_fixpix-

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.1b450.txt outimages="ptf10lya.1b450" rawpath="raw$" bias="bias" flat="ptf10lya.1b450.flat"

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.2b450.txt outimages="ptf10lya.2b450" rawpath="raw$" bias="bias" flat="ptf10lya.1b450.flat"

#determine wavelength calibration -- 1d and 2d
gswavelength ptf10lya.1b450.arc fl_inter+

#transform the CuAr spectrum, for checking that the transformation is OK
#output spectrum has prefix t
#gstransform ptf10lya.1b450.arc wavtran=ptf10lya.1b450.arc 

#transform the science spectrum
gstransform ptf10lya.1b450 wavtran=ptf10lya.1b450.arc 

#transform the science spectrum
gstransform ptf10lya.2b450 wavtran=ptf10lya.1b450.arc 

display tptf10lya.1b450.fits[sci,1] 1
display tptf10lya.2b450.fits[sci,1] 1

#sky subtraction: the sample is manually selected after inspection of the image
# output has an s prefixed to the front
gsskysub tptf10lya.1b450 long_sample="130:220,280:360" fl_inter-

gsskysub tptf10lya.2b450 long_sample="130:220,280:360" fl_inter-

#examine the sky subtracted spectrum
# ds9 stptf10lya.1b450.fits & # do this in a terminal
display stptf10lya.1b450[sci,1]  # or do this in iraf

#Extract the specctrum
gsextract stptf10lya.1b450 fl_inter+

gsextract stptf10lya.2b450 fl_inter+

gscalibrate estptf10lya.1b450  sfunc=sensb

gscalibrate estptf10lya.2b450  sfunc=sensb

oned # to load oned spectral reduction iraf package

splot cestptf10lya.1b450.fits  # just to check

wspectext cestptf10lya.1b450.fits[sci,1] ptf10lya.1blue.dat header-

wspectext cestptf10lya.2b450.fits[sci,1] ptf10lya.2blue.dat header-

######## the object spectrum R750 ###############

#normalize the flat field
gsflat @ptf10lya.r750.flat.txt ptf10lya.r750.flat order=23 rawpath="raw$" bias="bias" fl_inter-

#reduce CuAr.  Not flat fielded, gaps are not fixpixed
gsreduce @ptf10lya.r750.arc.txt outimages="ptf10lya.r750.arc" rawpath="raw$" fl_flat- bias="bias" fl_fixpix-

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.1r750.txt outimages="ptf10lya.1r750" rawpath="raw$" bias="bias" flat="ptf10lya.r750.flat"

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.2r750.txt outimages="ptf10lya.2r750" rawpath="raw$" bias="bias" flat="ptf10lya.r750.flat"

#determine wavelength calibration -- 1d and 2d
gswavelength ptf10lya.r750.arc fl_inter+

#transform the CuAr spectrum, for checking that the transformation is OK
#output spectrum has prefix t
#gstransform ptf10lya.r750.arc wavtran=ptf10lya.r750.arc 

#transform the science spectrum
gstransform ptf10lya.1r750 wavtran=ptf10lya.r750.arc 

#transform the science spectrum
gstransform ptf10lya.2r750 wavtran=ptf10lya.r750.arc 

display tptf10lya.1r750.fits[sci,1] 1
display tptf10lya.2r750.fits[sci,1] 2

#sky subtraction: the sample is manually selected after inspection of the image
# output has an s prefixed to the front
gsskysub tptf10lya.1r750 long_sample="150:250,280:380" fl_inter-

gsskysub tptf10lya.2r750 long_sample="150:250,280:380" fl_inter-

#examine the sky subtracted spectrum
# ds9 stptf10lya.1r750.fits & # do this in a terminal
display stptf10lya.1r750[sci,1]  # or do this in iraf

#Extract the specctrum
gsextract stptf10lya.1r750 fl_inter+

gsextract stptf10lya.2r750 fl_inter+

gscalibrate estptf10lya.1r750  sfunc=sensr

gscalibrate estptf10lya.2r750  sfunc=sensr

oned # to load oned spectral reduction iraf package

splot cestptf10lya.1r750.fits  # just to check

wspectext cestptf10lya.1r750.fits[sci,1] ptf10lya.1red.dat header-

wspectext cestptf10lya.2r750.fits[sci,1] ptf10lya.2red.dat header-


######## the object spectrum B455 ###############

#normalize the flat field
gsflat @ptf10lya.1b455.flat.txt ptf10lya.1b455.flat order=23 rawpath="raw$" bias="bias" fl_inter-

#reduce CuAr.  Not flat fielded, gaps are not fixpixed
gsreduce @ptf10lya.1b455.arc.txt outimages="ptf10lya.1b455.arc" rawpath="raw$" fl_flat- bias="bias" fl_fixpix-

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.1b455.txt outimages="ptf10lya.1b455" rawpath="raw$" bias="bias" flat="ptf10lya.1b455.flat"

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.2b455.txt outimages="ptf10lya.2b455" rawpath="raw$" bias="bias" flat="ptf10lya.1b455.flat"

#determine wavelength calibration -- 1d and 2d
gswavelength ptf10lya.1b455.arc fl_inter+

#transform the CuAr spectrum, for checking that the transformation is OK
#output spectrum has prefix t
#gstransform ptf10lya.1b455.arc wavtran=ptf10lya.1b455.arc 

#transform the science spectrum
gstransform ptf10lya.1b455 wavtran=ptf10lya.1b455.arc 

#transform the science spectrum
gstransform ptf10lya.2b455 wavtran=ptf10lya.1b455.arc 

display tptf10lya.1b455.fits[sci,1] 1
display tptf10lya.2b455.fits[sci,1] 1

#sky subtraction: the sample is manually selected after inspection of the image
# output has an s prefixed to the front
gsskysub tptf10lya.1b455 long_sample="50:150,205:305" fl_inter-

gsskysub tptf10lya.2b455 long_sample="50:150,205:305" fl_inter-

#examine the sky subtracted spectrum
# ds9 stptf10lya.1b455.fits & # do this in a terminal
display stptf10lya.1b455[sci,1]  # or do this in iraf
display stptf10lya.2b455[sci,1]  # or do this in iraf

#Extract the specctrum
gsextract stptf10lya.1b455 fl_inter+

gsextract stptf10lya.2b455 fl_inter+

gscalibrate estptf10lya.1b455  sfunc=sensb

gscalibrate estptf10lya.2b455  sfunc=sensb

oned # to load oned spectral reduction iraf package

splot cestptf10lya.1b455.fits  # just to check

wspectext cestptf10lya.1b455.fits[sci,1] ptf10lya.1blue455.dat header-

wspectext cestptf10lya.2b455.fits[sci,1] ptf10lya.2blue455.dat header-

######## the object spectrum R755 ###############

#normalize the flat field
gsflat @ptf10lya.r755.flat.txt ptf10lya.r755.flat order=23 rawpath="raw$" bias="bias" fl_inter-

#reduce CuAr.  Not flat fielded, gaps are not fixpixed
gsreduce @ptf10lya.r755.arc.txt outimages="ptf10lya.r755.arc" rawpath="raw$" fl_flat- bias="bias" fl_fixpix-

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.1r755.txt outimages="ptf10lya.1r755" rawpath="raw$" bias="bias" flat="ptf10lya.r755.flat"

#gsreduce subtracts bias, mosaics detectors, flat fields, gets approx wav. cal.
gsreduce @ptf10lya.2r755.txt outimages="ptf10lya.2r755" rawpath="raw$" bias="bias" flat="ptf10lya.r755.flat"

#determine wavelength calibration -- 1d and 2d
gswavelength ptf10lya.r755.arc fl_inter+

#transform the CuAr spectrum, for checking that the transformation is OK
#output spectrum has prefix t
#gstransform ptf10lya.r755.arc wavtran=ptf10lya.r755.arc 

#transform the science spectrum
gstransform ptf10lya.1r755 wavtran=ptf10lya.r755.arc 

#transform the science spectrum
gstransform ptf10lya.2r755 wavtran=ptf10lya.r755.arc 

display tptf10lya.1r755.fits[sci,1] 1
display tptf10lya.2r755.fits[sci,1] 2

#sky subtraction: the sample is manually selected after inspection of the image
# output has an s prefixed to the front
gsskysub tptf10lya.1r755 long_sample="80:180,220:320" fl_inter-

gsskysub tptf10lya.2r755 long_sample="80:180,220:320" fl_inter-

#examine the sky subtracted spectrum
# ds9 stptf10lya.1r755.fits & # do this in a terminal
display stptf10lya.1r755[sci,1]  # or do this in iraf

#Extract the specctrum
gsextract stptf10lya.1r755 fl_inter+

gsextract stptf10lya.2r755 fl_inter+

gscalibrate estptf10lya.1r755  sfunc=sensr

gscalibrate estptf10lya.2r755  sfunc=sensr

oned # to load oned spectral reduction iraf package

splot cestptf10lya.1r755.fits  # just to check

wspectext cestptf10lya.1r755.fits[sci,1] ptf10lya.1red755.dat header-

wspectext cestptf10lya.2r755.fits[sci,1] ptf10lya.2red755.dat header-



########### Make text files #############3
onedspec

wspectext cestptf10lya.1r750.fits[sci,1] ptf10lya.red.dat header-
wspectext cestptf10lya.1b450.fits[sci,1] ptf10lya.blue.dat header-

fixspec,'ptf10lya.blue.dat',redspec='ptf10lya.red.dat',bluerange=[3800,5500],redrange=[5501,9900],xra=[3700,10000],outfile='ptf10lya.dat'
