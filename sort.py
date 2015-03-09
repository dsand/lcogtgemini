'''
Created on Jan 14, 2015

@author: cmccully
'''
import os
from glob import glob
import shutil
import pyfits
import numpy as np
# gunzip all of the files
os.system('gunzip *.gz')
# Sort all of the gemini files into their respective directories so I can 
# reduce the data.
# make a bias directory:
if not os.path.exists('bias'):
    os.mkdir('bias')
# move all of the bias frames into the bias directory
biasims = glob('*bias.fits')
biasmjds = []
biasinstruments = []
for im in biasims:
    biasmjds.append(float(pyfits.getval(im, 'MJD-OBS', ext=1)))
    biasinstruments.append(pyfits.getval(im, 'INSTRUME', ext=0).strip())
    shutil.move(im, 'bias/')
    
biasmjds = np.array(biasmjds)
biasinstruments = np.array(biasinstruments)
biasims = np.array(biasims)

ims = glob('*.fits')
scimjds = []
scinames = []
scidates = []
sciinstruments = []
for im in ims:
# Get the science target names
    if pyfits.getval(im, 'OBSTYPE', ext=0) == 'OBJECT' and (pyfits.getval(im, 'OBSCLASS') == 'science' or pyfits.getval(im, 'OBSCLASS') == 'partnerCal') :
        obj = pyfits.getval(im, 'OBJECT')
        # remove any white space
        obj = obj.replace(' ', '_')
        obj = obj.lower()
        # for each science target make a new directory
        if not os.path.exists(obj):
            os.mkdir(obj)
        mjd = float(pyfits.getval(im, 'MJD-OBS', ext=1))
        scimjds.append(mjd)
        scinames.append(obj)
        scidate = pyfits.getval(im, 'DATE-OBS') 
        scidates.append(scidate)
        sciinstrument  = pyfits.getval(im, 'INSTRUME', ext = 0).strip()
        
        sciinstruments.append(sciinstrument)
        if not os.path.exists(obj + '/' + scidate):
            os.mkdir(obj + '/' + scidate)
        # Move the science data into the new directory
        shutil.move(im, obj + '/' + scidate + '/')

        # Find the closest bias frame and copy it into the science directory
        # get only the bias frames from the right site
        w = biasinstruments == sciinstrument
        ind = np.argmin(np.abs(mjd - biasmjds[w]))

        shutil.copy('bias/'+biasims[w][ind],obj + '/' + scidate + '/bias.fits')    
        
sciinstruments = np.array(sciinstruments)
scimjds = np.array(scimjds)
scinames = np.array(scinames)
scidates = np.array(scidates)

ims = glob('*.fits')
for im in ims:
    # Find the calibration files taken for these observations
    mjd = float(pyfits.getval(im, 'MJD-OBS', ext=1))
    instrument  = pyfits.getval(im, 'INSTRUME', ext = 0).strip()
    w = instrument == sciinstruments
    ind = np.argmin(np.abs(mjd - scimjds[w]))
    
    # move the calibration frames into their respective object folders
    shutil.move(im, scinames[w][ind]+'/'+scidates[w][ind]+'/')