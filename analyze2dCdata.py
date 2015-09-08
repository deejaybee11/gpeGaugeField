# -*- coding: utf-8 -*-
"""
Created on Thu Apr 09 13:01:42 2015

@author: Dylan Brown

The MIT License (MIT)

Copyright (c) 2015

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from tkinter import filedialog as tk
import os
import struct
from astropy.io import fits
from matplotlib import cm
import glob

plt.close("all")

root = tk.Tk()
root.withdraw()

folder = tk.askdirectory()

if not os.path.exists(folder+"/pngs"):
    os.makedirs(folder+"/pngs")
    

numfiles=0
for i in glob.glob(folder+"/fits/psi" + "*.fits"):
    numfiles +=1
    
for i in glob.glob(folder+"/pngs/" + "*.png"):
    os.remove(i)

filename = folder + "/groundState.fits"
fitsImage = fits.open(filename, mode='readonly')
image = fitsImage[0].data
fitsImage.close()

plt.ioff()
fig, ax = plt.subplots()
ax.imshow(image, cmap = cm.afmhot)
plt.savefig(folder+'/groundState.png',dpi = 250)
plt.close('all')
#
filename = folder + "/initPsi.fits"
fitsImage = fits.open(filename, mode='readonly')
image = fitsImage[0].data
fitsImage.close()

plt.ioff()
fig, ax = plt.subplots()
ax.imshow(image, cmap = cm.afmhot)
plt.savefig(folder+'/initPsi.png',dpi = 250)
plt.close('all')

filename = folder + "/energyX.fits"
fitsImage = fits.open(filename, mode='readonly')
image = fitsImage[0].data
fitsImage.close()

plt.ioff()
fig, ax = plt.subplots()
ax.imshow(image, cmap = cm.afmhot)
plt.savefig(folder+'/energyX.png',dpi = 250)
plt.close('all')

filename = folder + "/energyY.fits"
fitsImage = fits.open(filename, mode='readonly')
image = fitsImage[0].data
fitsImage.close()

plt.ioff()
fig, ax = plt.subplots()
ax.imshow(image, cmap = cm.afmhot)
plt.savefig(folder+'/energyY.png',dpi = 250)
plt.close('all')
#

##
for i in range(0,numfiles):

    print("Importing psi file number " + str(i))    
    
    filename = folder + "/fits/psi" + str(i) + ".fits"
    fitsImage = fits.open(filename,mode='readonly')
    image = (fitsImage[0].data)
    fitsImage.close()
    
    plt.ioff()
#for i in range(0,numfiles):
    
    fig, ax = plt.subplots()
    ax.imshow(image, cmap = cm.afmhot)
    
    if i<10:
        plt.savefig(folder+"/pngs"+'/Psi00'+str(i)+'.png',dpi = 250)
    elif i>=10:
        plt.savefig(folder+"/pngs"+'/Psi0'+str(i)+'.png',dpi = 250) 

    plt.close('all')
#for i in range(0,numfiles):
#    print("Importing phi file number " + str(i))    
#    
#    filename = folder + "/fits/phi" + str(i) + ".fits"
#    fitsImage = fits.open(filename,mode='readonly')
#    image2 = (fitsImage[0].data)
#    fitsImage.close()
#    
#    plt.ioff()
#
#    
#    fig, ax = plt.subplots()
#    ax.imshow(np.fft.fftshift(image2), cmap = cm.afmhot)
#    
#    if i<10:
#        plt.savefig(folder+"/pngs"+'/Phi00'+str(i)+'.png',dpi = 250)
#    elif i>=10:
#        plt.savefig(folder+"/pngs"+'/Phi0'+str(i)+'.png',dpi = 250) 
#
#    plt.close('all')
    
