# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:31:37 2018

@author: orchisama
"""
from __future__ import division
from pcmfile import *;
import quantize as qt
import window as win
import mdct as m
import psychoac as psy
import blockswitch as bs
import matplotlib.pyplot as plt


# create the audio file objects of the appropriate audioFile type
inFile= PCMFile("Castanets.wav")

# open input file and get its coding parameters
codingParams= inFile.OpenForReading()

# set additional coding parameters that are needed for encoding/decoding
codingParams.nSamplesPerBlock = 512
codingParams.nScaleBits = 4


#Before data-reading loop, create two arrays
priorBlock = np.zeros((codingParams.nChannels, codingParams.nSamplesPerBlock))
overlapAndAdd = np.zeros((codingParams.nChannels, codingParams.nSamplesPerBlock))
flag = False
sfBands = psy.ScaleFactorBands(psy.AssignMDCTLinesFromFreqLimits(codingParams.nSamplesPerBlock,
                                                             codingParams.sampleRate))
                                                                 
#store signal here                                                                 
signal = np.array([], dtype = float)
transients = []
PE = np.array([], dtype = float)

while True:
    data = inFile.ReadDataBlock(codingParams)
    nBits = 8
    if not data : break
    for iCh in range(codingParams.nChannels):
       
        #concatenate current block of new data samples with priorBlock
        newBlock = np.concatenate([priorBlock[iCh], data[iCh]])
        #save current block to pass on to detectTransient function
        dataBlock = np.copy(newBlock)
        #also keep appending it to longer signal chain
        signal = np.concatenate([signal, dataBlock])
        #save current block of new data as priorBlock for next pass
        priorBlock[iCh] = data[iCh]
        
        #window 2*nSamplesPerBlock set of samples
        newBlock = win.KBDWindow(newBlock)
        #MDCT
        MDCTdata = m.MDCT(newBlock,codingParams.nSamplesPerBlock,codingParams.nSamplesPerBlock)
        #find MDCT scale        
        maxLine = np.max(np.abs(newBlock))
        MDCTscale = qt.ScaleFactor(maxLine,codingParams.nScaleBits)
        
        #make some changes here to detect transients
        [tr, pe] = bs.detectTransient(sfBands, dataBlock, MDCTdata, MDCTscale, codingParams.sampleRate, 1)
        transients.append(tr)
        PE = np.append(PE, pe)
               
      

#find indices of all blocks that have transients
transient_blocks = np.where(transients)[0]
#transient position in samples
transient_pos = transient_blocks * codingParams.nSamplesPerBlock * 2.0
Fs = codingParams.sampleRate
#time vector
t = np.arange(0, np.size(signal)/Fs, 1.0/Fs)

plt.figure()
plt.plot(t,signal)
#plot transients as vertical lines
xcoords = transient_pos/Fs
for xc in xcoords:
    plt.axvline(x=xc, color = 'r', linestyle = '--')
plt.title('Castanets')
plt.xlabel('Time in seconds')
plt.xlim(xmin = 0, xmax = 26)
plt.ylabel('Amplitude')

# close the files
inFile.Close(codingParams)

#print(PE[transient_blocks])





