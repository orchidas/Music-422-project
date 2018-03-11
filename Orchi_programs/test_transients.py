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
inFile= PCMFile("audio/Castanets.wav")

# open input file and get its coding parameters
codingParams= inFile.OpenForReading()

# set additional coding parameters that are needed for encoding/decoding
codingParams.nSamplesPerBlock = 512
codingParams.nScaleBits = 4


#Before data-reading loop, create two arrays
priorBlock = np.zeros((codingParams.nChannels, codingParams.nSamplesPerBlock))
overlapAndAdd = np.zeros((codingParams.nChannels, codingParams.nSamplesPerBlock))
flag = False
sfBands = psy.ScaleFactorBands(psy.AssignMDCTLinesFromFreqLimits(codingParams.nSamplesPerBlock/2,
                                                             codingParams.sampleRate))
                                                                 
#store signal here                                                                 
signal = np.array([], dtype = float)
transients = [[] for i in range(codingParams.nChannels)]
percEntropy = [[] for i in range(codingParams.nChannels)]
E = [[] for i in range(codingParams.nChannels)]
thres_change = 100
thres_mag = 1000
count = 0

while True:
    data = inFile.ReadDataBlock(codingParams)
    nBits = 8
    if not data : break
    #also keep appending it to longer signal chain
    signal = np.concatenate([signal, data[iCh]])
    
    for iCh in range(codingParams.nChannels):
       
        newBlock = data[iCh]
#        #save current block to pass on to detectTransient function
#        dataBlock = newBlock       
#        #window 2*nSamplesPerBlock set of samples
#        newBlock = win.KBDWindow(newBlock)
#        #MDCT
#        MDCTdata = m.MDCT(newBlock,codingParams.nSamplesPerBlock/2,codingParams.nSamplesPerBlock/2)
#        #find MDCT scale        
#        maxLine = np.max(np.abs(newBlock))
#        MDCTscale = qt.ScaleFactor(maxLine,codingParams.nScaleBits)
#        
#        #find perceptual entropy for each channel
#        percEntropy[iCh].append(bs.detectTransient(sfBands, dataBlock, MDCTdata, MDCTscale, codingParams.sampleRate))
#    
#        #detect transients by detecting change in PE from block to block    
#        if((count > 0 and percEntropy[iCh][count] - percEntropy[iCh][count-1] >= thres_change)
#            or percEntropy[iCh][count] > thres_mag):
#            transients[iCh].append(True)
#        else:
#            transients[iCh].append(False)
            
        #using Prateek's FFT transient detector
        transients[iCh].append(bs.transient_detection(newBlock))
    
#    count += 1
                   

#find indices of all blocks that have transients for the first channel only
transient_blocks = np.where(transients[0])[0]
np.savetxt('transients.out', transient_blocks, delimiter=',') 
#transient position in samples
transient_pos = transient_blocks * codingParams.nSamplesPerBlock
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
plt.xlim(xmin = 0, xmax = 8)
plt.ylabel('Amplitude')

# close the files
inFile.Close(codingParams)






