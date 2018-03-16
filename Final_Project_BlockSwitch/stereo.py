

from window import*
from mdct import*
from psychoac import*

import numpy as np 

from numpy import*


# BMLD in Bark Scale 

def mldBark(z):
    
    a = 1.25
    offset = 2.5

    mldBark = np.zeros_like(z)

    mldBark = np.power(10.0, a * (1 - np.cos(np.pi * (np.minimum(z, 15.5)/15.5)) - offset))

    # normalize
    mldBark = mldBark

    return mldBark


# BMLD in Hz scale 
def mldBarkFrequency(f):
    
    a = 1.25
    offset = 2.5
    cutoff = 3000.

    mldBark = np.zeros_like(f)

    mldBark = np.power(10.0, a * (1 - np.cos(np.pi * (np.minimum(f, cutoff)/cutoff) ) - offset))

    # normalize
    mldBark = mldBark/np.amax(mldBark)

    return mldBark

# Peak calculation as used in homework 
def findpeaks(Xwdb, sampleRate, N):

    peaks = []
    freqs = []

    length = np.size(Xwdb)

    # find peaks and order from max amplitude to min
    for sample in range(1,length-1):

        if (abs(Xwdb[sample]) > abs(Xwdb[sample-1]) and abs(Xwdb[sample]) > abs(Xwdb[sample+1]) and 10.0* log10(abs(Xwdb[sample]))>-30.0):

            peaks = np.append(peaks,Xwdb[sample])
            freqs = np.append(freqs,sample)

    # peaks = peaks.astype(int)
    peaks = absolute(peaks).astype(int)
    freqsIndex = asarray(freqs).astype(int)

    # parabolic interpolation
    estimateFreqs = []
    estimateAmp = []

    for idx in range(0,len(freqs)):

        a = abs(Xwdb[int(freqs[idx])-1])
        b = abs(Xwdb[int(freqs[idx])])
        r = abs(Xwdb[int(freqs[idx])+1])
        p = (1/2)*(a-r)/(a+r-2*b)
        A = b-(a-r)*(p/4)
        estimateFreqs = append(estimateFreqs,(freqs[idx]+p)*(sampleRate/N))
        estimateAmp = append(estimateAmp,A)

    return estimateAmp, estimateFreqs, freqsIndex


# Function to estimate the base thresholds 
def calcBTHR(data, MDCTdata, MDCTscale, sampleRate, sfBands, noDrop):
   
    N = len(data)

    # noDrop = True 
    nMDCTLines = len(MDCTdata)
    X_fft = fft.fft(HanningWindow(data))[0:N/2]

    masked_intensity = zeros_like(MDCTdata)

    # shift by 0.5 samples
    MDCTFreqs = sampleRate / 2.0 / nMDCTLines * (arange(0, nMDCTLines) + 0.5)

    # threshold in quiet
    threshold_in_quiet = Intensity(Thresh(MDCTFreqs))

    # Using JOS code written for 320A
    estimated_peak_amplitudes, estimated_peak_frequencies, index = findpeaks(X_fft,sampleRate,N)

    BW = 3
    masker_spl = zeros(len(index), dtype=float64)
    num_peaks = len(index)

    # aggregate intensity across the peaks, create maskers, sum maskers
    for i in range(num_peaks):
        masker_spl[i] = SPL((8.0 / 3.0 * 4.0 / (N ** 2.0)) * sum(abs(X_fft[index[i] - BW:index[i] + BW])**2.0))
        maskers = Masker(estimated_peak_frequencies[i],masker_spl[i],True)
        if noDrop:
            maskers.drop=0
        masked_intensity += (maskers).vIntensityAtFreq(MDCTFreqs)

    masked_intensity = (masked_intensity + threshold_in_quiet)

    return  SPL(masked_intensity)


# Estimate the SMR's for the two channels 
def stereoSMR(stereoThreshold, mdctSPL, sfBands):
    

    # two channels
    numChannels = 2
    SMRs = []

    for channel in range(numChannels):

        SMRs.append([])
        # for each band calculate max SMR
        for band in range(sfBands.nBands):

            lower = sfBands.lowerLine[band]
            upper = sfBands.upperLine[band]+1

            # get mask for this band
            mask = stereoThreshold[channel][lower:upper]
            x = mdctSPL[channel][lower:upper]
            # SMR for whole band
            bandSMR = x - mask

            # max
            if len(bandSMR) == 0:

                SMRs[channel].append(-96.0)

            else:

                SMRs[channel].append(np.amax(bandSMR))

    return SMRs

# Correlation based on MPEG-1
def corellationLR( dataL , dataR , lowLine , upperLine ):

	corr = abs(sum(dataL[lowLine:upperLine]**2 - dataR[lowLine:upperLine]**2 )) < -0.8*abs(sum(dataL[lowLine:upperLine]**2 + dataR[lowLine:upperLine]**2))

	return corr

# Correllation based on coherence 
def coherenceLR( dataL , dataR , lowLine, upperLine , threshold):

    coherence = np.mean(( abs(dataL[lowLine:upperLine]*dataR[lowLine:upperLine])**2 ) / ( ( abs(dataL[lowLine:upperLine])**2 ) * abs( dataR[lowLine:upperLine] )**2 )) >= threshold

    return coherence


# Calculate masked thresholds based on Johnston and Ferrira's Sum and Difference Stereo Coding Paper
def stereoMaskThresholds(data, MDCTdata, MDCTscale, sampleRate, sfBands, LRMS, codingParams): # sendMS):
    

    ################ L/R SMR calculation ################

    # calculate MDCT SPL for L/R
    MDCT_Spl_L = SPL(4.*MDCTdata[0]**2)-(6.02*MDCTscale[0]) 
    MDCT_Spl_R = SPL(4.*MDCTdata[1]**2)-(6.02*MDCTscale[1]) 
    MDCT_Spl_LR = [MDCT_Spl_L, MDCT_Spl_R]

    

    dropLR = True
    dropMS = False
    BTHR_L = calcBTHR(data[0], MDCTdata[0], MDCTscale[0], sampleRate, sfBands, False) # Has drop
    BTHR_R = calcBTHR(data[1], MDCTdata[1], MDCTscale[1], sampleRate, sfBands, False)
    
    THR_LR = [BTHR_L, BTHR_R]

    ################ M/S SMR calculation ################

    # transform time domain L/R data into M/S
    data_MS = [(data[0] + data[1]) / 2.0, (data[0] - data[1]) / 2.0]
    # transform MDCT L/R data into M/S
    MDCT_data_MS = [(MDCTdata[0] + MDCTdata[1]) / 2.0, (MDCTdata[0] - MDCTdata[1]) / 2.0]

    # calculate MDCT SPL for M/S
    MDCT_Spl_M = SPL(4.*MDCT_data_MS[0]**2)-(6.02*MDCTscale[0]) 
    MDCT_Spl_S = SPL(4.*MDCT_data_MS[1]**2)-(6.02*MDCTscale[1]) 
    MDCT_Spl_MS = [MDCT_Spl_M, MDCT_Spl_S]

    # calculate basic thresholds for MS
    BTHR_M = calcBTHR(data_MS[0], MDCT_data_MS[0], MDCTscale[0], sampleRate, sfBands, False) # Has drop
    BTHR_S = calcBTHR(data_MS[1], MDCT_data_MS[0], MDCTscale[1], sampleRate, sfBands, False)
    BTHR_M_mldBark = calcBTHR(data_MS[0], MDCT_data_MS[0], MDCTscale[0], sampleRate, sfBands, True) # No drop
    BTHR_S_mldBark = calcBTHR(data_MS[1], MDCT_data_MS[0], MDCTscale[1], sampleRate, sfBands, True)
    
    BTHR_MS = [BTHR_M, BTHR_S]


    # MDCT freqs [Hz]
    MDCT_freqs = (((np.arange(len(MDCTdata[0])) + 0.5) / len(MDCTdata[0])) * (sampleRate / 2.0))


    # get mldBarks
    bark = Bark(MDCT_freqs)

    mldBark = mldBarkFrequency(MDCT_freqs)

    
    mldBark_M = BTHR_M_mldBark * mldBark
    mldBark_S = BTHR_S_mldBark * mldBark

    THR_MS = [np.maximum(BTHR_M, np.minimum(BTHR_S, mldBark_S)), np.maximum(BTHR_S, np.minimum(BTHR_M, mldBark_M))]

    # get max SMRs for L/R
    SMR_LR = stereoSMR(THR_LR, MDCT_Spl_LR, sfBands)

    # get max SMRs for M/S
    SMR_MS = stereoSMR(THR_MS, MDCT_Spl_MS, sfBands)


    SMR = np.zeros_like(SMR_MS)
    LRMSmdctLines = np.zeros_like(MDCT_Spl_MS)
    # band by band take M/S or L/R
    for channel in range(2):

        for line in range(sfBands.nBands):
            lowLine = sfBands.lowerLine[line]
            highLine = sfBands.upperLine[line] + 1

            if LRMS[line]:
                # take M/S SMR
                SMR[channel][line] = SMR_MS[channel][line]
                # take M/S lines
                LRMSmdctLines[channel][lowLine:highLine] = MDCT_data_MS[channel][lowLine:highLine]
            else:
                # take L/R SMR
                SMR[channel][line] = SMR_LR[channel][line]
                # take L/R lines
                LRMSmdctLines[channel][lowLine:highLine] = MDCTdata[channel][lowLine:highLine]

    return SMR,LRMSmdctLines
