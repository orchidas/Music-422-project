from __future__ import division
import numpy as np
from window import *
import matplotlib.pyplot as plt
#import psychoac_ as psy
from mdct import *
#import mdct_ as m

def SPL(intensity):
    """
    Returns the SPL corresponding to intensity (in units where 1 implies 96dB)
    
    """
    spl = 10*np.log10(intensity) + 120 - 24
    spl[spl < -30] = -30
    return spl

def Intensity(spl):
    """
    Returns the intensity (in units of the reference intensity level) for SPL spl
    
    """
    return np.power(10,(spl-96)/10.0)

def Thresh(f):
    """Returns the threshold in quiet measured in SPL at frequency f (in Hz)"""
    #set all frequencies below 10Hz to 10Hz
    f[f < 10] = 10
    #convert to kHz
    f = f/1000.0
    SPL_thresh = 3.64*np.power(f,-0.8) - 6.5*np.exp(-0.6*np.power((f-3.3),2)) + (10**-3)*np.power(f,4)

    return SPL_thresh 

def Bark(f):
    """Returns the bark-scale frequency for input frequency f (in Hz) """
    f = f/1000.0
    bark = 13*np.arctan(0.76*f) + 3.5*np.arctan(np.power(f/7.5,2))
    return bark # TO REPLACE WITH YOUR CODE

class Masker:
    """
    a Masker whose masking curve drops linearly in Bark beyond 0.5 Bark from the
    masker frequency
    """

    def __init__(self,f,SPL,isTonal=True):
        """
        initialized with the frequency and SPL of a masker and whether or not
        it is Tonal
        """
        self.f = f
        self.SPL = SPL
        self.isTonal = isTonal
        

    def IntensityAtFreq(self,freq):
        """The intensity of this masker at frequency freq"""
        
        return 0 # TO REPLACE WITH YOUR CODE

    def IntensityAtBark(self,z):
        """The intensity of this masker at Bark location z"""
        masker_z = Bark(self.f)
        Lm = self.SPL
        dz = z - masker_z
        
        #find spreading function
        if(dz < -0.5):
            spread_function = -27*(np.abs(dz)-0.5)
        elif(dz > 0.5):
            spread_function = (-27 + 0.367*np.max([Lm - 40,0]))*(np.abs(dz)- 0.5)
        else:
            spread_function = 0
        
        if(self.isTonal):
            delta = 16
        else:
            delta = 6
            
        SPL_at_bark = spread_function - delta
        
        return (Intensity(SPL_at_bark))

    def vIntensityAtBark(self,zVec):
        """The intensity of this masker at vector of Bark locations zVec"""
        
        masker_z = Bark(self.f)
        #downshift by delta
        if(self.isTonal):
            delta = 24
        else:
            delta = 6
            
        Lm = self.SPL  
        dz = zVec - masker_z
        spread_function = np.zeros(np.size(zVec))
        
        #find spreading 
        inds1 = np.where(dz < -0.5) 
        inds2 = np.where(dz > 0.5)
        spread_function[inds1] = -27*(np.abs(dz[inds1])-0.5)
        spread_function[inds2] = (-27 + 0.367*np.max([Lm - 40,0]))*(np.abs(dz[inds2])- 0.5)
   
       
        spread_function = spread_function + Lm - delta
        
        #convert SPL to intensity        
        return Intensity(spread_function)
                


# Default data for 25 scale factor bands based on the traditional 25 critical bands
cbFreqLimits = [100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,
                2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500]  

def AssignMDCTLinesFromFreqLimits(nMDCTLines, sampleRate, flimit = cbFreqLimits):
    """
    Assigns MDCT lines to scale factor bands for given sample rate and number
    of MDCT lines using predefined frequency band cutoffs passed as an array
    in flimit (in units of Hz). If flimit isn't passed it uses the traditional
    25 Zwicker & Fastl critical bands as scale factor bands.
    """
    
    MDCTLineLoc = (sampleRate/2.0)/nMDCTLines
    #MDCT locations are offset by 0.5 times a line spacing from DFT locations
    MDCTLineFreqs = MDCTLineLoc * np.arange(nMDCTLines) + 0.5*MDCTLineLoc
    flimit = np.concatenate([[0], flimit, [sampleRate/2.0]])
    nLinesPerBand = np.zeros(np.size(flimit)-1, dtype = int)
    
    for k in range(1,np.size(flimit)):
        nLinesPerBand[k-1] = np.size(MDCTLineFreqs[(MDCTLineFreqs>flimit[k-1]) & 
        (MDCTLineFreqs<flimit[k])])
            
    
    return nLinesPerBand # TO REPLACE WITH YOUR CODE

class ScaleFactorBands:
    """
    A set of scale factor bands (each of which will share a scale factor and a
    mantissa bit allocation) and associated MDCT line mappings.

    Instances know the number of bands nBands; the upper and lower limits for
    each band lowerLimit[i in range(nBands)], upperLimit[i in range(nBands)];
    and the number of lines in each band nLines[i in range(nBands)]
    """

    def __init__(self,nLines):
        """
        Assigns MDCT lines to scale factor bands based on a vector of the number
        of lines in each band
        """
        self.nLines = nLines
        self.nBands = 25
        self.upperLine = np.cumsum(nLines)-1
        self.lowerLine = self.upperLine - nLines + 1
        
        


def getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Return Masked Threshold evaluated at MDCT lines.

    Used by CalcSMR, but can also be called from outside this module, which may
    be helpful when debugging the bit allocation code.
    """
    #do fft
    N = np.size(data)
    X_fft = np.fft.fft(data)
    X_fft = X_fft[:N/2]
    
    #find tonal maskers        
    (peak_SPL, peak_freq) = findPeakSPLandFreq(X_fft, N, sampleRate, 'rect')
    nMDCTLines = N/2.0
    MDCTLineLoc = (sampleRate/2.0)/nMDCTLines
    MDCTLineFreqs = MDCTLineLoc * np.arange(nMDCTLines) + 0.5*MDCTLineLoc
    zVec = Bark(MDCTLineFreqs)
    
    
    npeaks = peak_SPL.size
    mask = list()
    intensity = list(np.array([]))
    
    for i in range(npeaks):
        mask.append(Masker(peak_freq[i], peak_SPL[i], True))
        intensity.append(mask[i].vIntensityAtBark(zVec))
        
     
    mask_thresh_SPL =  SPL(sum(intensity) + Intensity(Thresh(MDCTLineFreqs)))
        
    return mask_thresh_SPL # TO REPLACE WITH YOUR CODE


def CalcSMRs(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Set SMR for each critical band in sfBands.

    Arguments:
                data:       is an array of N time domain samples
                MDCTdata:   is an array of N/2 MDCT frequency lines for the data
                            in data which have been scaled up by a factor
                            of 2^MDCTscale
                MDCTscale:  is an overall scale factor for the set of MDCT
                            frequency lines
                sampleRate: is the sampling rate of the time domain samples
                sfBands:    points to information about which MDCT frequency lines
                            are in which scale factor band

    Returns:
                SMR[sfBands.nBands] is the maximum signal-to-mask ratio in each
                                    scale factor band

    Logic:
                Performs an FFT of data[N] and identifies tonal and noise maskers.
                Sums their masking curves with the hearing threshold at each MDCT
                frequency location to the calculate absolute threshold at those
                points. Then determines the maximum signal-to-mask ratio within
                each critical band and returns that result in the SMR[] array.
    """
    
    """The difference in level between a signal component and the masking
    threshold at a certain frequency is sometimes referred to as the signal to
    mask ratio, SMR"""
    
    #get combined masking curve at each MDCT frequency location
    masked_thres = getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands)
    #MDCT data may have been scaled up, we need to scale it back down
    MDCTdata = MDCTdata/(2**MDCTscale)
    SMR = np.zeros(sfBands.nBands)
    #MDCT with KBD window
    signal_SPL = SPL((2.0/findWindowPower('kbd', np.size(data)))*(np.abs(MDCTdata)**2))
    
    for k in range(sfBands.nBands):
        SMR[k] = np.max(signal_SPL[sfBands.lowerLine[k]:sfBands.upperLine[k]+1] -
        masked_thres[sfBands.lowerLine[k]:sfBands.upperLine[k]+1])
        
        
    return SMR        
        
    

def findPeaks(fft_data):
    
    """
    Finds peaks in data
    
    """
    N = np.size(fft_data)
    localPeakPos = np.array([], dtype = int)
    
    #find local minima
    for k in range(1,N-1):
        if(fft_data[k] > fft_data[k-1] and fft_data[k] > fft_data[k+1]):
            localPeakPos = np.append(localPeakPos, k)
    
    return localPeakPos
    
    
def findPeakSPLandFreq(X_fft, N, Fs, win):
    
    """Finds peak SPL and frequencies"""
    
    #find peaks
    
    peak_inten = 4/((N**2)*findWindowPower(win,N))*(np.abs(X_fft)**2)
    peakPos = findPeaks(peak_inten)
    npeaks = peakPos.size
    peak_freq = np.zeros(npeaks)
    peak_SPL = SPL(peak_inten[peakPos] + peak_inten[peakPos-1] + peak_inten[peakPos+1])
    freqs_fft = np.linspace(0,Fs/2,N/2+1)
    #omit last half samples
    freqs_fft = freqs_fft[:N/2]
     
    #intensity weighted average frequency for peak bin
    for l in range(npeaks):
        curPeak = peakPos[l]
        f = [freqs_fft[curPeak-1], freqs_fft[curPeak], freqs_fft[curPeak+1]]
        w = [peak_inten[curPeak-1], peak_inten[curPeak], peak_inten[curPeak+1]]
        peak_freq[l] = np.average(f, weights = w)

     
    return (peak_SPL, peak_freq)


def findWindowPower(win, N):
    """return window power based on type of window"""
    
    if(win == 'hann'):
        return 0.375
    elif(win == 'rect'):
        return 1
    elif(win == 'sine'):
        return 0.5
    elif(win == 'kbd'):
        return (1.0 / N) * np.sum(KBDWindow(np.ones(N))**2)
    
    
    
            

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":
    
    #Q1b    
    Fs = 48000
    amp = np.array([0.54,0.21,0.10,0.06,0.03,0.01])
    freqs = np.array([440,550,660,880,4400,8800], dtype = float)
    L = 2048
    x = np.zeros(L, dtype = float)
    n = np.arange(L)    
    
    for k in range(6):
        x = x + amp[k]*np.cos(2*np.pi*freqs[k]*n/Fs)
        
#    block = np.array([512,1024,2048], dtype = int)
#    plt.figure()
#    
#    for k in range(3):
#        
#         N = block[k]
#         x_win = HanningWindow(x[:N])
#         X_fft = np.fft.fft(x_win)
#         X_fft = X_fft[:N/2]
#         FFT_spl = SPL(4/((N**2)*findWindowPower('hann',N))*(np.abs(X_fft)**2))
#         freqs_fft = np.linspace(0,Fs/2,N/2+1)
#         #omit last half samples
#         freqs_fft = freqs_fft[:N/2]
#         
#         #find peaks
#         npeaks = 6             
#         (peak_SPL, peak_freq) = findPeakSPLandFreq(X_fft, N, Fs, 'hann')
#         print('Peak frequencies : ', peak_freq)
#         print('Peak SPL : ', peak_SPL)
#         
#         plt.subplot(3,1,(k+1))
#         plt.plot(np.log10(freqs_fft), FFT_spl)
#         plt.plot(np.log10(peak_freq), peak_SPL,'ro')
#         plt.title('N = ' + str(N))
#   
#    plt.xlabel('Log frequency in Hz')
#    plt.ylabel('Normalized SPL in dB')
    
    #------------------------------------------------------------------------
    
    #Q1c
#    N = 1024
#    x_win = HanningWindow(x[:N])
#    X_fft = np.fft.fft(x_win)
#    X_fft = X_fft[:N/2]
#    FFT_spl = SPL(4/((N**2)*findWindowPower('hann',N))*(np.abs(X_fft)**2))
#    freqs_fft = np.linspace(0,Fs/2,N/2+1)
#    #omit last half samples
#    freqs_fft = freqs_fft[:N/2]
#    
#    plt.figure()
#    plt.plot(np.log10(freqs_fft/1000.0), FFT_spl)
#    plt.plot(np.log10(freqs_fft/1000.0), Thresh(freqs_fft),'r')
#    plt.title('N = ' + str(N))
#    plt.xlabel('Log frequency in kHz')
#    plt.ylabel('Normalized SPL in dB')
    
    #------------------------------------------------------------------------
    
    #Q1d
#    freqs = np.array([0,100,200,300,400,510,630,770,920,1080,1270,1400,1720,2000,2320,2700,3150,
#                      3700,4400,5300,6400,7700,9500,12000,15500])
#    bark = Bark(freqs)
#    print(bark)
    
    #-----------------------------------------------------------------------
    
    #Q1e
    
#    N = 1024
#    x_win = HanningWindow(x[:N])
#    X_fft = np.fft.fft(x_win)
#    X_fft = X_fft[:N/2]
#    FFT_spl = SPL(4/((N**2)*findWindowPower('hann',N))*(np.abs(X_fft)**2))
#    freqs_fft = np.linspace(0,Fs/2,N/2+1)
#    #omit last half samples
#    freqs_fft = freqs_fft[:N/2]
#     
#    #find peaks
#    npeaks = 6         
#    (peak_SPL, peak_freq) = findPeakSPLandFreq(X_fft, N, Fs, 'hann')
#    print('Peak frequencies : ', peak_freq)
#    print('Peak SPL : ', peak_SPL)
#    
#    zVec = Bark(freqs_fft)
#    
#    plt.figure()
#    plt.plot(zVec, FFT_spl)
#    plt.plot(zVec, Thresh(freqs_fft))
#    
#    mask = list()
#    intensity = list(np.array([]))
#    
#    for i in range(npeaks):
#        mask.append(Masker(peak_freq[i], peak_SPL[i], True))
#        intensity.append(mask[i].vIntensityAtBark(zVec))
#        plt.plot(zVec, SPL(intensity[i]))
#        
#     
#    plt.plot(zVec, SPL(sum(intensity) + Intensity(Thresh(freqs_fft))), 'k.') 
#    plt.title('N = ' + str(N))
#    plt.ylim(ymin = -40, ymax = 200)
#    plt.xlabel('Barks')
#    plt.ylabel('Normalized SPL in dB')
##    
#    #------------------------------------------------------------------------
#    
#    #Q1f
#    
#    nMDCTLines = 512
#    sc = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(nMDCTLines,Fs))
#    MDCTLineLoc = (Fs/2.0)/nMDCTLines
#    MDCTLineFreqs = MDCTLineLoc * sc.lowerLine + 0.5*MDCTLineLoc
#    xposition = Bark(MDCTLineFreqs)
#    
#    for xc in xposition:
#        plt.axvline(x=xc, color='k', linestyle='--')
    
    #-----------------------------------------------------------------------
    
    #Q1g
    N = 1024
    data = x[:N]
    MDCTdata = MDCT(data, N/2.0, N/2.0)
    nMDCTLines = N/2.0
    
    sfBands_true = psy.ScaleFactorBands(psy.AssignMDCTLinesFromFreqLimits(nMDCTLines,Fs))
    sfBands = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(nMDCTLines,Fs))
    
    mask_thres_true = psy.getMaskedThreshold(data, MDCTdata, 0, Fs, sfBands_true)
    mask_thres = getMaskedThreshold(data, MDCTdata, 0, Fs, sfBands )
    
    nMDCTLines = N/2.0
    MDCTLineLoc = (Fs/2.0)/nMDCTLines
    MDCTLineFreqs = MDCTLineLoc * np.arange(nMDCTLines) + 0.5*MDCTLineLoc
    zVec = Bark(MDCTLineFreqs)
    
    SMR_true = psy.CalcSMRs(data, MDCTdata, 0, Fs, sfBands_true)
    SMR = CalcSMRs(data, MDCTdata, 0, Fs, sfBands)    
    
    plt.figure()
    plt.plot(zVec, mask_thres_true, 'b')
    plt.plot(zVec, mask_thres, 'g')
    plt.plot(zVec, SPL(4*(np.abs(MDCTdata)**2)), 'r')
    
    plt.ylim(ymin = -40, ymax = 100)
    plt.xlabel('Barks')
    plt.ylabel('Normalized SPL in dB')
    
    print(SMR_true)
    print(SMR)
#                

    
    

    pass # TO REPLACE WITH YOUR CODE
