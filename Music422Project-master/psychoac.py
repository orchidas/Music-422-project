import numpy as np
import matplotlib as plt
from window import *
from scipy import signal
import mdct as mdct


import psychoac_ as ps_true

def SPL(intensity):
    """
    Returns the SPL corresponding to intensity (in units where 1 implies 96dB)
    """
    #return np.zeros_like(intensity) # TO REPLACE WITH YOUR CODE


    I_0 = 10**(-12)
    spl = 10*np.log10(np.abs(intensity)/I_0) -24.
    spl[spl < -30.] = -30.
    #spl = np.max(-30.*np.ones((len(spl))),spl)
    return spl



def Intensity(spl):
    """
    Returns the intensity (in units of the reference intensity level) for SPL spl
    """
    #return np.zeros_like(spl) # TO REPLACE WITH YOUR CODE

    I_0 = 10**(-12)
    return I_0 * 10**((24 + spl)/10)





def Thresh(f):
    """Returns the threshold in quiet measured in SPL at frequency f (in Hz)"""

    #return np.zeros_like(f) # TO REPLACE WITH YOUR CODE

    f[f < 10] = 10.

    return 3.64*(f/1000.)**(-0.8) - 6.5*np.exp(-0.6*(f/1000. - 3.3)**2) + 10**(-3.) * (f/1000.)**4.






def Bark(f):
    """Returns the bark-scale frequency for input frequency f (in Hz) """
    #return np.zeros_like(f) # TO REPLACE WITH YOUR CODE

    return 13*np.arctan(0.76*f/1000.) + 3.5*np.arctan((f/7500.)**2)




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
        #pass # TO REPLACE WITH YOUR CODE
        
        if isTonal:
            self.delta = 16
        else:
            self.delta = 6
            
        self.z_masker = Bark(f)
        self.L_M = SPL 
        

    def IntensityAtFreq(self,freq):
        """The intensity of this masker at frequency freq"""
        return 0 # TO REPLACE WITH YOUR CODE

    def IntensityAtBark(self,z):
        """The intensity of this masker at Bark location z"""
        #return 0 # TO REPLACE WITH YOUR CODE
        
        dz = z - self.z_masker
        
        if dz >= -0.5 and dz <= 0.5:
            return 0.  + self.L_M -self.delta
        elif(dz <=0.5):
            return -27.*(np.abs(dz) - 0.5)  + self.L_M -self.delta
        else:
            return (-27. + 0.367*np.max([self.L_M - 40. , 0])) * (np.abs(dz) - 0.5) + self.L_M -self.delta
        
        

    def vIntensityAtBark(self,zVec):
        """The intensity of this masker at vector of Bark locations zVec"""
        #return np.zeros_like(zVec) # TO REPLACE WITH YOUR CODE
        
        
        
        dz = zVec - self.z_masker
        
        #a = np.where(dz>=-0.5 and dz<=0.5 ,  - self.delta + self.L_M ,0)
        #b = np.where(dz<-0.5, -27.*(np.abs(dz) - 0.5) - self.delta + self.L_M ,0)
        #c = np.where(dz>-0.5, (-27. + 0.367*np.max([self.L_M - 40. , 0])) * (np.abs(dz) - 0.5) - self.delta + self.L_M ,0)
        
        
        a = np.zeros((len(zVec),)) + self.L_M -self.delta
        b = np.where(dz<-0.5, -27.*(np.abs(dz) - 0.5)  ,0)
        c = np.where(dz>0.5, (-27. + 0.367*np.max([self.L_M - 40. , 0])) * (np.abs(dz) - 0.5)  ,0)
        
        return a + b + c
        
        
        


# Default data for 25 scale factor bands based on the traditional 25 critical bands
cbFreqLimits = [100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500, 24000]  # TO REPLACE WITH THE APPROPRIATE VALUES

def AssignMDCTLinesFromFreqLimits(nMDCTLines, sampleRate, flimit = cbFreqLimits):
    """
    Assigns MDCT lines to scale factor bands for given sample rate and number
    of MDCT lines using predefined frequency band cutoffs passed as an array
    in flimit (in units of Hz). If flimit isn't passed it uses the traditional
    25 Zwicker & Fastl critical bands as scale factor bands.
    """
    #return np.zeros(len(flimit)) # TO REPLACE WITH YOUR CODE
    
    spacing = sampleRate/(nMDCTLines*2.)
    freqOfBin = np.arange(0,sampleRate/2., spacing) + spacing/2.
    
    N = len(flimit)
    numBins = np.zeros((N,))
    
    for i in range(0, N):
        if i ==0:
            numBins[i] = ( (freqOfBin <= flimit[i])).sum()
        elif i<N-1:
            numBins[i] = ((flimit[i-1] < freqOfBin) & (freqOfBin <= flimit[i])).sum()
        elif i==N-1:
            numBins[i] = (freqOfBin > flimit[i-1]).sum()
    
    return numBins.astype(int)








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
        #pass # TO REPLACE WITH YOUR CODE
        self.nBands = len(nLines)
        
        self.lowerLine = np.zeros((self.nBands,))
        for i in range(1, self.nBands):
            self.lowerLine[i] = np.sum(nLines[0:i])
            
        self.upperLine = np.zeros((self.nBands,))
        for i in range(0, self.nBands):
            self.upperLine[i] = np.sum(nLines[0:i+1]) -1
            
        #self.nLines = nLines.astype(int)   
        self.nLines = nLines
        self.lowerLine = self.lowerLine.astype(int)
        self.upperLine = self.upperLine.astype(int)
            
            
            


def getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Return Masked Threshold evaluated at MDCT lines.

    Used by CalcSMR, but can also be called from outside this module, which may
    be helpful when debugging the bit allocation code.
    """
    #return np.zeros_like(0) # TO REPLACE WITH YOUR CODE
    
    N = len(data)

    # Fine peaks
    intensity = (4*np.abs(np.fft.fft(HanningWindow(data)))**2) / (np.mean(HanningWindow(np.ones((N,)))) * N**2)
    peaks = findPeaks(intensity)
    peakSPL_data = peakSPL(intensity,peaks)
    peakFreq_data = peakFreq(intensity,peaks, sampleRate)
    
    # Threshold in quiet
    
    #f = np.arange(0,sampleRate/2,1)
    f_data = np.linspace(0,sampleRate/2 + 1/sampleRate,N/2)
    thresh_data = Thresh(f_data)
    
    
    
    # Maskers
    
    numPeaks = len(peaks)
    
    maskers = []
    maskerIntensities = np.zeros((numPeaks,N/2))
    barkVecFreq = Bark(f_data)
    
    for i in range(0,numPeaks):
        
        maskers.append(Masker(peakFreq_data[i],peakSPL_data[i],isTonal=True))
        maskerIntensities[i,:] = maskers[i].vIntensityAtBark(barkVecFreq)
    
    
    allMaskingCurves = np.vstack([maskerIntensities,thresh_data])
    #sumMaskers = SPL(np.sum(Intensity(maskerIntensities),0))
    #allMaskingCurves = np.vstack([sumMaskers,thresh_data])
    
    combinedMasking = np.max(allMaskingCurves,0)
    
    return combinedMasking
    





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
    #return np.zeros_like(0) # TO REPLACE WITH YOUR CODE
    
    
    maskingCurve = getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands)
    
    
    MDCTscale = (2**(-MDCTscale)*MDCTdata)
    MDCTdata_SPL =  SPL(4*MDCTscale**2)
    #MDCTdata_SPL = 96 + 10*np.log10(2/(   np.mean(SineWindow(np.ones((N,)))**2)) * np.abs(2**(-MDCTscale)*MDCTdata)**2) 
    
    
    SMR = np.zeros((sfBands.nBands),)
    
    for i in range(0, sfBands.nBands):
        SMR[i] = np.max( MDCTdata_SPL[sfBands.lowerLine[i]:sfBands.upperLine[i]+1] - maskingCurve[sfBands.lowerLine[i]:sfBands.upperLine[i]+1])
        
    
    
    
    
    return SMR
    
    
    
    
    
    
    
    

#-----------------------------------------------------------------------------
    

# My helper functions
    
def findPeaks(data):
    N = len(data)
    peaks = []
    for u in range(1,N/2-1):
        if ((data[u]>data[u-1])&(data[u]>data[u+1])):
            peaks = np.append(peaks,u)
                     
    return peaks

def peakSPL(intensity, peaks):
    N = len(peaks)
    peakSPL = np.zeros(N,)
    for i in range(0,N):
        peakSPL[i] = 96 + 10*np.log10(intensity[int(peaks[i])-1] + intensity[int(peaks[i])] + intensity[int(peaks[i])+1]  )

    return peakSPL

def peakFreq(intensity, peaks,Fs):
    N = len(peaks)
    M = len(intensity)
    peakFreq = np.zeros(N,)
    for i in range(0,N):
        weight =  intensity[int(peaks[i])-1] +  intensity[int(peaks[i])]+ intensity[int(peaks[i])+1]
        a = intensity[int(peaks[i])-1]*(peaks[i]-1)
        b = intensity[int(peaks[i])]*peaks[i]
        c = intensity[int(peaks[i])+1]*(peaks[i]+1)
        peakFreq[i] = (( a + b + c  )/weight)*Fs/M

    return peakFreq





#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    pass # TO REPLACE WITH YOUR CODE
    
    
    testData = np.arange(0.001,1,0.01)*10**(-1)
    
    SPL_test = SPL(testData)
    SPL_true = ps_true.SPL(testData)
    SPL_err = np.sum(SPL_test-SPL_true)
    
    Intensity_test = Intensity(SPL_true)
    Intensity_true = ps_true.Intensity(SPL_true)
    Intensity_err = np.sum(Intensity_test-Intensity_true)
    
        
    
    
    
    
    
    
    
    # 1.b)
    Fs = 48000
    n = np.arange(0,Fs)

    x = 0.54 * np.cos(2*np.pi*440*n/Fs) + 0.21 * np.cos(2*np.pi*550*n/Fs) + 0.10 * np.cos(2*np.pi*660*n/Fs) + 0.06 * np.cos(2*np.pi*880*n/Fs) + 0.03 * np.cos(2*np.pi*4400*n/Fs) + 0.01 * np.cos(2*np.pi*8800*n/Fs)
    
    x_512 = x[0:512]
    x_1024 = x[0:1024]
    x_2048 = x[0:2048]
    
    f_512 = np.linspace(0,Fs,512)
    f_1024 = np.linspace(0,Fs,1024)
    f_2048 = np.linspace(0,Fs,2048)
    
    
    I_512 = (4*np.abs(np.fft.fft(HanningWindow(x_512)))**2) / (np.mean(HanningWindow(np.ones((512,))))* 512**2)
    I_1024 = (4*np.abs(np.fft.fft(HanningWindow(x_1024)))**2) / (np.mean(HanningWindow(np.ones((1024,)))) * 1024**2)
    I_2048 = (4*np.abs(np.fft.fft(HanningWindow(x_2048)))**2) / (np.mean(HanningWindow(np.ones((2048,)))) * 2048**2)
    
    
    
    # Find peak locations
    peaks_512 = findPeaks(I_512)
    peaks_1024 = findPeaks(I_1024)
    peaks_2048 = findPeaks(I_2048)
    
    peakSPL_512 = peakSPL(I_512,peaks_512)
    peakSPL_1024 = peakSPL(I_1024,peaks_1024)
    peakSPL_2048 = peakSPL(I_2048,peaks_2048)
    
    peakFreq_512 = peakFreq(I_512,peaks_512, Fs)
    peakFreq_1024 = peakFreq(I_1024,peaks_1024, Fs)
    peakFreq_2048 = peakFreq(I_2048,peaks_2048, Fs)
    
    
    
#    fig1 = plt.figure()
#    
#    plt.subplot(3,1,1)
#    plt.semilogx(f_512,SPL(I_512),'.-', label = 'N = 512')
#    axes = plt.gca()
#    axes.set_ylim([0,100])
#    axes.set_xlim([100,Fs/2])
#         
#    #plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1b) N = 512')
#    plt.grid(True)
#    plt.legend()
#    
#    plt.subplot(3,1,2)
#    plt.semilogx(f_1024,SPL(I_1024),'.-', label = 'N = 1024')
#    axes = plt.gca()
#    axes.set_ylim([0,100])
#    axes.set_xlim([100,Fs/2])
#         
#    #plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1b) N = 1024')
#    plt.grid(True)
#    plt.legend()
#    
#    plt.subplot(3,1,3)
#    plt.semilogx(f_2048,SPL(I_2048),'.-', label = 'N = 2048')
#    axes = plt.gca()
#    axes.set_ylim([0,100])
#    axes.set_xlim([100,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1b N = 2048)')
#    plt.grid(True)
#    plt.legend()
#         
#         
#    plt.savefig("Q1b.png")
    
    
    
    
    
    
    # 1.c)
    
    
    test_f = np.arange(0,Fs/2,1)
    
    f_thresh = Thresh(test_f)
    f_thresh_true = ps_true.Thresh(test_f)
    f_thresh_err = np.sum(f_thresh-f_thresh_true)
    
    
#    fig2 = plt.figure()
#
#    plt.semilogx(f_1024,SPL(I_1024),'.-', label = 'N = 1024')
#    plt.semilogx(test_f,f_thresh, '-', label = 'Thresh in Quiet')
#    axes = plt.gca()
#    axes.set_ylim([-10,100])
#    axes.set_xlim([10,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1c)')
#    plt.grid(True)
#    plt.legend()
#         
#         
#    plt.savefig("Q1c.png")




    # 1.d)
    
    f_l = np.array([0, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500])
    f_l = f_l.astype(float)

    bark = Bark(f_l)
    bark_true = ps_true.Bark(f_l)
    bark_err = np.sum(bark -bark_true)




    # 1.e)
    testMasker = Masker(2000.,100,isTonal=True)
    testMaskerIntensity = testMasker.IntensityAtBark(10)
    
    # Test vectorized
    barkVec = np.arange(0,25)
    barkVecFreq = Bark(f_1024)
    
    testMaskerIntensityVec = testMasker.vIntensityAtBark(barkVec)
    
    
    # Make masker objects
    numPeaks = 6
    
    maskers = []
    maskerIntensities = np.zeros((6,1024))
    
    for i in range(0,6):
        
        maskers.append(Masker(peakFreq_1024[i],peakSPL_1024[i],isTonal=True))
        maskerIntensities[i,:] = maskers[i].vIntensityAtBark(barkVecFreq)
    
    
    
    #f_center = [50, 150, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500, 10500, 13500, 18000]
    
#    fig3 = plt.figure()
#
#    plt.semilogx(f_1024,SPL(I_1024),'.-', label = 'N = 1024')
#    plt.semilogx(test_f,f_thresh, '-', label = 'Thresh in Quiet')
#    
#    plt.semilogx(f_1024 ,maskerIntensities[0,:],'-', label = 'peak 1')
#    plt.semilogx(f_1024 ,maskerIntensities[1,:],'-', label = 'peak 2')
#    plt.semilogx(f_1024 ,maskerIntensities[2,:],'-', label = 'peak 3')
#    plt.semilogx(f_1024 ,maskerIntensities[3,:],'-', label = 'peak 4')
#    plt.semilogx(f_1024 ,maskerIntensities[4,:],'-', label = 'peak 5')
#    plt.semilogx(f_1024 ,maskerIntensities[5,:],'-', label = 'peak 6')
#    
#    axes = plt.gca()
#    axes.set_ylim([-10,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.e) Maskers')
#    plt.grid(True)
#    plt.legend()
#         
#         
#    plt.savefig("Q1e_all.png")
    
    
    # Make plot for all thresholds combined
    thresh_1024 = Thresh(f_1024)
    allMaskingCurves = np.vstack([maskerIntensities,thresh_1024])
    combinedMasking = np.max(allMaskingCurves,0)
    
    
#    fig4 = plt.figure()
#    
#
#    plt.semilogx(f_1024,SPL(I_1024),'.-', label = 'N = 1024')
#    plt.semilogx(f_1024,combinedMasking, '-', label = 'Combined Threshold')
#
#    axes = plt.gca()
#    axes.set_ylim([-10,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.e) Combined')
#    plt.grid(True)
#    plt.legend()
#         
#         
#    plt.savefig("Q1e_combined.png")

    



# 1.f)
    nMDCTLines = 512.
    sampleRate = 48000.
    testAssignMDCTLines = AssignMDCTLinesFromFreqLimits(nMDCTLines, sampleRate)
    testAssignMDCTLines_true = ps_true.AssignMDCTLinesFromFreqLimits(nMDCTLines, sampleRate)
    testAssignMDCTLines_err = np.sum(testAssignMDCTLines_true - testAssignMDCTLines)
    
    
    # scale factor bands testing
    testScaleFactorBands = ScaleFactorBands(testAssignMDCTLines)
    testScaleFactorBands_true = ps_true.ScaleFactorBands(testAssignMDCTLines)
    
    testLowerLine = testScaleFactorBands.lowerLine 
    testLowerLine_true = testScaleFactorBands_true.lowerLine 
    testLowerLine_err = np.sum(testLowerLine_true - testLowerLine  )
    
    testUpperLine = testScaleFactorBands.upperLine 
    testUpperLine_true = testScaleFactorBands_true.upperLine 
    testUpperLine_err = np.sum(testUpperLine_true - testUpperLine  )
    
    
    xcoords = (testScaleFactorBands.upperLine)*sampleRate/nMDCTLines
    
#    fig4 = plt.figure()
#    
#    
#
#
#    plt.semilogx(f_1024,SPL(I_1024),'.-', label = 'N = 1024')
#    plt.semilogx(f_1024,combinedMasking, '-', label = 'Combined Threshold')
#    
#    for xc in xcoords:
#        plt.axvline(x=xc,color='k')
#
#    axes = plt.gca()
#    axes.set_ylim([-10,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.f) ')
#    plt.grid(True)
#    plt.legend()
#         
#         
#    plt.savefig("Q1f.png")
    
    
    
    
# 1.g)    

    # Testing getMaskedThreshold
    
    
    
#    combinedMasking = getMaskedThreshold(x_1024, 0, 0, Fs, 0)
#    combinedMasking_true = ps_true.getMaskedThreshold(x_1024, 0, 0, 48000, sfBands)
#    combinedMasking_err = np.sum(combinedMasking_true - combinedMasking)
    
#    fig5 = plt.figure()
#    
#    
#    plt.semilogx(f_1024,SPL(I_1024),'.-', label = 'N = 1024')
#    
#    plt.semilogx(f_1024[0:512],combinedMasking, '-', label = 'Combined Threshold')
#    
#
#    axes = plt.gca()
#    axes.set_ylim([-10,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.g getMaskedThreshold) ')
#    plt.grid(True)
#    plt.legend()
#         
#         
#    plt.savefig("Q1g.png")
    
    
    # Testing CalcSMRs
    
    data = x_1024
    N = 1024
    sampleRate = 48000
    MDCTscale = 0
    MDCTdata = (2**MDCTscale) * mdct.MDCT(SineWindow(data), N/2, N/2, isInverse=False)
    MDCTLines = AssignMDCTLinesFromFreqLimits(N/2, sampleRate)
    sfBands = ScaleFactorBands(MDCTLines)
    
    SMRs = CalcSMRs(data, MDCTdata, MDCTscale, sampleRate, sfBands)
    SMRs_true = ps_true.CalcSMRs(data, MDCTdata, MDCTscale, sampleRate, sfBands)
    SMRs_err = np.sum(SMRs-SMRs_true)
    
    
    combinedMasking = getMaskedThreshold(x_1024, 0, 0, Fs, sfBands)
    combinedMasking_true = ps_true.getMaskedThreshold(x_1024, 0, 0, 48000, sfBands)
    combinedMasking_err = np.sum(combinedMasking_true - combinedMasking)
    
    
    
    
    
    fig6 = plt.figure()
    
    
    plt.semilogx(f_1024,SPL(I_1024),'.-', label = 'N = 1024')
    
    plt.semilogx(f_1024[0:512],combinedMasking, '-', label = 'Combined Threshold')
    plt.semilogx(f_1024[0:512],combinedMasking_true, '-', label = 'Combined Threshold')
    

    axes = plt.gca()
    axes.set_ylim([-10,100])
    axes.set_xlim([50,Fs/2])
         
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('SPL (dB)')
    plt.title('1.g getMaskedThreshold) ')
    plt.grid(True)
    plt.legend() 
    
    
    
    
    
    
    
    
    
    
