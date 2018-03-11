from __future__ import division
import numpy as np
from window import *
import matplotlib.pyplot as plt
from psychoac import *
from mdct import *
import quantize as Q
import testbitalloc as tb


# Question 1.b)
def BitAllocUniform(bitBudget, maxMantBits, nBands, nLines, SMR=None):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are uniformely distributed for the mantissas.
    """
    
    #Bits is the number of bits per line per band
    N = np.sum(nLines)
    bits = np.zeros(nBands)
    bits = np.ones(nBands)*np.floor(bitBudget/N)
    
    #uniformly distribute remaining bits
    bitBudget = bitBudget - N*np.floor(bitBudget/N)
    
    while(bitBudget >= np.min(nLines)):
        for k in range(nBands):
            if(nLines[k] < bitBudget):
                bits[k] += 1
                bitBudget -= nLines[k]
        
    return bits.astype(np.int64) # TO REPLACE WITH YOUR VECTOR

def BitAllocConstSNR(bitBudget, maxMantBits, nBands, nLines, peakSPL):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep a constant
    quantization noise floor (assuming a noise floor 6 dB per bit below
    the peak SPL line in the scale factor band).
    
    """
    
    bits = np.zeros(nBands)
    SNR = peakSPL
    ind = 0
    
    #while(bitBudget >= np.min(nLines)):
    while(bitBudget >= np.min(nLines)):
        ind = np.argmax(SNR)
        if(bitBudget >= nLines[ind]):
            bits[ind] += 1
            bitBudget -= nLines[ind]
        SNR[ind] -= 6.02
        
    #get rid of all ones in bits array and reallocate them elsewhere
    found_ones = np.where(bits == 1)
    bits[found_ones] -= 1
    bitBudget += np.sum(nLines[found_ones])
    
    #make sure number of bits per band does not exceed maxMantBits
    found_max = np.where(bits > maxMantBits)
    extra_bits = bits[found_max] - maxMantBits
    bits[found_max] -= extra_bits
    bitBudget += np.sum(nLines[found_max]*extra_bits)
        
    while(bitBudget >= np.min(nLines)):
        ind = np.argmax(SNR)
        if(bitBudget >= nLines[ind] and bits[ind] > 0 and bits[ind] < maxMantBits):
            bits[ind] += 1
            bitBudget -= nLines[ind]
        SNR[ind] -= 6.02
            
    return bits.astype(np.int64) # TO REPLACE WITH YOUR VECTOR

def BitAllocConstMNR(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep the quantization
    noise floor a constant distance below (or above, if bit starved) the
    masked threshold curve (assuming a quantization noise floor 6 dB per
    bit below the peak SPL line in the scale factor band).
    """
    
    bits = np.zeros(nBands)
    
    while(bitBudget >= np.min(nLines)):    
        ind = np.argmax(SMR)
        if(bitBudget >= nLines[ind]):
            bits[ind] += 1
            bitBudget -= nLines[ind]
        SMR[ind] -= 6.02
        
    #get rid of all ones in bits array and reallocate them elsewhere
    found_ones = np.where(bits == 1)
    bits[found_ones] -= 1
    bitBudget += np.sum(nLines[found_ones])
    
    #make sure number of bits per band does not exceed maxMantBits
    found_max = np.where(bits > maxMantBits)
    extra_bits = bits[found_max] - maxMantBits
    bits[found_max] -= extra_bits
    bitBudget += np.sum(nLines[found_max]*extra_bits)
        
    while(bitBudget >= np.min(nLines)):
        ind = np.argmax(SMR)
        if(bitBudget >= nLines[ind] and bits[ind] > 0 and bits[ind] < maxMantBits):
            bits[ind] += 1
            bitBudget -= nLines[ind]
        SMR[ind] -= 6.02
        
    return bits.astype(np.int64) # TO REPLACE WITH YOUR VECTOR

# Question 1.c)
def BitAlloc(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Allocates bits to scale factor bands so as to flatten the NMR across the spectrum

       Arguments:
           bitBudget is total number of mantissa bits to allocate
           maxMantBits is max mantissa bits that can be allocated per line
           nBands is total number of scale factor bands
           nLines[nBands] is number of lines in each scale factor band
           SMR[nBands] is signal-to-mask ratio in each scale factor band

        Return:
            bits[nBands] is number of bits allocated to each scale factor band

        Logic:
           Maximizing SMR over blook gives optimization result that:
               R(i) = P/N + (1 bit/ 6 dB) * (SMR[i] - avgSMR)
           where P is the pool of bits for mantissas and N is number of bands
           This result needs to be adjusted if any R(i) goes below 2 (in which
           case we set R(i)=0) or if any R(i) goes above maxMantBits (in
           which case we set R(i)=maxMantBits).  (Note: 1 Mantissa bit is
           equivalent to 0 mantissa bits when you are using a midtread quantizer.)
           We will not bother to worry about slight variations in bit budget due
           rounding of the above equation to integer values of R(i).
    """
    #Find optimum bit budget allocation for each band
#    avgSMR = np.mean(SMR) 
#    R = np.zeros(nBands)
#    for k in range(nBands):
#        R[k] = np.round((bitBudget/nBands)/nLines[k] + (SMR[k] - avgSMR)/6.02)
#        if(R[k] < 2):
#            R[k] = 0
#        elif(R[k] > maxMantBits):
#            R[k] = maxMantBits
#        bitBudget -= nLines[k]*R[k]
#        
    #waterfilling works better
    R = BitAllocConstMNR(bitBudget, maxMantBits, nBands, nLines, SMR.copy())
        
    return R.astype(np.int64) # TO REPLACE WITH YOUR CODE

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
        
    plt.figure()
         
    N = 2048
    x_win = HanningWindow(x)
    X_fft = np.fft.fft(x_win)
    X_fft = X_fft[:N/2]
    data = SineWindow(x)
    MDCTdata = MDCT(data, N/2.0, N/2.0)
    MDCT_spl = np.maximum(96. + 10*np.log10(4.*MDCTdata**2),-30)
    freqs_mdct = np.linspace(0,Fs/2,N/2+1) + 0.5*(Fs/N)
    #omit last half samples
    freqs_mdct = freqs_mdct[:N/2]
     
    #find peaks
    npeaks = 6             
    (peak_SPL, peak_freq) = findPeakSPLandFreq(MDCTdata, N, Fs, 'hann')
    
    
    nBands = 25
    dataRate = 192000
    nMDCTLines = N/2.0
    bitBudget = (dataRate * nMDCTLines)/Fs - 2*(nBands*4) - 4
    sfBands = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(nMDCTLines,Fs))
    
    #uniform bit allocation
    bits_uniform = BitAllocUniform(bitBudget, 16, nBands, sfBands.nLines)
    
    #get highest SPL in each band
    peakSPL = np.zeros(nBands)
    for k in range(nBands):
        peakSPL[k] = np.max(MDCT_spl[sfBands.lowerLine[k]:sfBands.upperLine[k]+1])
    #bit allocation with constant quantization noise    
    bits_constSNR = BitAllocConstSNR(bitBudget, 16, nBands, sfBands.nLines, peakSPL.copy())
    
    #get SMR
    SMR = CalcSMRs(data, MDCTdata, 0, Fs, sfBands)
    #bit allocation with quantization noise floor constant noise floor below
    bits_constMNR = BitAllocConstMNR(bitBudget, 16, nBands, sfBands.nLines, SMR.copy())
    
    #get optimum bits in each band according to formula
    bits_opt = BitAlloc(bitBudget, 16, nBands, sfBands.nLines, SMR.copy())
    
    #plotting stuff
    zVec = Bark(freqs_mdct)    
    mask = list()
    intensity = list(np.array([]))
    
    for i in range(npeaks):
        mask.append(Masker(peak_freq[i], peak_SPL[i], True))
        intensity.append(mask[i].vIntensityAtBark(zVec))
        
        
    bits = [bits_uniform, bits_constSNR, bits_constMNR, bits_opt]
    freqs_bark = np.array([10**-9,100,200,300,400,510,630,770,920,1080,1270,1400,1720,2000,2320,2700,3150,
                      3700,4400,5300,6400,7700,9500,12000,15500])
    
    for i in range(4):
        noise_floor = peakSPL - 6.02 * bits[i]
        plt.figure(i+1)
        plt.semilogx(freqs_mdct, MDCT_spl)
        plt.semilogx(freqs_mdct, SPL(sum(intensity) + Intensity(Thresh(freqs_mdct))), 'k')
        plt.semilogx(freqs_bark, noise_floor)
        plt.ylim(ymin = np.min(noise_floor), ymax = 150)
        plt.xlim(xmin = 10, xmax = 2*10**4)
        plt.xlabel('Log frequency in Hz')
        plt.ylabel('Normalized SPL in dB')
    
    #test stuff
    print(tb.TestBitAlloc(bits_uniform ,bitBudget, 16, nBands, sfBands.nLines, SMR.copy(), peakSPL.copy()))    
    print(tb.TestBitAlloc(bits_constSNR ,bitBudget, 16, nBands, sfBands.nLines, SMR.copy(), peakSPL.copy()))
    print(tb.TestBitAlloc(bits_constMNR, bitBudget, 16, nBands, sfBands.nLines, SMR.copy(), peakSPL.copy()))

    pass # TO REPLACE WITH YOUR CODE
