import numpy as np
from audiofile import * # base class
import matplotlib as plt
from window import *
from scipy import signal
import mdct as mdct
import mdct_ as mdct_true
import psychoac as ps
import psychoac_ as ps_true

import testbitalloc

#import ScaleFactorBands, AssignMDCTLinesFromFreqLimits  # defines the grouping of MDCT lines into scale factor bands

import bitalloc_ as bitalloc_true

# Question 1.b)
def BitAllocUniform(bitBudget, maxMantBits, nBands, nLines, SMR=None):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are uniformely distributed for the mantissas.
    """
    #return zeros(nBands) # TO REPLACE WITH YOUR VECTOR
    
    
    numLines = np.sum(nLines)
    numBitsPerLine = int(np.floor(bitBudget/numLines))
    
    bitAlloc = (numBitsPerLine * np.ones((nBands))).astype(int)
    sampsLeft = bitBudget - numLines*numBitsPerLine
    for i in range(0,nBands):
        if (sampsLeft >nLines[i]):
            bitAlloc[i] = bitAlloc[i] +1
            sampsLeft = sampsLeft - nLines[i]
    
    
    return bitAlloc
    
    

def BitAllocConstSNR(bitBudget, maxMantBits, nBands, nLines, peakSPL):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep a constant
    quantization noise floor (assuming a noise floor 6 dB per bit below
    the peak SPL line in the scale factor band).
    """
    #return zeros(nBands) # TO REPLACE WITH YOUR VECTOR
    
    numLines = np.sum(nLines)
    
#    bitsNeededPerBand = np.ceil(peakSPL/6.02)
#    
#    #bitAlloc = np.zeros((numLines,))
#    
#    #for i in range(0, numLines):
#    #    noiseSPL = peakSPL[i] - 6
#    #    numBits = 
#    
#    bitsNeededPerBand = np.where(bitsNeededPerBand > 1., bitsNeededPerBand , 0)
#    bitsNeededPerBand = np.where(bitsNeededPerBand <= maxMantBits, bitsNeededPerBand , maxMantBits)
#
#    while( np.sum(bitsNeededPerBand*nLines) > bitBudget):
#        bitsNeededPerBand = bitsNeededPerBand -1.
#        #np.where(bitsNeededPerBand <1, 0)
#        bitsNeededPerBand = np.where(bitsNeededPerBand > 1, bitsNeededPerBand , 0)
#        
#    
#    numBitsLeft = bitBudget - np.sum(bitsNeededPerBand*nLines)
#    
#    while (numBitsLeft > np.min(nLines )):    
#        for i in range(0,nBands):
#            if(numBitsLeft >= nLines[i] and   bitsNeededPerBand[i] >1  and   bitsNeededPerBand[i] < maxMantBits ):
#                bitsNeededPerBand[i] = bitsNeededPerBand[i] + 1
#                numBitsLeft = numBitsLeft - nLines[i]
#    
#    return bitsNeededPerBand    


    SPL = peakSPL.copy()
    bitsNeededPerBand = np.zeros((nBands,))
    bitsUsed = 0
    
    while( bitBudget > bitsUsed):
        index = np.argmax(SPL)
        bitsNeededPerBand[index] = bitsNeededPerBand[index] + 1 
        bitsUsed = bitsUsed + nLines[index]
        SPL[index] = SPL[index] - 6.02

    bitsNeededPerBand = np.where(bitsNeededPerBand > 1, bitsNeededPerBand , 0)
    
    bitsNeededPerBand = np.where(bitsNeededPerBand < maxMantBits, bitsNeededPerBand , maxMantBits)
    return bitsNeededPerBand.astype(int)









def BitAllocConstMNR(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep the quantization
    noise floor a constant distance below (or above, if bit starved) the
    masked threshold curve (assuming a quantization noise floor 6 dB per
    bit below the peak SPL line in the scale factor band).
    """
    #return zeros(nBands) # TO REPLACE WITH YOUR VECTOR


    numLines = np.sum(nLines)
    
    #bitsNeededPerBand = np.ceil(SMR/6.02)
    
    #bitAlloc = np.zeros((numLines,))
    
    #for i in range(0, numLines):
    #    noiseSPL = peakSPL[i] - 6
    #    numBits = 
    
#    bitsNeededPerBand = np.where(bitsNeededPerBand > 1., bitsNeededPerBand , 0)
#    bitsNeededPerBand = np.where(bitsNeededPerBand <= maxMantBits, bitsNeededPerBand , maxMantBits)
#
#    while( np.sum(bitsNeededPerBand*nLines) > bitBudget):
#        bitsNeededPerBand = bitsNeededPerBand -1.
#        #np.where(bitsNeededPerBand <1, 0)
#        bitsNeededPerBand = np.where(bitsNeededPerBand > 1, bitsNeededPerBand , 0)
#        
#    
#    numBitsLeft = bitBudget - np.sum(bitsNeededPerBand*nLines)
#    
#    while (numBitsLeft > np.min(nLines )):    
#        for i in range(0,nBands):
#            if(numBitsLeft >= nLines[i] and   bitsNeededPerBand[i] >1  and   bitsNeededPerBand[i] < maxMantBits ):
#                bitsNeededPerBand[i] = bitsNeededPerBand[i] + 1
#                numBitsLeft = numBitsLeft - nLines[i]
#    
#    return bitsNeededPerBand 
    
    SMRwork = SMR.copy()
    bitsNeededPerBand = np.zeros((nBands,))
    bitsUsed = 0
    
    while( bitBudget > bitsUsed):
        index = np.argmax(SMRwork)
        bitsNeededPerBand[index] = bitsNeededPerBand[index] + 1 
        bitsUsed = bitsUsed + nLines[index]
        SMRwork[index] = SMRwork[index] - 6.02

    bitsNeededPerBand = np.where(bitsNeededPerBand > 1, bitsNeededPerBand , 0)
    bitsNeededPerBand = np.where(bitsNeededPerBand < maxMantBits, bitsNeededPerBand , maxMantBits)
    return bitsNeededPerBand.astype(int)









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
    #return zeros(nBands) # TO REPLACE WITH YOUR CODE


    avgSMR = np.mean(SMR)
    nBins = np.sum(nLines)
    SMRwork = SMR.copy()
    
    
    bitsNeededPerBand = np.round(bitBudget/nBins + (1./ 6.02 ) * (SMR - avgSMR))
    
    bitsNeededPerBand = np.where(bitsNeededPerBand > 1, bitsNeededPerBand , 0)
    bitsNeededPerBand = np.where(bitsNeededPerBand < maxMantBits, bitsNeededPerBand , maxMantBits)
    
    bitsUsed = np.sum(bitsNeededPerBand*nLines)
    
    while (bitsUsed <bitBudget):
        index = np.argmax(SMRwork)
        bitsNeededPerBand[index] = bitsNeededPerBand[index] + 1 
        bitsUsed = bitsUsed + nLines[index]
        SMRwork[index] = SMRwork[index] - 6.02
    
    
    bitsNeededPerBand = np.where(bitsNeededPerBand > 1, bitsNeededPerBand , 0)
    bitsNeededPerBand = np.where(bitsNeededPerBand < maxMantBits, bitsNeededPerBand , maxMantBits)
    
    return bitsNeededPerBand.astype(int)








#-----------------------------------------------------------------------------

# My helper functions
    
 # calculate peak SPL
def calcPeakSPLPerBand(MDCTdata, sfBands, nBands ):
    peakSPL = np.zeros((nBands,))
    
    for i in range(0, nBands):
        peakSPL[i] = np.max( ps.SPL(MDCTdata[sfBands.lowerLine[i] : sfBands.upperLine[i] ]) )
     
     

    return peakSPL


def makeNoiseFloor(SPL, sfBands, bits):
    
    nf = np.zeros((np.sum(sfBands.nLines),))
    
    for i in range(0, len(SPL)):
        nf[sfBands.lowerLine[i]: sfBands.upperLine[i]+1] = SPL[i] - 6.02 * bits[i]
        
    return nf    
    
    
    
    
    
    

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    #pass # TO REPLACE WITH YOUR CODE



    # 1.c) 
    Fs = 48000
    n = np.arange(0,Fs)
    N = 2048

    x = 0.54 * np.cos(2*np.pi*440*n/Fs) + 0.21 * np.cos(2*np.pi*550*n/Fs) + 0.10 * np.cos(2*np.pi*660*n/Fs) + 0.06 * np.cos(2*np.pi*880*n/Fs) + 0.03 * np.cos(2*np.pi*4400*n/Fs) + 0.01 * np.cos(2*np.pi*8800*n/Fs)
    
    x_2048 = x[0:N]
    MDCTscale = 0
    
    MDCTdata = (2**MDCTscale) * mdct.MDCT(SineWindow(x_2048), N/2, N/2, isInverse=False)
    MDCTdata_true = (2**MDCTscale) * mdct_true.MDCT(SineWindow(x_2048), N/2, N/2, isInverse=False)
    MDCTLines = ps.AssignMDCTLinesFromFreqLimits(N/2, Fs)
    sfBands = ps.ScaleFactorBands(MDCTLines)
    
    bitBudget = 2526     # calculated by hand for 128kb/s/ch
    maxMantBits = 2**4 - 1
    nBands = 25
    
    # test uniformly distribute 
    uniformBitAlloc = BitAllocUniform(bitBudget, maxMantBits, nBands, sfBands.nLines, SMR=None)
    uniformBitAlloc_true = bitalloc_true.BitAllocUniform(bitBudget, maxMantBits, nBands, sfBands.nLines, SMR=None)
    
    
    
    # test peakSPL, 6dB below
    peakSPL = calcPeakSPLPerBand(4*MDCTdata**2, sfBands, nBands )
    
    constSNRBitAlloc = BitAllocConstSNR(bitBudget, maxMantBits, nBands, sfBands.nLines, peakSPL)
    constSNRBitAlloc_true = bitalloc_true.BitAllocConstSNR(bitBudget, maxMantBits, nBands, sfBands.nLines, peakSPL)
    
    
    
    
    
    # Test conctant SMR
    SMR = ps_true.CalcSMRs(x_2048,  MDCTdata, MDCTscale, Fs, sfBands )
    combinedMasking = ps.getMaskedThreshold(x_2048, 0, 0, Fs, sfBands)
    
    constMNRBitAlloc = BitAllocConstMNR(bitBudget, maxMantBits, nBands, sfBands.nLines, SMR)
    constMNRBitAlloc_true = bitalloc_true.BitAllocConstMNR(bitBudget, maxMantBits, nBands, sfBands.nLines, SMR)
    
    
    
    
    
    
    
    
    
    # testing Bitallocations
    #testConstSNR = testbitalloc.TestBitAlloc(constSNRBitAlloc , bitBudget, maxMantBits, nBands, sfBands.nLines, SMR, peakSPL)
    
    
    
    
    
    
    # Make plots for 1.c)
    
    f_2048 = np.linspace(0,Fs,2048)
    MDCTspl = ps.SPL(4*MDCTdata**2)
    combinedMasking = ps.getMaskedThreshold(x_2048, 0, 0, Fs, sfBands)
    
    
    
    
    nfUniform = makeNoiseFloor(peakSPL, sfBands, uniformBitAlloc)
    nfConstSNR = makeNoiseFloor(peakSPL, sfBands, constSNRBitAlloc)
    nfConstNMR = makeNoiseFloor(peakSPL, sfBands, constMNRBitAlloc)
    
    
#    fig1 = plt.figure()
#    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
#    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
#    plt.semilogx(f_2048[0:1024], nfUniform, '-', label = 'bit allocation ')
#    
#    axes = plt.gca()
#    axes.set_ylim([-20,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.c) Uniformly Distributed ')
#    plt.grid(True)
#    plt.legend() 
#    plt.savefig("Q1c_uni.png")
#    
#    
#    
#    fig2 = plt.figure()
#    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
#    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
#    plt.semilogx(f_2048[0:1024], nfConstSNR, '-', label = 'bit allocation ')
#    
#    axes = plt.gca()
#    axes.set_ylim([-20,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.c) Conctant SNR ')
#    plt.grid(True)
#    plt.legend() 
#    plt.savefig("Q1c_SNR.png")
#    
#    
#    
#    fig3 = plt.figure()
#    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
#    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
#    plt.semilogx(f_2048[0:1024], nfConstNMR, '-', label = 'bit allocation ')
#    
#    axes = plt.gca()
#    axes.set_ylim([-20,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.c) Constant NMR ')
#    plt.grid(True)
#    plt.legend() 
#    plt.savefig("Q1c_NMR.png")
    
    
    
    
    
    
    
    
    
    
    # 1. d) 
    
    bitBudget_192 = 3892     # calculated by hand for 128kb/s/ch
    
    # test uniformly distribute 
    uniformBitAlloc_192 = BitAllocUniform(bitBudget_192, maxMantBits, nBands, sfBands.nLines, SMR=None)
    uniformBitAlloc_192_true = bitalloc_true.BitAllocUniform(bitBudget_192, maxMantBits, nBands, sfBands.nLines, SMR=None)
    
    
    
    # test peakSPL, 6dB below
    peakSPL_192 = calcPeakSPLPerBand(4*MDCTdata**2, sfBands, nBands )
    
    constSNRBitAlloc_192 = BitAllocConstSNR(bitBudget_192, maxMantBits, nBands, sfBands.nLines, peakSPL_192)
    constSNRBitAlloc_192_true = bitalloc_true.BitAllocConstSNR(bitBudget_192, maxMantBits, nBands, sfBands.nLines, peakSPL_192)
    
    
    
    
    
    # Test conctant SMR
    SMR_192 = ps_true.CalcSMRs(x_2048,  MDCTdata, MDCTscale, Fs, sfBands )
    combinedMasking_192 = ps.getMaskedThreshold(x_2048, 0, 0, Fs, sfBands)
    
    constMNRBitAlloc_192 = BitAllocConstMNR(bitBudget_192, maxMantBits, nBands, sfBands.nLines, SMR)
    constMNRBitAlloc_192_true = bitalloc_true.BitAllocConstMNR(bitBudget_192, maxMantBits, nBands, sfBands.nLines, SMR_192)
    
    
    
    
    
    
        # Make plots for 1.d)
    
    f_2048 = np.linspace(0,Fs,2048)
    MDCTspl = ps.SPL(4*MDCTdata**2)
    combinedMasking = ps.getMaskedThreshold(x_2048, 0, 0, Fs, sfBands)
    
    
    
    
    nfUniform_192 = makeNoiseFloor(peakSPL_192, sfBands, uniformBitAlloc_192)
    nfConstSNR_192 = makeNoiseFloor(peakSPL_192, sfBands, constSNRBitAlloc_192)
    nfConstNMR_192 = makeNoiseFloor(peakSPL_192, sfBands, constMNRBitAlloc_192)
    
    
#    fig4 = plt.figure()
#    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
#    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
#    plt.semilogx(f_2048[0:1024], nfUniform_192, '-', label = 'bit allocation ')
#    
#    axes = plt.gca()
#    axes.set_ylim([-20,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.d) Uniformly Distributed ')
#    plt.grid(True)
#    plt.legend() 
#    plt.savefig("Q1d_uni.png")
#    
#    
#    
#    fig5 = plt.figure()
#    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
#    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
#    plt.semilogx(f_2048[0:1024], nfConstSNR_192, '-', label = 'bit allocation ')
#    
#    axes = plt.gca()
#    axes.set_ylim([-20,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.d) Conctant SNR ')
#    plt.grid(True)
#    plt.legend() 
#    plt.savefig("Q1d_SNR.png")
#    
#    
#    
#    fig6 = plt.figure()
#    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
#    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
#    plt.semilogx(f_2048[0:1024], nfConstNMR_192, '-', label = 'bit allocation ')
#    
#    axes = plt.gca()
#    axes.set_ylim([-20,100])
#    axes.set_xlim([50,Fs/2])
#         
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('SPL (dB)')
#    plt.title('1.d) Constant NMR ')
#    plt.grid(True)
#    plt.legend() 
#    plt.savefig("Q1d_NMR.png")
    
    
    
    
    
    
    
    
    
    
    
    # 2.a)
    
    
    
    bitBudget_128 = 2526
    bitBudget_192 = 3892
    
    
    bitAlloc_2a_128 = BitAlloc(bitBudget_128, maxMantBits, nBands, sfBands.nLines, SMR)
    
    
    bitAlloc_2a_192 = BitAlloc(bitBudget_192, maxMantBits, nBands, sfBands.nLines, SMR)
    
    nfOpt_128 = makeNoiseFloor(peakSPL_192, sfBands, bitAlloc_2a_128)
    nfOpt_192 = makeNoiseFloor(peakSPL_192, sfBands, bitAlloc_2a_192)
    
    fig7 = plt.figure()
    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
    plt.semilogx(f_2048[0:1024], nfOpt_128, '-', label = 'bit allocation ')
    
    axes = plt.gca()
    axes.set_ylim([-20,100])
    axes.set_xlim([50,Fs/2])
         
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('SPL (dB)')
    plt.title('2.b) 128 Bits Optimal ')
    plt.grid(True)
    plt.legend() 
    plt.savefig("Q2a_128.png")
    
    
    fig8 = plt.figure()
    plt.semilogx(f_2048[0:1024], MDCTspl,'.-', label = 'MDCT SPL')
    plt.semilogx(f_2048[0:1024], combinedMasking, '-', label = 'Masked Threshold')
    plt.semilogx(f_2048[0:1024], nfOpt_192, '-', label = 'bit allocation ')
    
    axes = plt.gca()
    axes.set_ylim([-20,100])
    axes.set_xlim([50,Fs/2])
         
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('SPL (dB)')
    plt.title('2.b) 192 Bits Optimal ')
    plt.grid(True)
    plt.legend() 
    plt.savefig("Q2a_192.png")
    
    
    
    
    
    
    
    
    
    #IMDCTdata = (2**(-MDCTscale)) * mdct.MDCT(SineWindow(MDCTdata_quant, N/2, N/2, isInverse=true)



























