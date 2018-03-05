import numpy as np

# Question 1.b)
def BitAllocUniform(bitBudget, maxMantBits, nBands, nLines, SMR=None):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are uniformely distributed for the mantissas.
    """

    uniform = np.ones(np.size(nLines)) * np.floor(bitBudget/np.sum(nLines))

    k = 0
 
    while True:

        totalLevel = np.sum(uniform*nLines) + nLines[k]
            
        if totalLevel <=bitBudget:
                
                uniform[k] += 1
                k += 1
                k %= nBands
                
        else:

            return uniform 


# Or all the values in ana array and return that value 
def ifAllTrue(bands):
    
    OrAllVals = False
    
    for i in range(len(bands)):
        
        OrAllVals  = OrAllVals or bands[i]
        
    
    return OrAllVals
    
    
    
def BitAllocConstSNR(bitBudget, maxMantBits, nBands, nLines, peakSPL):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep a constant
    quantization noise floor (assuming a noise floor 6 dB per bit below
    the peak SPL line in the scale factor band).
    """
    
    constSNRbits = np.zeros(nBands , dtype = int)
    
    bitBudgetEachBand = [True]*nBands
    
    while ifAllTrue(bitBudgetEachBand):
        
        whichBandsHaveBudget = np.arange(nBands)[bitBudgetEachBand]

        
        truePeaks = (peakSPL - constSNRbits*6.02)[bitBudgetEachBand]

    
        maxIndex = np.argmax(truePeaks)

        
        if bitBudget - nLines[whichBandsHaveBudget[maxIndex]] < 0 :

            
            bitBudgetEachBand[whichBandsHaveBudget[maxIndex]] = False

            
        else : 
            
            constSNRbits[whichBandsHaveBudget[maxIndex]] += 1
            
            bitBudget -= nLines[whichBandsHaveBudget[maxIndex]]
            
            bitBudgetEachBand[whichBandsHaveBudget[maxIndex]] =  np.signbit(constSNRbits[whichBandsHaveBudget[maxIndex]] - maxMantBits)
            

    return constSNRbits
    

def BitAllocConstMNR(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep the quantization
    noise floor a constant distance below (or above, if bit starved) the
    masked threshold curve (assuming a quantization noise floor 6 dB per
    bit below the peak SPL line in the scale factor band).
    """

    
    constMNRbits = np.zeros(nBands , dtype = int)
    
    bitBudgetEachBand = [True]*nBands
    
    while ifAllTrue(bitBudgetEachBand):
        
        whichBandsHaveBudget = np.arange(nBands)[bitBudgetEachBand]
        
        truePeaks = (SMR - constMNRbits*6.02)[bitBudgetEachBand]
    
        maxIndex = np.argmax(truePeaks)
        
        if bitBudget - nLines[whichBandsHaveBudget[maxIndex]] < 0 :
            
            bitBudgetEachBand[whichBandsHaveBudget[maxIndex]] = False
            
        else : 
            
            constMNRbits[whichBandsHaveBudget[maxIndex]] += 1
            
            bitBudget -= nLines[whichBandsHaveBudget[maxIndex]]
            
            bitBudgetEachBand[whichBandsHaveBudget[maxIndex]] =  np.signbit(constMNRbits[whichBandsHaveBudget[maxIndex]] - maxMantBits)
            

    return constMNRbits


def BitAllocStereo(bitBudget,maxMantBits, nBands, nLines, SMR, LRorMS , bitReservoir):
    
    bits = np.zeros(nBands, dtype=int)
    valid = np.ones(nBands, dtype=bool)
    totalBits = int(bitBudget + bitReservoir)
    MS_Threshold = -5.
    LR_Threshold = -10.

    while valid.any():
        iMax = np.arange(nBands)[valid][np.argmax((SMR-bits*6.)[valid])]
        if LRorMS[iMax]: # MS
            if max(SMR-(bits - 1)*6.)<(MS_Threshold): valid[iMax] = False
        else: # LR
            if max(SMR-(bits - 1)*6.)<(LR_Threshold): valid[iMax] = False
        # print max(SMR-(bits-1)*6.)
        if (totalBits - nLines[iMax]) >=0:
            bits[iMax] += 1
            totalBits -= nLines[iMax]
            if bits[iMax] >= maxMantBits:
                valid[iMax] = False
        else:
            valid[iMax] = False


    totalBits+=sum(nLines[bits==1])
    bits[bits==1]=0

    bitDifference=totalBits - bitReservoir

    return bits,bitDifference

#Testing code
if __name__ == "__main__":


    pass