"""
codec.py -- The actual encode/decode functions for the perceptual audio codec

-----------------------------------------------------------------------
© 2009 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
"""

import numpy as np  # used for arrays

# used by Encode and Decode
import window as win # current window used for MDCT
# from window import KBDWindow as win # current window used for MDCT
from mdct import MDCT,IMDCT  # fast MDCT implementation (uses numpy FFT)
from quantize import *  # using vectorized versions (to use normal versions, uncomment lines 18,67 (30,79) below defining vMantissa and vDequantize)
# used only by Encode
from psychoac import CalcSMRs  # calculates SMRs for each scale factor band
from bitalloc_ import BitAlloc  #allocates bits to scale factor bands given SMRs


def Decode(scaleFactor,bitAlloc,mantissa,overallScaleFactor,codingParams):
    """Reconstitutes a single-channel block of encoded data into a block of
    signed-fraction data based on the parameters in a PACFile object"""

    rescaleLevel = 1.*(1<<overallScaleFactor)
    
    #determine correct size of halfN
    if codingParams.win_state == 0 :
        codingParams.a = codingParams.nMDCTLinesLong
        codingParams.b = codingParams.nMDCTLinesLong
        #f_alpha = l_alpha = 4.
    
    #next block is a short block
    elif codingParams.win_state == 1:
        codingParams.a = codingParams.nMDCTLinesLong
        codingParams.b = codingParams.nMDCTLinesShort
        #f_alpha = 4.
        #l_alpha = 6.
    
    #next block is a stop transition block
    elif codingParams.win_state == 2:
        codingParams.a = codingParams.nMDCTLinesShort
        codingParams.b = codingParams.nMDCTLinesShort
        #f_alpha = l_alpha = 6.
 
    #next block is a short block after a 
    elif codingParams.win_state == 3:
        codingParams.a = codingParams.nMDCTLinesShort
        codingParams.b = codingParams.nMDCTLinesLong
        #f_alpha = 6.
        #l_alpha = 4.
 
    else:
        raise ValueError('Unknown window state:' + str(codingParams.win_state))
        
        
    halfN = (codingParams.a + codingParams.b)/2
    N = 2*halfN
    # vectorizing the Dequantize function call
    #vDequantize = np.vectorize(Dequantize)

    # reconstitute the first halfN MDCT lines of this channel from the stored data
    mdctLine = np.zeros(halfN,dtype=np.float64)
    iMant = 0
    for iBand in range(codingParams.sfBands.nBands):
        nLines =codingParams.sfBands.nLines[iBand]
        if bitAlloc[iBand]:
            mdctLine[iMant:(iMant+nLines)]=vDequantize(scaleFactor[iBand], mantissa[iMant:(iMant+nLines)],codingParams.nScaleBits, bitAlloc[iBand])
        iMant += nLines
    mdctLine /= rescaleLevel  # put overall gain back to original level


    # IMDCT and window the data for this channel
    mdctData =  IMDCT(mdctLine, codingParams.a, codingParams.b)
    data = win.compose_kbd_window(mdctData, codingParams.a, codingParams.b, 4., 4.)# takes in halfN MDCT coeffs
    #data = win.compose_sine_window(mdctData, codingParams.a, codingParams.b)
        
    return data


def Encode(data,codingParams):
    """Encodes a multi-channel block of signed-fraction data based on the parameters in a PACFile object"""
    scaleFactor = []
    bitAlloc = []
    mantissa = []
    overallScaleFactor = []

    # loop over channels and separately encode each one
    for iCh in range(codingParams.nChannels):
        (s,b,m,o) = EncodeSingleChannel(data[iCh],codingParams)
        scaleFactor.append(s)
        bitAlloc.append(b)
        mantissa.append(m)
        overallScaleFactor.append(o)
    # return results bundled over channels
    return (scaleFactor,bitAlloc,mantissa,overallScaleFactor)


def EncodeSingleChannel(data,codingParams):
    """Encodes a single-channel block of signed-fraction data based on the parameters in a PACFile object"""

    # window data for side chain FFT and also window and compute MDCT
    timeSamples = data
    
    if codingParams.win_state == 0 :
        codingParams.a = codingParams.nMDCTLinesLong
        codingParams.b = codingParams.nMDCTLinesLong
        #f_alpha = l_alpha = 4.
    
    #next block is a short block
    elif codingParams.win_state == 1:
        codingParams.a = codingParams.nMDCTLinesLong
        codingParams.b = codingParams.nMDCTLinesShort
#        f_alpha = 4
#        l_alpha = 6.
    
    #next block is a stop transition block
    elif codingParams.win_state == 2:
        codingParams.a = codingParams.nMDCTLinesShort
        codingParams.b = codingParams.nMDCTLinesShort
#        f_alpha = l_alpha = 6.
 
    #next block is a short block after a 
    elif codingParams.win_state == 3:
        codingParams.a = codingParams.nMDCTLinesShort
        codingParams.b = codingParams.nMDCTLinesLong
#        f_alpha = 6.
#        l_alpha = 4.
 
    else:
        raise ValueError('Unknown window state:' + str(codingParams.win_state))
        
    # prepare various constants
    halfN = (codingParams.a + codingParams.b)/2
    N = 2*halfN
    nScaleBits = codingParams.nScaleBits
    maxMantBits = (1<<codingParams.nMantSizeBits)  # 1 isn't an allowed bit allocation so n size bits counts up to 2^n
    if maxMantBits>16: maxMantBits = 16  # to make sure we don't ever overflow mantissa holders
    sfBands = codingParams.sfBands
    
    # vectorizing the Mantissa function call
    #vMantissa = np.vectorize(Mantissa)


    # compute target mantissa bit budget for this block of halfN MDCT mantissas
    bitBudget = codingParams.targetBitsPerSample * halfN  # this is overall target bit rate
    bitBudget -=  nScaleBits*(sfBands.nBands +1)  # less scale factor bits (including overall scale factor)
    bitBudget -= codingParams.nMantSizeBits*sfBands.nBands  # less mantissa bit allocation bits

        
    mdctTimeSamples = win.compose_kbd_window(data, codingParams.a, codingParams.b, 4., 4.)
    #mdctTimeSamples = win.compose_sine_window(data, codingParams.a, codingParams.b)
    mdctLines = MDCT(mdctTimeSamples, codingParams.a, codingParams.b)

    # compute overall scale factor for this block and boost mdctLines using it
    maxLine = np.max( np.abs(mdctLines) )
    overallScale = ScaleFactor(maxLine,nScaleBits)  #leading zeroes don't depend on nMantBits
    mdctLines *= (1<<overallScale)


    # compute the mantissa bit allocations
    # compute SMRs in side chain FFT
    SMRs = CalcSMRs(timeSamples, mdctLines, overallScale, codingParams.sampleRate, sfBands)
    # perform bit allocation using SMR results
    bitAlloc = BitAlloc(bitBudget, maxMantBits, sfBands.nBands, sfBands.nLines, SMRs)

    # given the bit allocations, quantize the mdct lines in each band
    scaleFactor = np.empty(sfBands.nBands,dtype=np.int32)
    
    nMant=halfN
    for iBand in range(sfBands.nBands):
        if not bitAlloc[iBand]: nMant-= sfBands.nLines[iBand]  # account for mantissas not being transmitted
    mantissa=np.empty(nMant,dtype=np.int32)
    iMant=0
    for iBand in range(sfBands.nBands):
        lowLine = sfBands.lowerLine[iBand]
        highLine = sfBands.upperLine[iBand] + 1  # extra value is because slices don't include last value
        nLines= sfBands.nLines[iBand]
        
        #if there are no lines in the given critical band
        if(highLine - lowLine <= 0):
            scaleFactor[iBand] = 0
            mantissa[iMant:iMant+nLines] = 0
        else:
            mdctVals = mdctLines[lowLine:highLine]
            scaleLine = np.max(np.abs( mdctVals ) )
            scaleFactor[iBand] = ScaleFactor(scaleLine, nScaleBits, bitAlloc[iBand])
            
            if bitAlloc[iBand]:
                mantissa[iMant:iMant+nLines] = vMantissa(mdctVals,scaleFactor[iBand], nScaleBits, bitAlloc[iBand])
                iMant += nLines
    # end of loop over scale factor bands

    # return results
    return (scaleFactor, bitAlloc, mantissa, overallScale)



