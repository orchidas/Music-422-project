"""
codec.py -- The actual encode/decode functions for the perceptual audio codec

-----------------------------------------------------------------------
© 2009 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
"""

import numpy as np  # used for arrays

# used by Encode and Decode
#from window import SineWindow as win # current window used for MDCT
import window as win
# from window import KBDWindow as win # current window used for MDCT
from mdct import MDCT,IMDCT  # fast MDCT implementation (uses numpy FFT)
from quantize import *  # using vectorized versions (to use normal versions, uncomment lines 18,67 (30,79) below defining vMantissa and vDequantize)
# used only by Encode
from psychoac import CalcSMRs, AssignMDCTLinesFromFreqLimits  # calculates SMRs for each scale factor band
from bitalloc_ import BitAlloc  #allocates bits to scale factor bands given SMRs


def Decode(scaleFactor,bitAlloc,mantissa,overallScaleFactor,codingParams):
    """Reconstitutes a single-channel block of encoded data into a block of
    signed-fraction data based on the parameters in a PACFile object"""

    rescaleLevel = 1.*(1<<overallScaleFactor)
    halfN = codingParams.nMDCTLines
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
    #data = win( IMDCT(mdctLine, halfN, halfN) )  # takes in halfN MDCT coeffs
    half_long = 512
    half_short = 64
    
    #correct inverse MDCT based on window state
    if codingParams.win_state == 0:
        mdctLine = IMDCT(mdctLine, half_long, half_long)
        data = win.compose_kbd_window(mdctLine, half_long, half_long, 4., 4.)
    
    elif codingParams.win_state == 1:
        mdctLine = IMDCT(mdctLine, half_long, half_short)
        data = win.compose_kbd_window(mdctLine, half_long, half_short, 4., 4.)
       
    elif codingParams.win_state == 2:
        mdctLine = IMDCT(mdctLine, half_short, half_short)
        data = win.compose_kbd_window(mdctLine, half_short, half_short, 4., 4.)
        
    elif codingParams.win_state == 3:
        mdctLine = IMDCT(mdctLine, half_short, half_long)
        data = win.compose_kbd_window(mdctLine, half_short, half_long, 4., 4.)
    
    else:
        raise ValueError('Unknown window state:' + str(codingParams.win_state))

    # end loop over channels, return reconstituted time samples (pre-overlap-and-add)
    return data


def Encode(data,codingParams):
    """Encodes a multi-channel block of signed-fraction data based on the parameters in a PACFile object"""
    scaleFactor = []
    bitAlloc = []
    mantissa = []
    overallScaleFactor = []
    #window_state = []

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

    # prepare various constants
    halfN = codingParams.nMDCTLines
    N = 2*halfN
    nScaleBits = codingParams.nScaleBits
    maxMantBits = (1<<codingParams.nMantSizeBits)  # 1 isn't an allowed bit allocation so n size bits counts up to 2^n
    if maxMantBits>16: maxMantBits = 16  # to make sure we don't ever overflow mantissa holders
    sfBands = codingParams.sfBands
    # vectorizing the Mantissa function call
    #vMantissa = np.vectorize(Mantissa)
    
    # window data for side chain FFT and also window and compute MDCT based on right windpw state
    timeSamples = data
    half_long = 512
    half_short = 64
    
    if codingParams.win_state == 0:
        halfN = half_long
        mdctTimeSamples = win.compose_kbd_window(data, half_long, half_long, 4., 4.)
        mdctLines = MDCT(mdctTimeSamples, half_long, half_long)[:halfN]
 
    elif codingParams.win_state == 1:
        halfN = (half_long + half_short)/2
        mdctTimeSamples = win.compose_kbd_window(data, half_long, half_short, 4., 4.)
        mdctLines = MDCT(mdctTimeSamples, half_long, half_short)[:halfN]
 
    elif codingParams.win_state == 2:
        halfN = half_short
        mdctTimeSamples = win.compose_kbd_window(data, half_short, half_short, 4., 4.)
        mdctLines = MDCT(mdctTimeSamples, half_short, half_short)[:halfN]
 
    elif codingParams.win_state == 3:
        halfN = (half_long + half_short)/2
        mdctTimeSamples = win.compose_kbd_window(data, half_short, half_long, 4., 4.)
        mdctLines = MDCT(mdctTimeSamples, half_short, half_long)[:halfN]
 
    else:
        raise ValueError('Unknown window state:' + str(codingParams.win_state))
        
    #make sure you have the right sfBands values
    sfBands.nLines = AssignMDCTLinesFromFreqLimits(halfN, codingParams.sampleRate)
    sfBands.upperLine = np.cumsum(sfBands.nLines)-1
    sfBands.lowerLine = codingParams.sfBands.upperLine - sfBands.nLines + 1

    # compute target mantissa bit budget for this block of halfN MDCT mantissas
    bitBudget = codingParams.targetBitsPerSample * halfN  # this is overall target bit rate
    bitBudget -=  nScaleBits*(sfBands.nBands +1)  # less scale factor bits (including overall scale factor)
    bitBudget -= codingParams.nMantSizeBits*sfBands.nBands  # less mantissa bit allocation bits

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
            mdctVals = 0
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



