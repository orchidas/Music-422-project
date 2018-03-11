"""
codec.py -- The actual encode/decode functions for the perceptual audio codec

-----------------------------------------------------------------------
© 2009 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
"""



import numpy as np  # used for arrays
from stereo import*
# used by Encode and Decode
import window as win # current window used for MDCT
# from window import KBDWindow as win # current window used for MDCT
from mdct import MDCT,IMDCT  # fast MDCT implementation (uses numpy FFT)
from quantize import *  # using vectorized versions (to use normal versions, uncomment lines 18,67 below defining vMantissa and vDequantize)
# used only by Encode
from psychoac import CalcSMRs  # calculates SMRs for each scale factor band
from bitalloc import BitAllocStereo  #allocates bits to scale factor bands given SMRs
import matplotlib.pyplot as plt


def Decode(scaleFactor,bitAlloc,mantissa,overallScaleFactor,codingParams,LRorMS):
    """Reconstitutes a single-channel block of encoded data into a block of
    signed-fraction data based on the parameters in a PACFile object"""

    mdctDecoded=[]

    if codingParams.win_state == 0 :
        codingParams.a = codingParams.nMDCTLinesLong
        codingParams.b = codingParams.nMDCTLinesLong
        f_alpha = l_alpha = 4.
    
    #next block is a short block
    elif codingParams.win_state == 1:
        codingParams.a = codingParams.nMDCTLinesLong
        codingParams.b = codingParams.nMDCTLinesShort
        f_alpha = 4.
        l_alpha = 6.
    
    #next block is a stop transition block
    elif codingParams.win_state == 2:
        codingParams.a = codingParams.nMDCTLinesShort
        codingParams.b = codingParams.nMDCTLinesShort
        f_alpha = l_alpha = 6.
 
    #next block is a short block after a 
    elif codingParams.win_state == 3:
        codingParams.a = codingParams.nMDCTLinesShort
        codingParams.b = codingParams.nMDCTLinesLong
        f_alpha = 6.
        l_alpha = 4.
 
    else:
        raise ValueError('Unknown window state:' + str(codingParams.win_state))
        
    halfN = (codingParams.a + codingParams.b)/2
    N = 2*halfN


    for iCh in range(codingParams.nChannels):

        rescaleLevel = 1.*(1<<overallScaleFactor[iCh])
        # halfN = codingParams.nMDCTLines


        # reconstitute the first halfN MDCT lines of this channel from the stored data

        mdctDecoded.append(np.zeros(halfN,dtype=np.float64))
        iMant = 0


        for iBand in range(codingParams.sfBands.nBands):
            nLines =codingParams.sfBands.nLines[iBand]
            if bitAlloc[iCh][iBand]:
                mdctDecoded[iCh][iMant:(iMant+nLines)]=vDequantize(scaleFactor[iCh][iBand], mantissa[iCh][iMant:(iMant+nLines)],codingParams.nScaleBits, bitAlloc[iCh][iBand])
            iMant += nLines
        mdctDecoded[iCh] /= rescaleLevel  # put overall gain back to original level


    # recombine into L and R only
    mdctLineL=mdctDecoded[0]
    mdctLineR=mdctDecoded[1]


    for iBand in range(codingParams.sfBands.nBands):
        if LRorMS[iBand]:
            lowLine = codingParams.sfBands.lowerLine[iBand]
            highLine = codingParams.sfBands.upperLine[iBand] + 1  


            # Reconstruction, L=M-S and R=M+S
            mdctLineL[lowLine:highLine]=mdctDecoded[0][lowLine:highLine] - mdctDecoded[1][lowLine:highLine]
            mdctLineR[lowLine:highLine]=mdctDecoded[0][lowLine:highLine]+ mdctDecoded[1][lowLine:highLine]

    print mdctLineL


    # IMDCT and window the data for each channel
    dataL = win.compose_sine_window(IMDCT(mdctLineL, codingParams.a, codingParams.b),codingParams.a, codingParams.b)  # takes in halfN MDCT coeffs
    dataR = win.compose_sine_window(IMDCT(mdctLineR, codingParams.a, codingParams.b),codingParams.a, codingParams.b)

    # print dataL
    # print dataR

    data=dataL,dataR

    # end loop over channels, return reconstituted time samples (pre-overlap-and-add)
    return data


def Encode(data,codingParams, myhuffyman):


    scaleFactor = []
    bitAlloc = []
    mantissa = []
    overallScaleFactor = []
    tableID = []
    mantissaSignBits = []
    huffmanMantissa = []

    threshold = 0.75

    sfBands = codingParams.sfBands


    fftL = np.fft.fft(data[0])
    fftR = np.fft.fft(data[1])

    LRorMS = np.zeros(sfBands.nBands , dtype = 'int')

    count = 0

    coherence = False

    for band in range(sfBands.nBands):

        
        lowLine = sfBands.lowerLine[band]
        highLine = sfBands.upperLine[band] + 1

        if(count >= 6):
            if coherence:

                LRorMS[band] = coherenceLR( fftL, fftR , lowLine , highLine, threshold)

            else :

                LRorMS[band] = corellationLR( fftL, fftR , lowLine , highLine)

        count += 1


    print " L/R or M/S : " + str(LRorMS) 

    (scaleFactor,bitAlloc,mantissa,overallScaleFactor) = EncodeTwoChannels(data, codingParams,LRorMS , myhuffyman)


    for iCh in range(codingParams.nChannels):
        
        (signBits, unsignedMantissas) = removeMantissaSignBits(codingParams,mantissa[iCh],bitAlloc[iCh])

        (codeHuffman,tabID) = myhuffyman.encodeHuffman(codingParams,unsignedMantissas,bitAlloc[iCh])
        print "Table Being Used : ", tabID
        bitsAvailable = sum(bitAlloc[iCh]*codingParams.sfBands.nLines)
        huffmanBitsRequired = sum(len(huffman) for huffman in codeHuffman) + len(signBits) + codingParams.numTableBits
        myhuffyman.giveBits(bitsAvailable-huffmanBitsRequired)
       
        huffmanMantissa.append(codeHuffman)
        tableID.append(tabID)
        mantissaSignBits.append(signBits)

    return (scaleFactor,bitAlloc,overallScaleFactor,LRorMS, mantissaSignBits,huffmanMantissa,tableID)


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
#    vMantissa = np.vectorize(Mantissa)

    # compute target mantissa bit budget for this block of halfN MDCT mantissas
    bitBudget = codingParams.targetBitsPerSample * halfN  # this is overall target bit rate
    bitBudget -=  nScaleBits*(sfBands.nBands +1)  # less scale factor bits (including overall scale factor)
    bitBudget -= codingParams.nMantSizeBits*sfBands.nBands  # less mantissa bit allocation bits


    # window data for side chain FFT and also window and compute MDCT
    timeData = data
    mdctTimeData = win(data)
    mdctLines = MDCT(mdctTimeData, halfN, halfN)[:halfN]

    # compute overall scale factor for this block and boost mdctLines using it
    maxLine = np.max( np.abs(mdctLines) )
    overallScale = ScaleFactor(maxLine,nScaleBits)  #leading zeroes don't depend on nMantBits
    mdctLines *= (1<<overallScale)



    # compute the mantissa bit allocations
    # compute SMRs in side chain FFT
    SMRs = CalcSMRs(timeData, mdctLines, overallScale, codingParams.sampleRate, sfBands)
    # perform bit allocation using SMR results

    # global bitAlloc
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
        scaleLine = np.max(np.abs( mdctLines[lowLine:highLine] ) )
        scaleFactor[iBand] = ScaleFactor(scaleLine, nScaleBits, bitAlloc[iBand])
        if bitAlloc[iBand]:
            mantissa[iMant:iMant+nLines] = vMantissa(mdctLines[lowLine:highLine],scaleFactor[iBand], nScaleBits, bitAlloc[iBand])
            iMant += nLines
    # end of loop over scale factor bands

    # return results
    return (scaleFactor, bitAlloc, mantissa, overallScale)


def removeMantissaSignBits(codingParams,mantissa,bitAlloc):
    mantissaSignBits = []
    unsignedMantissas = []

    iMant=0
    for iBand in range(codingParams.sfBands.nBands):
        if bitAlloc[iBand]:
            for j in range(codingParams.sfBands.nLines[iBand]):
                streamFormat = '0' + str(bitAlloc[iBand]) + 'b'
                mantissaStringFormat = format(mantissa[iMant+j],streamFormat)
                mantissaSignBits.append(int(mantissaStringFormat[0],2))
                unsignedMantissas.append(int(mantissaStringFormat[1:],2))
            iMant += codingParams.sfBands.nLines[iBand]

    return (mantissaSignBits, unsignedMantissas)

def EncodeTwoChannels(data,codingParams,LRorMS , myhuffyman):

    # timeSamples = data
    
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
    if maxMantBits > 16: maxMantBits = 16  # to make sure we don't ever overflow mantissa holders
    sfBands = codingParams.sfBands

    # compute target mantissa bit budget for this block of halfN MDCT mantissas
    bitBudget = codingParams.targetBitsPerSample * halfN  # this is overall target bit rate
    bitBudget -= nScaleBits*(sfBands.nBands +1)  # less scale factor bits (including overall scale factor)
    bitBudget -= codingParams.nMantSizeBits*sfBands.nBands  # less mantissa bit allocation bits
    '''Subtract the bits needed for the table ID'''
    bitBudget -= codingParams.numTableBits
    '''Get extra bits saved from myhuffyman encoding to do bit alloc'''

    codingParams.bitReservoir += myhuffyman.snatchBits()

    timeData=[]
    mdctTimeData=[]
    mdctLines=[]
    maxLine=[]
    overallScale=[]

    for iCh in range(codingParams.nChannels):
        # window data for side chain FFT and also window and compute MDCT
        timeData.append(data[iCh])

        if (len(data[iCh]) != (codingParams.a + codingParams.b)):

            print "Length not equal "

            codingParams.a = codingParams.b = int(len(data[iCh])/2)

        mdctTimeData.append(win.compose_sine_window(data[iCh], codingParams.a, codingParams.b))
        mdctLines.append(MDCT(mdctTimeData[iCh], codingParams.a, codingParams.b)[:halfN])



        # compute overall scale factor for this block and boost mdctLines using it
        maxLine.append(np.max( np.abs(mdctLines[iCh]) ) )
        overallScale.append(ScaleFactor(maxLine[iCh],nScaleBits) ) #leading zeroes don't depend on nMantBits
        mdctLines[iCh] *= (1<<overallScale[iCh])

    # compute the mantissa bit allocations
    # compute SMRs in side chain FFT

    # print mdctLines
    (SMRs,LRorMSmdctLines) = stereoMaskThresholds(timeData, mdctLines, overallScale, codingParams.sampleRate, sfBands, LRorMS, codingParams)

    bitAlloc=[]
    scaleFactor=[]
    mantissa=[]

    # perform bit allocation using SMR results
    for iCh in range(codingParams.nChannels):
        ba,bitDifference=BitAllocStereo(bitBudget,maxMantBits, sfBands.nBands, sfBands.nLines, SMRs[iCh], LRorMS , codingParams.bitReservoir)
        bitAlloc.append(ba)
        codingParams.bitReservoir+=bitDifference

        print "Extra Bits Available : " + str(codingParams.bitReservoir)
        print "Bits in Bit Reservoir : " + str(myhuffyman.getBitReservoir())

        # given the bit allocations, quantize the mdct lines in each band
        scaleFactor.append(np.empty(sfBands.nBands,dtype=np.int32))
        nMant=halfN
        for iBand in range(sfBands.nBands):
            if not bitAlloc[iCh][iBand]: nMant-= sfBands.nLines[iBand]  # account for mantissas not being transmitted
        mantissa.append(np.empty(nMant,dtype=np.int32))
        iMant=0
        for iBand in range(sfBands.nBands):
            lowLine = sfBands.lowerLine[iBand]
            highLine = sfBands.upperLine[iBand] + 1  
            nLines= sfBands.nLines[iBand]
            print "Lowline is : " + str(lowLine)
            print "HighLine is :" + str(highLine)

            if(highLine - lowLine <= 0 ):
                scaleFactor[iCh][iBand] = 0
                mantissa[iCh][iMant:iMant+nLines] = 0 
                # scaleLine = np.array([scaleLine])

            else: 

                scaleLine = np.max(np.abs( LRorMSmdctLines[iCh][lowLine:highLine] ))

                scaleFactor[iCh][iBand] = ScaleFactor(scaleLine, nScaleBits, bitAlloc[iCh][iBand])
                if bitAlloc[iCh][iBand]:
                    mantissa[iCh][iMant:iMant+nLines] = vMantissa(LRorMSmdctLines[iCh][lowLine:highLine],scaleFactor[iCh][iBand], nScaleBits, bitAlloc[iCh][iBand])
                    iMant += nLines
    
    return (scaleFactor, bitAlloc, mantissa, overallScale)



