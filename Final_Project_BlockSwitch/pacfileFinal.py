from audiofile import * # base class
from bitpack import *  # class for packing data into an array of bytes where each item's number of bits is specified
import codec # module where the actual PAC coding functions reside(this module only specifies the PAC file format)
from psychoac import ScaleFactorBands, AssignMDCTLinesFromFreqLimits  # defines the grouping of MDCT lines into scale factor bands
import sys

import numpy as np  # to allow conversion of data blocks to numpy's array object
from Huffman import *
from blockswitch import *
MAX16BITS = 32767

class PACFile(AudioFile):
    """
    Handlers for a perceptually coded audio file I am encoding/decoding
    """

    # a file tag to recognize PAC coded files
    tag='PAC '

    def ReadFileHeader(self):
        """
        Reads the PAC file header from a just-opened PAC file and uses it to set
        object attributes.  File pointer ends at start of data portion.
        """
        # check file header tag to make sure it is the right kind of file
        tag=self.fp.read(4)
        if tag!=self.tag: raise "Tried to read a non-PAC file into a PACFile object"
        # use struct.unpack() to load up all the header data
        (sampleRate, nChannels, numSamples, nMDCTLines, nScaleBits, nMantSizeBits) \
                 = unpack('<LHLLHH',self.fp.read(calcsize('<LHLLHH')))
        nBands = unpack('<L',self.fp.read(calcsize('<L')))[0]
        nLines=  unpack('<'+str(nBands)+'H',self.fp.read(calcsize('<'+str(nBands)+'H')))
        sfBands=ScaleFactorBands(nLines)
        # load up a CodingParams object with the header data
        myParams=CodingParams()
        myParams.sampleRate = sampleRate
        myParams.nChannels = nChannels
        myParams.numSamples = numSamples
        myParams.nMDCTLines = myParams.nSamplesPerBlock = nMDCTLines
        myParams.nScaleBits = nScaleBits
        myParams.nMantSizeBits = nMantSizeBits
        # add in scale factor band information
        myParams.sfBands =sfBands

        myParams.nMDCTLinesLong = 512
        myParams.nMDCTLinesShort = 64
        myParams.nMDCTLinesTrans = (512 + 64)/2
       
        nLinesShort = AssignMDCTLinesFromFreqLimits(myParams.nMDCTLinesShort, myParams.sampleRate)
        myParams.sfBandsShort = ScaleFactorBands(nLinesShort)
        nLinesTrans = AssignMDCTLinesFromFreqLimits(myParams.nMDCTLinesTrans, myParams.sampleRate)
        myParams.sfBandsTrans = ScaleFactorBands(nLinesTrans)
        nLinesLong = AssignMDCTLinesFromFreqLimits(myParams.nMDCTLinesLong, myParams.sampleRate)
        myParams.sfBandsLong = ScaleFactorBands(nLinesLong)
        
        myParams.a = 512
        myParams.b = 512

        

        # start w/o all zeroes as data from prior block to overlap-and-add for output
        overlapAndAdd = []
        for iCh in range(nChannels): overlapAndAdd.append( np.zeros(nMDCTLines, dtype=np.float64) )
        myParams.overlapAndAdd=overlapAndAdd
        return myParams

    def ReadDataBlock(self, codingParams,myhuffyman):
        """
        Reads a block of coded data from a PACFile object that has already
        executed OpenForReading() and returns those samples as reconstituted
        signed-fraction data
        """
        # loop over channels (whose coded data are stored separately) and read in each data block
        data=[]
        
        for iCh in range(codingParams.nChannels):
            data.append(np.array([],dtype=np.float64))  # add location for this channel's data
            # read in string containing the number of bytes of data for this channel (but check if at end of file!)
            s=self.fp.read(calcsize("<L"))  # will be empty if at end of file
            if not s:
                # hit last block, see if final overlap and add needs returning, else return nothing
                if codingParams.overlapAndAdd:
                    overlapAndAdd=codingParams.overlapAndAdd
                    codingParams.overlapAndAdd=0  # setting it to zero so next pass will just return
                    return overlapAndAdd
                else:
                    return
            # not at end of file, get nBytes from the string we just read
            nBytes = unpack("<L",s)[0] # read it as a little-endian unsigned long
            # read the nBytes of data into a PackedBits object to unpack
            pb = PackedBits()
            pb.SetPackedData( self.fp.read(nBytes) ) # PackedBits function SetPackedData() converts strings to internally-held array of bytes
            if pb.nBytes < nBytes:  raise "Only read a partial block of coded PACFile data"

            codingParams.win_state = pb.ReadBits(2)
            
            #window state is going to be same for both channels, we should only do this once
            if(iCh == 0):

                if codingParams.win_state == 0:
                    codingParams.nMDCTLines = codingParams.nMDCTLinesLong
                    codingParams.a = 512
                    codingParams.sfBands = codingParams.sfBandsLong
                    
                
                elif codingParams.win_state == 1:
                    codingParams.nMDCTLines = codingParams.nMDCTLinesTrans
                    codingParams.a = 64
                    codingParams.sfBands = codingParams.sfBandsTrans
                   
                
                elif codingParams.win_state == 2:
                    codingParams.nMDCTLines = codingParams.nMDCTLinesShort
                    codingParams.b = 64
                    codingParams.sfBands = codingParams.sfBandsShort
                
                elif codingParams.win_state == 3:
                    codingParams.nMDCTLines = codingParams.nMDCTLinesTrans
                    codingParams.sfBands = codingParams.sfBandsTrans
                    codingParams.b = 512
                
                else:
                    raise ValueError('Invalid window state: ' + str(codingParams.win_state))
    
    
                overallScaleFactor=np.zeros((codingParams.nChannels),dtype='int')
                scaleFactor=np.zeros((codingParams.nChannels,codingParams.sfBands.nBands),dtype='int')
                mantissa=np.zeros((codingParams.nChannels,codingParams.nMDCTLines),dtype='int')
                bitAlloc=np.zeros((codingParams.nChannels,codingParams.sfBands.nBands),dtype='int')
                LRMS=np.zeros(codingParams.sfBands.nBands,dtype='int')


            # extract the data from the PackedBits object
            overallScaleFactor[iCh] = pb.ReadBits(codingParams.nScaleBits)  # overall scale factor
            '''Extract the table ID used to myhuffyman encode for this block'''
            codingParams.numTableBits = 4
            tableID = pb.ReadBits(codingParams.numTableBits)
            # scaleFactor=[]
            # bitAlloc=[]
            # mantissa[iCh]=np.zeros(codingParams.nMDCTLines,np.int32)  # start w/ all mantissas zero
            for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                ba = pb.ReadBits(codingParams.nMantSizeBits)
                if ba: ba+=1  # no bit allocation of 1 so ba of 2 and up stored as one less
                bitAlloc[iCh][iBand]=ba  # bit allocation for this band
                scaleFactor[iCh][iBand]=pb.ReadBits(codingParams.nScaleBits)  # scale factor for this band
                if bitAlloc[iCh][iBand]:
                    # if bits allocated, extract those mantissas and put in correct location in matnissa array
                    '''Get the sign bits back for this band'''
                    mantissaSignBits = np.empty(codingParams.sfBands.nLines[iBand],np.int32)
                    for j in range(codingParams.sfBands.nLines[iBand]):
                        mantissaSignBits[j] = pb.ReadBits(1)  # one bit for each sign
                    m=np.empty(codingParams.sfBands.nLines[iBand],np.int32)
                    for j in range(codingParams.sfBands.nLines[iBand]):
                        '''Place to add myhuffyman decoding'''
                        m[j]= myhuffyman.decodeHuffman(pb,tableID,bitAlloc[iCh][iBand])
                        '''Recover signed mantissa code using the information of bitAlloc to see where to put the sign bit'''
                        m[j] += mantissaSignBits[j] * (2**(bitAlloc[iCh][iBand]-1))


                    mantissa[iCh][codingParams.sfBands.lowerLine[iBand]:(codingParams.sfBands.upperLine[iBand]+1)] = m
            # done unpacking data (end loop over scale factor bands)
           
           
            #print("Channel :", iCh, "Scale Factor :", scaleFactor[iCh])

            # CUSTOM DATA:
            # unpack LRMS
            for iBand in range(codingParams.sfBands.nBands):
                LRMS[iBand]=pb.ReadBits(1)

        # recombine into L and R
        # (DECODE HERE) decode the unpacked data for both channels

        #print("Scale Factor :", scaleFactor)
        decodedData = self.Decode(scaleFactor,bitAlloc,mantissa,overallScaleFactor,codingParams,LRMS,myhuffyman)

        # print decodedData

        for iCh in range(codingParams.nChannels):
            # overlap-and-add first half, and append it to the data array (saving other half for next overlap-and-add)
            data[iCh] = np.concatenate( (data[iCh], np.add(codingParams.overlapAndAdd[iCh],decodedData[iCh][:codingParams.a]) ) )  # data[iCh] is overlap-and-added data
            codingParams.overlapAndAdd[iCh] = decodedData[iCh][codingParams.a:]  # save other half for next pass

        # end loop over channels, return signed-fraction samples for this block
        return data

    def WriteFileHeader(self,codingParams):
        """
        Writes the PAC file header for a just-opened PAC file and uses codingParams
        attributes for the header data.  File pointer ends at start of data portion.
        """
        # write a header tag
        self.fp.write(self.tag)
        # make sure that the number of samples in the file is a multiple of the
        # number of MDCT half-blocksize, otherwise zero pad as needed
        if not codingParams.numSamples%codingParams.nMDCTLines:
            codingParams.numSamples += (codingParams.nMDCTLines
                        - codingParams.numSamples%codingParams.nMDCTLines) # zero padding for partial final PCM block

        # # also add in the delay block for the second pass w/ the last half-block (JH: I don't think we need this, in fact it generates a click at the end)
        # codingParams.numSamples+= codingParams.nMDCTLines  # due to the delay in processing the first samples on both sides of the MDCT block

        # write the coded file attributes
        self.fp.write(pack('<LHLLHH',
            codingParams.sampleRate, codingParams.nChannels,
            codingParams.numSamples, codingParams.nMDCTLines,
            codingParams.nScaleBits, codingParams.nMantSizeBits  ))
        # create a ScaleFactorBand object to be used by the encoding process and write its info to header
        # print codingParams.nMDCTLines
        sfBands=ScaleFactorBands( AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLines,
                                                                codingParams.sampleRate)
                                )
        codingParams.sfBands=sfBands
        # print codingParams.sampleRate
        # print vars(sfBands)
        self.fp.write(pack('<L',sfBands.nBands))
        self.fp.write(pack('<'+str(sfBands.nBands)+'H',*(sfBands.nLines.tolist()) ))
        # self.fp.write(pack('<'+str(sfBands.nBands)+'H',*(sfBands.nLines) ))

        # start w/o all zeroes as prior block of unencoded data for other half of MDCT block
        priorBlock = []
        for iCh in range(codingParams.nChannels):
            priorBlock.append(np.zeros(codingParams.nMDCTLines,dtype=np.float64) )
        codingParams.priorBlock = priorBlock
        codingParams.bitReservoir=0 
        codingParams.curBlock=0

        codingParams.win_state = 0
        codingParams.a = 512
        codingParams.b = 512
        codingParams.nMDCTLinesLong = 512
        codingParams.nMDCTLinesShort = 64
        codingParams.nMDCTLinesTrans = (512 + 64)/2
        nLinesShort = AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLinesShort, codingParams.sampleRate)
        codingParams.sfBandsShort = ScaleFactorBands(nLinesShort)
        nLinesTrans = AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLinesTrans, codingParams.sampleRate)
        codingParams.sfBandsTrans = ScaleFactorBands(nLinesTrans)
        nLinesLong = AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLinesLong, codingParams.sampleRate)
        codingParams.sfBandsLong = ScaleFactorBands(nLinesLong)


        return

    def WriteDataBlock(self,data, codingParams,myhuffyman):
        """
        Writes a block of signed-fraction data to a PACFile object that has
        already executed OpenForWriting()"""

        # combine this block of multi-channel data w/ the prior block's to prepare for MDCTs twice as long
        fullBlockData=[]
        for iCh in range(codingParams.nChannels):
            fullBlockData.append( np.concatenate( ( codingParams.priorBlock[iCh], data[iCh]) ) )
            codingParams.priorBlock[iCh] = np.copy(data[iCh])  # current pass's data is next pass's prior block data

        # (ENCODE HERE) Encode the full block of multi=channel data

        if(codingParams.win_state == 0):
            codingParams.sfBands = codingParams.sfBandsLong
        elif(codingParams.win_state == 2):
            codingParams.sfBands = codingParams.sfBandsShort
        else:
            codingParams.sfBands = codingParams.sfBandsTrans

        (scaleFactor,bitAlloc,overallScaleFactor,LRMS,mantissaSignBits,mantissa,tableID) = self.Encode(fullBlockData,codingParams,myhuffyman)  # returns a tuple with all the block-specific info not in the file header
        #print("Scale factor", scaleFactor)
        
        # for each channel, write the data to the output file
        for iCh in range(codingParams.nChannels):
            
            # determine the size of this channel's data block and write it to the output file
            nBytes = codingParams.nScaleBits  # bits for overall scale factor
            '''Add bits needed for informaing which myhuffyman table used'''
            nBytes += codingParams.numTableBits
            iMant=0
            for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to get its bits
                nBytes += codingParams.nMantSizeBits+codingParams.nScaleBits    # mantissa bit allocation and scale factor for that sf band
                if bitAlloc[iCh][iBand]:
                    # if non-zero bit allocation for this band, add in bits for scale factor and each mantissa (0 bits means zero)
                    '''Allocate bits for mantissaSignBits'''
                    nBytes += codingParams.sfBands.nLines[iBand]
                    #nBytes += bitAlloc[iCh][iBand]*codingParams.sfBands.nLines[iBand]  # no bit alloc = 1 so actuall alloc is one higher
                    '''For myhuffyman encoding'''
                    for j in range(codingParams.sfBands.nLines[iBand]):
                        codeBitLength = len(mantissa[iCh][iMant+j])
                        nBytes += codeBitLength
                    iMant += codingParams.sfBands.nLines[iBand]  # add to mantissa offset if we passed mantissas for this band
                                # end computing bits needed for this channel's data
            # end computing bits needed for this channel's data

            # CUSTOM DATA:
            # add space for LRMS array
            nBytes += len(LRMS)

            nBytes += 2

            # now convert the bits to bytes (w/ extra one if spillover beyond byte boundary)
            if nBytes%BYTESIZE==0:  nBytes /= BYTESIZE
            else: nBytes = nBytes/BYTESIZE + 1
            self.fp.write(pack("<L",int(nBytes))) # stores size as a little-endian unsigned long

            # create a PackedBits object to hold the nBytes of data for this channel/block of coded data
            pb = PackedBits()
            pb.Size(nBytes)

            pb.WriteBits(codingParams.win_state,2)

            
            # now pack the nBytes of data into the PackedBits object
            pb.WriteBits(overallScaleFactor[iCh],codingParams.nScaleBits)  # overall scale factor
            '''Write myhuffyman table ID'''
            pb.WriteBits(tableID[iCh],codingParams.numTableBits)
            iMant=0  # index offset in mantissa array (because mantissas w/ zero bits are omitted)
            for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                ba = bitAlloc[iCh][iBand]
                if ba: ba-=1  # if non-zero, store as one less (since no bit allocation of 1 bits/mantissa)
                pb.WriteBits(ba,codingParams.nMantSizeBits)  # bit allocation for this band (written as one less if non-zero)
                pb.WriteBits(scaleFactor[iCh][iBand],codingParams.nScaleBits)  # scale factor for this band (if bit allocation non-zero)
                if bitAlloc[iCh][iBand]:
                    '''Write bits for mantissaSignBits'''
                    for j in range(codingParams.sfBands.nLines[iBand]):
                        pb.WriteBits(mantissaSignBits[iCh][iMant+j],1) # one bit for each sign
                    for j in range(codingParams.sfBands.nLines[iBand]):
                        #pb.WriteBits(mantissa[iCh][iMant+j],bitAlloc[iCh][iBand])     # mantissas for this band (if bit allocation non-zero) and bit alloc <>1 so is 1 higher than the number
                        '''For myhuffyman encoding'''
                        codeBitLength = len(mantissa[iCh][iMant+j])

                        # print "*********" + str(mantissa[iCh][iMant+j])
                        pb.WriteBits(int(mantissa[iCh][iMant+j],2),codeBitLength)
                    iMant += codingParams.sfBands.nLines[iBand]  # add to mantissa offset if we passed mantissas for this band
        # done packing (end loop over scale factor bands)            # done packing (end loop over scale factor bands)

            # CUSTOM DATA:
            # pack LRMS array
            for iBand in range(codingParams.sfBands.nBands):
                pb.WriteBits(LRMS[iBand],1)

            # finally, write the data in this channel's PackedBits object to the output file
            self.fp.write(pb.GetPackedData())
            
        # end loop over channels, done writing coded data for all channels
        return

    def Close(self,codingParams):
        """
        Flushes the last data block through the encoding process (if encoding)
        and closes the audio file
        """
        # determine if encoding or encoding and, if encoding, do last block
        if self.fp.mode == "wb":  # we are writing to the PACFile, must be encode
            # we are writing the coded file -- pass a block of zeros to move last data block to other side of MDCT block

            #fixing last block of data that's to be written
            if codingParams.win_state == 0 or codingParams.win_state == 3:
                nMDCTLines = codingParams.nMDCTLinesLong
            else:
                nMDCTLines = codingParams.nMDCTLinesShort

            data = []
            for iCh in range(codingParams.nChannels): data.append( np.zeros(nMDCTLines, dtype=np.float) )                 
            self.WriteDataBlock(data, codingParams,myhuffyman)
            
        self.fp.close()

    def Encode(self,data,codingParams,myhuffyman):
        """
        Encodes multichannel audio data and returns a tuple containing
        the scale factors, mantissa bit allocations, quantized mantissas,
        and the overall scale factor for each channel.
        """
        #Passes encoding logic to the Encode function defined in the codec module
        return codec.Encode(data,codingParams,myhuffyman)

    def Decode(self,scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams,LRMS,myhuffyman):
        """
        Decodes a single audio channel of data based on the values of its scale factors,
        bit allocations, quantized mantissas, and overall scale factor.
        """
        #Passes decoding logic to the Decode function defined in the codec module
        return codec.Decode(scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams,LRMS)





if __name__=="__main__":

    import sys
    import time
    from pcmfile import * # to get access to WAV file handling

    myhuffyman = Huffman()
    W = WindowState()
    nSubBlocks = 8
    #some variables I want to keep tab of 
    readSignal = np.array([], dtype = float)
    writeSignal = np.array([], dtype = float)
    transients = []
    count = 1
    countBlockEncode = 0
    countBlockDecode = 0

    input_filename = "../audio/Castanets.wav"
    coded_filename = "../audio/castanets_BS_HUFFMAN_STEREO_reducedBits.pac"
    output_filename = "../audio/castanets_BS_HUFFMAN_STEREO_reducedBits.wav"


    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
        coded_filename = sys.argv[1][:-4] + ".pac"
        output_filename = sys.argv[1][:-4] + "_decoded.wav"


    print "\nRunning the PAC coder ({} -> {} -> {}):".format(input_filename, coded_filename, output_filename)
    elapsed = time.time()

    for Direction in ("Encode","Decode"):
#    for Direction in ("Decode"):

        # create the audio file objects
        if Direction == "Encode":
            print "\n\tEncoding PCM file ({}) ...".format(input_filename),
            inFile= PCMFile(input_filename)
            outFile = PACFile(coded_filename)
        else: # "Decode"

            print "Saved "+str(myhuffyman.getBitReservoir())+" Bits!"
            print "\n\tDecoding PAC file ({}) ...".format(coded_filename),
            inFile = PACFile(coded_filename)
            outFile= PCMFile(output_filename)
            

        # open input file
        codingParams=inFile.OpenForReading()  # (includes reading header)

        # pass parameters to the output file
        if Direction == "Encode":
            # set additional parameters that are needed for PAC file
            # (beyond those set by the PCM file on open)
            codingParams.nMDCTLines = 512
            codingParams.nScaleBits = 4
            codingParams.nMantSizeBits = 4
            codingParams.targetBitsPerSample = 4
            # tell the PCM file how large the block size is
            codingParams.nSamplesPerBlock = codingParams.nMDCTLines
            codingParams.numTableBits = 4
            codingParams.win_state = 0
        else: # "Decode"
            # set PCM parameters (the rest is same as set by PAC file on open)
            codingParams.bitsPerSample = 16
            codingParams.numTableBits = 4
        # only difference is in setting up the output file parameters


        # open the output file
        outFile.OpenForWriting(codingParams) # (includes writing header)

        # Read the input file and pass its data to the output file to be written
        firstBlock = True  # when de-coding, we won't write the first block to the PCM file. This flag signifies that

        curBlock = [[] for x in range(codingParams.nChannels)]
        transBlock = [[] for x in range(codingParams.nChannels)]
        
        while True:

            if Direction == "Encode":
                data=inFile.ReadDataBlock(codingParams)
            else:
                data=inFile.ReadDataBlock(codingParams,myhuffyman)

            if not data: break  # we hit the end of the input file

            # don't write the first PCM block (it corresponds to the half-block delay introduced by the MDCT)
            if firstBlock and Direction == "Decode":
                firstBlock = False
                continue    

#            if Direction == "Encode" :
#                # readSignal = np.concatenate([readSignal,data[0]])
#                countBlockEncode += 1
#            elif Direction == "Decode" :
#                # writeSignal = np.concatenate([writeSignal, data[0]])
#                countBlockDecode += 1
                
            
            #while encoding
            if Direction == "Encode" :
                
                #detect transient in block
                newBlock = np.copy(data[0])
                is_transient = transient_detection(newBlock)
                #is_transient = False
                codingParams.win_state = W.nextBuffer(is_transient)
                
                if(is_transient):
                    transients.append(count)
                count += 1
            
                #print('Window state :' , codingParams.win_state)
        
                #depending on window state, send the right amount of data to be encoded            
            
                if codingParams.win_state == 0:
                    sys.stdout.write("_ ")  # just to signal how far we've gotten to user
                elif codingParams.win_state == 1:
                    sys.stdout.write("/ ")
                elif codingParams.win_state == 2:
                    sys.stdout.write("^ ")
                else:
                    sys.stdout.write("\\ ")
                
                #send to WriteDataBlock for encoding
                if(codingParams.win_state == 0 or codingParams.win_state == 3):
                    
                        outFile.WriteDataBlock(data, codingParams,myhuffyman)
                        
                elif(codingParams.win_state == 1):
                        
                        #encode start transition window
                        for iCh in range(codingParams.nChannels):
                            #first 64 samples
                            curBlock[iCh] = data[iCh][:codingParams.nMDCTLinesShort]
                        
                        outFile.WriteDataBlock(curBlock, codingParams,myhuffyman)
                        
                        #now code the short blocks
                        codingParams.win_state = 2
                        W.setState(2)
                        
                        for k in range(1,nSubBlocks):
                            for iCh in range(codingParams.nChannels):
                                transBlock[iCh] = data[iCh][np.arange(codingParams.nMDCTLinesShort) + k*codingParams.nMDCTLinesShort]
                            
                            outFile.WriteDataBlock(transBlock, codingParams,myhuffyman) 
                            
                elif(codingParams.win_state == 2):
                    
                    #print('Two consecutive transients')
                    #print "shape of Tranition Block : " + str(np.shape(transBlock))
                    
                    for k in range(nSubBlocks):
                        for iCh in range(codingParams.nChannels):
                            transBlock[iCh] = data[iCh][np.arange(codingParams.nMDCTLinesShort) + k*codingParams.nMDCTLinesShort]
                            
                        outFile.WriteDataBlock(transBlock, codingParams,myhuffyman) 
                            
                else :
                    raise ValueError('Invalid window state: ' + str(codingParams.win_state))
            
            #while decoding
            elif Direction == "Decode":
                
                if codingParams.win_state == 0:
                    sys.stdout.write("_ ")  # just to signal how far we've gotten to user
                elif codingParams.win_state == 1:
                    sys.stdout.write("/ ")
                elif codingParams.win_state == 2:
                    sys.stdout.write("^ ")
                else:
                    sys.stdout.write("\\ ")
                    
                outFile.WriteDataBlock(data, codingParams)
                
        
            #sys.stdout.write(".")  # just to signal how far we've gotten to user
         
            sys.stdout.flush()
        # end loop over reading/writing the blocks

        # close the files
        inFile.Close(codingParams)
        outFile.Close(codingParams)
    # end of loop over Encode/Decode

    elapsed = time.time()-elapsed
    print "\nDone with Encode/Decode test\n"
    print elapsed ," seconds elapsed"
    # plt.plot(codec.whichbuffer/codec.numTime)

    # plt.title('Average Bits Allocated for the whole audio ( nMDCTLines = 1024 )')
    # plt.xlabel('Frequency Bands')
    # plt.ylabel('Number of Bits Allocated')

    # plt.show()
