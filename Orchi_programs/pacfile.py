"""
pacfile.py -- Defines a PACFile class to handle reading and writing audio
data to an audio file holding data compressed using an MDCT-based perceptual audio
coding algorithm.  The MDCT lines of each audio channel are grouped into bands,
each sharing a single scaleFactor and bit allocation that are used to block-
floating point quantize those lines.  This class is a subclass of AudioFile.

-----------------------------------------------------------------------
© 2009 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------

See the documentation of the AudioFile class for general use of the AudioFile
class.

Notes on reading and decoding PAC files:

    The OpenFileForReading() function returns a CodedParams object containing:

        nChannels = the number of audio channels
        sampleRate = the sample rate of the audio samples
        numSamples = the total number of samples in the file for each channel
        nMDCTLines = half the MDCT block size (block switching not supported)
        nSamplesPerBlock = MDCTLines (but a name that PCM files look for)
        nScaleBits = the number of bits storing scale factors
        nMantSizeBits = the number of bits storing mantissa bit allocations
        sfBands = a ScaleFactorBands object
        overlapAndAdd = decoded data from the prior block (initially all zeros)

    The returned ScaleFactorBands object, sfBands, contains an allocation of
    the MDCT lines into groups that share a single scale factor and mantissa bit
    allocation.  sfBands has the following attributes available:

        nBands = the total number of scale factor bands
        nLines[iBand] = the number of MDCT lines in scale factor band iBand
        lowerLine[iBand] = the first MDCT line in scale factor band iBand
        upperLine[iBand] = the last MDCT line in scale factor band iBand


Notes on encoding and writing PAC files:

    When writing to a PACFile the CodingParams object passed to OpenForWriting()
    should have the following attributes set:

        nChannels = the number of audio channels
        sampleRate = the sample rate of the audio samples
        numSamples = the total number of samples in the file for each channel
        nMDCTLines = half the MDCT block size (format does not support block switching)
        nSamplesPerBlock = MDCTLines (but a name that PCM files look for)
        nScaleBits = the number of bits storing scale factors
        nMantSizeBits = the number of bits storing mantissa bit allocations
        targetBitsPerSample = the target encoding bit rate in units of bits per sample

    The first three attributes (nChannels, sampleRate, and numSamples) are
    typically added by the original data source (e.g. a PCMFile object) but
    numSamples may need to be extended to account for the MDCT coding delay of
    nMDCTLines and any zero-padding done in the final data block

    OpenForWriting() will add the following attributes to be used during the encoding
    process carried out in WriteDataBlock():

        sfBands = a ScaleFactorBands object
        priorBlock = the prior block of audio data (initially all zeros)

    The passed ScaleFactorBands object, sfBands, contains an allocation of
    the MDCT lines into groups that share a single scale factor and mantissa bit
    allocation.  sfBands has the following attributes available:

        nBands = the total number of scale factor bands
        nLines[iBand] = the number of MDCT lines in scale factor band iBand
        lowerLine[iBand] = the first MDCT line in scale factor band iBand
        upperLine[iBand] = the last MDCT line in scale factor band iBand

Description of the PAC File Format:

    Header:

        tag                 4 byte file tag equal to "PAC "
        sampleRate          little-endian unsigned long ("<L" format in struct)
        nChannels           little-endian unsigned short("<H" format in struct)
        numSamples          little-endian unsigned long ("<L" format in struct)
        nMDCTLines          little-endian unsigned long ("<L" format in struct)
        nScaleBits          little-endian unsigned short("<H" format in struct)
        nMantSizeBits       little-endian unsigned short("<H" format in struct)
        nSFBands            little-endian unsigned long ("<L" format in struct)
        for iBand in range(nSFBands):
            nLines[iBand]   little-endian unsigned short("<H" format in struct)

    Each Data Block:  (reads data blocks until end of file hit)

        for iCh in range(nChannels):
            nBytes          little-endian unsigned long ("<L" format in struct)
            as bits packed into an array of nBytes bytes:
                overallScale[iCh]                       nScaleBits bits
                for iBand in range(nSFBands):
                    scaleFactor[iCh][iBand]             nScaleBits bits
                    bitAlloc[iCh][iBand]                nMantSizeBits bits
                    if bitAlloc[iCh][iBand]:
                        for m in nLines[iBand]:
                            mantissa[iCh][iBand][m]     bitAlloc[iCh][iBand]+1 bits
                <extra custom data bits as long as space is included in nBytes>

"""
from __future__ import division
from audiofile import * # base class
from bitpack import *  # class for packing data into an array of bytes where each item's number of bits is specified
import codec    # module where the actual PAC coding functions reside(this module only specifies the PAC file format)
from psychoac import ScaleFactorBands, AssignMDCTLinesFromFreqLimits  # defines the grouping of MDCT lines into scale factor bands
import sys
from blockswitch import *
import window as win
import mdct as m
import quantize as qt

import numpy as np  # to allow conversion of data blocks to numpy's array object
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
        # start w/o all zeroes as data from prior block to overlap-and-add for output
        overlapAndAdd = []
        for iCh in range(nChannels): overlapAndAdd.append( np.zeros(nMDCTLines, dtype=np.float64) )
        myParams.overlapAndAdd=overlapAndAdd
        return myParams


    def ReadDataBlock(self, codingParams):
        """
        Reads a block of coded data from a PACFile object that has already
        executed OpenForReading() and returns those samples as reconstituted
        signed-fraction data
        
        How can I print from within this function? - it's not returning the right size of data
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
                  
            # CUSTOM DATA:
            # < now can unpack any custom data passed in the nBytes of data >
            #unpack 2 bits of data for determining window state
            codingParams.win_state = pb.ReadBits(2)
            
            #determine block size and overlap add facttors given window state
            if codingParams.win_state == 0:
                codingParams.nMDCTLines = 512
                olaSamplesCur = olaSamplesStore = 512
            
            elif codingParams.win_state == 1:
                codingParams.nMDCTLines = (512 + 64) / 2
                olaSamplesCur = 512
                olaSamplesStore = 512
            
            elif codingParams.win_state == 2:
                codingParams.nMDCTLines = 64
                olaSamplesCur = 64
                olaSamplesStore = 64
            
            elif codingParams.win_state == 3:
                codingParams.nMDCTLines = (512 + 64) / 2
                olaSamplesCur = 64
                olaSamplesStore = 64
            
            else:
                raise ValueError('Invalid window state: ' + str(codingParams.win_state))
                
            #make sure you have the right sfBands values
            codingParams.sfBands.nLines = AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLines, codingParams.sampleRate)
            codingParams.sfBands.upperLine = np.cumsum(codingParams.sfBands.nLines)-1
            codingParams.sfBands.lowerLine = codingParams.sfBands.upperLine - codingParams.sfBands.nLines + 1
            
            
            # extract the data from the PackedBits object
            overallScaleFactor = pb.ReadBits(codingParams.nScaleBits)  # overall scale factor
            scaleFactor=[]
            bitAlloc=[]
            mantissa=np.zeros(int(codingParams.nMDCTLines),np.int32)  # start w/ all mantissas zero
            for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                ba = pb.ReadBits(codingParams.nMantSizeBits)
                if ba: ba+=1  # no bit allocation of 1 so ba of 2 and up stored as one less
                bitAlloc.append(ba)  # bit allocation for this band
                scaleFactor.append(pb.ReadBits(codingParams.nScaleBits))  # scale factor for this band
                if bitAlloc[iBand]:
                    # if bits allocated, extract those mantissas and put in correct location in matnissa array
                    m=np.empty(codingParams.sfBands.nLines[iBand],np.int32)
                    for j in range(codingParams.sfBands.nLines[iBand]):
                        m[j]=pb.ReadBits(bitAlloc[iBand])     # mantissas for this band (if bit allocation non-zero) and bit alloc <>1 so encoded as 1 lower than actual allocation
                    mantissa[codingParams.sfBands.lowerLine[iBand]:(codingParams.sfBands.upperLine[iBand]+1)] = m
            # done unpacking data (end loop over scale factor bands)


            # (DECODE HERE) decode the unpacked data for this channel, overlap-and-add first half, and append it to the data array (saving other half for next overlap-and-add)
            decodedData = self.Decode(scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams)

            

            
            #for debugging
#            if(codingParams.win_state != 0):
#                print("For debugging :", codingParams.win_state, np.size(decodedData), np.size(codingParams.overlapAndAdd[iCh]))            
            
            #take care of overlap add with different block sizes
            data[iCh] = np.concatenate( (data[iCh],np.add(codingParams.overlapAndAdd[iCh],decodedData[:olaSamplesCur]) ) )  # data[iCh] is overlap-and-added data
            codingParams.overlapAndAdd[iCh] = decodedData[olaSamplesStore:]  # save other half for next pass


            
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
        sfBands=ScaleFactorBands( AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLines,
                                                                codingParams.sampleRate))
        codingParams.sfBands=sfBands
        self.fp.write(pack('<L',sfBands.nBands))
        self.fp.write(pack('<'+str(sfBands.nBands)+'H',*(sfBands.nLines.tolist()) ))
        # start w/o all zeroes as prior block of unencoded data for other half of MDCT block
        priorBlock = []
        for iCh in range(codingParams.nChannels):
            priorBlock.append(np.zeros(codingParams.nMDCTLines,dtype=np.float64) )
        codingParams.priorBlock = priorBlock
        
        #create new members of codingParams here
        codingParams.win_state = 0
        
        return


    def WriteDataBlock(self,data, codingParams):
        """
        Writes a block of signed-fraction data to a PACFile object that has
        already executed OpenForWriting()"""

        # combine this block of multi-channel data w/ the prior block's to prepare for MDCTs twice as long
        fullBlockData=[]
        
        # if status for data[iCh] is start or short
        
        #find out block length and number of subblocks for each window state
        if codingParams.win_state == 0:
            numSubBlocks = nDivide = 1
            nLines = 512
        elif codingParams.win_state == 1:
            numSubBlocks = 1
            nDivide = 1
            nLines = (512 + 64)/2
        elif codingParams.win_state == 2
        :
            numSubBlocks = nDivide = 8
            nLines = 64
        elif codingParams.win_state == 3:
            numSubBlocks = nDivide = 1
            nLines = (512 + 64)/2
        else: 
            raise ValueError('Invalid window state: ' + str(codingParams.win_state))
            
        #make sure you have the right sfBands values depeding on window being used         
        codingParams.sfBands.nLines = AssignMDCTLinesFromFreqLimits(nLines, codingParams.sampleRate)
        codingParams.sfBands.upperLine = np.cumsum(codingParams.sfBands.nLines)-1
        codingParams.sfBands.lowerLine = codingParams.sfBands.upperLine - codingParams.sfBands.nLines + 1
       
        # store a copy of data in a temp array
        tempData = np.copy(data)
        # loop through the subblocks and process them (assuming same number of subblocks for all channels)
        for iBlk in range(numSubBlocks):
            fullBlockData=[]
            for iCh in range(codingParams.nChannels):
                b = int(codingParams.nMDCTLines/nDivide)
                fullBlockData.append( 
                        np.concatenate( 
                            ( codingParams.priorBlock[iCh], 
                                tempData[iCh][np.arange(b) + iBlk*b] ) 
                                      ) 
                                    )
                # current pass's data is next pass's prior block data
                codingParams.priorBlock[iCh] = tempData[iCh][np.arange(b) + iBlk*b]
                
            # (ENCODE HERE) Encode the full block of multi=channel data

            # print " LeData Being Encoded : " + str(fullBlockData)
            (scaleFactor,bitAlloc,mantissa, overallScaleFactor) = self.Encode(fullBlockData,codingParams)  # returns a tuple with all the block-specific info not in the file header


            # for each channel, write the data to the output file
            for iCh in range(codingParams.nChannels):
    
                # determine the size of this channel's data block and write it to the output file
                nBytes = codingParams.nScaleBits  # bits for overall scale factor
                for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to get its bits
                    nBytes += codingParams.nMantSizeBits+codingParams.nScaleBits    # mantissa bit allocation and scale factor for that sf band
                    if bitAlloc[iCh][iBand]:
                        # if non-zero bit allocation for this band, add in bits for scale factor and each mantissa (0 bits means zero)
                        nBytes += bitAlloc[iCh][iBand]*codingParams.sfBands.nLines[iBand]  # no bit alloc = 1 so actuall alloc is one higher
                # end computing bits needed for this channel's data
    
                # CUSTOM DATA:
                # < now can add space for custom data, if desired>
                #add space for window_state
                nBytes += 2
    
                # now convert the bits to bytes (w/ extra one if spillover beyond byte boundary)
                if nBytes%BYTESIZE==0:  nBytes /= BYTESIZE
                else: nBytes = nBytes/BYTESIZE + 1
                self.fp.write(pack("<L",int(nBytes))) # stores size as a little-endian unsigned long
    
                # create a PackedBits object to hold the nBytes of data for this channel/block of coded data
                pb = PackedBits()
                pb.Size(nBytes)
    
                # CUSTOM DATA:
                # < now can add in custom data if space allocated in nBytes above>
                #to determine window state
                pb.WriteBits(codingParams.win_state, 2)
                
                # now pack the nBytes of data into the PackedBits object
                pb.WriteBits(overallScaleFactor[iCh],codingParams.nScaleBits)  # overall scale factor
                iMant=0  # index offset in mantissa array (because mantissas w/ zero bits are omitted)
                for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                    ba = bitAlloc[iCh][iBand]
                    if ba: ba-=1  # if non-zero, store as one less (since no bit allocation of 1 bits/mantissa)
                    pb.WriteBits(ba,codingParams.nMantSizeBits)  # bit allocation for this band (written as one less if non-zero)
                    pb.WriteBits(scaleFactor[iCh][iBand],codingParams.nScaleBits)  # scale factor for this band (if bit allocation non-zero)
                    if bitAlloc[iCh][iBand]:
                        for j in range(codingParams.sfBands.nLines[iBand]):
                            pb.WriteBits(mantissa[iCh][iMant+j],bitAlloc[iCh][iBand])     # mantissas for this band (if bit allocation non-zero) and bit alloc <>1 so is 1 higher than the number
                        iMant += codingParams.sfBands.nLines[iBand]  # add to mantissa offset if we passed mantissas for this band
                # done packing (end loop over scale factor bands)
    
                # finally, write the data in this channel's PackedBits object to the output file
                self.fp.write(pb.GetPackedData())
                # end loop over channels, done writing coded data for all channels
        
        #end loop over subblocks
        return

    def Close(self,codingParams):
        """
        Flushes the last data block through the encoding process (if encoding)
        and closes the audio file
        """
        # determine if encoding or encoding and, if encoding, do last block
        if self.fp.mode == "wb":  # we are writing to the PACFile, must be encode
            # we are writing the coded file -- pass a block of zeros to move last data block to other side of MDCT block
            
            data = [ np.zeros(codingParams.nMDCTLines,dtype=np.float),
                     np.zeros(codingParams.nMDCTLines,dtype=np.float) ]
            self.WriteDataBlock(data, codingParams)
        self.fp.close()


    def Encode(self,data,codingParams):
        """
        Encodes multichannel audio data and returns a tuple containing
        the scale factors, mantissa bit allocations, quantized mantissas,
        and the overall scale factor for each channel.
        """
        #Passes encoding logic to the Encode function defined in the codec module
        return codec.Encode(data,codingParams)

    def Decode(self,scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams):
        """
        Decodes a single audio channel of data based on the values of its scale factors,
        bit allocations, quantized mantissas, and overall scale factor.
        """
        #Passes decoding logic to the Decode function defined in the codec module
        return codec.Decode(scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams)








#-----------------------------------------------------------------------------

# Testing the full PAC coder (needs a file called "input.wav" in the code directory)
if __name__=="__main__":

    import sys
    import time
    from pcmfile import * # to get access to WAV file handling

    input_filename = "Castanets.wav"
    coded_filename = "coded.pac"
    output_filename = "Castanets_output.wav"

    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
        coded_filename = sys.argv[1][:-4] + ".pac"
        output_filename = sys.argv[1][:-4] + "_decoded.wav"


    print "\nRunning the PAC coder ({} -> {} -> {}):".format(input_filename, coded_filename, output_filename)
    elapsed = time.time()

    for Direction in ("Encode", "Decode"):
  #  for Direction in ("Encode"):

        # create the audio file objects
        if Direction == "Encode":
            print "\n\tEncoding PCM file ({}) ...".format(input_filename),
            inFile= PCMFile(input_filename)
            outFile = PACFile(coded_filename)
        else: # "Decode"
            print "\n\tDecoding PAC file ({}) ...".format(coded_filename),
            inFile = PACFile(coded_filename)
            outFile= PCMFile(output_filename)
        # only difference is file names and type of AudioFile object

        # open input file
        codingParams=inFile.OpenForReading()  # (includes reading header)

        # pass parameters to the output file
        if Direction == "Encode":
            # set additional parameters that are needed for PAC file
            # (beyond those set by the PCM file on open)
            codingParams.nMDCTLines = 512
            codingParams.nScaleBits = 4
            codingParams.nMantSizeBits = 4
            #codingParams.targetBitsPerSample = 5.34/2
            #set target bits per sample for given datarate 
            dataRate = 705600
            nBitsPerBlock = dataRate/codingParams.sampleRate * codingParams.nMDCTLines
            codingParams.targetBitsPerSample = (nBitsPerBlock - 2*(4*25) - 4)/(codingParams.nMDCTLines)
            # tell the PCM file how large the block size is
            codingParams.nSamplesPerBlock = codingParams.nMDCTLines
        else: # "Decode"
            # set PCM parameters (the rest is same as set by PAC file on open)
            codingParams.bitsPerSample = 16
        # only difference is in setting up the output file parameters


        # open the output file
        outFile.OpenForWriting(codingParams) # (includes writing header)

        # Read the input file and pass its data to the output file to be written
        firstBlock = True  # when de-coding, we won't write the first block to the PCM file. This flag signifies that
        W = WindowState()
        #first block won't have transients since we are zero padding anyway
        codingParams.win_state = W.state
        curBlock = inFile.ReadDataBlock(codingParams)
        
#        #need a separate sfbands object to do this
#        transient_sfbands = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLines/2, codingParams.sampleRate))
#        percEntropy = [[] for i in range(codingParams.nChannels)]
#        thres_change = 200
#        thres_mag = 1000
#        count = 0
        
        #start reading/writing data
        while True:
            
            count = 0
            nextBlock = inFile.ReadDataBlock(codingParams)  
            #reached end of file            
            if not curBlock: break  # we hit the end of the input file
        
            
            # don't write the first PCM block (it corresponds to the half-block delay introduced by the MDCT)
            if firstBlock and Direction == "Decode":
                
                firstBlock = False
                #do transient detection while encoding only

                curBlock[0] = np.copy(nextBlock[0])
                curBlock[1] = np.copy(nextBlock[1])

                continue
            
            #do block switching

            if nextBlock and Direction == "Encode":
                
                #find transient for each channel of data
                 for iCh in range(codingParams.nChannels):
                     
                     #newBlock = np.concatenate([curBlock[iCh], nextBlock[iCh]])
                     newBlock = nextBlock[iCh]
                     
#                     transient detection with perceptual entropy
#                     #save current block to pass on to detectTransient function
#                     dataBlock = newBlock
#                     #window 2*nSamplesPerBlock set of samples
#                     newBlock = win.KBDWindow(newBlock)
#                     #MDCT
#                     MDCTdata = m.MDCT(newBlock,codingParams.nSamplesPerBlock/2,codingParams.nSamplesPerBlock/2)
#                     #find out MDCT scale
#                     maxLine = np.max(np.abs(newBlock))
#                     MDCTscale = qt.ScaleFactor(maxLine,codingParams.nScaleBits)
#        
#                     #find perceptual entropy for each channel
#                     percEntropy[iCh].append(detectTransient(transient_sfbands, dataBlock, MDCTdata, MDCTscale, codingParams.sampleRate))
#    
#                     #detect transients by detecting change in PE from block to block    
#                     if((count > 0 and percEntropy[iCh][count] - percEntropy[iCh][count-1] >= thres_change)
#                        or percEntropy[iCh][count] > thres_mag):
#                        is_transient = True
#                     else:
#                        is_transient = False
                     
                 #transient detection using Prateek's FFT Based method
                 is_transient = transient_detection(newBlock)
                 
                 count += 1
                
                 #determine window state
                 W.state = W.nextBuffer(is_transient)
                 codingParams.win_state = W.state
            
            #write the current block with the correct window_state
            outFile.WriteDataBlock(curBlock,codingParams)
            curBlock = nextBlock
            
            
            if codingParams.win_state == 0:
                sys.stdout.write("_ ")  # just to signal how far we've gotten to user
            elif codingParams.win_state == 1:
                sys.stdout.write("/ ")
            elif codingParams.win_state == 2:
                sys.stdout.write("^ ")
            else:
                sys.stdout.write("\\ ")
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
