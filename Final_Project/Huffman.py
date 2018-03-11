# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from collections import deque
import cPickle as pickle


class Histogram:
    
    # Stores the occurence probability of each symbol
    probability = dict()
    # Stores the occurence frequency of each symbol
    statistics = dict()
    # Stores the Huffman Nodes to be pop out (prioritized based on freq)
    queue = deque()
    
    #---PUBLIC METHODS---#
    def __init__(self):
        self.LOW_FREQ = 10
        self.ESCAPE_CODE = -1
        pass

    
    def getMatchScore(self,blockHistogram):
        # use probability to compare and give score
        histogramProb = self.probability
        # Sort the HistogramProb based on the mantissa code value for alignment
        sortedHistogramProb = OrderedDict(sorted(histogramProb.items(), key=lambda t: t[0]))
        sortedHistogramProb = np.array(sortedHistogramProb.values())
         # Sort the block HistogramProb based on the mantissa code value for alignment
        sortedBlockHistogramProb = OrderedDict(sorted(blockHistogram.probability.items(),key=lambda t: t[0]))
        sortedBlockHistogramProb = np.array(sortedBlockHistogramProb.values())
        #print "Block prob: ",sortedBlockHistogramProb
        difference = sum((sortedHistogramProb - sortedBlockHistogramProb)**2)
        return (3.0-difference)
        
    
    def generateStatistics(self,mantissaCode):
        # Count the occurence freq for each manitssaCode
        for i in range(len(mantissaCode)):
            if mantissaCode[i] in self.statistics:
                self.statistics[mantissaCode[i]] = self.statistics[mantissaCode[i]] + 1
            else:
                self.statistics[mantissaCode[i]] = 1
        totalCount = sum(self.statistics.values())
        # Transform statistics to probabilities for finding histogram matches
        for key, value in self.statistics.items():
            self.probability[key] = value/float(totalCount)

    
    def makeHuffmanNodeQueue(self):
        # sort the statistics from low freq to high freq
        sortedStatistics = OrderedDict(sorted(self.statistics.items(), key=lambda t: t[1]))
        sortedStatistics = sortedStatistics.items()
        escapeFreq = 0
        for pair in sortedStatistics:
            # pair[0] -> mantissa code
            # pair[1] -> frequency
            # Skip very low frequency symbols and replace with escape code
            if pair[1] < self.LOW_FREQ:
                escapeFreq += 1
            else:
                newNode = HuffmanNode(pair[0],pair[1])
                self.queue.append(newNode)
        # Add the node for escape code
        newNode = HuffmanNode(self.ESCAPE_CODE ,escapeFreq)
        self.queue.append(newNode)
        self.queue = deque(sorted(self.queue, key=lambda t: t.freq))


    def appendToHuffmanQueue(self,huffmanNode):
        self.queue.append(huffmanNode)
        self.queue = deque(sorted(self.queue, key=lambda t: t.freq))

    
    def getNextPair(self):
        if len(self.queue)==1:
            return (self.queue.popleft(),None)
        firstNode = self.queue.popleft()
        secondNode = self.queue.popleft()
        return (firstNode,secondNode)

#------------------------------------------------------------------------#
class HuffmanTable:
    
   
    def __init__(self,manitssaCodeToHuffmanCode):
        self.encodingTable = manitssaCodeToHuffmanCode
        self.decodingTable = dict()
        for key, value in manitssaCodeToHuffmanCode.items():
            self.decodingTable[value] = key

#------------------------------------------------------------------------#

class Huffman:
    
    #---PUBLIC METHODS---#
    def __init__(self):
        with open('huffmanTables/huffmanTables.pickle', 'rb') as handle:
            self.huffmanTables = pickle.load(handle)
        with open('huffmanTables/histograms.pickle', 'rb') as handle:
            self.histograms = pickle.load(handle)
        self.ESCAPE_CODE = -1
        self.bitsInReservoir = 0
    
    
    def encodeHuffman(self,codingParams,mantissaCode,bitAlloc):
        bestHuffmanCodedMantissa = []

        bestTableID = 1
        
        for ID in self.huffmanTables.keys():
            huffmanCodedMantissa = []
            huffmanTable = self.huffmanTables[ID]
            iMant=0
            for iBand in range(codingParams.sfBands.nBands):
                if bitAlloc[iBand]:
                    for j in range(codingParams.sfBands.nLines[iBand]):
                        code = mantissaCode[iMant+j]
                        if code in huffmanTable.encodingTable.keys():
                            huffmanCodedMantissa.append(huffmanTable.encodingTable[code])
                        # If mantissa code not presented in table, use escape code and code the origin mantissa
                        else:
                            form = '0' + str(bitAlloc[iBand]) + 'b'
                            mantString = format(code,form)
                            huffmanCodedMantissa.append(huffmanTable.encodingTable[self.ESCAPE_CODE]+mantString)
                    iMant += codingParams.sfBands.nLines[iBand]
            totalHuffmanBitLength = sum(len(huff) for huff in huffmanCodedMantissa)
            if  ID == 1:
                totalBitLength = totalHuffmanBitLength
                bestHuffmanCodedMantissa = huffmanCodedMantissa
                bestTableID = ID
            if  totalHuffmanBitLength < totalBitLength:
                totalBitLength = totalHuffmanBitLength
                bestHuffmanCodedMantissa = huffmanCodedMantissa
                bestTableID = ID
        return (bestHuffmanCodedMantissa,bestTableID)
    
   
    def decodeHuffman(self,bitReader,tableID,bitAlloc):
        huffmanTable = self.huffmanTables[tableID]
        huffmanCode = ""
        huffmanCodes = huffmanTable.decodingTable.keys()
        mantissa = self.decodeHuffmanHelper(bitReader,huffmanTable,huffmanCodes,huffmanCode)
        if mantissa == self.ESCAPE_CODE:
            return bitReader.ReadBits(bitAlloc)
        return mantissa
    
    
    def decodeHuffmanHelper(self,bitReader,huffmanTable,huffmanCodes,huffmanCode):
        while True:
            if huffmanCode in huffmanCodes:
                return huffmanTable.decodingTable[huffmanCode]
            if(bitReader.ReadBits(1)==0):
                huffmanCode += repr(0)
            else:
                huffmanCode += repr(1)

    
    def giveBits(self,numBits):
        self.bitsInReservoir += numBits

    
    def snatchBits(self):
        bitReserve = 0

        bitsInReservoirThreshold = 10
        bitDepsoitTransferPercentage = 1.0/5
        if self.bitsInReservoir > bitsInReservoirThreshold:
            bitReserve = self.bitsInReservoir*bitDepsoitTransferPercentage
            self.bitsInReservoir -= bitReserve
        elif self.bitsInReservoir < 0:
            bitReserve = self.bitsInReservoir
            self.bitsInReservoir = 0
        return bitReserve

    def getBitReservoir(self):
        return self.bitsInReservoir

 
"""---TESTING!----"""
if __name__ == "__main__":



    pass
