from Huffman import*

class HuffmanNode:
    freq = None
    mantissaCode = None
    zero = None #left Node
    one = None  #right Node
    def __init__(self,mantissaCode=None,freq=None,left=None,right=None):
        self.zero = left
        self.one = right
        self.freq = freq
        self.mantissaCode = mantissaCode

#------------------------------------------------------------------------#

class HuffmanTrainer:
    
    histogram = Histogram()
    huffmanCodeTable = dict()
    root = HuffmanNode()
    
    #---PUBLIC METHODS---#
    #----------------------------------------------------------#
    #. HuffmanTrainer constructor
    #. This function takes in the ID of the Huffman table to be
    #. constructed and then init the HuffmanTrainer object.
    #. @Param:  table ID
    #. @Return: void
    #----------------------------------------------------------#
    def __init__(self,tableID):
        self.tableID = tableID
        self.histogram = Histogram()

    #----------------------------------------------------------#
    #. countFreq(mantissaCode)
    #. This function takes in a list of mantissa code (*without
    #. sign bit*) and count the frequency of the occurence of
    #. paricular mantissa code.
    #. @Param:  a list of mantissa code (without sign bit)
    #. @Return: void
    #----------------------------------------------------------#
    def countFreq(self,mantissaCode):
        self.histogram.generateStatistics(mantissaCode)
    
    #----------------------------------------------------------#
    #. constructHuffmanTable()
    #. This function construct a Huffman code table based on the
    #. statistics it collects so far and store the table with
    #. corresponding ID in order for future decode lookup.
    #. @Param:  void
    #. @Return: void
    #----------------------------------------------------------#
    def constructHuffmanTable(self):
        self.histogram.makeHuffmanNodeQueue()
        self.__buildEncodingTree()
        self.__buildEncodingTable()
        # Save the huffman code table to file to read in for decode
        with open('huffmanTables.pickle', 'rb') as handle:
            huffmanTables = pickle.load(handle)
            huffmanTables[self.tableID] = HuffmanTable(self.huffmanCodeTable)
            with open('huffmanTables.pickle', 'wb') as handle:
                pickle.dump(huffmanTables, handle)
        # Save the corresponding histogram to file to read in for decode
        with open('histograms.pickle', 'rb') as handle:
            histograms = pickle.load(handle)
            histograms[self.tableID] = self.histogram
            with open('histograms.pickle', 'wb') as handle:
                pickle.dump(histograms, handle)

    #---PRIVATE METHODS---#
    #----------------------------------------------------------#
    #. __buildEncodingTree()
    #. Private method to build the huffman encoding tree and sets
    #. the root huffman node for building encoding table.
    #. @Param:  void
    #. @Return: void
    #----------------------------------------------------------#
    def __buildEncodingTree(self):
        while True:
            (firstNode,secondNode) = self.histogram.getNextPair()
            if secondNode == None:
                self.root = firstNode
                break
            joinedNode = HuffmanNode(None,firstNode.freq+secondNode.freq,firstNode,secondNode)
            self.histogram.appendToHuffmanQueue(joinedNode)

    #----------------------------------------------------------#
    #. __buildEncodingTableHelper()
    #. Private recursive helper method to build the huffman
    #. encoding table.
    #. @Param:  void
    #. @Return: void
    #----------------------------------------------------------#
    def __buildEncodingTableHelper(self,node,huffmanCode):
        if node.mantissaCode != None:
            # Need to transform huffman code from string to int
            self.huffmanCodeTable[node.mantissaCode] = huffmanCode
            return
        self.__buildEncodingTableHelper(node.zero,huffmanCode+"0")
        self.__buildEncodingTableHelper(node.one,huffmanCode+"1")

    #----------------------------------------------------------#
    #. __buildEncodingTable()
    #. Private method to build the huffman encoding table.
    #. @Param:  void
    #. @Return: void
    #----------------------------------------------------------#
    def __buildEncodingTable(self):
        huffmanCode = ""
        self.__buildEncodingTableHelper(self.root,huffmanCode)

#------------------------------------------------------------------------#