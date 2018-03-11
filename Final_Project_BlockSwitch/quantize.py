"""
quantize.py -- routines to quantize and dequantize floating point values
between -1.0 and 1.0 ("signed fractions")
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
# import quantize_

### Problem 1.a.i ###
def QuantizeUniform(aNum,nBits):
    """
    Uniformly quantize signed fraction aNum with nBits
    """
    #Notes:
    #The overload level of the quantizer should be 1.0

    aQuantizedNum = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE


    ### YOUR CODE STARTS HERE ###

    if nBits > 0 : 

        maxVal = (1 << nBits - 1) - 1 

        codeGain = (1 << nBits) - 1.0 

        if aNum < 0 :

            numSign = 1 

        else : 

            numSign = 0 

        if abs(aNum) >=1:

            aQuantizedNum = maxVal

        else:

            aQuantizedNum = int((codeGain*abs(aNum) + 1.0)/ 2.0)

        if numSign:

            aQuantizedNum += (1 << nBits - 1)

    ### YOUR CODE ENDS HERE ###

    return aQuantizedNum

### Problem 1.a.i ###
def DequantizeUniform(aQuantizedNum,nBits):
    """
    Uniformly dequantizes nBits-long number aQuantizedNum into a signed fraction
    """


    aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###

    if nBits > 0 :

        codeGain = (1 << nBits) - 1.0

        numSign = 0

        if aQuantizedNum & (1 << nBits - 1) == (1 << nBits - 1): 

            aQuantizedNum -= (1 << nBits - 1)

            numSign = 1


        aNum = 2.0*aQuantizedNum / codeGain

        if numSign:

            aNum *= -numSign


    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.a.ii ###
def vQuantizeUniform(aNumVec, nBits):
    """
    Uniformly quantize vector aNumberVec of signed fractions with nBits
    """

    aQuantizedNumVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE
    aVec = aNumVec.copy()

    #Notes:
    #Make sure to vectorize properly your function as specified in the homework instructions

    ### YOUR CODE STARTS HERE ###

    if nBits > 0 :

        maxVal = (1 << nBits - 1) - 1 

        codeGain = (1 << nBits) - 1.0

        aQuantizedNumVec[abs(aVec) >= 1] = maxVal

        aQuantizedNumVec[abs(aVec) < 1] = ((abs(aVec)[abs(aVec) < 1.0]*codeGain + 1.0)/2.0).astype(int)

        signArray = np.zeros(len(aVec),dtype = bool)

        signArray[np.sign(aVec) == -1] = True 

        aQuantizedNumVec[signArray] += (1 << nBits - 1)


    ### YOUR CODE ENDS HERE ###

    return aQuantizedNumVec

### Problem 1.a.ii ###
def vDequantizeUniform(aQuantizedNumVec, nBits):
    """
    Uniformly dequantizes vector of nBits-long numbers aQuantizedNumVec into vector of  signed fractions
    """

    aNumVec = np.zeros_like(aQuantizedNumVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    aQuantizedVec = aQuantizedNumVec.copy()

    ### YOUR CODE STARTS HERE ###

    if nBits > 0 : 

        codeGain = (1 << nBits) - 1.0

        signArray = np.zeros(len(aQuantizedVec), dtype = bool)

        signArray[(aQuantizedVec & (1 << nBits - 1)) == (1 << nBits - 1)] = True 

        aQuantizedVec[signArray] -=  (1 << nBits - 1)

        aNumVec = 2.0*aQuantizedVec/codeGain

        aNumVec[signArray] = - aNumVec[signArray]

    ### YOUR CODE ENDS HERE ###

    return aNumVec

### Problem 1.b ###
def ScaleFactor(aNum, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point scale factor for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """
    #Notes:
    #The scale factor should be the number of leading zeros

    scale = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE


    ### YOUR CODE STARTS HERE ###

    if nScaleBits > 0 and nMantBits > 0: 

        nBits = (1 << nScaleBits) - 1 + nMantBits

        quantizedCode = QuantizeUniform(abs(aNum),nBits)

        quantizedCode <<= 1  # Remove the sign bit from the code 

        while scale < ( (2**nScaleBits) - 1) and ((quantizedCode & (1 << (nBits - 1))) == 0) :

            scale += 1
            quantizedCode <<= 1


    ### YOUR CODE ENDS HERE ###

    return scale

### Problem 1.b ###
def MantissaFP(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    mantissa = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###

    if nScaleBits > 0 and nMantBits > 0 :

        nBits = (1 << nScaleBits) - 1 + nMantBits

        quantizedCode = QuantizeUniform(abs(aNum),nBits)

        signMantissa = 1 << nMantBits - 1

        scaleBit = 1 << nBits - 1
        
        quantizedCode <<= (scale + 1)

        if scale < (2**nScaleBits - 1):

            
            quantizedCode -= scaleBit
            quantizedCode <<= 1
            
        quantizedCode >>= (nBits + 1 - nMantBits)

        if(aNum < 0):

            quantizedCode += signMantissa

        mantissa = quantizedCode


    ### YOUR CODE ENDS HERE ###

    return mantissa


### Problem 1.b ###
def DequantizeFP(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for floating-point scale and mantissa given specified scale and mantissa bits
    """

    aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###

    if nScaleBits > 0 and nMantBits > 0 :

        nBits = (1 << nScaleBits) - 1 + nMantBits

        signMantissa = 1 << nMantBits - 1

        signNumber = 1 << nBits - 1

        codeSign = False

        if signMantissa & mantissa:

            codeSign = True 

            mantissa -= signMantissa


        if scale != (2**nScaleBits - 1):

            mantissa += 1 << (nMantBits - 1)


        if scale < (2**nScaleBits - 1) - 1 : 

            mantissa = (mantissa << 1) + 1

            mantissa <<= 2**nScaleBits - 3 - scale


        if codeSign:

            code = mantissa + signNumber

        else : 

            code = mantissa


    aNum = DequantizeUniform(code, nBits)

    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.c.i ###
def Mantissa(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the block floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    mantissa = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###

    if nScaleBits > 0 and nMantBits > 0 :

        nBits = (1 << nScaleBits) - 1 + nMantBits
        signMantissa = 1 << nMantBits - 1 
        signNumber = 1 << nBits - 1

        code = QuantizeUniform(abs(aNum),nBits)
        code <<=scale + 1
        code >>=(nBits + 1 - nMantBits)

        if aNum < 0 :

            code += signMantissa

        mantissa = code 

    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.c.i ###
def Dequantize(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for block floating-point scale and mantissa given specified scale and mantissa bits
    """

    aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###

    if nScaleBits > 0 and nMantBits > 0 :

        nBits = (1 << nScaleBits) - 1 + nMantBits
        signMantissa = 1 << nMantBits - 1 
        signNumber = 1 << nBits - 1

        codeSign = False

        if mantissa & signMantissa :

            codeSign = True

            mantissa -= signMantissa


        mantissa = mantissa << (2**nScaleBits - 1 - scale)

        if scale < (2**nScaleBits - 1) and mantissa != 0:

            code = mantissa + (1 << (2**nScaleBits - scale - 2))

        if codeSign:

            code += signNumber



    aNum = DequantizeUniform(code, nBits)



    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.c.ii ###
def vMantissa(aNumVec, scale, nScaleBits=3, nMantBits=5):
    """
    Return a vector of block floating-point mantissas for a vector of  signed fractions aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    mantissaVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###

    if nScaleBits > 0 and nMantBits > 0:

        nBits = (1 << nScaleBits) - 1 + nMantBits
        signMantissa = 1 << nMantBits - 1 
        signNumber = 1 << nBits - 1

        signArray = np.zeros(len(aNumVec), dtype = bool)
        signArray[np.sign(aNumVec) == -1] = True

        code = vQuantizeUniform(abs(aNumVec),nBits)


        code <<= (scale + 1)
        code >>= nBits + 1 - nMantBits

        code[signArray] += signMantissa

        mantissaVec = code




    ### YOUR CODE ENDS HERE ###

    return mantissaVec

### Problem 1.c.ii ###
def vDequantize(scale, mantissaVec, nScaleBits=3, nMantBits=5):
    """
    Returns a vector of  signed fractions for block floating-point scale and vector of block floating-point mantissas given specified scale and mantissa bits
    """

    aNumVec = np.zeros_like(mantissaVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    mantissa = mantissaVec.copy()

    ### YOUR CODE STARTS HERE ###

    if nScaleBits > 0 and nMantBits > 0 :


        nBits = (1 << nScaleBits) - 1 + nMantBits
        signMantissa = 1 << nMantBits - 1 
        signNumber = 1 << nBits - 1

        negativeIndices = (mantissa & signMantissa == signMantissa)

        mantissa[negativeIndices] -= signMantissa

        code = mantissa << (2**nScaleBits - 1) - scale

        if scale < (2**nScaleBits - 1):

            code[mantissa > 0] += 1 << (2**nScaleBits - 2 - scale)


        code[negativeIndices] += signNumber


        aNumVec = vDequantizeUniform(code,nBits)

    ### YOUR CODE ENDS HERE ###

    return aNumVec

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###
    pass
    # pass # THIS DOES NOTHING

    # number = np.array([0.4,-0.3,-0.9])

    # number_1 =  0.6

    # code_quant = QuantizeUniform(number_1,8)
    # code_quant_1 = quantize_.QuantizeUniform(number_1,8)

    # dequant = DequantizeUniform(code_quant,8)
    # dequant_1 = quantize_.DequantizeUniform(code_quant_1,8)

    # print "Uniform Dequantized numbers from my solution and actual solution"

    # print dequant
    # print dequant_1

    # codes = vQuantizeUniform(number,8)
    # codes_1 = quantize_.vQuantizeUniform(number,8)

    # num = vDequantizeUniform(codes,8)
    # num_1 = quantize_.vDequantizeUniform(codes_1,8)

    # print "Dequantized vectors from my solution and actual solution"
    # print num 
    # print num_1

    # # print()

    # print "Scale Factors from my solution and actual solutions"

    # scale = ScaleFactor(number_1,3,5)
    # scale_1 = quantize_.ScaleFactor(number_1,3,5)

    # print scale 
    # print scale_1

    
    # mantissa_1 = vMantissa(number,scale,3,5)
    # mantissa_2 = quantize_.vMantissa(number,scale_1,3,5)

    # print "BFP Mantissa from my solution and actual solutions"
    # print mantissa_1 
    # print mantissa_2

    
    # dequantized = vDequantize(scale,mantissa_1,3,5)
    # dequantized_1 = quantize_.vDequantize(scale_1,mantissa_2,3,5)

    # print "BFP Dequantized from my solution and actual solutions"
    # print dequantized 
    # print dequantized_1

 


    ### YOUR TESTING CODE ENDS HERE ###

