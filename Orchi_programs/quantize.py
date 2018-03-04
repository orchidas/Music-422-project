"""
quantize.py -- routines to quantize and dequantize floating point values
between -1.0 and 1.0 ("signed fractions")
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

from __future__ import division
import numpy as np



### Problem 1.a.i ###
def QuantizeUniform(aNum,nBits):
    """
    Uniformly quantize signed fraction aNum with nBits
    """
    #Notes:
    #The overload level of the quantizer should be 1.0


    ### YOUR CODE STARTS HERE ###
    

    if (np.abs(aNum) >= 1):
    	code = np.int(2**(nBits - 1)-1)
    else :
    	code = np.int(np.floor(((2**nBits - 1) * np.abs(aNum) + 1)/2))

    #determine sign bit
    if(aNum >= 0):
    	aQuantizedNum = code
    else :
    	aQuantizedNum = code + 2**(nBits-1)


    ### YOUR CODE ENDS HERE ###

    return aQuantizedNum

    

### Problem 1.a.i ###
def DequantizeUniform(aQuantizedNum,nBits):
    """
    Uniformly dequantizes nBits-long number aQuantizedNum into a signed fraction
    """

    #aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###

    #determine sign bit
    s = aQuantizedNum >> (nBits-1)
    if (s == 0) :
    	sign = 1
    else:
    	sign = -1

    #determine code
    code = aQuantizedNum & ((1 << nBits-1)-1)

    number = (2*code)/(2**nBits - 1)

    aNum = sign*np.abs(number)
    ### YOUR CODE ENDS HERE ###

    return aNum



### Problem 1.a.ii ###
def vQuantizeUniform(aNumVec, nBits):
    """
    Uniformly quantize vector aNumberVec of signed fractions with nBits
    """

    aQuantizedNumVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    #Notes:
    #Make sure to vectorize properly your function as specified in the homework instructions

    ### YOUR CODE STARTS HERE ###

    s = np.zeros(aNumVec.size, dtype = int)
    code = np.zeros(aNumVec.size, dtype = int)

    #find sign bits for vector of inputs
    s[np.where(aNumVec < 0)] = 1

    #find codes
    overloadPos = np.where(np.abs(aNumVec) >= 1)
    normalPos = np.where(np.abs(aNumVec) < 1)

    if(np.size(overloadPos) > 0):
        code[overloadPos] = 2**(nBits-1)-1

    code[normalPos] = np.floor(((2**nBits - 1) * np.abs(aNumVec[normalPos]) + 1)/2)

    #final quantized vector 
    aQuantizedNumVec = code + s*(2**(nBits-1))



    ### YOUR CODE ENDS HERE ###

    return aQuantizedNumVec



### Problem 1.a.ii ###
def vDequantizeUniform(aQuantizedNumVec, nBits):
    """
    Uniformly dequantizes vector of nBits-long numbers aQuantizedNumVec into vector of  signed fractions
    """

    aNumVec = np.zeros_like(aQuantizedNumVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    sign = np.zeros(aQuantizedNumVec.size, dtype = int)
    code = np.zeros(aQuantizedNumVec.size, dtype = int)
    number = np.zeros(aQuantizedNumVec.size, dtype = float)

    #find sign bits
    s = aQuantizedNumVec >> (nBits-1)
    sign[np.where(s == 0)] = 1
    sign[np.where(s == 1)] = -1

    #find code
    code = aQuantizedNumVec & ((1 << nBits-1)-1)

    number = (2*code)/(2**nBits - 1)

    #final dequantized number
    aNumVec = np.multiply(sign,np.abs(number))


    ### YOUR CODE ENDS HERE ###

    return aNumVec



### Problem 1.b ###
def ScaleFactor(aNum, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point scale factor for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """
    #Notes:
    #The scale factor should be the number of leading zeros

    ### YOUR CODE STARTS HERE ###

    #step 1 - quantize number as R bit code
    R = 2**(nScaleBits)-1+nMantBits
    aQuantizedNum = QuantizeUniform(aNum, R)
    #ignore sign bit 
    code = aQuantizedNum & ((1 << R-1)-1)
    
    #step 2 - count number of leading zeros in quantized number
    k = -1
    msb = 0
    while (msb == 0 and k < R-1):
        msb = (1 << (R-2)) & code
        code = code << 1
        k = k + 1
    scale = k
        
    
#    if(code != 0):
#        scale = (R-1) - int(np.ceil(np.log2(code)))
#    else:
#        scale = (R-1)
    
    #scale can't be greater than 2**(Rs-1)
    if(scale > 2**(nScaleBits)-1):
        scale = 2**(nScaleBits)-1
        

    ### YOUR CODE ENDS HERE ###

    return scale

### Problem 1.b ###
def MantissaFP(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    ### YOUR CODE STARTS HERE ###

    #step 1 - quantize number as R bit code

    R = (1 << nScaleBits)-1+nMantBits
    aQuantizedNum = QuantizeUniform(abs(aNum), R)
    scale_max = 2**(nScaleBits)-1


    #step 2
    #method 1 - works for test case, not for problem 2

#    if(scale == scale_max):
#        #scale + 1 is because of sign bit
#        mantissa = aQuantizedNum & ((1 << (nMantBits - 1))-1)
#        
#    else:
#         #scale + 2 is because of sign bit plus one after leading zeros
#         nShifts =  (R-1)-(scale+1)
#         mantissa = aQuantizedNum & ((1 << nShifts) - 1)
#         if(mantissa > 2**(nMantBits-1)-1):
#             mantissa = mantissa >> (nShifts - (nMantBits-1))
             
    #method 2 - similar to what prateek did
             
    #mask scale factor zeros and sign bit
    aQuantizedNum = aQuantizedNum << (scale + 1)
    
    if(scale < scale_max):
        #mask leading one
        aQuantizedNum = aQuantizedNum - (1 << (R-1))
        #left shift to get rid of leading zero created by mask
        aQuantizedNum = aQuantizedNum << 1
    
    #right shift to get rid of trailing zeros
    mantissa = aQuantizedNum >> (R+1-nMantBits)
    
    if(aNum < 0):
        mantissa += 1 << nMantBits - 1
        

    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.b ###
def DequantizeFP(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for floating-point scale and mantissa given specified scale and mantissa bits
    """

    #aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE
    
    R = (1 << nScaleBits)-1+nMantBits
    #max value of scale
    scale_max = 2**(nScaleBits)-1
    #get sign bit
    s = mantissa >> (nMantBits - 1)
    #mantissa code
    mant_code = mantissa & ((1 << (nMantBits-1))-1)
    

    ### YOUR CODE STARTS HERE ###
    if(scale == scale_max):
        #stays the same
        mant_code = mant_code
    elif (scale == scale_max - 1):
        #add a leading one
        mant_code = mant_code + (2**(nMantBits-1))
    else:
        #add a leading one
        mant_code = mant_code + (2**(nMantBits-1))
        #add trailing zeros
        nzeros = (R-1) - (scale + 1 + nMantBits - 1)
        mant_code = mant_code << nzeros
        #add a one after abcd
        mant_code = mant_code + 2**(nzeros-1)
        
    aNum = DequantizeUniform(mant_code + s*(2**(R-1)), R)

    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.c.i ###
def Mantissa(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the block floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    ### YOUR CODE STARTS HERE ###

    R = (1 << nScaleBits) - 1 + nMantBits
    aQuantizedNum = QuantizeUniform(abs(aNum), R)

    #step 2 - cannot assume one after leading zeros
    #method 1 

#    #scale + 2 gets replaced with scale + 1 because there's no 1 after leading zeros
#    mantissa = aQuantizedNum & ((1 << (R-(scale+1)))-1)
#    if(mantissa > 2**(nMantBits-1)-1):
#        #one more right shift is needed because we are not ignoring the 1 after leading zeros
#        mantissa = mantissa >> (R - (scale+1) -nMantBits + 1)

    #method 2 - easier method
    #mask scale factor zeros and sign bit
    aQuantizedNum = aQuantizedNum << (scale + 1)
    
    #right shift to get rid of trailing zeros
    mantissa = aQuantizedNum >> (R+1-nMantBits)

    
    #add sign bit to first bit of mantissa
    if(aNum < 0):
        mantissa = mantissa + (1<<(nMantBits-1)) 

    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.c.i ###
def Dequantize(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for block floating-point scale and mantissa given specified scale and mantissa bits
    """

    aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    R = 2**(nScaleBits)-1+nMantBits
    #max value of scale
    scale_max = 2**(nScaleBits)-1
    #get sign bit
    s = mantissa >> (nMantBits - 1)
    #mantissa code
    mant_code = mantissa & ((1 << (nMantBits-1))-1)
    

    if(scale == scale_max or mant_code == 0):
        #stays the same
        mant_code = mant_code
    else:
        #add trailing zeros - no minus 1 after nMantBits because there's no leading one
        nzeros = R - (scale + 1 + nMantBits-1)
        mant_code = mant_code << nzeros
        #add a one after abcd
        mant_code = mant_code + 2**(nzeros-1)
        
    aNum = DequantizeUniform(mant_code + s*(2**(R-1)), R)
    
    

    ### YOUR CODE ENDS HERE ###

    return aNum
    


### Problem 1.c.ii ###
#there is a problem in this function - argghhh - FIXED IT!!
def vMantissa(aNumVec, scale, nScaleBits=3, nMantBits=5):
    """
    Return a vector of block floating-point mantissas for a vector of  signed fractions aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    mantissaVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
#    R = 2**(nScaleBits)-1+nMantBits
#    aQuantizedVec = vQuantizeUniform(aNumVec, R)
#    #find sign bit separately
#    s = aQuantizedVec >> (R-1)
#
#    #step 2 - cannot assume one after leading zeros
#
#    #scale + 2 gets replaced with scale + 1 because there's no 1 after leading zeros
#    mantissaVec = aQuantizedVec & ((1 << (R-(scale+1)))-1)
#    
#    ind = np.where(mantissaVec > 2**(nMantBits-1)-1)
#    #one more right shift is needed because we are not ignoring the 1 after leading zeros
#    mantissaVec[ind] = mantissaVec[ind] >> (R - (scale+1) -nMantBits + 1)
#    
#    #add sign bit to first bit of mantissa
#    mantissaVec = mantissaVec + s*(2**(nMantBits-1)) 
    
    #this is the way it's done in the non-vectorized function    
    R = (1 << nScaleBits) - 1 + nMantBits
    aQuantizedVec = vQuantizeUniform(abs(aNumVec),R)
    aQuantizedVec = aQuantizedVec << (scale+1)
    mantissaVec = aQuantizedVec >> (R+1-nMantBits)
    inds = np.where(aNumVec < 0)
    mantissaVec[inds] += (1<<(nMantBits-1)) 
    

    ### YOUR CODE ENDS HERE ###

    return mantissaVec


### Problem 1.c.ii ###
def vDequantize(scale, mantissaVec, nScaleBits=3, nMantBits=5):
    """
    Returns a vector of  signed fractions for block floating-point scale and vector of block floating-point mantissas given specified scale and mantissa bits
    """

    aNumVec = np.zeros_like(mantissaVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    R = 2**(nScaleBits)-1+nMantBits
    #max value of scale
    scale_max = 2**(nScaleBits)-1
    #get sign bit
    s = mantissaVec >> (nMantBits - 1)
    #mantissa code
    mant_code = mantissaVec & ((1 << (nMantBits-1))-1)
    
  
    if(scale == scale_max):
        #stays the same
        mant_code = mant_code
    else:
        #find non zero values of mantissa code
        ind = np.where(mant_code != 0)
        #add trailing zeros - no minus 1 after nMantBits because there's no leading one
        nzeros = R - (scale + 1 + nMantBits-1)
        mant_code[ind] = mant_code[ind] << nzeros
        #add a one after abcd
        mant_code[ind] = mant_code[ind] + 2**(nzeros-1)
        
    aNumVec = vDequantizeUniform(mant_code + s*(2**(R-1)), R)

    ### YOUR CODE ENDS HERE ###

    return aNumVec

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###



    pass # THIS DOES NOTHING

    ### YOUR TESTING CODE ENDS HERE ###

