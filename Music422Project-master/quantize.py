"""
quantize.py -- routines to quantize and dequantize floating point values
between -1.0 and 1.0 ("signed fractions")
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np

### Problem 1.a.i ###
def QuantizeUniform(aNum,nBits):
    """
    Uniformly quantize signed fraction aNum with nBits
    """
    #Notes:
    #The overload level of the quantizer should be 1.0

    #aQuantizedNum = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    if(aNum >= 0):
        s = 0
    else:
        s = 1
        
    if(np.abs(aNum)>=1):
        code = 2**(nBits-1) - 1
    else:
        code = np.floor(((2**(nBits)-1)*np.abs(aNum)+1)/2)
        
        
    if(s ==0):
       return int(code)
    else:
       return int(code + 2**(nBits-1))
            

    ### YOUR CODE ENDS HERE ###

    #return 

### Problem 1.a.i ###
def DequantizeUniform(aQuantizedNum,nBits):
    """
    Uniformly dequantizes nBits-long number aQuantizedNum into a signed fraction
    """

    aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    if (aQuantizedNum>>(nBits-1)):
        sign = -1
    else:
        sign = 1
    
    if (sign == -1):
        aNum =  sign*float(2*np.abs(aQuantizedNum - 2**(nBits-1)))/(2**nBits -1)    
    else:
        aNum =  sign*float(2*np.abs(aQuantizedNum))/(2**nBits -1) 
    


    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.a.ii ###
def vQuantizeUniform(aNumVec, nBits):
    """
    Uniformly quantize vector aNumberVec of signed fractions with nBits
    """

    #aQuantizedNumVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    #Notes:
    #Make sure to vectorize properly your function as specified in the homework instructions

    ### YOUR CODE STARTS HERE ###
    
    s = np.signbit(aNumVec)
    
    aQuantizedNumVec = s*2**(nBits-1) +  np.where(np.abs(aNumVec.copy()) >= 1 ,2**(nBits-1) - 1 ,np.floor(((2**nBits-1)*np.abs(aNumVec.copy())+1)/2))
    

    ### YOUR CODE ENDS HERE ###

    return aQuantizedNumVec.astype(int)

### Problem 1.a.ii ###
def vDequantizeUniform(aQuantizedNumVec, nBits):
    """
    Uniformly dequantizes vector of nBits-long numbers aQuantizedNumVec into vector of  signed fractions
    """

    #aNumVec = np.zeros_like(aQuantizedNumVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
   
    
    s = np.where((aQuantizedNumVec.copy()).astype(float)>=2**(nBits-1),-1,1)
    s_zer = np.where((aQuantizedNumVec.copy()).astype(float)>=2**(nBits-1),1,0)
    
    aNumVec = s.astype(float)*2*np.abs((aQuantizedNumVec.copy()).astype(float)-s_zer.astype(float)*2**(nBits-1))/(2**nBits - 1)
    #aNumVec = s*2*np.abs((aQuantizedNumVec.copy()).astype(float))/(2**nBits - 1)
    
    

    ### YOUR CODE ENDS HERE ###

    return aNumVec

### Problem 1.b ###
def ScaleFactor(aNum, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point scale factor for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """
    #Notes:
    #The scale factor should be the number of leading zeros

    #scale = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    #R = nScaleBits+nMantBits
    tempR = 2**nScaleBits - 1 + nMantBits
    
    aNumQuant = QuantizeUniform(aNum,tempR)
    scale = 0
    
    
    if (aNumQuant<= 2**(tempR-1)):
    
        for i in range(0, 2**nScaleBits - 1):
            if (aNumQuant<2**(tempR - i - 2)):
                scale = scale + 1
    else:
        for i in range(0, 2**nScaleBits - 1):
            if (aNumQuant - 2** (tempR-1) <2**(tempR - i - 2)):
                scale = scale + 1
    
    
    
    

    ### YOUR CODE ENDS HERE ###

    return scale

### Problem 1.b ###
def MantissaFP(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """
    
    #mantissa = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE




    ### YOUR CODE STARTS HERE ###
    
    tempR = 2**nScaleBits - 1 + nMantBits
    
    aNumQuant = QuantizeUniform(aNum,tempR)
    
    if (scale == 2**nScaleBits - 1):
        mantissa = np.signbit(aNum)*2**(nMantBits-1) + ((aNumQuant-np.signbit(aNum)*2**(tempR-1))>>(tempR - scale -nMantBits))
    else:
        mantissa = np.signbit(aNum)*2**(nMantBits-1) + ((aNumQuant-np.signbit(aNum)*2**(tempR-1))>>(tempR - scale -nMantBits -1)) - 2**(nMantBits-1)
    
    
    
    
    
    
    
    

    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.b ###
def DequantizeFP(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for floating-point scale and mantissa given specified scale and mantissa bits
    """

    #aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    R = 2**nScaleBits - 1 + nMantBits
    
    if(mantissa>=2**(nMantBits-1)):
        s = 1
    else:
        s = 0
        
    if (scale==2**(nScaleBits)-1):
        num = s*2**(R-1) + (mantissa - s*2**(nMantBits-1))
    else:
        num = s*2**(R-1) + ((mantissa - s*2**(nMantBits-1)) + 2**(nMantBits-1) <<(2**nScaleBits- 2 - scale)) + (1<<(2**(nScaleBits)- 3 - scale))
        
    
        
    
    aNum = DequantizeUniform(num,R)
    
    

    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.c.i ###
def Mantissa(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the block floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    #mantissa = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    



    ### YOUR CODE STARTS HERE ###
    
    tempR = 2**nScaleBits - 1 + nMantBits
    
    aNumQuant = QuantizeUniform(aNum,tempR)
    
    if (scale == 2**nScaleBits - 1):
        mantissa = np.signbit(aNum)*2**(nMantBits-1) + ((aNumQuant-np.signbit(aNum)*2**(tempR-1))>>(tempR - scale -nMantBits ))
    else:
        mantissa = np.signbit(aNum)*2**(nMantBits-1) + ((aNumQuant-np.signbit(aNum)*2**(tempR-1))>>(tempR - scale -nMantBits )) 
    
    
    
    

    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.c.i ###
def Dequantize(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for block floating-point scale and mantissa given specified scale and mantissa bits
    """

    #aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    R = 2**nScaleBits - 1 + nMantBits
    
    if(mantissa>=2**(nMantBits-1)):
        s = 1
    else:
        s = 0
        
    if (scale==2**(nScaleBits)-1):
        #num = s*2**(R-1) + (mantissa - s*2**(nMantBits-1))
        num = s*2**(R-1) + (mantissa - s*2**(nMantBits-1))
    elif(mantissa - s*2**(nMantBits-1) ==0):
            num = s*2**(R-1)
    else:
        #num = s*2**(R-1) + ((mantissa - s*2**(nMantBits-1)) + 2**(nMantBits-1) <<(2**nScaleBits- 2 - scale)) + (1<<(2**(nScaleBits)- 3 - scale))
        num = s*2**(R-1) + ((mantissa - s*2**(nMantBits-1))  <<(2**nScaleBits- 1 - scale)) + (1<<(2**(nScaleBits)- 2 - scale))
        
    
        
    
    aNum = DequantizeUniform(num,R)
    
    
    
    
    

    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.c.ii ###
def vMantissa(aNumVec, scale, nScaleBits=3, nMantBits=5):
    """
    Return a vector of block floating-point mantissas for a vector of  signed fractions aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    #mantissaVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    
    tempR = 2**nScaleBits - 1 + nMantBits
    
    aNumQuant = vQuantizeUniform(aNumVec,tempR)



    if (scale == 2**nScaleBits - 1):
        mantissaVec = np.signbit(aNumVec.copy())*2**(nMantBits-1) + ((aNumQuant.copy()-np.signbit(aNumVec.copy())*2**(tempR-1))>>(tempR - scale -nMantBits ))
    else:
        mantissaVec = np.signbit(aNumVec.copy())*2**(nMantBits-1) + ((aNumQuant-np.signbit(aNumVec.copy())*2**(tempR-1))>>(tempR - scale -nMantBits )) 
    
    
    
    

    ### YOUR CODE ENDS HERE ###

    return mantissaVec

### Problem 1.c.ii ###
def vDequantize(scale, mantissaVec, nScaleBits=3, nMantBits=5):
    """
    Returns a vector of  signed fractions for block floating-point scale and vector of block floating-point mantissas given specified scale and mantissa bits
    """

    #aNumVec = np.zeros_like(mantissaVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    
    
    
    
    R = 2**nScaleBits - 1 + nMantBits
    
    s = np.where(mantissaVec>=2**(nMantBits-1),1,0)

        
    if (scale==2**(nScaleBits)-1):
        #num = s*2**(R-1) + (mantissa - s*2**(nMantBits-1))
        num = s*2**(R-1) + (mantissaVec - s*2**(nMantBits-1))
    else:
        num = np.where(mantissaVec - s*2**(nMantBits-1) ==0,s*2**(R-1),s*2**(R-1) + ((mantissaVec - s*2**(nMantBits-1))  <<(2**nScaleBits- 1 - scale)) + (1<<(2**(nScaleBits)- 2 - scale)))
        
        
        
    #elif(mantissaVec.copy() - s*2**(nMantBits-1) ==0):
    #        num = s*2**(R-1)
    #else:
        #num = s*2**(R-1) + ((mantissaVec.copy() - s*2**(nMantBits-1))  <<(2**nScaleBits- 1 - scale)) + (1<<(2**(nScaleBits)- 2 - scale))
        
    
    
    aNumVec = vDequantizeUniform(num,R)
    
    
    

    ### YOUR CODE ENDS HERE ###

    return aNumVec

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    pass # THIS DOES NOTHING

    ### YOUR TESTING CODE ENDS HERE ###

