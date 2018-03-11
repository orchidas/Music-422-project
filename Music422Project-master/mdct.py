"""
- mdct.py -- Computes reasonably fast MDCT/IMDCT using numpy FFT/IFFT
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
import mdct_ as mdct_true
import matplotlib.pyplot as plt
import time

# Clean up plots
plt.close('all')

### Problem 1.a ###
def MDCTslow(data, a, b, isInverse=False):
    """
    Slow MDCT algorithm for window length a+b following pp. 130 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    (and where 2/N factor is included in forward transform instead of inverse)
    a: left half-window length
    b: right half-window length
    """

    ### YOUR CODE STARTS HERE ###

    
    N = float(a+b)
    n_0 = (float(b)+1.0)/2.0
    


    if(isInverse ==False):
        data_altered = np.zeros(int(N/2),)
        n = np.arange(0.,N)
        for k in range(0,int(N/2)):
            data_altered[k] = (2/N)*np.sum(data*np.cos(((2*np.pi)/N)*(n+n_0)*(float(k)+0.5)   ))
        
        
        return data_altered
        
        
        
        
        
        
    elif(isInverse ==True):
        data_altered = np.zeros(int(N),)
        k = np.arange(0.,N/2)
        for n in range(0,int(N)):
            data_altered[n] = 2*np.sum(data*np.cos((2*np.pi/N)*(n+n_0)*(k+0.5)))
        
        return data_altered
    
    else:       
         return np.zeros_like( data ) # CHANGE THIS   





    ### YOUR CODE ENDS HERE ###

### Problem 1.c ###
def MDCT(data, a, b, isInverse=False):
    """
    Fast MDCT algorithm for window length a+b following pp. 141-143 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    (and where 2/N factor is included in forward transform instead of inverse)
    a: left half-window length
    b: right half-window length
    """

    ### YOUR CODE STARTS HERE ###

    #return np.zeros_like( data ) # CHANGE THIS
    
    
    N = float(a+b)
    n_0 = (float(b)+1.0)/2.0
    k = np.arange(0.,N/2)
    k_imdct = np.arange(0.,N)
    n = np.arange(0.,N)
    
    
    


    if(isInverse==False):
        preTwiddle = np.exp(-1j*np.pi*n/N)
        postTwiddle = np.exp(-1j*2*np.pi*n_0 * (k+0.5)/N)
        return (2/N)*np.real(postTwiddle*np.fft.fft(data*preTwiddle)[0:int(N/2)])
        
        
        
        
    elif(isInverse==True):
        preTwiddle = np.exp(1j*2*np.pi*k_imdct*n_0/N)
        postTwiddle = np.exp(1j*np.pi*(n+n_0)/N)
        #return N*np.real(postTwiddle*np.fft.ifft(preTwiddle*np.concatenate((data, -np.flip(data,0) ))))
        return N*np.real(postTwiddle*np.fft.ifft(preTwiddle*np.concatenate((data, -data[::-1] ))))
        








    ### YOUR CODE ENDS HERE ###

def IMDCT(data,a,b):

    ### YOUR CODE STARTS HERE ###

    #return np.zeros_like( data ) # CHANGE THIS



    return MDCT(data, a, b, isInverse=True)
    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    pass # THIS DOES NOTHING
    
    
    
    # 1. b)
    
    x = np.array([3,3,3,3,2,1,0,-1,-2,-3,-4,-4])
    
    firstBlock = np.array([0,0,0,0])
    lastBlock = np.array([0,0,0,0])
    
    output = np.zeros((12,))
    
    tempLast4 = np.zeros((4,))
    
    for i in range(0,4):
        if(i == 0):
            temp = np.concatenate((firstBlock,x[0:4]))
            temp = 0.5 * MDCTslow(MDCTslow(temp, 4,4, False), 4, 4, True)
            #output[0:4] = temp[0:4] + firstBlock
            tempLast4 = temp[4:8]
        elif(i == 3):
            temp = np.concatenate((x[(i-1)*4:i*4],lastBlock))
            temp = 0.5 * MDCTslow(MDCTslow(temp, 4,4, False), 4, 4, True)
            output[(i-1)*4: (i-1)*4+4] = temp[0:4] + tempLast4
            
        else:
            #temp = np.concatenate((x[i*4:i*4+4],x[(i-1)*4:i*4]))
            temp = x[(i-1)*4:(i-1)*4+8]
            temp = 0.5 * MDCTslow(MDCTslow(temp, 4,4, False), 4, 4, True)
            output[(i-1)*4:(i-1)*4+4] = temp[0:4] + tempLast4
            tempLast4 = temp[4:8]
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # My other testing
    # Make a test signal
    testSignal = np.arange(0,1,1./2048.)
    #testSignal = np.arange(0,1,0.2)
    halfBlock = 1024
    
    # Test the mdct against solution
    mdct_test_true = mdct_true.MDCTslow(testSignal, halfBlock,halfBlock, False)
    mdct_test = MDCTslow(testSignal, halfBlock,halfBlock, False)
    
    mdct_err = np.sum(mdct_test_true - mdct_test)
    
    # Test the Imdct against solution
    Imdct_test_true = mdct_true.MDCTslow(mdct_test_true, halfBlock,halfBlock, True)
    Imdct_test = MDCTslow(mdct_test_true, halfBlock,halfBlock, True)
    
    Imdct_err = np.sum(Imdct_test_true -Imdct_test)

    
    #plt.plot(mdct_test)
    #plt.plot(mdct_test_true)
    #plt.show()
    
    
    
    #Testing fast MDCT
    
    mdctFast_test_true = mdct_true.MDCT(testSignal, halfBlock,halfBlock, False)
    mdctFast_test = MDCT(testSignal, halfBlock,halfBlock, False)
    
    mdctFast_err = np.sum(mdctFast_test_true - mdctFast_test)
    
    # Test the Imdct against solution
    ImdctFast_test_true = mdct_true.IMDCT(mdctFast_test_true, halfBlock,halfBlock)
    ImdctFast_test = IMDCT(mdctFast_test, halfBlock,halfBlock)
    
    ImdctFast_err = np.sum(ImdctFast_test_true -ImdctFast_test)
    
    
    
    
    
    # Timing for 1c
    # slow
    t1 = time.time()
    mdct_test = MDCTslow(testSignal, halfBlock,halfBlock, False)
    Imdct_test = MDCTslow(mdct_test_true, halfBlock,halfBlock, True)
    t2 = time.time()
    
    time_slow = t2-t1
    
    # fast
    t3 = time.time()
    mdctFast_test = MDCT(testSignal, halfBlock,halfBlock, False)
    ImdctFast_test = IMDCT(mdctFast_test, halfBlock,halfBlock)
    t4 = time.time()
    
    time_fft = t4-t3
    
    speedup = time_slow/time_fft
    
    
    
    
    
    
    

    ### YOUR TESTING CODE ENDS HERE ###

