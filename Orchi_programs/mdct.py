"""
- mdct.py -- Computes reasonably fast MDCT/IMDCT using numpy FFT/IFFT
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
import time 
import window as win
import matplotlib.pyplot as plt


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
    N = a+b
    n = np.arange(0,N)
    n0 = (b+1.0)/2
    temp = np.copy(data)    
    
    if(isInverse is False):
    #compute MDCT
        for k in range(0,N/2):
            data[k] = 2.0/N * np.sum(temp * np.cos(2*np.pi/N * (n + n0)*(k + 0.5)))
       
        temp = np.copy(data[:N/2])
        k = np.arange(0,N/2)
    
    else:   
    #compute IMDCT 
        for n in range(0,N):
            data[n] = 2.0*np.sum(temp * np.cos(2*np.pi/N*(n + n0)*(k + 0.5)))
        
        
    #return np.zeros_like( data ) # CHANGE THIS
    return data

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
    N = a+b
    n0 = (b+1.0)/2.0
    n = np.arange(0,N)
    k = np.arange(0,N/2)
    
    if(isInverse is False):
        
        pre_twiddle = np.exp(-1j*np.pi*n/N)
        data = pre_twiddle * data
        fft_data = np.fft.fft(data)
        post_twiddle = np.exp(-1j*2*np.pi/N*n0*(k+0.5))*(2.0/N)
        data = np.real(post_twiddle*fft_data[:N/2])
    
    else:
        data = IMDCT(data, a, b)
    
    return data

    ### YOUR CODE ENDS HERE ###

def IMDCT(data,a,b):

    ### YOUR CODE STARTS HERE ###
    N = a+b
    n0 = (b+1.0)/2.0
    n = np.arange(0,N)
    k = np.arange(0,N)
    
    #a[::-1] flips the data
    conc_data = np.concatenate([data, -data[::-1]])    
    pre_twiddle = np.exp(1j * 2*np.pi/N * n0 * k)
    data = pre_twiddle * conc_data
    ifft_data =  np.fft.ifft(data)
    post_twiddle = np.exp(1j*np.pi/N*(n+n0))*N
    data = np.real(post_twiddle * ifft_data)
            

    return data
    # CHANGE THIS
    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    #Q1b
    x = np.array([3,3,3,3,2,1,0,-1,-2,-3,-4,-4], dtype = float)
    x_prev = np.zeros(4, dtype = float)
   
    
    for i in range(4):
        
        if(i == 3):
            x_cur = np.zeros(4)
        else:
            x_cur = x[i*4:(i+1)*4]
            
        y = np.append(x_prev, x_cur)
        
        #do MDCT and inverse MDCT
        
        #slow MDCT        
        #res_cur = MDCTslow(y, 4, 4)
        
        #fast MDCT/IMDCT
        res_cur = MDCT(y,4,4,False)
        res_cur = IMDCT(res_cur,4,4)
        
        res_cur = res_cur/2.0
        
        if(i == 0):
            res = np.zeros(4, dtype = float)
        else:
            res = np.append(res, (res_prev+res_cur[:4]))
        
        res_prev = res_cur[4:]
        x_prev = x_cur
        
    print(res)

############################################################################
    
    #Q1c     
    #test time for MDCT and fast MDCT
#    x = np.random.randn(2048)
#    temp = np.copy(x)
#    
#    t_slow = time.time()
#    res = MDCTslow(x,1024,1024)
#    t_slow = time.time() - t_slow
#    
#    t_fast = time.time()
#    res = MDCT(temp,1024,1024,True)
#    t_fast = time.time() - t_fast
#    
#    print(t_slow/t_fast)

###########################################################################

    #Q1f
    N = 1024
    n = np.arange(N)
    fs = 44100.0
    x = np.cos(2*np.pi*3000*n/fs)
    x1 = np.copy(x)
    x2 = np.copy(x)
        
    x_sine = win.SineWindow(x)
    X_sine_fft = np.fft.fft(x_sine)
    #we put 2/N in the forward transform, so we need to remove a factor of 4/N^2
    sine_FFT_spl = 96 + 10*np.log10(4.0/((N**2)*0.5)*(np.abs(X_sine_fft)**2))
    X_sine_MDCT = MDCT(x_sine,512,512,False)
    sine_MDCT_spl = 96 + 10*np.log10(4.0*(np.abs(X_sine_MDCT)**2))
    
    x_hann = win.HanningWindow(x1)
    X_hann_fft = np.fft.fft(x_hann)
    hann_FFT_spl = 96 + 10*np.log10(4.0/((N**2)*0.375)*(np.abs(X_hann_fft)**2))
    
    freqs_fft = np.linspace(0,fs/2,N/2+1)
    #omit last half of spectrum because it is symmetric
    freqs_fft = freqs_fft[:N/2]
    freqs_mdct = np.linspace(0,fs/2,N/2+1)
    freqs_mdct = freqs_mdct[:N/2]
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(n, x_sine)
    plt.subplot(2,1,2)
    plt.plot(n, x_hann)
    
    ax = plt.figure()
    plt.plot(freqs_fft,sine_FFT_spl[:N/2],'r', label = 'Sine window FFT')
    plt.plot(freqs_mdct, sine_MDCT_spl,'b', label = 'Sine window MDCT')
    plt.plot(freqs_fft, hann_FFT_spl[:N/2],'g', label = 'Hann window FFT')    
    plt.xlabel('Frequency in Hz')
    plt.ylabel('SPL in dB')
    plt.legend()
    plt.show()

    
    
    pass # THIS DOES NOTHING

    ### YOUR TESTING CODE ENDS HERE ###

