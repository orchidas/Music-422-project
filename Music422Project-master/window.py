"""
window.py -- Defines functions to window an array of data samples
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
import matplotlib.pyplot as plt
import mdct 


# Clean up plots
plt.close('all')

### Problem 1.d ###
def SineWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray sine-windowed
    Sine window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###

    #return np.zeros_like(dataSampleArray) # CHANGE THIS

    N = len(dataSampleArray)
    n = np.arange(0,N)
    return np.sin(np.pi*(n+0.5)/N)*dataSampleArray





    ### YOUR CODE ENDS HERE ###


def HanningWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray Hanning-windowed
    Hann window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###

    #return np.zeros_like(dataSampleArray) # CHANGE THIS

    N = len(dataSampleArray)
    n = np.arange(0,N)
    return 0.5*(1 - np.cos(2*np.pi*(n+0.5)/N))*dataSampleArray




    ### YOUR CODE ENDS HERE ###


### Problem 1.d - OPTIONAL ###
def KBDWindow(dataSampleArray,alpha=4.):
    """
    Returns a copy of the dataSampleArray KBD-windowed
    KBD window is defined following pp. 108-109 and pp. 117-118 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###

    return np.zeros_like(dataSampleArray) # CHANGE THIS
    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    #pass # THIS DOES NOTHING
    
    
     #testSignal = np.arange(0,1,1./2048.)
     testSignal = np.ones((1000))
     
     sinWindow_test = SineWindow(testSignal)
     HannWindow_test = HanningWindow(testSignal)
    
    
     #plt.plot(HannWindow_test)
     #plt.plot()
     #plt.show()
     
     
     
     
     # 1.f)
     N = 1024
     n = np.arange(0,N)
     fs = 44100
     x = np.cos(2*np.pi*3000*n/fs)
     
     
     sinFFT = np.fft.fft(SineWindow(x))
     HanningFFT = np.fft.fft(HanningWindow(x))
     sinMDCT = mdct.MDCT(SineWindow(x),512,512)
     
     sinFFT_SPL = 96 + 10*np.log10(4/(N**2 *   np.mean(SineWindow(np.ones((N,)))**2)) * np.abs(sinFFT)**2) 
     
     HanningFFT_SPL = 96 + 10*np.log10(4/(N**2 *   np.mean(HanningWindow(np.ones((N,)))**2)) * np.abs(HanningFFT)**2) 
     
     sinMDCT_SPL = 96 + 10*np.log10(2/(   np.mean(SineWindow(np.ones((N,)))**2)) * np.abs(sinMDCT)**2) 
     
     f_fft = np.linspace(0,fs,N)
     f_MDCT = np.linspace(0,fs/2,N/2)
     
     plt.plot(f_fft,sinFFT_SPL, label = 'Sine FFT')
     plt.plot(f_fft,HanningFFT_SPL, label = 'Hann FFT')
     plt.plot(f_MDCT,sinMDCT_SPL, label = 'Sine MDCT')
     
     axes = plt.gca()
     axes.set_xlim([0,fs/2])
     
     plt.xlabel('Frequency (Hz)')
     plt.ylabel('SNR (dB)')
     plt.title('1f)')
     plt.grid(True)
     plt.legend()
     
     
     plt.savefig("Q1f.png")
     
     
     
    
    
    
    
    
    

    ### YOUR TESTING CODE ENDS HERE ###

