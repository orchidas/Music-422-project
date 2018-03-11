"""
window.py -- Defines functions to window an array of data samples
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np

import matplotlib.pyplot as plt

from mdct import*

### Problem 1.d ###
def SineWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray sine-windowed
    Sine window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###


    N = len(dataSampleArray)

    theta = np.pi*(np.arange(N) + 0.5) / N 

    return dataSampleArray*np.sin(theta)
    

    ### YOUR CODE ENDS HERE ###


def HanningWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray Hanning-windowed
    Hann window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###

    N = len(dataSampleArray)

    theta = 2*np.pi*(np.arange(N) + 0.5) / N 

    return dataSampleArray*0.5*(1.0 - np.cos(theta))

    ### YOUR CODE ENDS HERE ###


### Problem 1.d - OPTIONAL ###
def KBDWindow(dataSampleArray,alpha=4.):
    """
    Returns a copy of the dataSampleArray KBD-windowed
    KBD window is defined following pp. 108-109 and pp. 117-118 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###


    N = len(dataSampleArray)

    kbdArg = (4.0 * np.arange(N/2.0 + 1) / N - 1.0 )

    besselOut = np.i0(alpha*np.pi*np.sqrt(1.0 - kbdArg**2 ))/np.i0(alpha*np.pi)

    kbdWind = np.zeros(N, dtype = float )

    kbdWind[: N/2] = np.cumsum(besselOut[0:len(besselOut) - 1]**2)/np.sum(besselOut**2)

    #kbdWind [ N/2 : ] = np.flip(kbdWind[: N/2],0)
    kbdWind [ N/2 : ] = np.flipud(kbdWind[: N/2])

    kbdWind = np.sqrt(kbdWind)

    return dataSampleArray*kbdWind


def compose_kbd_window(dataSampleArray, left, right, left_alpha=4., right_alpha=4.):
    """ Compose a hybrid Kaiser-Bessel Derived window for block-switched MDCT
    windows. Parameters left, right control the size of the window segments,
    while the alpha parameters tune the frequency selectivity vs. rolloff. """

    # Make sure that left + right is the size of the window provided
    if left + right != len(dataSampleArray):
        msg = 'Signal size, {} , must match the composed size left+right: {}' \
               .format(str(len(dataSampleArray)), str(left+right))
        raise ValueError(msg)
    
    # Create a window for size left
    a_ones = np.ones(2*left)
    a_window = KBDWindow(a_ones, alpha=left_alpha)[:left]

    # Create a window for size right
    b_ones = np.ones(2*right)
    b_window = KBDWindow(b_ones, alpha=right_alpha)[right:]
    
    #for plotting only    
#    if(left == 64 or right == 64):
#        plt.figure
#        plt.plot(np.concatenate([a_window, b_window]))
#        plt.show()
        
    return dataSampleArray * np.concatenate([a_window, b_window])
    
    ### YOUR CODE ENDS HERE ###
    
    
def compose_sine_window(dataSampleArray, left, right):
    """ Compose a hybrid Kaiser-Bessel Derived window for block-switched MDCT
    windows. Parameters left, right control the size of the window segments,
    while the alpha parameters tune the frequency selectivity vs. rolloff. """

    # Make sure that left + right is the size of the window provided
    if left + right != len(dataSampleArray):
        msg = 'Signal size, {} , must match the composed size left+right: {}' \
               .format(str(len(dataSampleArray)), str(left+right))
        raise ValueError(msg)
    
    # Create a window for size left
    a_ones = np.ones(2*left)
    a_window = SineWindow(a_ones)[:left]

    # Create a window for size right
    b_ones = np.ones(2*right)
    b_window = SineWindow(b_ones)[right:]
    w = np.concatenate([a_window, b_window])
            
    return (dataSampleArray * w)



#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###



        N = 1024
        Fs = 44100
        f = 3000.0
        n = np.arange(N)
        x = np.cos(2 * np.pi * f * n / Fs)
        
        avgSine = 0.5
        avgKBD = 1.0 / N * np.sum(KBDWindow(np.ones(N))**2)
        avgHanning = 0.375

        signalSineWindowed = SineWindow(x)
        signalHanningWindowed = HanningWindow(x)
        signalKBDWindowed = KBDWindow(x)

        SPLSine = 96.0 + 10.0 * np.log10(4.0 / (avgSine * N ** 2) * np.abs(np.fft.fft(signalSineWindowed)[ : N / 2]) ** 2)

        SPLHanning = 96.0 + 10.0 * np.log10(4.0 / (avgHanning * N ** 2) * np.abs(np.fft.fft(signalHanningWindowed)[ : N / 2]) ** 2)

        SPLKBD = 96.0 + 10.0 * np.log10(4.0 / (avgKBD * N ** 2) * np.abs(np.fft.fft(signalKBDWindowed)[ : N / 2]) ** 2)

        MDCTSine = MDCT(signalSineWindowed, N / 2, N / 2)
        MDCTKBD = MDCT(signalKBDWindowed, N / 2, N / 2)
        SPLSine_MDCT = 96.0 + 10.0 * np.log10(2.0 / avgSine * MDCTSine ** 2)
        SPLmdctKBD_MDCT = 96.0 + 10.0 * np.log10(2.0 / avgKBD * MDCTKBD ** 2)

        freqs= np.arange(N / 2) * Fs / N
  
        plt.plot(freqs, SPLSine, freqs, SPLHanning, freqs + (Fs / (2.0*N)), SPLSine_MDCT, freqs + (Fs / (2.0*N)), SPLmdctKBD_MDCT)
        plt.ylabel('SPL in dB')
        plt.xlim(0, Fs / 2)
        plt.xlabel('Frequency (Hz)')
        plt.title('Comparison of FFT and MDCT With Different Windows')
        plt.legend(('FFT + Sine','FFT + Hanning','MDCT + Sine','MDCT + KBD(alpha = 4)'))

        plt.show()


    ### YOUR TESTING CODE ENDS HERE ###

