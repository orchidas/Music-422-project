"""
- mdct.py -- Computes reasonably fast MDCT/IMDCT using numpy FFT/IFFT
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np

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

    NhalfLength = (a + b)/2

    kFreq = np.arange(0,NhalfLength)
    nTime = np.arange(0,2*NhalfLength)


    N0 = (1.0 + b)/2.0


    if isInverse:
        

        slowIMDCT = np.zeros(2*NhalfLength , dtype = float)

        for i in nTime:

            slowIMDCT[i] = np.sum(data*np.cos(2 * np.pi * (kFreq + 0.5) * (i + N0) / (2*NhalfLength)))

        slowIMDCT *= 2

        return slowIMDCT


    else : 
        

        slowMDCT = np.zeros(NhalfLength , dtype = float)

        for i in kFreq:

            slowMDCT[i] = np.sum(data*np.cos(2 * np.pi * (i + 0.5) * (nTime + N0)/(2*NhalfLength)))

        slowMDCT *= 1.0/NhalfLength


        return slowMDCT

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


    NhalfLength = (a + b)/2


    kFreq = np.arange(0,NhalfLength)
    nTime = np.arange(0,2*NhalfLength)


    N0 = (1.0 + b)/2.0


    if isInverse : 


        preTwiddle = np.exp((np.pi * 2j * N0/(2*NhalfLength))* np.arange(2*NhalfLength, dtype = float))

        postTwiddle = (2*NhalfLength) * np.exp ( (1j* np.pi / (2*NhalfLength)) * np.linspace(N0 , 2*NhalfLength + N0 - 1 , 2*NhalfLength))

        #newData = np.concatenate((data, -np.flip(data,0)))
        newData = np.concatenate((data, -np.flipud(data)))

        fastIMDCT = np.real(postTwiddle * np.fft.ifft(preTwiddle * newData))


        return fastIMDCT


    else : 

        preTwiddle = data*np.exp((np.pi * (-1j)/(2*NhalfLength))* np.arange(2*NhalfLength, dtype = float))

        postTwiddle = (1.0/NhalfLength)*np.exp (( -2j * np.pi * N0 / (2*NhalfLength)) * np.linspace(0.5 , NhalfLength - 0.5 , NhalfLength))

        preTwiddlePositiveFFT = np.fft.fft(preTwiddle)[ : NhalfLength]

        fastMDCT = np.real(postTwiddle * preTwiddlePositiveFFT)

        return fastMDCT



    ### YOUR CODE ENDS HERE ###

def IMDCT(data,a,b):

    ### YOUR CODE STARTS HERE ###

    return MDCT(data , a , b , isInverse = True) 

    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

            x = np.array([3, 3, 3, 3, 2 , 1 , 0, -1, -2, -3, -4, -4], dtype = float)

            a = 4
            b = a

            newLen = a + b + len(x)

            xout = np.zeros(0, dtype = float)
            priorBlock = np.zeros(a, dtype = float)
            for n in range(newLen / a - 1):
                newData = np.concatenate((np.zeros(a), x, np.zeros(b)))[a * n:b * n + (a+b)]/np.sqrt(2.0)
                MDCTout = MDCTslow(newData, a, b)
                IMDCToutput = MDCTslow(MDCTout, a, b,True)/np.sqrt(2.0)
                xout = np.append(xout, IMDCToutput[:b] + priorBlock)
                priorBlock = IMDCToutput[a:]

            print 'Input :', str(x)
            print 'Output :', str(xout)
        
    ### YOUR TESTING CODE ENDS HERE ###

