import numpy as np
from window import *

from mdct import *

from quantize import*

import matplotlib as mat, matplotlib.pyplot as plt

def SPL(intensity):
    """
    Returns the SPL corresponding to intensity (in units where 1 implies 96dB)
    """
    return np.maximum(-30,96.0 + 10*np.log10(intensity + np.finfo(np.float32).eps))

def Intensity(spl):
    """
    Returns the intensity (in units of the reference intensity level) for SPL spl
    """
    return 10.0**((spl-96.0)/10.0)

def Thresh(f):
    """Returns the threshold in quiet measured in SPL at frequency f (in Hz)"""

    f = np.maximum(10.0, f)

    quietSPL = 0

    quietSPL += 3.64*(f/1000.0)**(-0.8)

    quietSPL += -6.5*np.exp(-0.6*(f/1000.0 - 3.3)**2)

    quietSPL += (10.0**-3)*(f/1000.0)**4.0

    return quietSPL

def Bark(f):
    """Returns the bark-scale frequency for input frequency f (in Hz) """

    barkScale = 0 

    barkScale += 13.0*np.arctan(0.00076*f)

    barkScale += 3.5*np.arctan((f/7500.0)**2)

    return barkScale

class Masker:
    """
    a Masker whose masking curve drops linearly in Bark beyond 0.5 Bark from the
    masker frequency
    """

    def __init__(self,f,SPL,isTonal=True):
        """
        initialized with the frequency and SPL of a masker and whether or not
        it is Tonal
        """

        self.spl = SPL
        self.freq = f

        if isTonal:

            self.maskerDrop = 16.0 # in dB

        else : 

            self.maskerDrop = 6.0

    def vIntensityAtFreq(self,fVec):
        """Vectorized intensity of this masker at frequencies in fVec"""

        return self.vIntensityAtBark(Bark(fVec))

        
    def IntensityAtFreq(self,freq):
        """The intensity of this masker at frequency freq"""

        Z = Bark(freq)
        return self.IntensityAtBark(Z)

    def IntensityAtBark(self,z):
        """The intensity of this masker at Bark location z"""

        maskerZ = Bark(self.freq)

        maskeddB = 0

        maskeddB += self.spl - self.maskerDrop

        if abs(maskerZ - z) > 0.5 :

            if maskerZ > z:

                maskeddB += -27.0*(maskerZ - (z + 0.5))

            else:

                maskeddB += (-27.0 + 0.367*np.maximum(0.0, self.spl - 40.0))*(z - (maskerZ + 0.5))


        return  Intenstiy(maskeddB)

    def vIntensityAtBark(self,zVec):
        """The intensity of this masker at vector of Bark locations zVec"""

        maskerZ = Bark(self.freq)
        vmaskeddB = np.ones(len(zVec))*(self.spl - self.maskerDrop)

        ifmaskerZbigger = (maskerZ - zVec > 0.5)

        vmaskeddB[ifmaskerZbigger] += -27.0*(maskerZ - (zVec[ifmaskerZbigger] + 0.5))

        ifmaskerZsmaller = (maskerZ - zVec < -0.5)

        vmaskeddB[ifmaskerZsmaller] += (-27.0 + 0.367*np.maximum(0.0,self.spl - 40.0))*(zVec[ifmaskerZsmaller] - (maskerZ + 0.5))


        return Intensity(np.array(vmaskeddB))


# Default data for 25 scale factor bands based on the traditional 25 critical bands
cbFreqLimits = [0,100.0, 200.0, 300.0, 400.0, 510.0, 630.0, 770.0, 920.0, 1080.0, 1270.0,1480.0, 1720.0, 2000.0, 2320.0, 2700.0, 3150.0, 3700.0, 4400.0, 5300.0,
                6400.0, 7700.0, 9500.0, 12000.0, 15500.0]



def AssignMDCTLinesFromFreqLimits(nMDCTLines, sampleRate, flimit = cbFreqLimits):
    """
    Assigns MDCT lines to scale factor bands for given sample rate and number
    of MDCT lines using predefined frequency band cutofsampleRate passed as an array
    in flimit (in units of Hz). If flimit isn't passed it uses the traditional
    25 Zwicker & Fastl critical bands as scale factor bands.
    """

    lines = []

    lowerLine =  -1

    for line in range(len(flimit)):

        if flimit[line] <= (sampleRate/(2.0*nMDCTLines))*(nMDCTLines - 0.5):

            upperLine = int(flimit[line]/(sampleRate/(2.0*nMDCTLines)) - 0.5)

            lines.append(upperLine - lowerLine)

            lowerLine = upperLine

        else : 

            lines.append(nMDCTLines - lowerLine - 1)

            break



    return lines # TO REPLACE WITH YOUR CODE

class ScaleFactorBands:
    """
    A set of scale factor bands (each of which will share a scale factor and a
    mantissa bit allocation) and associated MDCT line mappings.

    Instances know the number of bands nBands; the upper and lower limits for
    each band lowerLimit[i in range(nBands)], upperLimit[i in range(nBands)];
    and the number of lines in each band nLines[i in range(nBands)]
    """

    def __init__(self,nLines):
        """
        Assigns MDCT lines to scale factor bands based on a vector of the number
        of lines in each band
        """

        self.lowerLine = np.zeros(len(nLines), dtype = int)
        self.nLines = np.array(nLines, dtype=int)
        self.upperLine = np.zeros(len(nLines) , dtype = int)
        self.upperLine[0] = nLines[0] - 1

        for i in range(1,len(nLines)):

            self.upperLine[i] = self.upperLine[i-1] + nLines[i]
            self.lowerLine[i] = self.upperLine[i-1] + 1

        self.nBands = len(nLines)



def getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Return Masked Threshold evaluated at MDCT lines.

    Used by CalcSMR, but can also be called from outside this module, which may
    be helpful when debugging the bit allocation code.
    """

    numBands = sfBands.nBands

    XhalfFFT = np.fft.fft(HanningWindow(data))[ : len(data)/2]

    XIntensity = (8.0/3.0)*(4.0/(len(data)**2))*(np.abs(XhalfFFT)**2)

    spl = SPL(XIntensity)

    maskersPerPeak = []

    Z = Bark((float(sampleRate)/len(data))*np.linspace(0.5 , (len(data) + 1)/2 , len(data)/2))

    masked = np.zeros(len(data)/2)

    for i in range(1,len(data)/2 - 1):

        if(XIntensity[i] > XIntensity[i-1] and XIntensity[i] > XIntensity[i + 1]):

            IAgg = XIntensity[i-1] + XIntensity[i] + XIntensity[i + 1]

            fInterp = (float(sampleRate)/len(data))*np.sum(np.array([i - 1, i , i + 1])*np.array([XIntensity[i-1],XIntensity[i],XIntensity[i + 1]]))/IAgg

            peakSPL = SPL(IAgg)

            if peakSPL > Thresh(fInterp):

                maskersPerPeak.append(Masker(fInterp , peakSPL, isTonal = True))

            XIntensity[i-1] = 0 
            XIntensity[i] = 0 
            XIntensity[i+1] = 0 


    for bands in range(sfBands.nBands):

            IAgg = 0.0 
            fInterp = 0.0

            for lines in range(sfBands.lowerLine[bands] , sfBands.upperLine[bands] + 1):

                IAgg += XIntensity[lines]

                fInterp += lines*XIntensity[lines]*(float(sampleRate)/len(data))


            fInterp /= IAgg
            ActSPL = SPL(IAgg)


            if ActSPL > Thresh(fInterp):


                maskersPerPeak.append(Masker(fInterp , ActSPL, isTonal = False))


    for i in range(len(maskersPerPeak)):

        masked += maskersPerPeak[i].vIntensityAtBark(Z)


    return masked + Intensity(Thresh((float(sampleRate)/len(data))*(np.linspace(0.5 , (len(data) + 1)/2 , len(data)/2))))


def CalcSMRs(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Set SMR for each critical band in sfBands.

    Arguments:
                data:       is an array of N time domain samples
                MDCTdata:   is an array of N/2 MDCT frequency lines for the data
                            in data which have been scaled up by a factor
                            of 2^MDCTscale
                MDCTscale:  is an overall scale factor for the set of MDCT
                            frequency lines
                sampleRate: is the sampling rate of the time domain samples
                sfBands:    points to information about which MDCT frequency lines
                            are in which scale factor band

    Returns:
                SMR[sfBands.nBands] is the maximum signal-to-mask ratio in each
                                    scale factor band

    Logic:
                Performs an FFT of data[N] and identifies tonal and noise maskers.
                Sums their masking curves with the hearing threshold at each MDCT
                frequency location to the calculate absolute threshold at those
                points. Then determines the maximum signal-to-mask ratio within
                each critical band and returns that result in the SMR[] array.
    """

    MDCTscaled = MDCTdata/(2**MDCTscale)

    MDCTIntensity = 4*(MDCTscaled)**2

    MDCTSPL = SPL(MDCTIntensity)

    masked = getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands)

    SMR = np.zeros(sfBands.nBands)

    for i in range(sfBands.nBands):


        SMR[i] = np.max(MDCTSPL[sfBands.lowerLine[i] : sfBands.upperLine[i] + 1] - masked[sfBands.lowerLine[i]:sfBands.upperLine[i] + 1])


    return SMR

#-----------------------------------------------------------------------------


#Testing code
if __name__ == "__main__":


    pass

    