# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 16:19:46 2018

@author: orchisama
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import psychoac as psy


def calculatePerceptualEntropy(sfBands, energy, SMR):
    """
    Calculates perceptual entropy according to the formula given in Pg 294
    """
    
    pe = 0.0
    thres = 10**-9
    
    for k in range(sfBands.nBands):
        #get total signal energy in band
        energy_band = np.sum(energy[sfBands.lowerLine[k]:sfBands.upperLine[k]+1])
        #make sure there are no numerical errors due to division by 0
        if(SMR[k] < thres):
            pe += 0.0
        else:
            pe += sfBands.nLines[k] * np.log2(1 + np.sqrt(energy_band/SMR[k]))
        
    return pe
    
    
def detectTransient(sfBands, data, MDCTdata, MDCTscale, sampleRate):
    """
    Function to detect if given block has a transient based on perceptual entropy
    sfBands - scale Factor Bands object
    MDCTdata - MDCT data in a block
    win - 0 or 1 (1 bit), type of window - sine/KBD
    """
    
#    if (win == 0):
#        window = 'sine'
#    else:
    window = 'kbd'
    
    #play around with the threshold value
    thres = 2500
    signal_inten = (2.0/psy.findWindowPower(window, np.size(data)))*(np.abs(MDCTdata)**2)
    SMR = psy.CalcSMRs(data, MDCTdata, MDCTscale, sampleRate, sfBands)
    #convert SMR to intensity for perceptual entropy calculation
    PE = calculatePerceptualEntropy(sfBands, signal_inten, psy.Intensity(SMR))
    
    if(PE > thres):
        return (True,PE)
    else:
        return (False, PE)
    
#-----------------------------------------------------------------------------
    
class WindowState(object):

    def __init__(self):
        self.state = 0

    def nextBuffer(self, is_transient):
        """ External method to advance the state machine's state. """
        if is_transient:
            return self.transient()
        else:
            return self.no_transient()

    def transient(self):
        """ Internal method to transition state based on onset presence """
        #if next buffer is detected as a transient
        if self.state == 0:
            #long window to transition window
            self.state = 1
        elif self.state == 1:
            #transition window to short window
            self.state = 2
        elif self.state == 2:
            #keep working with short windows
            self.state = 2
        elif self.state == 3:
            #short window to transition window
            self.state = 1
        return self.state

    def no_transient(self):
        """ Internal method to transition state when no onset is present. """
        #if next buffer is not a transient
        if self.state == 0:
            #if not a transient, keep using long windows
            self.state = 0
        elif self.state == 1:
            #if next buffer is not a transient and we are using a transition window,
            #keep using short windows
            self.state = 2
        elif self.state == 2:
            #if using short windows, switch to stop transition window
            self.state = 3
        elif self.state == 3:
            #after short transition window, switch back to long window
            self.state = 0

    return self.state

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":
    
    
    #does nothing
    pass

