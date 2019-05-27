#!/usr/bin/ipython
# -*-Coding:Utf-8 -*

# This program is free software; you can redistribute and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# This is a modified version of the GaussSum code to generate the IR Spectra

import math
import numpy

colors = "bgrcmyk"

def lorentzian(x,peak,height,width):
    """The lorentzian curve.

    f(x) = a/(1+a)

    where a is FWHM**2/4
    """
    a = width**2./4.
    return float(height)*a/( (peak-x)**2 + a )

class Spectrum(object):
    """Convolutes and stores spectrum data.

    Usage:
     Spectrum(start,end,numpts,peaks,width,formula)

    where
     peaks is [(pos,height),...]
     formula is a function such as gaussianpeak or delta
    

    >>> t = Spectrum(0,50,11,[[(10,1),(30,0.9),(35,1)]],5,delta)
    >>> t.spectrum
    array([[ 0.        ],
           [ 1.        ],
           [ 1.        ],
           [ 1.        ],
           [ 0.        ],
           [ 0.89999998],
           [ 1.89999998],
           [ 1.89999998],
           [ 1.        ],
           [ 0.        ],
           [ 0.        ]],'d')
    """
    def __init__(self,start,end,numpts,peaks,width,formula):
        self.start = start
        self.end = end
        self.numpts = numpts
        self.peaks = peaks
        self.width = width
        self.formula = formula

        # len(peaks) is the number of spectra in this object
        self.spectrum = numpy.zeros( (numpts,len(peaks)),"d")
        self.xvalues = numpy.arange(numpts)*float(end-start)/(numpts-1) + start
        for i in range(numpts):
            x = self.xvalues[i]
            for spectrumno in range(len(peaks)):
                for (pos,height) in peaks[spectrumno]:
                    self.spectrum[i,spectrumno] = self.spectrum[i,spectrumno] + formula(x,pos,height,width)



def activity_to_intensity(activity, frequency, excitation, temperature):
    """Convert Raman acitivity to Raman intensity according to
    Krishnakumar et al, J. Mol. Struct., 2004, 702, 9."""

    excitecm = 1 / (1e-7 * excitation)
    f = 1e-13
    above = f * (excitecm - frequency)**4 * activity
    exponential = -6.626068e-34 * 299792458 * frequency / (1.3806503e-23 * temperature)
    below = frequency * (1 - math.exp(exponential))
    return above / below

def get_scaling_factors(filename, scale):
    """Read in scaling factors from an existing output file.

    Note: Scale is prepopulated with the general scaling factor
    """
    inputfile = open(filename, "r")
  #  inputfile = "TD.log"
    line = inputfile.readline()
    line = inputfile.readline()
    i = 0
    line = inputfile.readline().split('\t')
    while len(line) > 6: # Read in the individual scaling factors
        sf = line[-2]
        if sf != '':
            scale[i] = float(sf)
        i += 1
        line = inputfile.readline().split('\t')
    inputfile.close()
    return scale


def Vibfreq(json,logfilename,start,end,numpts,FWHM,typeofscale,scalefactor,excitation,temperature):
    
    def dofreqs(name,act):

        filename=name+"Spectrum.txt"
        freq = json["results"]["freq"]["vibrational_freq"].copy() # Copy so that it won't be changed by the routine
        if (json["results"]["freq"]["vibrational_sym"]):
            vibsyms = json["results"]["freq"]["vibrational_sym"]
        else:
            vibsyms = ['?'] * len(freq)

        # Handle the scaling of the frequencies
        scale = [scalefactor] * len(freq)
        if typeofscale == "Gen":
            general = True

        for i in range(len(freq)): # Scale the freqs
            freq[i] = freq[i]*scale[i]

        # Convolute the spectrum
        #print(start,end,numpts,list(zip(freq,act)),FWHM,lorentzian)
        spectrum = Spectrum(start,end,numpts,[list(zip(freq,act))],FWHM,lorentzian)
        if name == "Raman":
            intensity = [activity_to_intensity(activity, frequency, excitation, temperature)
                         for activity, frequency in zip(act, freq)]
            spectrum_intensity = Spectrum(start,end,numpts,[list(zip(freq, intensity))],
                                           FWHM,lorentzian)

        outputfile = open(filename,"w")
 ##       screen.write("Writing scaled spectrum to "+filename+"\n")
        outputfile.write("Spectrum\t%s\t\tNormal Modes\n" % ["", "\t"][name=="Raman"])
        outputfile.write("Freq (cm-1)\t%s act\t%s\tMode\tLabel\tFreq (cm-1)\t%s act\t" % (name, ["", "Intensity\t"][name=="Raman"],name))
        outputfile.write("%sScaling factors\tUnscaled freq\n" % ["", "Intensity\t"][name=="Raman"])
        width = end-start
        for x in range(0,numpts):
            if spectrum.spectrum[x,0]<1e-20:
                spectrum.spectrum[x,0] = 0.
            realx = width*(x+1)/numpts+start
            outputfile.write(str(realx)+"\t"+str(spectrum.spectrum[x,0]))
            if name == "Raman":
                outputfile.write("\t%f" % spectrum_intensity.spectrum[x,0])
            if x<len(freq): # Write the activities (assumes more pts to plot than freqs - fix this)
                outputfile.write("\t\t"+str(x+1)+"\t"+vibsyms[x]+"\t"+str(freq[x])+"\t"+str(act[x]))
                if name == "Raman":
                    outputfile.write("\t%f" % intensity[x])
                outputfile.write("\t"+str(scale[x])+"\t"+str(json["results"]["freq"]["vibrational_freq"][x]))
            outputfile.write("\n")
        outputfile.close()

#        if root:
#
#            if general==True:
#                title = "%s spectrum scaled by %f" % (name, scalefactor)
#            else:
#                title = "%s spectrum scaled by individual scaling factors" % name
#
#            if name == "IR":
#                g = MPLPlot()
#                g.setlabels("Frequency (cm$^{-1}$)", "IR activity")
#                g.subplot.invert_xaxis()
#                g.subplot.invert_yaxis()
#
#                g.data(zip(spectrum.xvalues,spectrum.spectrum[:,0]), title=title, lines=True)
#                g.subplot.legend(loc="lower right", prop={'size':8})
#                DisplayPlot(root, g, "%s Spectrum" % name)
#
#            if name == "Raman":
#                for type, spec in [("activity", spectrum),
#                                   ("intensity", spectrum_intensity)]:
#                    g = MPLPlot()
#                    g.setlabels("Frequency (cm$^{-1}$)", "Raman %s" % type)
#                    g.data(zip(spec.xvalues,spec.spectrum[:,0]), lines=True, title=title)
#                    g.subplot.legend(loc="upper right", prop={'size':8})
#                    DisplayPlot(root, g, "%s Spectrum" % name)


    ############## START OF MAIN FUNCTION #############

  ##  screen.write("Starting to analyse the vibrational frequencies\n")

    # Create a new output folder if necessary (returns the location of it in any case)
   ## gaussdir=folder(screen,logfilename)

#    gaussdir = "."
    vib_int = []
    vib_raman = []
    try:
        vib_int= json["results"]["freq"]["vibrational_int"]
    except KeyError:
        vib_int = []
    try:
        vib_raman = json["results"]["freq"]["vibrational_raman"]
    except KeyError :
        vib_raman = []
  
    if (len(vib_int)) > 0 :
 #   if (json["results"]["freq"]["vibrational_int"]):
        dofreqs("IR",vib_int)
    else :
        print("vibrational intensities not found.")


    if (len(vib_raman)) > 0 :
  #  if (json["results"]["freq"]["vibrational_raman"]):
        dofreqs("Raman",vib_raman)
    else :
        print("vibrational Raman not found.")

##  screen.write("Finished\n")

#return True

    ############## END OF MAIN FUNCTION #############
    
    
    
