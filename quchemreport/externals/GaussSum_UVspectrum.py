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

# This is a modified version of the GaussSum code to generate the UV Spectra


import math
import numpy
from cclib.parser.utils import convertor

colors = "bgrcmyk"

class GaussianSpectrum(object):
    """An optimised version of Spectrum for convoluting gaussian curves.

    Usage:
     GaussianSpectrum(start,end,numpts,peaks,width)
    where
     peaks -- ( [List of peaks],[ [list of heights],[list of heights],..] )
    """
    def __init__(self,start,end,numpts,peaks,width):
        self.start = start
        self.end = end
        self.numpts = numpts
        self.peaks = peaks[0]
        self.heights = peaks[1]
        self.width = width

        # make heights a local variable as it's accessed in the inner loop
        heights = self.heights

        # len(heights) is the number of spectra in this object
        data = []
        self.xvalues = numpy.arange(self.numpts)*float(self.end-self.start)/(self.numpts-1) + self.start
        A = -2.7726/self.width**2
        for x in self.xvalues:
            tot = [0]*len(self.heights) # The total for each spectrum for this x value

            for peakno in range(len(self.peaks)): # For each peak
                pos = self.peaks[peakno]
                exponent = math.exp(A*(pos-x)**2)
                for spectrumno in range(len(heights)):
                    tot[spectrumno] += heights[spectrumno][peakno]*exponent

            data.append(tot)

        self.spectrum = numpy.swapaxes(numpy.array(data),0,1)
    
def ET(logfilename,et_energies,et_oscs,et_rotats,
       start,end,numpts,FWHM):
           
#######################################################
#                 UV-Visible section                  #
#######################################################

    t = GaussianSpectrum(start,end,numpts,
                         ( et_energies,[[x*2.174e8/FWHM for x in et_oscs]] ),
                         FWHM)

    outputfile=open("UVSpectrum.txt","w")
    outputfile.write("Energy (cm-1)\tWavelength (nm)\tAbs\t<--UV Spectrum\tUV-Vis transitions-->\tEnergy (cm-1)\tWavelength (nm)\tOsc. strength\n")
    
    width=end-start
    for x in range(numpts):
        realx=width*x/numpts+start
        outputfile.write(str(realx)+"\t"+str(convertor(realx,"cm-1","nm"))+"\t"+str(t.spectrum[0,x]))
        if x<len(et_energies): # Write the oscillator strengths out also
            outputfile.write("\t\t\t"+str(et_energies[x])+"\t"+str(convertor(et_energies[x],"cm-1","nm")) \
                             +"\t"+str(et_oscs[x]))
        outputfile.write("\n")
    outputfile.close()
        
#        g = MPLPlot()
#        xvalues_nm = [convertor(x,"cm-1","nm") for x in t.xvalues]
#        g.setlabels("Wavelength (nm)", r"$\epsilon$")
#        g.data(zip(xvalues_nm,t.spectrum[0,:]))
#        energies_nm = [convertor(x,"cm-1","nm") for x in logfile.etenergies]
#        oscdata = [(x, y) for x, y in zip(energies_nm, logfile.etoscs) if start < x < end]
#        g.data(oscdata, vlines=True, y2axis="Oscillator strength")
#        
#        g.subplot.set_xlim(left=start, right=end)
        #DisplayPlot(root, g, "UV-Vis Spectrum")

    if len(et_rotats) > 0 :
#######################################################
# Circular dichroism section from here to end of main #
#######################################################

        # Equation 8 in Stephens, Harada, Chirality, 2010, 22, 229.
        # This uses Delta, the half width at 1/e height.
        # We use Sigma, the full width at 1/e height.
        #   Delta = Sigma / 2

        endwaveno = convertor(start,"nm","cm-1")
        startwaveno = convertor(end,"nm","cm-1")
        sigma = convertor(FWHM, "eV", "cm-1")
        Delta = sigma / 2.
        
        prefactor = 1.0 / (2.296e-39 * math.sqrt(math.pi) * Delta)
        peakmax = []
        for i in range(len(et_rotats)):
            peakmax.append(prefactor * et_rotats[i] *
                           et_energies[i] * 1e-40)
        # FWHM is sqrt(ln2) times sigma
        real_FWHM = math.sqrt(math.log(2)) * sigma
        t = GaussianSpectrum(startwaveno,endwaveno,numpts,
                             ( et_energies,[peakmax] ),
                             real_FWHM)
        
                # Write CDSpectrum.txt containing info on the CD spectrum
  
        outputfile=open("CDSpectrum.txt","w")
        outputfile.write("Energy (cm-1)\tWidth (nm)\tAbs\t<--CD Spectrum\tStates-->\tWavelength (nm)\tEnergy (cm-1)\tR(length)\n")
        
        
#        def arrayData():
#            if os.path.isfile("CDSpectrum.txt"):
#                
#                                                return [t] 
        
        width = endwaveno - startwaveno
        for x in range(numpts):
            # Need to use float, otherwise it's an integer
            # realx = float(width)*x/numpts+start
            realx = width * x / numpts + startwaveno
            outputfile.write("%f\t%f\t%f" % (realx, 1.0e7/realx, t.spectrum[0,x]))
            if x < len(et_energies): # Write the R values out also
                outputfile.write( "\t\t\t%f\t%f\t%f" % (convertor(et_energies[x], "cm-1", "nm"),
                                                        et_energies[x],et_rotats[x]) )
            outputfile.write("\n")
        outputfile.close()
        
#        g = MPLPlot()
#        g.setlabels("Wavelength (nm)", r"$\epsilon$")
#        
#        xvalues_nm = [convertor(x,"cm-1","nm") for x in t.xvalues]
#        g.data(zip(xvalues_nm,t.spectrum[0,:]))
#        energies_nm = [convertor(x,"cm-1","nm") for x in logfile.etenergies]
#        oscdata = [(x, y) for x, y in zip(energies_nm, logfile.etrotats) if start < x < end]
#        g.data(oscdata, vlines=True, y2axis="R (length) / $10^{-40}$")
#        
#        g.subplot.set_xlim(left=start, right=end)
#        g.secondaxis.set_ylim(bottom=0)
    
        #######################################################
        #                      End of main                    #
        ####################################################### 
        
        
        
