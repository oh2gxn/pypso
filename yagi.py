#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# @file
# @brief Code for evaluating Yagi-Uda antennas with NEC2
# @author Janne Toivola, OH2GXN, 2015

import os
import sys
import numpy
import argparse
import re
import tempfile
import subprocess


class Yagi:
    """Represents a Yagi-Uda antenna and interfaces the NEC2 simulator.
    NOTE: There are three different representations involved:
    - CSV file:   absolute, non-negative element positions and lengths
    - PSO vector: incremental element positions, and absolute lengths
    - NEC file:   absolute real coordinates of element ends (and much more)
    """

    def __init__(self, dimensions, diam=0.01, cond=2.4938e7):
        """Constructs a Yagi antenna based on the given dimensions.
        Typically, these are directly from build instructions.

        Arguments:
        dimensions -- a numpy matrix with a row for each element,
        and 2 to 4 columns:
	- element position on the boom (in meters),
        - element length (in m, centered / symmetric),
	- optional element diameter (in m, default 10mm), and
        - optional element conductivity (in S, default Al).

	The first element is typically a reflector, the second one is
        the driven element: at least two elements are required.

        diam -- default element diameter, if dimensions.shape[1]<3

        cond -- default element conductivity, if dimensions.shape[1]<4
        """

        # NOTE: Conductivity of Al is 2.4938e7 Siemens
        
        # Check matrix size
        # NOTE: Does not allow a plain dipole (vector),
        #       since the reflector is required
        shape = dimensions.shape
        if (len(shape) != 2):
            raise RuntimeError("Invalid matrix of antenna dimensions: %s" %
                               repr(dimensions))

        # Copy values
        matrix = numpy.zeros((shape[0],4), numpy.float)
        if (shape[1] < 4):
            matrix[:,3] = cond # default conductivity
            if (shape[1] < 3):
                matrix[:,2] = diam # default diameter
            else:
                matrix[:,2] = dimensions[:,2]
        else:
            matrix[:,2:4] = dimensions[:,2:4]
        matrix[:,1] = dimensions[:,1] # lengths
        
        # Eliminate offset
        # NOTE: internally, reflector position == 0,
        # but NEC files etc. have center of gravity (CoG) considered
        Xmin = numpy.min(dimensions[:,0])
        matrix[:,0] = dimensions[:,0] - Xmin
        
        self.dimensions = matrix

        # FIXME: None of the below can be loaded from a CSV file matrix?
        # These are "additional details" that don't fit the neat PSO world.
        # Have them packed on a single row at the beginning of CSV???
        self.gap = 0.0 # meters between element center and boom center
        self.boom = numpy.zeros(3, numpy.float) # len, diam, conductivity
        self.pole = numpy.zeros(3, numpy.float) # len, diam, conductivity
        self.radiator = 2 # driven by the second element (1st is reflector)
        self.frequencies = numpy.array([51.510, 0.020, 51.510]) # 6m call freq.


    @staticmethod
    def loadCSV(filename):
        return Yagi( numpy.loadtxt( filename, delimiter=',' ))
        
    def saveCSV(self, filename):
        numpy.savetxt( filename, self.dimensions, delimiter=',' )
        

    def copy(self):
        """Creates a full copy of this object."""
        cp = Yagi(self.dimensions)
        cp.gap = self.gap
        cp.boom = numpy.copy(self.boom)
        cp.pole = numpy.copy(self.pole)
        cp.directions  = numpy.copy(self.directions)
        cp.frequencies = numpy.copy(self.frequencies)
        return cp


    def setBoom(self, length, diameter, gap=0.01, cond=2.4938e7):
        """Adds the specified boom to the design.

        Arguments:        
        length   -- how far the forward end is from the reflector (in m)
        diameter -- pipe diameter (in m), zero == no boom in simulation
        gap  -- insulation between the boom and elements (in m)
        cond -- conductivity of the material (in S)
        """
        self.gap = gap # insulation between elements and boom
        self.boom = numpy.array([length,diameter,cond], numpy.float)
        # TODO: offset etc.?

    def updateBoom(self, diameter=None, cond=None, gap=None):
        """Adds or updates a boom of suitable length to the design

        Arguments:
        diameter -- pipe diameter (in m), 0 == no boom, None == previous
        cond -- conductivity (in S), None == previous
        gap -- center-to-center distance between the elements and boom (in m)
        """
        if diameter is None:
            d = self.boom[1]
        else:
            d = diameter
        if cond is None:
            c = self.boom[2]
        else:
            c = cond
        if gap is None:
            g = self.gap
        else:
            g = gap
        # automatic length:
        length = (numpy.max(self.dimensions[:,0]) -
                  numpy.min(self.dimensions[:,0]))
        if d > 0.0:
            self.gap  = g
            self.boom = numpy.array([length,d,c], numpy.float)
        else:
            self.gap  = 0.0
            self.boom = numpy.zeros(3, numpy.float)

    def setPole(self, height, diameter, cond=2.4938e7):
        """Adds a supporting vertical pole.

        Arguments:
        height -- level of yagi boom from the ground in m, 0 => boom at ground
        diameter -- pipe diameter in m, 0 => boom at ground
        cond -- conductivity in S
        """
        # TODO: a typical yagi is not mounted at the very top
        if (height > 0.0 and diameter > 0.0):
            self.pole = numpy.array([height,diameter,cond], numpy.float)
        else:
            self.pole = numpy.zeros(3, numpy.float)

    def setRadiator(self, index):
        """Selects which one of E elements acts as the driven element.
        Arguments:
        index -- an integer in [1,E], typically 2 (since 1st is reflector)
        """
        if 0 < index and index <= self.dimensions.shape[0]:
            self.radiator = index

    def setFrequencies(self, frlist):
        """Sets the frequency range used for NEC. It has to be uniform and
        contiguous.

        Arguments:
        frlist -- list of 1 to 3 values: minimum, [[increment], maximum]
        """
        if frlist is None:
            return
        minimum = frlist[0]
        if len(frlist) > 1:
            if len(frlist) > 2:
                increment = frlist[1]
                maximum   = frlist[2]
            else:
                maximum   = frlist[1]
                increment = maximum - minimum
        else:
            increment = 0.020
            maximum   = minimum
        self.frequencies = numpy.array([minimum, increment, maximum])


    def toVector(self):
        """Represents the length & position of each element as a vector.
        Used for PSO or other algorithms updating the free parameters.
        Positions are incremental, so any unsorted vector of non-negative
        floats is valid.
        """
        E = self.dimensions.shape[0]
        matrix = numpy.copy(self.dimensions[:,0:2])
        matrix[1:E,0] = numpy.diff(matrix[:,0]) # increments, not absolute
        # TODO: maybe (len0, pos0, len1, pos1,...) is not good?
        return matrix.flatten()

    def fromVector(self, dimensions):
        """Creates a new yagi object according to a given vector.
        Useful for initializing PSO with a set of randomized yagis.
        """
        newYagi = self.copy()
        newYagi.updateFromVector(dimensions)
        return newYagi

    def updateFromVector(self, newDimensions):
        """Set the free parameters according to updated values.
        NOTE: a setter method – does not create a new yagi object.
        """
        E = self.dimensions.shape[0]
        matrix = newDimensions.reshape((E,2))
        self.dimensions[:,0] = numpy.cumsum(matrix[:,0]) # absolute, not incr.
        self.dimensions[:,1] = matrix[:,1]
        self.updateBoom()
        

    def fprintNEC(self, stream):
        """Prints a NEC2 compatible description of the antenna.
        NOTE: The issue is that NEC2 considers also the input signal.
        """
        
        # Comment lines and Comment End
        W = self.dimensions # shorter name for wires
        E = W.shape[0]
        stream.write("CM %d-el yagi: " % E)
        stream.write("NEC file generated by yagi.py by OH2GXN\n")
        stream.write("CE\n")
        
        # Geometry
        S = 11 # number of wire segments in simulation, TODO: dynamic?
        Z = 1.0 # default height, NOTE: elements need Z > 0
        if self.pole[1] > 0.0:
            Z = self.pole[0] # top of the mast
        if self.boom[1] > 0.0:
            # approx. CoG when boom diameter >> element diameter
            # TODO: better approximation?
            CoG = 0.5*self.boom[0]
            Zelem = Z + 0.5*self.boom[1] + self.gap + 0.5*numpy.max(W[:,2])
        else:
            # approx. CoG with uniform element material
            # TODO: better approximation?
            CoG = numpy.dot(W[:,0],W[:,1])/numpy.sum(W[:,1])
            Zelem = Z
                    
        # Elements: tags [1:E]
        for e in range(0,E):
            # NOTE: element centers at the same plane Z=Zelem
            # Wire Tag Segs X1 Y1 Z1 X2 Y2 Z2 Radius
            stream.write("GW %d %d %g %g %g %g %g %g %g\n" %
                         ((e+1), S,
                          W[e,0]-CoG, -0.5*W[e,1], Zelem,
                          W[e,0]-CoG,  0.5*W[e,1], Zelem,
                          W[e,2] * 0.5))
        # Boom at Y = 0.0, tag = E+1
        if self.boom[1] > 0.0:
            stream.write("GW %d %d %g 0.0 %g %g 0.0 %g %g\n" %
                         ((E+1), S,
                         -CoG, Z,
                         self.boom[0]-CoG, Z,
                         self.boom[1] * 0.5))
        # Mast at (X,Y) = (0.0,0.0)
        if self.pole[1] > 0.0:
            stream.write("GW %d %d 0.0 0.0 0.0 0.0 0.0 %g %g\n" %
                         ((E+2), S, self.pole[0], self.pole[1] * 0.5))
        stream.write("GE 1\n") # Geometry End

        # Ground:
        # - type -1: free space
        # - type 0: finite ground with a reflection coefficient approximation
        # - type 1: perfectly conducting ground
        # - type 2: finite ground using the Sommerfeld-Norton method
        gntype = -1
        dielectric   = 13    # TODO: is this a reasonable value?
        conductivity = 0.005 # Siemens
        stream.write("GN %d 0 0 0 %g %g 0 0 0 0\n" %
                     (gntype, dielectric, conductivity))

        # Excitation with voltage source (driven element = tag 2)
        # FIXME: YU7EF 7L6+8L4 dual band yagi is driven by the 4th element!
        V = (10.0, 0.0) # Volts, no phase considerations, TODO: free parameter
        stream.write("EX 0 %d %d 0 %g %g\n" % (self.radiator, (S//2)+1, 
                                               V[0], V[1]))

        # Wire losses
        for e in range(0,E):
            # Loading Type Tag StartSeg EndSeg Conductivity
            stream.write("LD 5 %d 1 %d %g\n" % ((e+1), S, W[e,3]))
        # NOTE: possibility for trap coils etc. with Type 4 (R,L,C not Cond)
        if self.boom[1] > 0.0:
            stream.write("LD 5 %d 1 %d %g\n" % ((E+1), S, self.boom[2]))
        if self.pole[1] > 0.0:
            stream.write("LD 5 %d 1 %d %g\n" % ((E+2), S, self.pole[2]))
        
        # TODO: Transmission Line: TL Tag1 Seg1 Tag2 Seg2 Z0 Len/VF 0 0 0 0
        # where the mechanical cable Length is divided with the Velocity Factor
        # e.g. (abs(W[2,0]-CoG) + Z) / 0.66, but needs to be excited from
        # an additional wire element Tag1 at the bottom of the mast...
        # and the match+balun for the driven element???

        # Frequencies
        # TODO: adjustable frequencies according to evaluation criteria
        # IARU Region 1 Band Plan for 6 meters:
        # - SSB: 50.100 - 50.200 MHz
        # - FM simplex: 51.410 - 51.590 MHz, 51.510 MHz calling freq.
        #frange = numpy.array([51.410, 0.020, 51.590]) # [min, inc, max] MHz
        #frange = numpy.array([51.510, 0.020, 51.510]) # [min, inc, max] MHz
        frange = self.frequencies
        F = int((frange[2] - frange[0])/frange[1] + 1) # increment > 0 !
        stream.write("FR 0 %d 0 0 %g %g\n" % (F, frange[0], frange[1]))

        # Report: which directions included in the simulation
        # FIXME: bigger increments, 1.0 deg => 38MB of text @ output
        # TODO: adjustable report according to criteria
        azimuth  = numpy.array([0.0, 1.0, 360.0]) # [min, inc, max] deg
        Rz = int((azimuth[2] - azimuth[0])/azimuth[1]) # +1?
        altitude = numpy.array([0.0, 1.0, 90.0]) # NOTE: deg from zenith!
        Ry = int((altitude[2] - altitude[0])/altitude[1] + 1) # +1
        stream.write("RP 0 %d %d 1000 %g %g %g %g\n" %
                     (Ry, Rz, altitude[0], azimuth[0], altitude[1], azimuth[1]))

        stream.write("EN\n") # The End


    def evaluate(self, criterion):
        """Runs NEC2 and evaluates the given design criterion."""
        # TODO: criterion coupled with the NEC FR and RP cards?
        # Possible hack:
        # - collect all frequency range expressions of type "[a:i:b,c:j:d]"
        # - for each in the union of all frequencies: run NEC, collect the stats
        # - evaluate the criterion without further sanitization?

        # TODO: proper FR card and other details!
        
        # Run nec2c, which is a translation of the original NEC2
        # NOTE: it uses ordinary input & output files ALWAYS,
        # unless hacked to use stdout for better composability?
        tmpnec = tempfile.NamedTemporaryFile(suffix='.nec')
        basename = os.path.splitext(tmpnec.name)[0]
        if len(basename) < 2:
            basename = tmpnec.name # splitext failed, just use X.nec.out
        tmpout = '%s.out' % basename
        necfile = open(tmpnec.name, 'w')
        self.fprintNEC( necfile )
        #self.fprintNEC( tmpnec.file ) # NOTE: tmpnec.file does not exist
        necfile.close()
        try:
            sys.stderr.write("Running: nec2c -i %s -o %s\n" %
                             (tmpnec.name,tmpout))
            stdout = subprocess.check_output([ "nec2c", 
                                               "-i", tmpnec.name,
                                               "-o", tmpout ])
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Error: Failed to run nec2c.\n") # print e too?
            stdout = ""
        except OSError as e:
            sys.stderr.write("Error: Could not find nec2c.\n")
            stdout = ""            
        #sys.stderr.write("Output: stdout='%s';\n" % stdout)
            
        try:
            with open(tmpout, 'r') as results:
                for line in results:
                    if '=' in line:
                        sys.stdout.write("NEC2: %s" % line) # DEBUG
                        
        # Example output from NEC2C, for each FR:
        # ...
        # ---------- POWER BUDGET ---------
        # INPUT POWER   =  5.4174E-02 Watts
        # RADIATED POWER=  4.8132E-02 Watts
        # STRUCTURE LOSS=  6.0415E-03 Watts
        # NETWORK LOSS  =  0.0000E+00 Watts
        # EFFICIENCY    =   88.85 Percent
        # ...

        # TODO: extract the relevant expressions into suitable arrays, like
        #       efficiency[51410:20:51590], maybe sparse or lazy ???
        except IOError:
            sys.stderr.write("Error: Could not read NEC results.\n")

        tmpnec.close() # or .file.close() ?
        try:
            os.remove(tmpout)
        except (OSError, IOError):
            pass

        # just some test expression:
        weight = numpy.sum( numpy.multiply( self.dimensions[:,0], 
                                            self.dimensions[:,1] ))
        return eval(criterion)


#####
# For using this stuff "stand-alone" from command line:
if __name__ == '__main__':
    # Parse command line arguments
    note='Run NEC2 for a given Yagi antenna.'
    parser = argparse.ArgumentParser(description=note)
    parser.add_argument('file', metavar='elements.csv',
        help='file with element length,position,[diameter,[conductivity]]')
    parser.add_argument('-b', '--boom', metavar='length,diameter',
        help='optional boom length and diameter')
    parser.add_argument('-f', '--frequency', metavar='min:inc:max',
        help='optional operating frequency or frequency range in MHz')
    parser.add_argument('-g', '--gap', metavar='distance',
        help='optional insulating gap between element center and boom center')
    parser.add_argument('-p', '--pole', metavar='height,diameter',
        help='optional antenna mast height and diameter')
    parser.add_argument('-r', '--radiator', metavar='index',
        help='optional driven element index (1...E, default=2)')
    parser.add_argument('-o', '--output', metavar='yagi.nec',
        help='optional NEC output file')
    # TODO: others
    args = parser.parse_args()

    # Build the antenna
    antenna = Yagi.loadCSV( args.file )
    gap = 0.0
    if args.gap is not None:
        gap = float(args.gap)
    if args.boom is not None:
        length,diameter = map(float, args.boom.split(','))
        antenna.setBoom(length, diameter, gap)
        # TODO: optional length
    if args.frequency is not None:
        frlist = map(float, args.frequency.split(':'))
        antenna.setFrequencies(frlist)
    if args.pole is not None:
        length,diameter = map(float, args.pole.split(','))
        antenna.setPole(length, diameter)
    if args.radiator is not None:
        antenna.setRadiator(int(args.radiator))
    # TODO: other parameters?

    # Print the NEC stuff
    fid = sys.stdout
    if args.output is not None:
        fid = open(args.output, 'w') # TODO: 'wb'?
        # TODO: temp file, then run NEC2
    antenna.fprintNEC(fid)
    fid.close()

    # TODO: run NEC
    #sys.stderr.write("SWR: %f\n" % antenna.evaluate("SWR"))
    
