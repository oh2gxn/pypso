#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# @file
# @brief Code for evaluating Yagi-Uda antennas with NEC2
# @author Janne Toivola, OH2GXN, 2015

import sys
import numpy
import argparse


class Yagi:
    """Represents a Yagi-Uda antenna and interfaces NEC2 program."""

    def __init__(self, dimensions, diam=0.01, cond=2.4938e7):
        """Constructs a Yagi antenna based on the given dimensions.

        Arguments:
        dimensions -- a numpy matrix with a row for each element,
        and 2 to 4 columns:
	    - element position (in meters) on the boom,
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
        self.gap = 0.0
        self.boom = numpy.zeros(3, numpy.float)
        self.pole = numpy.zeros(3, numpy.float)
        self.directions = None  # TODO
        self.frequencies = None # TODO

    # TODO: __init__(filename) for loading from file..?
        

    def copy(self):
        """Creates a full copy of this object."""
        cp = yagi(self.dimensions)
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
        if gap is not None:
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



    def toVector(self):
        """Represents the length & position of each element as a vector.
        Used for PSO or other algorithms updating the free parameters.
        """
        # TODO: maybe (len0, pos0, len1, pos1,...) is not good?
        return self.dimensions[:,0:2].flatten()

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
        self.dimensions[:,0:2] = newDimensions.reshape((E,2))
        

    def fprintNEC(self, stream):
        """Prints a NEC2 compatible description of the antenna."""
        
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
        V = (10.0, 0.0) # Volts, no phase considerations, TODO: free parameter
        stream.write("EX 0 %d %d 0 %g %g\n" % (2, (S//2)+1, V[0], V[1]))

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
        # IARU Region 1 Band Plan:
        # - SSB: 50.100 - 50.200 MHz
        # - FM simplex: 51.410 - 51.590 MHz, 51.510 MHz calling freq.
        #frange = numpy.array([51.410, 0.020, 51.590]) # [min, inc, max] MHz
        frange = numpy.array([51.510, 0.020, 51.510]) # [min, inc, max] MHz
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

        # Example output from NEC2C, for each FR:
        # ...
        # ---------- POWER BUDGET ---------
        # INPUT POWER   =  5.4174E-02 Watts
        # RADIATED POWER=  4.8132E-02 Watts
        # STRUCTURE LOSS=  6.0415E-03 Watts
        # NETWORK LOSS  =  0.0000E+00 Watts
        # EFFICIENCY    =   88.85 Percent
        # ...
        
        return 0.0


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
    parser.add_argument('-p', '--pole', metavar='height,diameter',
        help='optional antenna mast height and diameter')
    parser.add_argument('-o', '--output', metavar='yagi.nec',
        help='optional NEC output file')
    # TODO: others
    args = parser.parse_args()

    # Build the antenna
    antenna = Yagi( numpy.loadtxt( args.file, delimiter=',' ))
    if args.boom is not None:
        len,dia = map(float, args.boom.split(','))
        antenna.setBoom(len,dia)
        # TODO: optional length
    if args.pole is not None:
        len,dia = map(float, args.pole.split(','))
        antenna.setPole(len,dia)    
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
    
