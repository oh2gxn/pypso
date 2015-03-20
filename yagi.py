# Filename : yagi.py
# Author: Janne Toivola, OH2GXN
# $Id$

import sys
import numpy
import argparse


class Yagi:
    '''An object representing a Yagi-Uda antenna and interfacing NEC2 command
    line program'''

    def __init__(self, dimensions, diam=0.01, cond=2.4938e7):
        '''Constructs an E-element Yagi antenna based on the given dimensions
	   (where E > 1). The parameter is E-by-2, E-by-3, or E-by-4 matrix:
           [pos0,len0,diam0,cond0; pos1,len1,diam1,cond1...], where
	   - pos  :: yagi element position (in meters) on the boom,
           - len  :: element length (in m, symmetric),
	   - diam :: optional element diameter (in m), and
           - cond :: optional element conductivity (in S).
	   The first element is typically the reflector, the second one is
           the driven element. If the matrix contains only 2 columns, the
           given default element diameter is used for all elements.'''

        # NOTE: Conductivity of Al is 2.4938e7 Siemens
        
	# check matrix size, NOTE: does not allow a plain dipole (vector)
        shape = dimensions.shape
        if (len(shape) != 2):
            raise RuntimeError("Invalid matrix of antenna dimensions: %s" %
                               repr(dimensions))
        matrix = numpy.zeros((shape[0],4), numpy.float)
        matrix[:,0:2] = dimensions[:,0:2] # position,length
        if (shape[1] < 4):
            matrix[:,3] = cond # default conductivity
            if (shape[1] < 3):
                matrix[:,2] = diam # default diameter
            else:
                matrix[:,2] = dimensions[:,2]
        else:
            matrix[:,2:4] = dimensions[:,2:4]
        # TODO: eliminate offset (reflector position == 0)
        self.dimensions = matrix


    def setBoom(self, length, diameter, cond=2.4938e7):
        '''Adds the specified boom to the design.'''
	self.boom = numpy.array([length,diameter,cond], numpy.float)


    def setPole(self, height, diameter, cond=2.4938e7):
        '''Adds a supporting vertical pole.'''
       	self.pole = numpy.array([height,diameter,cond], numpy.float)


    def toVector(self):
        '''Represents the free parameters (element len,pos) as a vector.'''
        return self.dimensions[:,0:2].flatten()


    def readNEC(self, stream):
        '''Reads other simulation details, like EX & FR, from file'''
        # TODO


    def fprintNEC(self, stream):
        '''Prints a NEC2 compatible description of the antenna.'''
        stream.write("CM %d-el yagi\n" % self.dimensions.shape[0])
        stream.write("CE\n")
        # TODO


    def evaluate(self, criterion):
        '''Runs NEC2 and evaluates the given design criterion.'''
        # TODO: coupled with the NEC RP card?
        return 0.0



def addParameters(parser):
    '''Command line parameters for stand-alone NEC runs.'''
    parser.add_argument('file', metavar='file.csv',
        help='file with element length,position,[diameter,[conductivity]]')
    # TODO: others


if __name__ == '__main__':
    # parse command line arguments
    note='Run NEC2 for a given Yagi antenna.'
    parser = argparse.ArgumentParser(description=note)
    addParameters(parser)
    args = parser.parse_args()

    # build the antenna
    antenna = yagi( numpy.loadtxt( args.file ))
    # TODO: other parameters

    # print the NEC stuff
    antenna.fprintNEC(sys.stdout)

    # run NEC
    sys.stderr.write("SWR: %f\n" % antenna.evaluate("SWR"))
    
