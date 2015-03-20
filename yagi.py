# Filename : yagi.py
# Author: Janne Toivola, OH2GXN
# $Id$

import sys
import numpy
import argparse


class Yagi:
    '''An object representing a Yagi-Uda antenna and interfacing NEC2 command line program'''

    def __init__(self, matrix, diam=10.0):
        '''Constructs an E-element Yagi antenna based on the given dimensions (where E > 1).
	   Matrix has E-by-2 or E-by-3 elements: [pos0, len0, diam0; pos1, len1, diam1...], where 
	   - pos  :: yagi element position (in mm) on the boom,
           - len  :: yagi element half length (in mm, symmetric),
	   - diam :: optional yagi element diameter (in mm).
	   The first element is typically the reflector, the second one is the driven element.
           If the matrix contains only 2 columns, the given default element diameter is used for all elements.'''

	# TODO: check matrix size?

        self.dimensions = matrix


    def setBoom(self, length, diameter):
        '''Adds the specified boom to the design.'''
	self.boom = numpy.array([length,diameter])


    def setHeight(self, height):
        '''Adds a supporting vertical pole and conductive ground.'''
        # TODO

    def read

    def fprintNEC(self, fid):
        '''Prints a NEC2 compatible description of the antenna.'''
        # TODO


    def evaluate(self, criterion):
        '''Runs NEC2 and evaluates the given design criterion.'''
        # TODO
        return 0.0



def addParameters(parser):
    '''Command line parameters for stand-alone NEC runs.'''
    parser.add_argument('positions', metavar='pos0,pos1,...',
                        help='element positions')
    parser.add_argument('lengths', metavar='len0,len1,...',
                        help='element lengths')
    parser.add_argument('diameters', metavar='diam0,diam1,...',
                        help='element diameters')
    # TODO: others


if __name__ == '__main__':
    sys.stderr.write('NOTE: This code was not meant to be run stand-alone.\n')
    parser = argparse.ArgumentParser(description='Run NEC2 for a given Yagi antenna.')
    addParameters(parser)
    args = parser.parse_args()
    # TODO: ant = Yagi(args)
