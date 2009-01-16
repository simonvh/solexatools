#!/usr/bin/env python
import sys
import re
import os
from numpy import max
from optparse import OptionParser
from solexatools import peak_stats
from solexatools.track import SimpleTrack

parser = OptionParser()
parser.add_option("-p", "--peakfile", dest="peakfile", help="Peaks in (fixedStep) Wiggle/bed format", metavar="FILE")
parser.add_option("-d", "--datafile", dest="datafile", help="Data in (fixedStep) Wiggle format", metavar="FILE")
parser.add_option("-f", "--format", dest="format", help="Output format: all|number|max|maxfeature|catch", metavar="F", default="all")
parser.add_option("-z", "--zeroes", dest="zeroes", help="Pprint zeroes", action="store_true", default=False)

(options, args) = parser.parse_args()

if not options.peakfile or not options.datafile:
	parser.print_help()
	sys.exit(0)

peakfile = options.peakfile
datafile = options.datafile
format = options.format

peaks = SimpleTrack(peakfile)
data = SimpleTrack(datafile)

formatter = {
	"all": all_formatter,
	"number": number_formatter,
	"max" : max_formatter,
	"maxfeature": maxfeature_formatter,
	"catch": catch_formatter,
}

result = peakstats(peaks, data, formatter[format])

if format == "catch":
	print "Catch format currently not completely working!"

for row in result:
	print row
