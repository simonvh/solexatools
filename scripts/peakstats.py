#!/usr/bin/env python
import sys
from optparse import OptionParser
from solexatools import peak_stats
from solexatools.track import SimpleTrack
from os.path import basename,splitext

parser = OptionParser()
parser.add_option("-p", "--peakfile", dest="peakfile", help="Peaks in (fixedStep) Wiggle/bed format", metavar="FILE")
parser.add_option("-d", "--datafile", dest="datafile", help="Data in (fixedStep) Wiggle format", metavar="FILE")
parser.add_option("-f", "--format", dest="format", help="Output format: all|number|max|mean|maxfeature|catch", metavar="F", default="all")
parser.add_option("-z", "--zeroes", dest="zeroes", help="Print zeroes", action="store_true", default=False)

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
	"all": peak_stats.all_formatter,
	"number": peak_stats.number_formatter,
	"max" : peak_stats.max_formatter,
	"mean" : peak_stats.mean_formatter,
	"maxfeature": peak_stats.maxfeature_formatter,
	"window": peak_stats.bin_formatter,
	"catch": peak_stats.catch_formatter,
}

options = {"bins":10}

result = peak_stats.peak_stats(peaks, data, formatter[format], options)

if format == "catch":
	name = splitext(basename(datafile))[0]
	print "## %s" % name
	print "## %s" % name

for row in result:
	print row
