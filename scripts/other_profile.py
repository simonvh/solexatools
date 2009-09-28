#!/usr/bin/env python
import sys
from optparse import OptionParser
from solexatools import peak_stats
from solexatools.track import SimpleTrack

parser = OptionParser()
parser.add_option("-p", "--peakfile", dest="peakfile", help="Peaks in (fixedStep) Wiggle/bed format", metavar="FILE")
parser.add_option("-d", "--datafile", dest="datafile", help="Data in (fixedStep) Wiggle format", metavar="FILE")
parser.add_option("-f", "--format", dest="format", help="Output format: all|number|max|mean|maxfeature|catch", metavar="F", default="all")
parser.add_option("-b", "--binsize", dest="binsize", help="binsize", metavar="F", default=500, type=int)
parser.add_option("-w", "--window", dest="window", help="window", metavar="F", default=20000, type=int)

(options, args) = parser.parse_args()

if not options.peakfile or not options.datafile:
	parser.print_help()
	sys.exit(0)

peakfile = options.peakfile
datafile = options.datafile
format = options.format
window = options.window
binsize = options.binsize

peaks = SimpleTrack(peakfile)
data = SimpleTrack(datafile)

formatter = {
	"all": peak_stats.all_formatter,
	"number": peak_stats.number_formatter,
	"max" : peak_stats.max_formatter,
	"mean" : peak_stats.mean_formatter,
	"maxfeature": peak_stats.maxfeature_formatter,
	"catch": peak_stats.catch_formatter,
}

result = peak_stats.binned_peak_stats(peaks, data, window, binsize, formatter[format])

if format == "catch":
	print "Catch format currently not completely working!"

for bin,val in result:
	print "%s\t%0.2e" % (bin,val)
