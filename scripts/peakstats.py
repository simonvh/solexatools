#!/usr/bin/env python
import sys
from optparse import OptionParser
from solexatools import peak_stats
from solexatools.track import SimpleTrack
from os.path import basename,splitext
import os

DEFAULT_BINS = 10
VERSION = "2.1"
usage = "Usage: %prog -p BEDFILE -d DATAFILE [options]\nVersion: " + str(VERSION)

parser = OptionParser(usage=usage, version="%prog " + str(VERSION))
parser.add_option("-p", "--peakfile", dest="peakfile", help="Peaks in BED format", metavar="FILE")
parser.add_option("-d", "--datafile", dest="datafile", help="Data in BAM, BED or WIG format", metavar="FILE")
parser.add_option("-f", "--format", dest="format", help="Output format: all|number|max|mean|sum|maxfeature|length|window", metavar="F", default="all")
parser.add_option("-z", "--zeroes", dest="zeroes", help="Print zeroes", action="store_true", default=False)
parser.add_option("-b", "--bins", dest="bins", help="Number of bins (only when format option 'window' is used)", type="int", default=DEFAULT_BINS)
parser.add_option("--rpkm", dest="rpkm", help="Report RPKM (only with BAM datafile and format option 'window')", default=False, action="store_true")
parser.add_option("--remove_dup", dest="remove_dup", help="Remove duplicates (only with BAM datafile and format option 'window')", default=False, action="store_true")
parser.add_option("--unique", dest="unique", help="Only uniquely mapped reads (only with BAM datafile and format option 'window')", default=False, action="store_true")

(options, args) = parser.parse_args()

if not options.peakfile or not options.datafile:
	parser.print_help()
	sys.exit(1)

peakfile = options.peakfile
datafile = options.datafile
format = options.format

peaks = SimpleTrack(peakfile)

formatter = {
	"all": peak_stats.all_formatter,
	"number": peak_stats.number_formatter,
	"max" : peak_stats.max_formatter,
	"mean" : peak_stats.mean_formatter,
	"maxfeature": peak_stats.maxfeature_formatter,
	"sum": peak_stats.sum_formatter,
	"window": peak_stats.bin_formatter,
	"length": peak_stats.length_formatter,
	"catch": peak_stats.catch_formatter,
}

formatter_options = {"bins":options.bins}


result = []
if datafile.endswith("bam"):
	if not os.path.exists(datafile + ".bai"):
		print "Please provide a sorted and indexed bam file"
		sys.exit(1)
	if options.format =="window":
		result = peak_stats.bam_binned_peak_stats(peaks, datafile, options.bins, options.rpkm, options.remove_dup, options.unique)
	else:
		result = peak_stats.bam_peak_stats(peaks, datafile, formatter[format], formatter_options)

else:
	data = SimpleTrack(datafile)
	result = peak_stats.peak_stats(peaks, data, formatter[format], formatter_options)

if options.format == "window" and options.bins == DEFAULT_BINS:
	sys.stderr.write("Using default of %d bins, optionally specify a different bin number with the -b option\n" % options.bins)

if format == "catch":
	name = splitext(basename(datafile))[0]
	print "## %s" % name
	print "## %s" % name

for row in result:
	if row:
		print row
