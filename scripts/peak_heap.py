#!/usr/bin/env python
#
# peak_heap.py
#
# Author: Simon van Heeringen <s.vanheeringen@ncmls.ru.nl
#
import sys
import os
from optparse import OptionParser
from solexatools.peak_heap import *

EXTS=[".bed", ".bam"]

parser = OptionParser()
parser.add_option("-p", "--peakfile(s)", dest="peakfile", help="Peaks in wiggle/bed format, seperated by commas", metavar="FILES")
parser.add_option("-d", "--datafile(s)", dest="datafile", help="Datafiles in bed format, seperated by commas", metavar="FILES")
parser.add_option("-m", "--nomerge", dest="merge", action="store_false", help="Don't merge overlapping peaks", default=True)

(options, args) = parser.parse_args()

if not options.peakfile or not options.datafile:
	parser.print_help()
	sys.exit(0)

peakfiles = [x.strip() for x in options.peakfile.split(",")]

datafiles = {}
for fname in options.datafile.split(","):
	root,ext = os.path.splitext(fname)
	if not ext in EXTS:
		print "Only bed or bam files are supported!"
		sys.exit(1)
	else:
		datafiles[os.path.basename(root)] = fname

result = peak_heap(peakfiles, datafiles, options.merge)
samples = sorted(datafiles.keys())
header = "peak\t" + "\t".join(samples)
print header
for peak,counts in result.items():
	s = peak
	for sample in samples:
		s += "\t%s" % counts[sample]
	print s


