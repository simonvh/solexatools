#!/usr/bin/env python
import sys
from datetime import datetime
from optparse import OptionParser
from solexatools.assign_peaks import *
from solexatools.track import *

parser = OptionParser()
parser.add_option("-p", "--peakfile", dest="peakfile", help="Peaks in Wiggle/bed format", metavar="FILE")
parser.add_option("-t", "--targetfile", dest="targetfile", help="Targets in 6-column bed format (including strand!!)", metavar="FILE")
parser.add_option("-e", "--extend", dest="extend", help="Extend targets by N bases", metavar="N", default=5000)
parser.add_option("-u", "--upstream", dest="upstream", help="Maximum upstream length", metavar="N", default=100000)
parser.add_option("-o", "--overview", action="store_true", dest="overview", help="Printe overview", default=False)
(options, args) = parser.parse_args()

genefile = options.targetfile
bedfile = options.peakfile
extend = int(options.extend)
upstream = int(options.upstream)
overview = options.overview

if not genefile or not bedfile:
	parser.print_help()
	sys.exit(0)

print "#assign_peaks.py %s" % (datetime.today().strftime("%d-%m-%Y %H:%M"))
print "#peaks: %s" % bedfile
print "#targets: %s" % genefile
print "#extension: %s maximum upstream: %s" % (extend, upstream) 

matrix = assign_peaks(genefile, bedfile, upstream, extend)
if overview:
	names = ["prox_upstream","downstream","in_gene","upstream"]
	str = ""
	for name in names:
		str += "\t%s" % name
	print str
	for gene,targets in matrix.items():
		cat = {"prox_upstream":[], "downstream":[], "in_gene":[], "upstream":[]}
		for (chr, middle, target, val) in targets:
			if val < 0:
				cat["downstream"].append(target)
			elif val == 0:
				cat["in_gene"].append(target)
			elif val <= extend:
				cat["prox_upstream"].append(target)
			elif val > extend:
				cat["upstream"].append(target)
			else:
				raise ShouldNotHappenError

		str = gene
		for name in names:
			str += "\t"  + (",".join(cat[name]))
		print str
else:
	for gene, target in matrix.items():
		for (chr, middle, target, val) in target:
			print "%s\t%s\t%s\t%s\t%s" % (gene, chr, middle, target, val)
