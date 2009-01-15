#!/usr/bin/env python
import sys
import re
import os
from numpy import max
from optparse import OptionParser
from SolexaTools.solexatools import SimpleTrack


parser = OptionParser()
parser.add_option("-p", "--peakfile", dest="peakfile", help="Peaks in (fixedStep) Wiggle/bed format", metavar="FILE")
parser.add_option("-d", "--datafile", dest="datafile", help="Data in (fixedStep) Wiggle format", metavar="FILE")
parser.add_option("-f", "--format", dest="format", help="Output format: all|number|max|catch", metavar="F", default="all")
parser.add_option("-z", "--zeroes", dest="zeroes", help="Pprint zeroes", action="store_true", default=False)

(options, args) = parser.parse_args()

if not options.peakfile or not options.datafile:
	parser.print_help()
	sys.exit(0)

peakfile = options.peakfile
datafile = options.datafile
format = options.format
catch_spacing = 10

peaks = SimpleTrack(peakfile)
data = SimpleTrack(datafile)

peak_feature = peaks.get_next_feature()
data_feature = None

if format == "catch":
	bla = datafile.replace(".wig", "")
	print "## %s" % bla

prev_data_seq = ""
while peak_feature:
	overlap = []

	while data_feature and peak_feature and (data_feature[0] > peak_feature[0]):
		if format == "catch":
			print "# %s:%s-%s" % peak_feature[0:3]
			print "%s\t%s\t%s\t%s" % (peak_feature[0], peak_feature[1], peak_feature[1] + catch_spacing, 0)		
			print "%s\t%s\t%s\t%s" % (peak_feature[0], peak_feature[2] - catch_spacing, peak_feature[2], 0)	

		peak_feature = peaks.get_next_feature()
	
	data_feature = data.get_next_feature()
	
	if peak_feature:
		while (data_feature and ((data_feature[0] < peak_feature[0]) or ((data_feature[0] == peak_feature[0]) and (data_feature[2] < peak_feature[1])))):
			data_feature = data.get_next_feature()

		while (data_feature and (data_feature[1] <= peak_feature[2] and data_feature[0] == peak_feature[0])):
			overlap.append(data_feature)
			data_feature = data.get_next_feature()

		#for row in fraction(overlap, 0.5):
		#	print "\t".join(map(str,row))
	
	
		if format == "number":
			if len(overlap) > 1:
				# Number
				print "\t".join(map(str, list(peak_feature[0:3]) +  [len(overlap)]))
			elif options.zeroes:
				print "\t".join(map(str, list(peak_feature[0:3]) +  [0]))
		elif format == "max":
			if len(overlap) > 1:
				# Maximum value
				print "\t".join(map(str, list(peak_feature[0:3]) +  [max_val(overlap)]))
			elif options.zeroes:
				print "\t".join(map(str, list(peak_feature[0:3]) +  [0]))

		elif format == "maxfeature":
			if len(overlap) > 1:
				print "\t".join(map(str, list(max_feature(overlap))))
		elif format == "all" or format == "catch":
			if format == "catch":
				print "# %s:%s-%s" % peak_feature[0:3]
			if len(overlap) < 2:
				print "%s\t%s\t%s\t%s" % (peak_feature[0], peak_feature[1], peak_feature[1] + catch_spacing, 0)		
				print "%s\t%s\t%s\t%s" % (peak_feature[0], peak_feature[2] - catch_spacing, peak_feature[2], 0)	
			else:
				if overlap[0][0] < peak_feature[1] + catch_spacing:
					print "%s\t%s\t%s\t%s" % (peak_feature[0], peak_feature[1], peak_feature[1] + catch_spacing, 0)		
				last = peak_feature[1]
				for f in overlap:
					if f[1] > last + catch_spacing:
						print "%s\t%s\t%s\t%s" % (peak_feature[0], last, last + catch_spacing, 0)		
					if f[1] >= last + catch_spacing * 3:
						print "%s\t%s\t%s\t%s" % (peak_feature[0], f[1] - catch_spacing, f[1], 0)		
				
					print "\t".join(map(str, f))
					last = f[2]
			
				if overlap[-1][2] <= peak_feature[2] - catch_spacing * 3:
					print "%s\t%s\t%s\t%s" % (peak_feature[0], overlap[-1][2],overlap[-1][2] + catch_spacing, 0)		
			
				if overlap[-1][2] < peak_feature[2]:
					print "%s\t%s\t%s\t%s" % (peak_feature[0], peak_feature[2] - catch_spacing, peak_feature[2], 0)		

	peak_feature = peaks.get_next_feature()

