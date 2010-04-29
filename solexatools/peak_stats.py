#!/usr/bin/env python
import sys
import re
import os
from numpy import max,mean,sum
from solexatools.track import SimpleTrack

def max_val(features):
	test = map(lambda x: x[3], features)
	return max(test)

def mean_val(features):
	test = map(lambda x: x[3], features)
	return mean(test)

def max_feature(features):
	features.sort(cmp=lambda x,y: cmp(x[3],y[3]))
	return features[-1]

def fraction(features, f):
	m = max_val(features)
	return filter(lambda x:x[3] >= f * m, features)

def max_formatter(peak, overlap, options={}):
	if overlap:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], max_val(overlap))
	else:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], 0)

def sum_formatter(peak, overlap, options={}):
	if overlap:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], sum([x[3] for x in overlap]))
	else:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], 0)

def mean_formatter(peak, overlap, options={}):
	if overlap:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], mean_val(overlap))
	else:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], 0)

def maxfeature_formatter(peak, overlap, options={}):
	if overlap:
		return "\t".join([str(x) for x in max_feature(overlap)])
	else:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], 0)
		
def number_formatter(peak, overlap, options={}):
    return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], len(overlap))

def tuple_number_formatter(peak, overlap, options={}):
    return (peak[0], peak[1], peak[2], len(overlap))

def all_formatter(peak, overlap, options={}):
	return "\n".join("%s\t%s\t%s\t%s" % x[:4] for x in overlap)

def tuple_all_formatter(peak, overlap, options={}):
	if overlap:
		return ["%s\t%s\t%s\t%s" % x[:4] for x in overlap]
	else:
		return []

def peak_and_tuple_all_formatter(peak, overlap, options={}):
	if overlap:
		return peak, ["%s\t%s\t%s\t%s" % x[:4] for x in overlap]
	else:
		return peak, []

def bin_formatter(peak, overlap, options={"bins":10}):
	nr_bins = options["bins"]
	l = (peak[2] - peak[1])/ float(nr_bins)
	bins = [0] * nr_bins
	

	for feature in overlap:
		
		#print "feature", feature
		m = feature[2]
		if feature[2] > peak[2]:
			m = peak[2]

		for i in range(int((feature[1] - peak[1])  / l ),int((m - peak[1] - 1) / l + 1)):
			if i >= 0:
				#print "BIN %i +1," % i, bins[99]
				bins[i] += 1	
	if len(peak) >= 5 and peak[4] and peak[4] == "-":
		bins = bins[::-1]
		
	return "%s\t%s\t%s\t" % (peak[0], peak[1], peak[2]) + "\t".join([str(x) for x in bins])

def dist_to_center_formatter(peak, overlap, options={}):
	dist = []
	center = peak[1] + (peak[2] - peak[1]) / 2
	for feature in overlap:
		dist.append((abs(feature[1] + (feature[2] - feature[1]) / 2 - center), feature[3]))
	return peak, dist
		


CATCH_SPACING = 10
def catch_formatter(peak, overlap):
	rstr = "# %s:%s-%s" % peak[0:3]
	#print peak, overlap
	if overlap:
		# left boundary
		first = overlap[0]
		#print "FIRST!", first
		if first[1] > peak[1] + CATCH_SPACING:
			rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[1] + CATCH_SPACING, 0)])
		else:
			rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[1], first[2], first[3])])
		
		last = first
		for f in overlap[1:-1]:
			if f[1] > last[2] + CATCH_SPACING:
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], last[2], last[2] + CATCH_SPACING, 0)
			if f[1] >= last[2] + CATCH_SPACING * 3:
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[1] - CATCH_SPACING, f[1], 0)
			rstr += "\n" + "\t".join(map(str, f[:4]))
			last = f

		if len(overlap) > 0 :
			f = overlap[-1]
			if f[2] > peak[2]:
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[1], peak[2], f[3])
			elif f[2] <= peak[2] - CATCH_SPACING:
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[1], f[2], f[3])
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[2], f[2] + CATCH_SPACING, 0)
			elif f[2] < peak[2]:
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[1], f[2], f[3])
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[2], peak[2], 0)
	
	else:
		rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[1] + CATCH_SPACING, 0)])
		rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[2] - CATCH_SPACING, peak[2], 0)])	
	return rstr


def peak_stats(peak_track, data_track, formatter=number_formatter, formatter_options={}, zeroes=True):
	ret = []
	peak_feature = peak_track.get_next_feature()
	#print "p1:", peak_feature	
	data_feature = None
	data_feature = data_track.get_next_feature()
	#print "d1:", data_feature	

	prev_data_seq = ""
	while peak_feature:
		#print "peak:", peak_feature
		overlap = []
		#print "step 1"	
		while data_feature and peak_feature and (data_feature[0] > peak_feature[0]):
			#print data_feature, peak_feature
			if zeroes:	
				ret.append(formatter(peak_feature, [], formatter_options))
			peak_feature = peak_track.get_next_feature()
			#print "p2:", peak_feature	
	
		if peak_feature:
			while (data_feature and ((data_feature[0] < peak_feature[0]) or ((data_feature[0] == peak_feature[0]) and (data_feature[2] < peak_feature[1])))):
#				print "Neee", data_feature
				data_feature = data_track.get_next_feature()
#				print "d2:", data_feature	
	
			while (data_feature and (data_feature[1] <= peak_feature[2] and data_feature[0] == peak_feature[0])):
#				print "Hier moet ik zijn", data_feature
				#print "Adding %s" % str(data_feature)
				overlap.append(data_feature)
				data_feature = data_track.get_next_feature()
				#print "d3:", data_feature	

			if len(overlap) > 0:
				ret.append(formatter(peak_feature, overlap, formatter_options))
			else:	
				if zeroes:
					ret.append(formatter(peak_feature, [], formatter_options))
				#sys.stderr.write("NO OVERLAP: %s\t%s\t%s\n" % peak_feature[:3])
		#print "step 2"	
		
		#print "peak before", peak_feature
		peak_feature = peak_track.get_next_feature()
		#print "peak after", peak_feature
		
		#print "step 3"	
		
		#print "p3:", peak_feature	
		#print data_feature
		data_feature = data_track.get_previous_feature()
		#print "step 4"	
		while data_feature and peak_feature and ((data_feature[2] >= peak_feature[1] and data_feature[0] == peak_feature[0]) or data_feature[0] > peak_feature[0]):
			#print data_feature, "-", peak_feature
			#print "GO BACK peak", peak_feature, " data", data_feature
			data_feature = data_track.get_previous_feature()
		#print "step 5"	
		
		#print "This is where we are: %s" % str(data_feature)
		if data_feature:
			data_feature = data_track.get_previous_feature()
		
		#print "step 3"	
		#print "This is where we are: %s" % str(data_feature)
		if not data_feature:
			data_feature = data_track.get_next_feature()
		#print "This is where we are: %s" % str(data_feature)
	return ret

def binned_peak_stats(peak_track, data_track, window, binsize, formatter=number_formatter, zeroes=True):
	ret = []
	peak_feature = peak_track.get_next_feature()
	#print "p1:", peak_feature	
	data_feature = None
	data_feature = data_track.get_next_feature()
	#print "d1:", data_feature	

	prev_data_seq = ""
	overlap = {}
	for bin in range(-window, window, binsize):
		overlap[bin] = []
		
	while peak_feature:
		chr = peak_feature[0]
		start = (peak_feature[1] + peak_feature[2]) / 2 + bin
		while data_feature and ((data_feature[2] >= start and data_feature[0] == chr) or data_feature[0] > chr):
			data_feature = data_track.get_previous_feature()
		if not data_feature:
			data_feature = data_track.get_next_feature()
		for bin in range(-window, window, binsize):
			start = (peak_feature[1] + peak_feature[2]) / 2 + bin
			end = start + binsize - 1 

			while (data_feature and ((data_feature[0] < chr) or ((data_feature[0] == chr) and (data_feature[2] < start)))):
				data_feature = data_track.get_next_feature()
	
			#if not data_feature or data_feature[1] > end or data_feature[0] != chr:
			#	overlap[bin].append(0)
			#else:
			while (data_feature and (data_feature[1] <= end and data_feature[0] == chr)):
				
				if  peak_feature[4] == "+":
					overlap[bin].append(data_feature[3])
					#print data_feature, "in bin", bin
				else:
					overlap[-bin - binsize].append(data_feature[3])
					#print data_feature, "in bin", -bin -binsize
					
				data_feature = data_track.get_next_feature()

			
		peak_feature = peak_track.get_next_feature()
	ret = []
	for bin in range(-window, window, binsize):
		#print overlap[bin]
		ret.append([bin, mean(overlap[bin])])
	return ret
