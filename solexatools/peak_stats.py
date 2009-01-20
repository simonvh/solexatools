#!/usr/bin/env python
import sys
import re
import os
from numpy import max
from solexatools.track import SimpleTrack

def max_val(features):
	test = map(lambda x: x[3], features)
	return max(test)

def max_feature(features):
	features.sort(cmp=lambda x,y: cmp(x[3],y[3]))
	return features[-1]

def fraction(features, f):
	m = max_val(features)
	return filter(lambda x:x[3] >= f * m, features)

def max_formatter(peak, overlap):
	if overlap:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], max_val(overlap))
	else:
		return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], 0)

def maxfeature_formatter(peak, overlap):
	return "\t".join([str(x) for x in max_feature(overlap)])

def number_formatter(peak, overlap):
    return "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[2], len(overlap))

def all_formatter(peak, overlap):
	return "\n".join("%s\t%s\t%s\t%s" % x for x in overlap)

CATCH_SPACING = 10
def catch_formatter(peak, overlap):
	rstr = "# %s:%s-%s" % peak[0:3]
	#print peak, overlap
	if overlap:
		# left boundary
		first = overlap.pop()
		#print "FIRST!", first
		if first[1] > peak[1] + CATCH_SPACING:
			rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[1] + CATCH_SPACING, 0)])
		else:
			rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[1], first[2], first[3])])
		
		last = first[2]
		for f in overlap[:-1]:
			if f[1] > last + CATCH_SPACING:
				rstr += "\n%s\t%s\t%s\t%s" % (peake[0], last, last + CATCH_SPACING, 0)
			if f[1] >= last + CATCH_SPACING * 3:
				rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[1] - CATCH_SPACING, f[1], 0)
			rstr += "\n" + "\t".join(map(str, f))
			last = f

		if len(overlap) > 0 :
			f = overlap[-1]
			rstr += "\n" + "\t".join(map(str, f))
			if f[2] > peak[2] - CATCH_SPACING:
				 rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[1], peak[2], f[3])
			else: 
				if f[2] <= peak[2] - CATCH_SPACING * 3:
					rstr += "\n%s\t%s\t%s\t%s" % (peak[0], f[2], f[2] + CATCH_SPACING, 0)
				if f[2] < peak[2]:
					rstr += "\n%s\t%s\t%s\t%s" % (peak[0], peak[2] - CATCH_SPACING, peak[2], 0)
	
	else:
		rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[1], peak[1] + CATCH_SPACING, 0)])
		rstr = "\n".join([rstr, "%s\t%s\t%s\t%s" % (peak[0], peak[2] - CATCH_SPACING, peak[2], 0)])	
	return rstr



def peak_stats(peak_track, data_track, formatter=number_formatter, zeroes=True):
	ret = []
	peak_feature = peak_track.get_next_feature()
	#print "p1:", peak_feature	
	data_feature = None
	data_feature = data_track.get_next_feature()
	#print "d1:", data_feature	

	prev_data_seq = ""
	while peak_feature:
		overlap = []

		while data_feature and peak_feature and (data_feature[0] > peak_feature[0]):
			if zeroes:	
				ret.append(formatter(peak_feature, []))
			peak_feature = peak_track.get_next_feature()
			#print "p2:", peak_feature	
	
		if peak_feature:
			while (data_feature and ((data_feature[0] < peak_feature[0]) or ((data_feature[0] == peak_feature[0]) and (data_feature[2] < peak_feature[1])))):
				data_feature = data_track.get_next_feature()
				#print "d2:", data_feature	
	
			while (data_feature and (data_feature[1] <= peak_feature[2] and data_feature[0] == peak_feature[0])):
				overlap.append(data_feature)
				data_feature = data_track.get_next_feature()
				#print "d3:", data_feature	

			ret.append(formatter(peak_feature, overlap))
		peak_feature = peak_track.get_next_feature()
		#print "p3:", peak_feature	

	return ret

