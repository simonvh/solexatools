#!/usr/bin/env python
import sys
import re
import os
from numpy import *
from datetime import datetime
from tempfile import NamedTemporaryFile
from optparse import OptionParser


def get_minimum_distance(a, b):
	# assign_peaks
	print a,b
	if not len(a) or not len(b):
		return [], [] 
	
	ax, bx = ix_(a, b)
	dist = bx - ax
	print "dits:", dist
	in_range = any(dist >= 0, 0)
	print "in_range:", in_range
	dist[dist < 0] = max(dist.flatten()) + 10 
	print dist
	close_idx = argmin(dist, 1)[in_range]
	print "LA:", close_idx,  dist[:,close_idx].min(axis=1)
	if len(close_idx):
		return close_idx, dist[:,close_idx].min(axis=1)
	else:
		return [], []

def min_gt_zero(a, c=0):
	# assign_peaks
	if len(a[a > 0]):
		m = min(a[a > 0])
		if m < c:
			m = c
		return a > m
	else:
		return ones(shape(a), dtype=bool)

	
def get_overlap(targets, genes, strand, extend=0):
	# assign_peaks
	# a and b are lists of start and end coordinates
	# coordinates of b will be extend by extend
	if not len(targets) or not len(genes):
		return [],[]
	
	# Define matrices
	(t_end, g_start) = ix_(targets[:,1], genes[:,0])
	(t_start, g_end) = ix_(targets[:,0], genes[:,1])

	# Get difference between gene_start and target_end
	dist_start = g_start - t_end
	# Get difference between target_start and gene_end
	dist_end = t_start - g_end

	overlap = dist_start <= 0
	dist_start[apply_along_axis(min_gt_zero, 1, dist_start, extend)] = -(extend + 1)
	dist_start[overlap] = 0
	
	overlap =  dist_end <= 0
	dist_end[apply_along_axis(min_gt_zero, 1, dist_end, extend)] = -(extend + 1)
	dist_end[overlap] = 0
	
	overlap = logical_not(logical_and(dist_start == 0, dist_end == 0))
	dist_start[logical_and(dist_start == 0, overlap)] = -(extend + 1)
	dist_end[logical_and(dist_end == 0, overlap)] = -(extend + 1)
	
	c = logical_and(dist_start < -extend, dist_end <= extend)
	c = logical_and(c, dist_end > 0)
	dist_start[c] = -dist_end[c]

	c = logical_and(dist_end < -extend, dist_start <= extend)
	c = logical_and(c, dist_start > 0)
	dist_end[c] = -dist_start[c]
	
	return where(vstack(len(targets) * [strand]), dist_start, dist_end)

def assign_peaks(genefile, bedfile, upstream, extend):
	upstream = int(upstream)
	extend = int(extend)
	
	extra = 0

	strand_map = {"+":True, "-":False, 1:True, -1:False, "1":True, "-1":False}
	str = {}
	pos = {}
	nme = {}

	for line in open(genefile):
		if not line.startswith("track"):
			(chr, start, end, name, score, strand) = line[:-1].split("\t")[:6]
			if not(pos.has_key(chr)):
				pos[chr] = []
				str[chr] = []
				nme[chr] = []
			strand= strand_map[strand]
	
			if strand:
			# + strand
				pos[chr].append([int(start) + extra, int(end)])
			else:
			# - strand
				pos[chr].append([int(start), int(end) - extra])
			str[chr].append(strand)
			nme[chr].append(name)

	for chr in pos.keys():
		pos[chr] = array(pos[chr])
		str[chr] = array(str[chr])
		nme[chr] = array(nme[chr])

	targets = {}
	target_names = {}
	for line in open(bedfile):
		if not line.startswith("track"):
			vals =  line[:-1].split("\t")
			name = ""
			if len(vals) > 3:
				name = vals[3]
			(chr, start, end) = vals[:3] 
			if not name:
				name = "%s:%s-%s" % (chr, start, end)
			if not targets.has_key(chr):
				targets[chr] = []
				target_names[chr] = []
			targets[chr].append([int(start),int(end)])
			target_names[chr].append(name)
		
	matrix = {}
	for chr in targets.keys():
		if pos.has_key(chr):
			t_all = array(targets[chr])
			t_middle = (t_all[:,0] + t_all[:,1]) / 2
			#print t_middle
			tn = array(target_names[chr])
	
			plus = pos[chr][:,0][str[chr]]
			minus = pos[chr][:,1][logical_not(str[chr])]
			
	
			for target_id,row in enumerate(get_overlap(t_all, pos[chr], str[chr], extend)):
				for gene_id, val in enumerate(row):
					if val >= -extend and val < upstream:
						matrix.setdefault(nme[chr][gene_id], []).append( [chr, t_middle[target_id], tn[target_id], val])
	return matrix

class SimpleTrack:
	def __init__(self, file, sorted=None):
		self.BUFSIZE = 10000
		self.fixed = self.is_fixedstep(file)
		if not(self.fixed) and not(sorted):
			self.fh = self.open_sorted_file(file)
		else:
			self.fh = open(file)
		self.chrom = None
		self.start = None
		self.step = None
		self.span = None
		self.p = re.compile(r'chrom=(\w+)\s+start=(\d+)\s+step=(\d+)\s+span=(\d+)')
		self.lines = self.fh.readlines(self.BUFSIZE)
		self.index = 0

	def open_sorted_file(self, file):
		temp = NamedTemporaryFile()
		tempname = temp.name
		os.system("sort -k1,1 -k2g,2 %s > %s" % (file, tempname))
		return temp

	def readline(self):
		#make use of a buffer to speed up the reading of a file
		if len(self.lines) and self.index < len(self.lines):
			self.index += 1
			return self.lines[self.index - 1]
		else:
			self.lines = self.fh.readlines(self.BUFSIZE)
			if not len(self.lines):
				return None
			self.index = 1
			return self.lines[0]

	def get_next_feature(self):
		line = self.readline()
		while (line and (line[0] == "#" or line.startswith("track"))):
			line = self.readline()
		if not(line):
			return None

		if self.fixed:
			m = self.p.search(line)
			if m:
				self.chrom = m.group(1)
				self.start = int(m.group(2))
				self.step = int(m.group(3))
				self.span = int(m.group(4))
				line = self.readline()
			value = float(line[:-1])
			retval = (self.chrom, self.start, self.start + self.span, value)
			self.start += self.step
			return retval
		else:
			vals = line[:-1].split("\t")
			(start, end) = map(int, vals[1:3])
			value = None
			if len(vals) == 4:
				try:
					value = float(vals[3])
				except:
					value = vals[3]
					
			return (vals[0], start, end, value)

	def is_fixedstep(self,file):
		f = open(file)
		for i in range(3):
			line = f.readline()
			if line.startswith("fixedStep"):
				f.close()
				return True
		f.close()
		return False

def max_val(features):
	test = map(lambda x: x[3], features)
	return max(test)

def max_feature(features):
	features.sort(cmp=lambda x,y: cmp(x[3],y[3]))
	return features[-1]

def fraction(features, f):
	m = max_val(features)
	return filter(lambda x:x[3] >= f * m, features)

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

