#!/usr/bin/env python
from numpy import *

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
