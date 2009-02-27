#!/usr/bin/env python
import re
import os
from tempfile import NamedTemporaryFile
from numpy import *

class SimpleTrackIterator:
	def __init__(self):
		self.c = 0

	def next(self):
		self.c += 1

		return self.c


class SimpleTrack:
	def __init__(self, file, mem=False, sorted=None):
		self.BUFSIZE = 100000
		self.format = self.guess_format(file)
		self.fixed = self.is_fixedstep(file)
		if not(self.fixed) and not(sorted):
			self.fh = self.open_sorted_file(file)
		else:
			self.fh = open(file)
		self.chrom = None
		self.start = None
		self.step = None
		self.span = None
		self.eof = False
		self.p = re.compile(r'chrom=(\w+)\s+start=(\d+)\s+step=(\d+)\s+span=(\d+)')
		self.line_index = []
		self.line_index.append(self.fh.tell())
		if mem:
			self.lines = self.fh.readlines()
			self.number_of_lines = len(self.lines)
		else:
			self.readlines()
		self.index = 0

	def guess_format(self, file):
		if file[-3:] in ["gff"]:
			return "gff"
		if self.is_fixedstep(file):
			return "fixed"
		f = open(file)
		line = f.readline()
		if line.find("track") != -1:
			if line.find("wiggle_0") != -1:
				return "wiggle_0"
			else:
				return self.guess_line_format(f.readline())
		else:
			return self.guess_line_format(line)

	def guess_line_format(self, line):
		vals = line.strip().split()
		if len(vals) == 3:
			return "bed3"
		if len(vals) == 4:
			return "bed4"
		if len(vals) == 6:
			return "bed6"
		if len(vals) == 9:
			try:
				int(vals[1]) + int(vals[2])
				return "bed9"
			except:
				return "gff"

		return "unknown"

	def open_sorted_file(self, file):
		temp = NamedTemporaryFile()
		tempname = temp.name
		os.system("sort -k1,1 -k2g,2 %s > %s" % (file, tempname))
		return temp

	def readlines(self):
		self.lines = self.fh.readlines(self.BUFSIZE)
		self.number_of_lines = len(self.lines)

	def readline(self):
		#make use of a buffer to speed up the reading of a file
		if self.lines and self.index < len(self.lines):
			self.index += 1
			return self.lines[self.index - 1]
		else:
			ind = self.fh.tell()
			self.readlines()
			if not self.lines:
				self.eof = True
				return None
			self.line_index.append(ind)
			self.index = 1
			return self.lines[0]

	def index_min_1(self):
		self.index -= 1

		if self.index < 0 :
			if len(self.line_index) == 1:
				self.index += 1
				return None
			del self.line_index[-1]
			self.fh.seek(self.line_index[-1])
			self.readlines()
			self.index = self.number_of_lines - 1


	def get_previous_feature(self):
		if self.index == 1 and len(self.line_index) == 1:
			return None
	
		if self.fixed and self.index == 2 and len(self.line_index) == 1:
			return None

		self.index_min_1()
		self.index_min_1()
		
		if self.fixed:
			if self.lines[self.index].startswith("fixed"):
				self.index_min_1()
			target = self.index
			while not self.lines[self.index].startswith("fixed"):
				self.index_min_1()
			self.index_min_1()
			
			for i in range(target - self.index - 1):
				self.get_next_feature()
		
		if self.eof:
			self.fh.seek(self.line_index[-1])
			self.readlines()
			self.index = self.number_of_lines - 1
			self.eof = False
		return self.get_next_feature()
	
	def get_next_feature(self):
		line = self.readline()
		#print line
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
			start, end, value, strand = None, None, None, None
			if self.format == "gff":
				(start, end) = map(int, vals[3:5])
				value = vals[8]
				strand = vals[6]
			else:
				(start, end) = map(int, vals[1:3])
				if len(vals) >= 4:
					try:
						value = float(vals[3])
					except:
						value = vals[3]
					if len(vals) >= 5:
						strand = vals[5]	
			return (vals[0], start, end, value, strand)

	def is_fixedstep(self,file):
		f = open(file)
		for i in range(3):
			line = f.readline()
			if line.startswith("fixedStep"):
				f.close()
				return True
		f.close()
		return False

class ExonTrack:
	def __init__(self, file):
			""" Reads a GTF file of exons"""	
			self.exons = {}
			self.chr = {}
			self.strand = {}
			
			# gene_id "NM_001005484"; transcript_id "NM_001005484"; 
			p = re.compile(r'(gene|transcript_id) "([^"]+)"')
			
			for line in open(file):
				vals = line[:-1].split("\t")
				name = p.search(vals[8]).group(2)

				self.chr[name] = vals[0]
				self.strand[name] = vals[6]
				if self.exons.has_key(name):
					self.exons[name].append((int(vals[3]), int(vals[4])))
				else:
					self.exons[name] = [(int(vals[3]), int(vals[4]))]
			
			# sort all exons based on start position
			for name in self.exons.keys():
				self.exons[name].sort(lambda x,y: x[0] - y[0])
	
	def get_exons(self):
		t = []
		for name in self.exons.keys():
			for exon in self.exons[name]:
				t.append((self.chr[name], exon[0], exon[1], name))
		return t

class TrackStats:
	def __init__(self, file):
		s = SimpleTrack(file)
		
		values = []
		f = s.get_next_feature()
		while f:
			values.append(f[3])
			f = s.get_next_feature()

		values = array(values, dtype=float)
		self.mean = mean(values)
		self.median = median(values)
		self.min= min(values)
		self.max = max(values)


if __name__ == "__main__":
	t = ExonTrack("test.gtf")
	for name in t.exons.keys():
		print name, t.exons[name]
	print t.get_exons()
