#!/usr/bin/env python
import re
import os
from tempfile import NamedTemporaryFile
from numpy import *

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
			if len(vals) >= 4:
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

class TrackStats():
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



