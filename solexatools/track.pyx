# cython: profile=True
import re
import os
import sys
from tempfile import NamedTemporaryFile
from numpy import *

from libc.string cimport strncmp
from libc.string cimport strtok
from libc.string cimport strlen
from libc.stdio cimport sprintf
from libc.stdio cimport EOF

cdef extern from "stdio.h":
	ctypedef void* FILE
	FILE* fopen(char* filename, char* mode)
	char *fgets(char *S, int count, FILE *stream)
	char *fputs(char *S, FILE *stream)
	int fclose(FILE *stream)

cdef extern from "stdlib.h":
	int atoi(char*)

# pysam for BAM/SAM support 
sam_support = True
try:
	import pysam
except ImportError:
	sam_support = False

class Track():
	
	def __init__(self, fname, sorted=False):
		pass

	def __iter__(self):
		return self

	def next(self):
		pass
		# if last 
		# raise StopIteration

#cdef class SamTrack:
#	cdef FILE* fh
#	
#	def __init__(self, fname, sorted=False):
#			
#		if sorted:
#			fh = fopen(<char *> fname, "r")
#		else:
#			tempname = self._sorted_sam_file(fname)
#			fh = fopen(<char *> tempname, "r")
#	
#	def __iter__(self):
#		return self
#		
#	def next(self):
#		cdef int start
#		cdef int flag 
#		cdef char* comment = "@"
#		cdef int i
#		cdef int unmapped_flag = 0x0004
#		cdef int strand_flag = 0x0010
#		cdef int length 
#		cdef char* p
#		cdef char* chrom
#		cdef FILE* fin
#		cdef FILE* fout
#		cdef char line[500]
#		cdef char out_line[500]
#
#		while (fgets(line, 1000, iself.fh) != NULL):
#			if strncmp(line, comment, 1):
#				p = strtok(line, "\t")
#				flag = atoi(strtok(NULL, "\t"))
#				if not (flag & unmapped_flag):
#					chrom = strtok(NULL, "\t")
#					start = atoi(strtok(NULL, "\t")) - 1
#					for i in range(5):
#						strtok(NULL, "\t")
#					length = strlen(strtok(NULL, "\t"))
#					
#					if flag & strand_flag:
#						return(chrom, start, start + length, 0, "-")
#					else:
#						return(chrom, start, start + length, 0, "+")
#		if (fgets(line, 1000, self.fh) == NULL):
#			raise StopIteration
#		
#	def _sorted_sam_file(self, fname):
#		temp = NamedTemporaryFile()
#		tempname = temp.name
#		os.system("sort -S 4G -k3,3 -k4g,4 %s > %s" % (fname, tempname))
#		return tempname


class SimpleTrack:
	def __init__(self, file, mem=False, sorted=False):
		self.BUFSIZE = 100000
		self.fh = None
		
		self.format = self.guess_format(file)
		#sys.stderr.write("FORMAT %s %s\n" % (file, self.format))
		
		if self.format == "unknown":
			print "Unknown format: %s" % file
			sys.exit(1)

		if self.format == "bam":
			self.fixed = False
			self.fh = self.open_sorted_bam_file(file)
		elif self.format == "sam":
			self.fixed = False
			self.fh = self.open_sorted_sam_file(file)
		else:
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
		self.p = re.compile(r'chrom=(\w+)\s+start=(\d+)\s+step=(\d+)(\s+span=(\d+))?')
		self.line_index = []
		self.lines = []
		self.line_index.append(self.fh.tell())
		if mem:
			self.lines = self.fh.readlines()
			self.number_of_lines = len(self.lines)
		else:
			self.readlines()
		self.index = 0

	def guess_format(self, file):
		if file[-3:].lower() in ["bam", "sam"]:
			if not sam_support:
				print "SAM/BAM support not available, please install pysam!"
				sys.exit(1)
			return  file[-3:].lower()
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
				try:
					int(vals[3]) + int(vals[4])
					return "gff"
				except:
					pass
		return "unknown"

	def open_sorted_file(self, file):
		temp = NamedTemporaryFile()
		tempname = temp.name
		os.system("sort -k1,1 -k2g,2 %s > %s" % (file, tempname))
		return temp

	def open_sorted_sam_file(self, fname):
		temp = NamedTemporaryFile()
		cdef char* tempname = temp.name
		cdef int start
		cdef int flag 
		cdef char* comment = "@"
		cdef int i
		cdef int unmapped_flag = 0x0004
		cdef int strand_flag = 0x0010
		cdef int length 
		cdef char* p
		cdef char* chrom
		cdef FILE* fin
		cdef FILE* fout
		cdef char line[500]
		cdef char out_line[500]

		fin = fopen(<char *> fname, "r")
		fout = fopen(<char *> temp.name, "w")
		while (fgets(line, 1000, fin) != NULL):
			if strncmp(line, comment, 1):
				p = strtok(line, "\t")
				flag = atoi(strtok(NULL, "\t"))
				if not (flag & unmapped_flag):
					chrom = strtok(NULL, "\t")
					start = atoi(strtok(NULL, "\t")) - 1
					for i in range(5):
						strtok(NULL, "\t")
					length = strlen(strtok(NULL, "\t"))
					
					if flag & strand_flag:
						sprintf(out_line, "%s\t%d\t%d\t0\t0\t-\n", chrom, start, start + length)
					else:
						sprintf(out_line, "%s\t%d\t%d\t0\t0\t+\n", chrom, start, start + length)
					fputs(out_line, fout)
		
		fclose(fin)
		fclose(fout)

		temp2 = NamedTemporaryFile()
		tempname2 = temp2.name
		os.system("sort -S 4G -k1,1 -k2g,2 %s > %s" % (tempname, tempname2))
		temp.close()
		return temp2

	def open_sorted_bam_file(self, file):
		temp = NamedTemporaryFile()
		tempname = temp.name
		pysam.sort("-n", file, tempname)
		temp2 = NamedTemporaryFile()
		sam = pysam.Samfile(tempname)
		for read in sam:
			if not read.is_unmapped:
				chrom = sam.getrname(read.rname)
				start = read.pos
				end = read.qend
				if end != read.qlen:
					raise Exception, "Read does not map completely, don't know what to do with this!"
				strand = "+"
				if read.is_reverse:
					strand = "-"
				temp2.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start, start + end, 0, 0, strand))
		temp2.seek(0)
		temp.close()
		return temp2

	def readlines(self):
		lines = self.fh.readlines(self.BUFSIZE)
		if lines:
			self.lines = lines 
			self.number_of_lines = len(self.lines)
		else:
			self.eof = True
		

	def readline(self):
		#make use of a buffer to speed up the reading of a file
		if self.lines and self.index < len(self.lines):
			self.index += 1
			return self.lines[self.index - 1]
		else:
			ind = self.fh.tell()
			self.readlines()
			if self.eof:
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
		#print "index", self.index
		if self.index <= 1 and len(self.line_index) <= 1:
			self.index_min_1()
			return None
	
		if self.fixed and self.index == 2 and len(self.line_index) == 1:
			return None

		self.index_min_1()
		if not self.eof:
			self.index_min_1()
		#print "index", self.index
		#print "len", len(self.lines)
		if self.fixed:
			#print "Fixed!"
			if self.lines[self.index].startswith("fixed"):
				#print "min 1"
				self.index_min_1()
			target = self.index
			
			if self.index == 0:
				return None
			while not self.lines[self.index].startswith("fixed"):
				#print "while - 1, index", self.index
				self.index_min_1()
			#print "extra - 1"
			self.index_min_1()
			
			for i in range(target - self.index - 2):
				self.get_next_feature()
		
		#if self.eof:
		#	self.fh.seek(self.line_index[-1])
		#	self.readlines()
		#	self.index = self.number_of_lines - 1
		self.eof = False
		return self.get_next_feature()
	
	def get_next_feature(self):
		line = self.readline()
		
		while (line and (line[0] == "#" or line.startswith("track") or line == "\n")):
			line = self.readline()
		if not(line):
			return None
		if self.fixed:
			m = self.p.search(line)
			if m:
				self.chrom = m.group(1)
				self.start = int(m.group(2))
				self.step = int(m.group(3))
				self.span = self.step
				if m.group(4):
					self.span = int(m.group(5))
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
					if len(vals) > 5:
						strand = vals[5]	
			#print self.index,  (vals[0], start, end, value, strand)
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
