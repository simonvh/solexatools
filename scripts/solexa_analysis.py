#!/usr/bin/env python
import os
import re
import sys 
import ConfigParser
from SolexaTools.solexatools import assign_peaks
from datetime import *
BUFSIZE = 10000

def to_bool(x):
	return str(x).strip().lower() not in ['f', 'false', '0', '', 'no']

class SolexaExperiment:
	def __init__(self, params):
		self.params = params
		self.name = params["name"]
		self.elandfiles = []
		for elandfile in [params[x] for x in params.keys() if x.startswith("elandfile")]:
			self.elandfiles.append(elandfile)

		# Bed file params
		self.bed_file = ""
		self.extend_to = int(self.params["bed_extend_to"])
		self.read_length = int(self.params["read_length"])
		self.mm = self.params["bed_use_mismatch"]
		if re.compile(r'[^012]').search(self.mm):
			print "Invalid value %s for parameter bed_use_mismatch" % self.mm
			sys.exit()
		
		# Wiggle
		self.wig_create = to_bool(self.params["wig_create"])
		self.wig_extend_to = int(self.params["wig_extend_to"])
		self.wig_color = self.params["wig_color"]
		if self.wig_color == "":
			self.wig_color = 1

		# Normalization
		self.do_normalization = to_bool(self.params["do_normalization"])
			
		# Macs
		self.peak_files = []
		self.macs_call_peaks = to_bool(self.params["macs_call_peaks"])
		if self.macs_call_peaks:
			self.macs_bandwidth = self.params["macs_bandwidth"].split(",")
			self.macs_pvalue = self.params["macs_pvalue"].split(",")


	def make_bed_file(self, bed_dir, overwrite=True):
		self.bed_file = os.path.join(bed_dir, "%s_%s.bed" % (self.name, self.extend_to)) 
		p = re.compile(r'U[%s]' % self.mm)
		strandmap = {"F":"+", "R":"-"}
		if os.path.exists(self.bed_file) and not overwrite:
			return
		out = open(self.bed_file, "w")
		seen = {}
		
		for file in self.elandfiles:
			
			f = open(file)
			lines = f.readlines(BUFSIZE)
			while lines:
				for line in lines:
					vals = line[:-1].split("\t")
					if p.match(vals[2]):
						chr = vals[6].split('/')[-1].replace(".fa", "")
						strand = strandmap[vals[8]]
						start = int(vals[7])
						end = start + self.read_length
						if strand == "+":
							end += self.extend_to - self.read_length
						else:
							start -= (self.extend_to - self.read_length)

						if not seen.has_key("%s:%s-%s" % (chr, start, strand)):
							out.write("%s\t%s\t%s\t%s\t0\t%s\n" % (chr, start, end, vals[2], strand))
							seen["%s:%s-%s" % (chr, start, strand)] = 1
				lines = f.readlines(BUFSIZE)
			f.close()
			
	def call_peaks(self, peak_dir):
		cwd = os.getcwd()
		os.chdir(peak_dir)
		for pvalue in self.macs_pvalue:
			for bandwidth in self.macs_bandwidth:
				name = "%s_%s_p%s_bw%s" % (self.name, self.read_length, pvalue, bandwidth)
				command = "macs -t %s --name %s format BED --tsize %s --pvalue %s --bw %s" % (self.bed_file, name, self.read_length, pvalue, bandwidth)
			        os.system(command)
				self.peak_files.append(os.path.join(peak_dir, "%s_peaks.bed" % name))
		os.chdir(cwd)
	
	def create_wiggle(self, sum_bed, wig_dir):
		
		self.wig_file = os.path.join(wig_dir, "%s_%s.wig" % (self.name, self.wig_extend_to)) 
		if self.bed_file == "":
			print "No bedfile found. Make bedfiles first!"
			sys.exit()

		command = "%s -i %s -c %s -e %s > %s" % (sum_bed, self.bed_file, self.wig_color, self.wig_extend_to, self.wig_file)
		print command
		os.system(command)

class SolexaAnalysis:
	def __init__(self, config):
		self.experiments = []
		self.do_directory_settings(config)
		self.do_experiment_settings(config)
		self.do_program_settings(config)
		
		if config.has_section("assignpeaks"):
			self.do_assignpeaks_settings(config)
		
	def do_program_settings(self, config):
		self.normtags = config.get("programs", "normtags")
		self.sum_bed = config.get("programs", "sum_bed")
		config.remove_section("programs")

	def do_directory_settings(self, config):
		out_dir = config.get("directories", "output_dir")
		# if empty use current directory
		if out_dir == "":
			out_dir = os.getcwd()

		# expand out_dir to full path
		out_dir = os.path.abspath(out_dir)
		
		# check if directory exists
		if not os.path.exists(out_dir):
			print "Output directory %s does not exist!" % out_dir
			sys.exit()

		# set directories 
		self.out_dir = out_dir
		self.eland_dir = os.path.join(out_dir, config.get("directories", "eland_dir"))
		self.wig_dir = os.path.join(out_dir, config.get("directories", "wig_dir"))
		self.peak_dir = os.path.join(out_dir, config.get("directories", "peak_dir"))
		self.bed_dir = os.path.join(out_dir, config.get("directories", "bed_dir"))

	def do_experiment_settings(self, config):
		experiments =  [x for x in config.sections() if x.startswith("experiment")]

		if len(experiments) == 0:
			print "No experiment section(s) specified!"
			sys.exit()
		
		for section in experiments:
			# Expand elandfiles is necessary
			params = dict(config.items("default_settings")) 
			for (k,v) in config.items(section):
				if k.startswith("elandfile"):
					if os.path.exists(v):
						params[k] = v
					elif os.path.exists(os.path.join(self.eland_dir, v)):
						params[k] = os.path.join(self.eland_dir, v)
					else:
						print "Cannot find eland file %s" % v
						sys.exit()
				else:
					params[k] = v
					
			print params
			self.add_experiment(SolexaExperiment(params))
	
	def do_assignpeaks_settings(self, config):
		self.do_assignpeaks = False
		if config.has_option("assignpeaks", "do_assignpeaks") and (config.getboolean("assignpeaks", "do_assignpeaks") == True):
			self.do_assignpeaks = True
			self.genefile = config.get("assignpeaks", "genefile")
			self.upstream = config.getint("assignpeaks", "upstream")
			self.extend = config.getint("assignpeaks", "extend")
		config.remove_section("assignpeaks")

	def add_experiment(self, exp):
		self.experiments.append(exp)

	def make_bed_files(self, overwrite=True):
		if not os.path.exists(self.bed_dir):
			os.mkdir(self.bed_dir)
		for exp in self.experiments:
			exp.make_bed_file(self.bed_dir, overwrite)
	
	def normalize(self, run_command=True):
		bedfiles = []
		for exp in self.experiments:
			if exp.do_normalization:
				bedfiles.append(exp.bed_file)
		print "Normalization: %s" % (str(bedfiles))
		
		if len(bedfiles) > 0:
			first = bedfiles[0]
			command = "%s --input %s %s" % (self.normtags, first, " ".join(bedfiles[1:]))
			if run_command:
				print command
				os.system(command)

		for exp in self.experiments:
			if exp.do_normalization:
				exp.bed_file = exp.bed_file.replace(".bed", "_norm.bed")
	
	def call_peaks(self):
		if not os.path.exists(self.peak_dir):
			os.mkdir(self.peak_dir)
		
		for exp in self.experiments:
			if exp.macs_call_peaks:
				exp.call_peaks(self.peak_dir)

	def create_wiggle(self):
		if not os.path.exists(self.wig_dir):
			os.mkdir(self.wig_dir)
		
		for exp in self.experiments:
			exp.create_wiggle(self.sum_bed, self.wig_dir)

	def assign_peaks(self):

		if self.do_assignpeaks:
			if not os.path.exists(self.genefile):
				print "Genefile %s not found! Cannot assign peaks to genes." % self.genefile
				sys.exit()
			for exp in self.experiments:
				for peakfile in exp.peak_files:
					overlap_file = peakfile.replace(".bed",  "_" + os.path.split(self.genefile)[-1])
					f = open(overlap_file, "w")
					matrix = assign_peaks(self.genefile, peakfile, self.upstream, self.extend)
					f.write("#assign_peaks.py %s\n" % (datetime.today().strftime("%d-%m-%Y %H:%M")))
					f.write("#peaks: %s\n" % peakfile)
					f.write("#targets: %s\n" % self.genefile)
					f.write("#extension: %s maximum upstream: %s\n" % (self.extend, self.upstream))
					for gene, target in matrix.items():
						for (chr, middle, target, val) in target:
							f.write("%s\t%s\t%s\t%s\t%s\n" % (gene, chr, middle, target, val))
					f.close()
	def run(self):
		s.make_bed_files()
		s.normalize()
		s.call_peaks()
		s.create_wiggle()
		s.assign_peaks()
		

if not os.path.exists("my_solexa.cfg"):
	print "No configuration file found!"
	print "This script is still under development by Simon although quite usable already"
	print "It will do some basic analyses:"
	print "\t* combination of lanes"
	print "\t* normalization"
	print "\t* peak calling"
	print "\t* wiggle file creaion"
	print "\t* assign peaks to genes"
	print "Multiple settings are possible, defined by a configuration file"
	print "Ask Simon for an example config"
	sys.exit()

config = ConfigParser.RawConfigParser()
config.read('my_solexa.cfg')

s = SolexaAnalysis(config)
s.run()

