#!/usr/bin/env python
import sys
from optparse import OptionParser
from tempfile import NamedTemporaryFile
from solexatools.track import *
from solexatools import peak_stats
from math import log

DEFAULT_UP = 250
DEFAULT_DOWN = 500

parser = OptionParser()
parser.add_option("-p", "--polII data", dest="polfile", help="PolII reads (BED format)", metavar="FILE")
parser.add_option("-g", "--genefile", dest="genefile", help="Genes in 6-BED, GTF or GFF format", metavar="FILE")
parser.add_option("-u", "--upstream", dest="upstream", 
		help="Promoter definition: upstream (default: %s)" % DEFAULT_UP, 
		metavar="FILE", default=DEFAULT_UP)
parser.add_option("-d", "--downstream", dest="downstream", 
		help="Promoter definition: downstream (default: %s)" % DEFAULT_DOWN, 
		metavar="FILE", default=DEFAULT_DOWN)
#parser.add_option("-z", "--zeroes", dest="zeroes", help="Print zeroes", action="store_true", default=False)

sys.stderr.write("NOT FINAL!! STILL UNDER DEVELOPMENT\n")
(options, args) = parser.parse_args()

if not options.polfile or not options.genefile:
	parser.print_help()
	sys.exit(0)

def genes_to_promoters_and_exons(file, up, down):
	exon_data = {}
	chr = {}
	strand = {}
	for line in open(file):
		vals = line[:-1].split("\t")
		name = vals[8].split(";")[0].split(" ")[1]

		chr[name] = vals[0]
		strand[name] = vals[6]
		if exon_data.has_key(name):
			exon_data[name].append((int(vals[3]), int(vals[4])))
		else:
			exon_data[name] = [(int(vals[3]), int(vals[4]))]

	exons = []
	location = {}
	promoter = []
	
	for gene, pos in exon_data.items():
		pos.sort(lambda x,y: x[0] - y[0])
		location[gene] = "%s:%s-%s" % (chr[gene], pos[0][0], pos[-1][1])
		if strand[gene] == "+":
			if pos[0][0] - up > 0:
				promoter.append([chr[gene], pos[0][0] - up, pos[0][0] + down, gene])
			else: 
				promoter.append([chr[gene], 0, pos[0][0] + down, gene])
			
			# Whole gene body

			if (pos[-1][1] - pos[0][0]) - down > 0:
				exons.append([chr[gene], pos[0][0] + down, pos[-1][1], gene])
				
			# EXONS ONLY
			#
			#if (pos[0][1] - pos[0][0]) - down > 0:
			#	exons.append([chr[gene], pos[0][0] + down, pos[0][1], gene])
			#for p in pos[1:]:
			#	if p[1] > pos[0][0] + down:
			#		if p[0] > pos[0][0] + down:
			#			exons.append([chr[gene], p[0] , p[1], gene])
			#		else:
			#			exons.append([chr[gene], pos[0][0] + down , p[1], gene])
						
		else:		
			if pos[-1][1] - down > 0:
				promoter.append([chr[gene], pos[-1][1] - down, pos[-1][1] + up, gene])
			else:
				promoter.append([chr[gene], 0, pos[-1][1] + up, gene])
			
			if pos[-1][1] - pos[0][0] > down:
				exons.append([chr[gene], pos[0][0], pos[-1][1] - down, gene])

			
			#if (pos[-1][1] - pos[-1][0]) - down > 0:
			#	exons.append([chr[gene], pos[-1][0], pos[-1][1] - down, gene])
			#for p in pos[:-1]:
			#	if p[0] < pos[-1][1] - down:
			#		if p[1] < pos[-1][1] - down:
			#			exons.append([chr[gene], p[0] , p[1], gene])
			#		else: 
			#			exons.append([chr[gene], p[0] , pos[-1][1] - down, gene])

	return promoter, exons, location

upstream = int(options.upstream)
downstream = int(options.downstream)
promoter, exons, location = genes_to_promoters_and_exons(options.genefile, upstream, downstream)
tmp = NamedTemporaryFile()
tmpname = tmp.name
f = open(tmpname, "w")
seen = {}
for row in promoter:
	p = "%s:%s-%s" % (row[0], row[1], row[2])
	s = "%s\t%s\t%s\t%s\n" % (row[0], row[1], row[2], row[3])
	if not seen.has_key(p):
		f.write(s)
		seen[p] = row[3]
f.close()

#for k,v in seen.items():
#	print "PROM\t%s\t%s" % (k,v)
prom_track = SimpleTrack(tmpname)
pol_track = SimpleTrack(options.polfile)
result = peak_stats.peak_stats(prom_track, pol_track, peak_stats.number_formatter)

prom_result = {}
for row in result:
	chr,start,end,num = row.split("\t")
	prom_result[seen["%s:%s-%s" % (chr, start,end)]] = int(num)


#GENES
tmp = NamedTemporaryFile()
tmpname = tmp.name
f = open(tmpname, "w")
seen = {}
for row in exons:
	p = "%s:%s-%s" % (row[0], row[1], row[2])
	s = "%s\t%s\t%s\t%s\n" % (row[0], row[1], row[2], row[3])
	if not seen.has_key(p):
		f.write(s)
		seen[p] = row[3]
f.close()
#for k,v in seen.items():
#	print "EXON\t%s\t%s" % (k,v)

exon_track = SimpleTrack(tmpname)
pol_track = SimpleTrack(options.polfile)
result = peak_stats.peak_stats(exon_track, pol_track, peak_stats.number_formatter)

gene_result = {}
gene_length = {}
for row in result:
	chr, start, end, num = row.split("\t")
	gene = seen["%s:%s-%s" % (chr, start,end)]
	chr,start,end,num = row.split("\t")
	if gene_result.has_key(gene):
		gene_result[gene] += int(num)
	else:
		gene_result[gene] = int(num)
	if gene_length.has_key(gene):
		gene_length[gene] += int(end) - int(start)
	else:
		gene_length[gene] = int(end) - int(start)


print "location\tgene_name\t# of tags / total_length\t# of tags in promoter\t# of tags gene body, norm to promoter\tlog2(prom/genebody)"
for gene in prom_result.keys():
	prom = float(prom_result[gene])
	body = 0
	body_norm = 0
	ratio = "NA"
	length = upstream
	if gene_result.has_key(gene):
		body = float(gene_result[gene])
		body_norm = body / float(gene_length[gene]) * (upstream + downstream)
		ratio = log((prom + 1)/ (body_norm + 1))/log(2)
		length = gene_length[gene] + upstream	
	print "%s\t%s\t%s\t%s\t%s\t%s" % (location[gene], gene, (prom + body) /  length , prom, body_norm, ratio )

