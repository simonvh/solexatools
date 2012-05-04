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

sys.stderr.write("WARNING: this script is not exhaustively tested!\n")
(options, args) = parser.parse_args()

if not options.polfile or not options.genefile:
	parser.print_help()
	sys.exit(0)

def genes_to_promoters_and_exons(file, up, down):
	exon_data = {}
	chr = {}
	strand = {}
	type = "BED"
	if file.endswith(".gff"):
		type = "GFF"
	line = ""
	f = open(file)
	while not line or line.startswith("browser") or line.startswith("track"):
		line = f.readline()
	vals = line.strip().split("\t")
	try:
		int(vals[3]), int(vals[4])
		type = "GFF"
	except:
		pass
	strand_col, start_col, end_col = 5, 1, 2
	if type == "GFF":
		strand_col, start_col, end_col = 6, 3, 4

	line_num = 1
	while line:
		
		vals = line[:-1].split("\t")
		if len(vals) <= strand_col:
			sys.stderr.write("ERROR: Strand column not present at line %s in file %s\n" % (line_num,file))
			sys.exit(1)
		
		if type == "BED":
			name = vals[3]
		else:
			try:
				name = vals[8].split(";")[0].split(" ")[1]
			except:
				name = vals[8]

		name = "%s_|_%s_%s" % (name, vals[0], vals[start_col])
		
		chr[name] = vals[0]
		
		strand[name] = vals[strand_col]
		if exon_data.has_key(name):
			exon_data[name].append((int(vals[start_col]), int(vals[end_col])))
		else:
			exon_data[name] = [(int(vals[start_col]), int(vals[end_col]))]

		line = f.readline()
		line_num += 1
	f.close()
	genes = []
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
			
			if (pos[-1][1] - pos[0][0]) - down > 0:
				genes.append([chr[gene], pos[0][0] + down, pos[-1][1], gene])
						
		else:		
			if pos[-1][1] - down > 0:
				promoter.append([chr[gene], pos[-1][1] - down, pos[-1][1] + up, gene])
			else:
				promoter.append([chr[gene], 0, pos[-1][1] + up, gene])
			
			if pos[-1][1] - pos[0][0] > down:
				genes.append([chr[gene], pos[0][0], pos[-1][1] - down, gene])

	return promoter, genes, location

upstream = int(options.upstream)
downstream = int(options.downstream)
promoter, genes, location = genes_to_promoters_and_exons(options.genefile, upstream, downstream)
tmp = NamedTemporaryFile()
tmpname = tmp.name
f = open(tmpname, "w")
seen = {}

genes = dict([(x[3], x[:3]) for x in genes])
promoter = dict([(x[3], x[:3]) for x in promoter])

for name,(chr,start,end) in promoter.items():
	p = "%s\t%s\t%s" % (chr,start,end)
	if not seen.has_key(p):
		f.write("%s\n" % p)
		seen[p] = 1
f.close()

prom_track = SimpleTrack(tmpname)
pol_track = SimpleTrack(options.polfile)
result = peak_stats.peak_stats(prom_track, pol_track, peak_stats.number_formatter)

prom_result = {}
for row in result:
	chr,start,end,num = row.split("\t")
	prom_result["%s\t%s\t%s" % (chr,start,end)] = float(num)

prom_result = dict([(n, prom_result["%s\t%s\t%s" % (chr,start,end)]) for n,(chr,start,end) in promoter.items()])

#print prom_result

#GENES
tmp = NamedTemporaryFile()
tmpname = tmp.name
f = open(tmpname, "w")
seen = {}
for name,(chr,start,end) in genes.items():
	p = "%s\t%s\t%s" % (chr,start,end)
	if not seen.has_key(p):
		f.write("%s\n" % p)
		seen[p] = 1
f.close()

exon_track = SimpleTrack(tmpname)
pol_track = SimpleTrack(options.polfile)
result = peak_stats.peak_stats(exon_track, pol_track, peak_stats.number_formatter)

gene_result = {}
gene_length = {}
for row in result:
	chr, start, end, num = row.split("\t")
	loc = "%s\t%s\t%s" % (chr,start,end)
	gene_result[loc] = float(num)

gene_result = dict([(n, gene_result["%s\t%s\t%s" % (chr,start,end)]) for n,(chr,start,end) in genes.items()])
gene_length = dict([(n, end - start) for n,(chr,start,end) in genes.items()])

#print gene_result

print "location\tgene_name\t# of tags / total_length\t# of tags in promoter\t# of tags gene body, norm to promoter\tlog2(prom/genebody)\tlength body/prom"
for gene in prom_result.keys():
	prom = prom_result[gene]
	body = 0
	body_norm = 0
	ratio = "NA"
	length = upstream
	rel_length = "NA"
	if gene_result.has_key(gene):
		body = gene_result[gene]
		body_norm = body / float(gene_length[gene]) * (upstream + downstream)
		if body_norm > 0 and prom > 0:
			ratio = log((prom)/ (body_norm))/log(2)
		else:
			body_norm = "NA"
		
		length = gene_length[gene] + upstream	+ downstream
		rel_length = gene_length[gene]  / float(upstream + downstream)
	name = gene.split("_|_")[0]
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (location[gene], name, (prom + body) /  length , prom, body_norm, ratio, rel_length )

