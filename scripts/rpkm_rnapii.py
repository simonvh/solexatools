#!/usr/bin/env python
import sys
from optparse import OptionParser
from solexatools import peak_stats
from solexatools.track import SimpleTrack
from os.path import basename,splitext
from subprocess import Popen,PIPE

parser = OptionParser()
parser.add_option("-f", "--featurefile", dest="featurefile", help="Features (genes) in  BED format", metavar="FILE")
parser.add_option("-d", "--datafile", dest="datafile", help="Data in BED format", metavar="FILE")

(options, args) = parser.parse_args()

if not options.featurefile or not options.datafile:
	parser.print_help()
	sys.exit(0)

featurefile = options.featurefile
datafile = options.datafile

p = Popen("cat %s | wc -l" % datafile, shell=True, stdout=PIPE)
total_reads = float(p.communicate()[0].strip())

ids = {}
features = SimpleTrack(featurefile)
f = features.get_next_feature()
while f:
	(chrom, start, end, val, strand) = f
	ids["%s:%s-%s" % (chrom, start, end)] = val
	f = features.get_next_feature()
	
features = SimpleTrack(featurefile)
data = SimpleTrack(datafile)
formatter = peak_stats.number_formatter
result = peak_stats.peak_stats(features, data, formatter)

print "ID\tchrom\tstart\tend\tnumber of reads\tRPKM"
for row in result:
	(chrom,start,end,num_reads) = row.split("\t")
	feature = "%s:%s-%s" % (chrom, start, end)
	start, end, num_reads = float(start), float(end), float(num_reads)
	rpkm = 1e9 * num_reads / (total_reads * (end - start))
	print "%s\t%s\t%d\t%d\t%d\t%s" % (ids[feature], chrom, start, end, num_reads, rpkm)


