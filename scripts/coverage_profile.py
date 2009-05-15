#!/usr/bin/env python
from numpy import *
import sys
from optparse import OptionParser,TitledHelpFormatter
from solexatools.track import *

VERSION = "1.0"

DEFAULT_POSITION = "middle"

usage = "usage: %prog -f <FILE> -d <FILE> -b <SIZE> -w <LENGTH> [-p POS]"
parser = OptionParser(version=VERSION, usage=usage, formatter=TitledHelpFormatter(max_help_position=40, short_first=1))
parser.add_option("-f", "--featurefile", dest="feature_file", help="GFF formatted file of features (NO BED format!)", metavar="FILE")
parser.add_option("-d", "--datafile", dest="data_file", help="BED formatted datafile (NO WIG or GFF format!)", metavar="FILE")
parser.add_option("-w", "--window", dest="window", help="Window length", metavar="LENGTH")
parser.add_option("-b", "--bin", dest="bin", help="Bin size", metavar="SIZE")
parser.add_option("-p", "--position", dest="position", help="Position of feature to consider: start, end or middle (strand is taken into account, default is middle)", metavar="POS", default=DEFAULT_POSITION)

(options, args) = parser.parse_args()

if not options.feature_file or not options.data_file or not options.window or not options.bin:
	parser.print_help()
	sys.exit(0)

data_file = options.data_file
feature_file = options.feature_file
bin = int(options.bin)
window = int(options.window)
position = options.position	

# Speed up reading of large files
BUFSIZE = 10000000

# Dict to store all features per chromosome
pos = {}

# Read file of data points (only start & end are used)
t = SimpleTrack(data_file)
f = t.get_next_feature()
while f:
	if not(pos.has_key(f[0])):
		pos[f[0]] = []
	pos[f[0]].append([f[1],f[2] + 1])
	f = t.get_next_feature()

# Transform lists to Numpy arrays
for chr in pos.keys():
	pos[chr] = array(pos[chr])

# Array containing all the bins
bins = arange(-window, window, bin)

# Read the whole file at once
f = open(feature_file)
lines = f.readlines()
f.close()

# Array to store results
result = array(zeros((len(lines), len(bins))), dtype=int32)

c = 0
for line in lines:
	vals = line[:-1].split("\t")
	strand = vals[6]
	if strand == "-1":
		strand = "-"
	chr = vals[0]
	middle = 0
	if (position == "start" and strand == "-") or (position == "end" and strand == "+"):
		middle = int(vals[4])
	elif (position == "end" and strand == "-") or (position == "start" and strand == "+"):
		middle = int(vals[3])
	else:
		middle = (int(vals[4]) + int(vals[3])) / 2
	
	# Get all data points within window range of the feature to minimize amount of calculations (and use less memory)
	if not(pos.has_key(chr)):
		in_range = []
	else:
		in_range = pos[chr][pos[chr][:,1] >= middle - window]
	if len(in_range) > 0:
		in_range = in_range[in_range[:,0] <= middle + window]

	if len(in_range) > 0:
		# Define matrices
		(f_start, b_start) = ix_(in_range[:,0], bins + middle)
		(f_end, b_end) = ix_(in_range[:,1], bins + middle + bin)
	
		# Calculate overlap
		overlap = maximum(f_end, b_end) -  minimum(f_start, b_start) - abs(b_start - f_start) - abs(b_end - f_end)
		overlap[overlap < 0] = 0
	
		# Take care of minus strands
		if vals[6] == "-" or vals[6] == "-1":
			result[c] = overlap.sum(0)[::-1]
		else:
			result[c] = overlap.sum(0)
		c += 1
	else:
		result[c] = array(zeros((len(bins))))

pos = -window

for i in range(len(result[0])):
	d = result[:,i]
	#print "%s\t%0.2f\t%0.2f\t%0.2f\t%s" % (pos, mean(d), std(d), median(d), len(d))
	print "%s\t%0.2e" % (pos, mean(d))
	pos += bin

#for mean_val, median_val in zip(result.mean(0), median(result)):
#	print "%s\t%si\t%s" % (pos, mean_val, median_val)
#	pos += bin


