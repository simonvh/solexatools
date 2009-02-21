import sys
from struct import *
from time import sleep

data = []
chr = []
c = 0.0
for line in open(sys.argv[1]):
	vals = line.strip().split("\t")
	if vals[0] not in chr:
		chr.append(vals[0])
	data.append(pack("III", chr.index(vals[0]), int(vals[1]), int(vals[2])))
	c += 1.0
	#data.append((vals[0], int(vals[1]), int(vals[2]), c))
	
print "Done"
while 1:
	sleep(1)
