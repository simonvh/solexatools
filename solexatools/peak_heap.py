""" 
peak_heap

Author: Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>

TODO:
- better name
- get rid of mergebed dependency

"""
from solexatools.track import SimpleTrack
from solexatools import peak_stats
from tempfile import NamedTemporaryFile
import os

MERGE = "mergebed.pl"

def _make_one_peak_file(peaks, filename, merge=True):
	temp = NamedTemporaryFile()
	tempname = temp.name
	temp.close()

	# check if all files exist
	for peakfile in peaks:
		if not os.path.exists(peakfile):
			raise IOError,  "%s does not exist!" % peakfile
	
	if merge:
		# check if we can run mergebed.pl
		if os.system("%s > /dev/null 2>&1" % MERGE):
			raise Exceptionr,  "%s is requires to run peak_heap" % MERGE

	all_files = " ".join(peaks)
	
	command = "cat %s |grep -v track> %s"
	ret = os.system(command % (all_files, tempname))
	if ret:
		raise Exception,  "cat didn't work"

	if merge:
		command = "%s -i %s > %s  2>/dev/null" % (MERGE, tempname, filename)
		ret = os.system(command)
		if ret:
			raise Exception, "Something went wrong running mergebed.pl"
	else:
		os.rename(tempname, filename)

def peak_heap(peaks, data, merge=True, rpkm=False, rmdup=True, rmrepeats=True):
	""" peak_heap(peaks, data)
			peaks: list of peakfiles in BED format
			data: dictionary containging name, read files in BED format items
		
		Throws all peaks on one big heap, and then for each sample determine number of reads per peak

		returns: dictionary: for each sample name a file with number of reads per peak
	"""

	result = {}
	
	if len(data.keys()) == 0:
		raise Exception, "No samples given"

	# First make one big peak file
	peakfile = NamedTemporaryFile().name
	_make_one_peak_file(peaks, peakfile, merge)
	
	count = {}
	peak_track = SimpleTrack(peakfile)
	f = peak_track.get_next_feature()
	while f:
		count["%s:%s-%s" % (f[0],f[1],f[2])] = {}
		f = peak_track.get_next_feature()
	
	# Now run peakstats for all samples
	for sample, fname in data.items():
		#print sample, file
		peak_track = SimpleTrack(peakfile)
		
		if fname.endswith("bam"):
			result = peak_stats.bam_binned_peak_stats(peak_track, fname, 1, rpkm, rmdup, rmrepeats)
		else:
			data_track = SimpleTrack(fname)
			result = peak_stats.peak_stats(peak_track, data_track, peak_stats.number_formatter)
		for row in result:
			vals = row.strip().split("\t")
			count["%s:%s-%s" % (vals[0], vals[1], vals[2])][sample] = int(vals[3])

	return count


