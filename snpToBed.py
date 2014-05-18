#!/usr/bin/python
import fileinput
from collections import OrderedDict


from optparse import OptionParser
parser = OptionParser()

#parser.add_option("-v", "--vcf", dest="vcffile", help="input FILE in .vcf format", metavar="FILE")
parser.add_option("-i", "--input", dest="inputfile", help="input FILE in .bed format", metavar="FILE")
parser.add_option("-o", "--output", dest="outputfile", help="output FILE in .bed format", metavar="FILE")
parser.add_option("-n", "--normalize",
                  action="store_true", dest="normalize", default=False,
                  help="normalize snp counter by sequence length")
parser.add_option("-s", "--min-size", dest="min_size", type = "int")
(options, args) = parser.parse_args()

vcf_file_name = 'cheetah_SNP_filtered.vcf'
introns_file_name = 'cheetah_exonss.bed'

if options.inputfile: introns_file_name = options.inputfile

output_file_name = 'top100_variable_exons.bed'
if options.outputfile: output_file_name = options.outputfile

#if options.vcffile: vcf_file_name = options.vcffile
if options.inputfile: exons_file_name = options.inputfile

#vcf_file = open(vcf_file_name, 'r')

# basic window object definition
class WindowInfo(object):
	def __init__(self, chromosome, window_start, window_end, name, cheetah_no = 7):
		self.chromosome = chromosome
		self.window_start = long(window_start)
		self.window_end = long(window_end)
		self.intron_name = name
		self.counter = 0
		self.cheetah_no = cheetah_no
		self.cheetah = [ 0 for x in range(cheetah_no)]
		self.line_infos = OrderedDict()
		self.sequence = ''
	#
	def process_line(self, line):
		if long(self.window_start) > long(line[1]) - 1 or long(self.window_end) < long(line[1]) - 1:
			#print "{0}, {1}, {2}".format(self.window_start, self.window_end, line[1])
			return
		self.counter += 1 #the same as before
		for x in range(self.cheetah_no):
			if line[- self.cheetah_no + x].startswith('0/1'):
				self.cheetah[x] += 1      
		self.line_infos[long(line[1])] = line
	def normalized_counter(self):
		return self.counter*1.0/(long(self.window_end) - long(self.window_start) + 1)
	def to_str(self):
		return "{0}\t{1}\t{2}".format(self.chromosome, self.window_start, self.window_end)  
# ----------
# main code
# ----------   
introns_file = open(introns_file_name, 'r') 
#load introns
intron_holder = OrderedDict()
for st in introns_file:
	intron_line = st.split()
	if (not intron_line[0]  in intron_holder.keys()): intron_holder[intron_line[0]] = OrderedDict()
	intron_holder[intron_line[0]][intron_line[1]] = WindowInfo(*intron_line[:4])

introns_file.close

#
WindowStart = 0
WindowEnd = 0
chromosome_name = ''

vcf_file = open(vcf_file_name, 'r')

list_line = list()
#current_chromosome_data
current_chrom_windows = OrderedDict()
current_window = list()
for i in vcf_file:
	if i.startswith('#'):
		continue
	list_line = i.strip('\n').split()
	if chromosome_name != list_line[0]:
		chromosome_name = list_line[0]
		current_chrom_windows = iter(intron_holder[list_line[0]].items())
		current_window = next(current_chrom_windows)
		WindowStart = current_window[1].window_start
		WindowEnd = current_window[1].window_end
		
	if long(list_line[1]) - 1 <= long(WindowEnd):
		current_window[1].process_line(list_line)
	else:
		try:
			while long(list_line[1]) - 1 > long(WindowEnd):
				current_window = next(current_chrom_windows)
				WindowStart = current_window[1].window_start
				WindowEnd = current_window[1].window_end
			current_window[1].process_line(list_line)
		except:
			continue
vcf_file.close()

min_size = 2000
if options.min_size: min_size = options.min_size

CounterList = [item for chrom_introns in intron_holder.values() for item in chrom_introns.values() if long(long(item.window_end) - long(item.window_start) + 1) >= long(min_size)]
print CounterList[:2]
SortedList = sorted(CounterList, key = lambda window: window.normalized_counter() if options.normalize else window.counter, reverse = True)
SortedByCheetah = [sorted(CounterList, key = lambda window: window.cheetah[i], reverse = True) for i in range(7)]

print len(set([x.counter for x in SortedList]))

#
##1. output data for top n most variable
outputAmount = 100
filtered_list = SortedList[:outputAmount]
g = open(output_file_name, 'w')
for x in filtered_list:
	g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(x.chromosome, x.window_start, x.window_end, x.intron_name, x.counter, x.normalized_counter() if options.normalize else x.counter))
g.close()
