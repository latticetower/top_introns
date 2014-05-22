#!/usr/bin/python
import fileinput
from collections import OrderedDict


from optparse import OptionParser
parser = OptionParser()

parser.add_option("-v", "--vcf", dest="vcf_file", help="input FILE in .vcf format", metavar="FILE")
parser.add_option("-i", "--input", dest="input_file", help="input FILE in .bed format", metavar="FILE")
parser.add_option("-o", "--output", dest="output_file", help="output FILE in .bed format", metavar="FILE")
parser.add_option("", "--normalize",
                  action="store_true", dest="normalize", default=False,
                  help="normalize snp counter by sequence length")

parser.add_option("-n", "", dest="top_n", type="int", default=100)

parser.add_option("-s", "--min-size", dest="min_size", type = "int", default=1)
(options, args) = parser.parse_args()

vcf_file_name = 'cheetah_SNP_filtered.vcf'
windows_file_name = 'cheetah_exonss.bed'

output_file_name = 'top100_variable_exons.bed'

if options.output_file: output_file_name = options.output_file

if options.vcf_file: vcf_file_name = options.vcf_file
if options.input_file: windows_file_name = options.input_file

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
windows_file = open(windows_file_name, 'r')
#load introns
windows_container = OrderedDict()
for st in windows_file:
	window_line = st.split()
	if (not window_line[0]  in windows_container.keys()): windows_container[window_line[0]] = OrderedDict()
	window_container[window_line[0]][window_line[1]] = WindowInfo(*window_line[:4])

windows_file.close

#
window_start = 0
window_end = 0
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
		current_chrom_windows = iter(window_container[list_line[0]].items())
		current_window = next(current_chrom_windows)
		window_start = current_window[1].window_start
		window_end = current_window[1].window_end

	if long(list_line[1]) - 1 <= long(WindowEnd):
		current_window[1].process_line(list_line)
	else:
		try:
			while long(list_line[1]) - 1 > long(WindowEnd):
				current_window = next(current_chrom_windows)
				window_start = current_window[1].window_start
				window_end = current_window[1].window_end
			current_window[1].process_line(list_line)
		except:
			continue
vcf_file.close()

min_size = 2000
if options.min_size: min_size = options.min_size

counter_list = [item for chrom_windows in window_container.values() for item in chrom_windows.values() if long(long(item.window_end) - long(item.window_start) + 1) >= long(min_size)]
sorted_list = sorted(counter_list, key = lambda window: window.normalized_counter() if options.normalize else window.counter, reverse = True)
sorted_by_cheetah = [sorted(counter_list, key = lambda window: window.cheetah[i], reverse = True) for i in range(7)]
#TODO: compute number of species from header line in .bed file
print len(set([x.counter for x in sorted_list]))

#
##1. output data for top n most variable
output_amount = options.top_n
filtered_list = sorted_list[: output_amount]
g = open(output_file_name, 'w')
for x in filtered_list:
	g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(x.chromosome, x.window_start, x.window_end, x.intron_name, x.counter, x.normalized_counter() if options.normalize else x.counter))
g.close()
