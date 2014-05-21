import fileinput
from collections import OrderedDict

fasta_file_name = 'cheetah_genome.fasta'#'test_input.fasta'
output_file_name = 'cheetah_exons.fasta' #'test.fasta'#
bed_file_name = 'top100_variable_exons.bed'#'top_exons.bed'#

from optparse import OptionParser
parser = OptionParser()

#parser.add_option("-v", "--vcf", dest="vcffile", help="input FILE in .vcf format", metavar="FILE")
parser.add_option("-i", "--input", dest="inputfile", help="input FILE in .bed format", metavar="FILE")
parser.add_option("-o", "--output", dest="outputfile", help="output file in fasta format", metavar="FILE")
parser.add_option("-f", "--fasta", dest="referencefile", help = "input FILE in .fasta format", metavar = "FILE")

(options, args) = parser.parse_args()

#if options.vcffile: vcf_file_name = options.vcffile
if options.inputfile: bed_file_name = options.inputfile
if options.outputfile: output_file_name = options.outputfile
if options.referencefile: fasta_file_name = options.referencefile


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
			return
		self.counter += 1 #the same as before
		for x in range(self.cheetah_no):
			if line[- self.cheetah_no + x].startswith('0/1'):
				self.cheetah[x] += 1      
		self.line_infos[long(line[1])] = line
	#
	def to_str(self):
		return "{0}\t{1}\t{2}".format(self.chromosome, self.window_start, self.window_end)  
	#
	def isoverlap(self, start, end):
		return (self.window_start <= start and start <= self.window_end) or (self.window_start <= end and end <= self.window_end)
	def print_to_nexus(self, output_folder_name):
		output_file_name = "{0}{1}.nexus".format(output_folder_name, self.intron_name)
		g = open(output_file_name, 'w')
		g.write("#NEXUS\n")
		g.write("Begin data;\n")
		g.write("Dimentions ntax = 8, nchar = {0};".format(len(self.sequence)))
		g.write("Format datatype=dna symbols=\"ACTG\" missing=? gap=-;")
		g.write("Matrix")
		g.write("Total {0}".format(self.sequence))
		g.close()#for i in range(7):
	def print_to_fasta(self):
		return ">{0}\t{1}\n{2}\n".format(self.chromosome, self.intron_name, self.sequence)
# ----------


def load_fasta(fasta_file_name):
	chromosome_name = ''
	buffer = ''
	offset = 0 # zero-based offset for current string  
	with open(fasta_file_name) as f:
		for line in f:
			if line.startswith('>'):
				buffer = ''
				offset = 0
				chromosome_name = line.strip('\n')[1:]
				print "Processing new chromosome {0}".format(chromosome_name)
			else:
				buffer = line.strip('\n')
				yield(chromosome_name, buffer, offset)
				offset += len(buffer)
		# 
	#
	
output_amount = 100
filtered_list = OrderedDict()
with open(bed_file_name) as f:
	for line in f:
		line_data = line.strip('\n').split()
		if (not line_data[0]  in filtered_list.keys()): filtered_list[line_data[0]] = OrderedDict()
		filtered_list[line_data[0]][line_data[1]] = WindowInfo(*line_data[:4])

current_chromosome = ''
chromosome_windows = list()
current_window = 0
output_file = open(output_file_name, 'w')

for chromosome in load_fasta(fasta_file_name):
	if current_chromosome != chromosome[0]: 
		current_chromosome = chromosome[0]
		if (not current_chromosome in filtered_list.keys()):
			current_chromosome = ''
			continue
		#if current_chromosome=='A2': break
		print "{0}".format(chromosome)
		chromosome_windows = iter(sorted(filtered_list[current_chromosome].items(), key= lambda w: long(w[0])))
		current_window = next(chromosome_windows)
	#break
	offset = long(chromosome[2])
	offset_end = long(offset + len(chromosome[1]) - 1)
	if long(offset_end) < long(current_window[1].window_start): 
		continue
	if long(current_window[1].window_start) < long(offset) and long(current_window[1].window_end) > long(offset_end):
		current_window[1].sequence += chromosome[1]
		continue
	if long(offset) <= long(current_window[1].window_end) and long(current_window[1].window_end) <= long(offset_end) and long(current_window[1].window_start < offset):
		current_window[1].sequence += chromosome[1][ : current_window[1].window_end - offset + 1]
		output_file.write(current_window[1].print_to_fasta())      
		try:
			#while long(current_window[1].window_end) < long(offset): 
			current_window = next(chromosome_windows)
		except:
			continue#
	try:
		while long(current_window[1].window_start) >= long(offset) and long(current_window[1].window_end) <= long(offset_end):
			current_window[1].sequence = chromosome[1][current_window[1].window_start - offset : current_window[1].window_end - offset + 1]
			output_file.write(current_window[1].print_to_fasta())
			#while long(current_window[1].window_end) < long(offset): 
			current_window = next(chromosome_windows)
	except:
		continue#
	if long(offset) <= long(current_window[1].window_start) and long(current_window[1].window_start) <= long(offset_end) and long(current_window[1].window_end) > long(offset_end):
		current_window[1].sequence = chromosome[1][current_window[1].window_start - offset : ]
		continue
	if long(offset) > long(current_window[1].window_end):
		try:
			while long(current_window[1].window_end) < long(offset): 
				current_window = next(chromosome_windows)
		except:
			continue
	#break
	
output_file.close()