#!/usr/bin/python
import fileinput
import os.path
from collections import OrderedDict

from optparse import OptionParser
parser = OptionParser()

#parser.add_option("-v", "--vcf", dest="vcffile", help="input FILE in .vcf format", metavar="FILE")
parser.add_option("-i", "--input", dest="inputfile", help="input FILE in .bed format", metavar="FILE")
parser.add_option("-o", "--output", dest="outputfolder", help="output folder for .nexus files", metavar="FOLDER")
parser.add_option("-f", "--fasta", dest="referencefile", help = "input FILE in .fasta format", metavar = "FILE")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

vcf_file_name = 'cheetah_SNP_filtered.vcf'
exons_file_name = 'top100_variable_exons.bed'
output_folder_name = 'results/'
fasta_file_name = 'cheetah_exons.fasta'

#if options.vcffile: vcf_file_name = options.vcffile
if options.inputfile: exons_file_name = options.inputfile
if options.outputfolder: output_folder_name = options.outputfolder
if options.referencefile: fasta_file_name = options.referencefile

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
	#
	def to_str(self):
		return "{0}\t{1}\t{2}".format(self.chromosome, self.window_start, self.window_end)  
	#
	def getIUPAC(self, letters):
		if len(letters) == 1:
			return letters[0]
		if len(set(letters) & set(['A', 'C', 'T', 'G'])) == 4:
			return 'N'
		if len(set(letters) & set(['C', 'T', 'G'])) == 3:
			return 'B'
		if len(set(letters) & set(['A', 'T', 'G'])) == 3:
			return 'D'
		if len(set(letters) & set(['A', 'T', 'C'])) == 3:
			return 'H'
		if len(set(letters) & set(['A', 'C', 'G'])) == 3:
			return 'V'
		if len(set(letters) & set(['A', 'G'])) == 2:
			return 'R'
		if len(set(letters) & set(['C', 'T'])) == 2:
			return 'Y'
		if len(set(letters) & set(['A', 'C'])) == 2:
			return 'M'
		if len(set(letters) & set(['T', 'G'])) == 2:
			return 'K'
		if len(set(letters) & set(['T', 'A'])) == 2:
			return 'W'
		if len(set(letters) & set(['C', 'G'])) == 2:
			return 'S'
	def getFromIUPAC(self, letter):  
		if letter == 'R': return ['A', 'G']
		if letter == 'Y': return ['C', 'T']
		if letter == 'M': return ['A', 'C']
		if letter == 'K': return ['T', 'G']
		if letter == 'W': return ['T', 'A']
		if letter == 'S': return ['C', 'G']
		if letter == 'B': return ['C', 'T', 'G']
		if letter == 'D': return ['A', 'T', 'G']
		if letter == 'H': return ['A', 'T', 'C']
		if letter == 'V': return ['A', 'C', 'G']
		if letter == 'N': return ['A', 'C', 'T', 'G']
		return [letter]
	#
	def sequence_for_cheetah(self, i):
		buffer = ''.join(self.sequence)
		for x in self.line_infos.items():
			#print x[1]
			[chromosome, pos, _ , ref, alt] = x[1][:5]
			index = long(pos) - 1#convert to zero-based
			index -= long(self.window_start)
			info = x[1][-7 + i]
			letters = list() 
			if info.startswith('0/0'):
				letters = self.getFromIUPAC(ref)
			else: 
				if info.startswith('1/1'):
					letters = self.getFromIUPAC(alt)
				else:
					letters = self.getFromIUPAC(alt) + self.getFromIUPAC(ref)
			if options.verbose: print "got letters for {3} - {4} - {5}: {0} {1}, {2}".format(letters, self.getIUPAC(letters), x, index, self.window_start, i)
			buffer = buffer[:index] + self.getIUPAC(letters) + buffer[index + 1:]
		if options.verbose: print "---\n result {0} string: {0}\n".format(self.chromosome, buffer)
		return buffer
	#
	def isoverlap(self, start, end):
		return (self.window_start <= start and start <= self.window_end) or (self.window_start <= end and end <= self.window_end)
	def print_to_nexus(self, output_folder_name):
		output_file_name = "{1}.nexus".format(output_folder_name, self.intron_name)
		dir_path = os.path.join(output_folder_name)
		if not os.path.isdir(dir_path): os.makedirs(dir_path)
		g = open(os.path.join(dir_path, output_file_name), 'w')
		g.write("#NEXUS\n")
		g.write("[ Title Phylogenetic Analysis]\nbegin data;\n")
		g.write("   dimensions ntax=7 nchar={0};\n".format(len(self.sequence)))
		g.write("   format missing=? gap=- matchchar=. datatype=dna interleave=yes;\n")
		g.write("   matrix\n\n")
		g.write("[!Domain=Data;]\n")
		for i in range(7):
			g.write("Cheetah{0} {1}\n".format(i + 1, self.sequence_for_cheetah(i)))
		g.write(";\nEnd;")
		g.close()#for i in range(7):
# ----------
# main code
# ----------   
exons_file = open(exons_file_name, 'r') 
#load introns
exon_holder = OrderedDict()
for st in exons_file:
	exon_line = st.split()
	if (not exon_line[0]  in exon_holder.keys()): exon_holder[exon_line[0]] = OrderedDict()
	exon_holder[exon_line[0]][exon_line[3]] = WindowInfo(*exon_line[:4])

exons_file.close

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
		if not chromosome_name in exon_holder.keys():
			print "\n\n--not found {0}\n".format(chromosome_name)
			continue
		current_chrom_windows = iter(sorted(exon_holder[list_line[0]].items(), key = lambda w: w[1].window_start))
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

#
def load_fasta(fasta_file_name):
	chromosome_name = ''
	exon_name = ''
	with open(fasta_file_name) as f:
		for line in f:
			if line.startswith('>'):
				[chromosome_name, exon_name] = line[1:].strip('\n').split()[:2]
				print "Loading new chromosome {0} {1}".format(chromosome_name, exon_name)
			else:
				buffer = line.strip('\n')
				yield(chromosome_name, exon_name, buffer)
				#print len(buffer)
		# 
	#

##1. output data for top n most variable
fasta_lines = load_fasta(fasta_file_name)
for chromosome in fasta_lines:
	#print "chrom info: {0}".format(chromosome)
	exon_holder[chromosome[0]][chromosome[1]].sequence = chromosome[2]
	exon_holder[chromosome[0]][chromosome[1]].print_to_nexus(output_folder_name)