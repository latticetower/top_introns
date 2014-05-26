import fileinput
from collections import OrderedDict
from window_info import WindowInfo
from optparse import OptionParser

class WindowToFastaConverter(object):
	def __init__(self, output_amount):
		self.output_amount = output_amount
		self.current_chromosome = ""
		
	def load_fasta(self, fasta_file_name):
		chromosome_name = ''
		buffer = ''
		offset = 0 # zero-based offset for current string
		with open(fasta_file_name) as f:
			for line in f:
				if line.startswith('>'):
					buffer = ''
					offset = 0
					chromosome_name = line.strip('\n')[1:]
					#if options.verbose: print("Processing new chromosome {0}".format(chromosome_name))
				else:
					buffer = line.strip('\n')
					yield(chromosome_name, buffer, offset)
					offset += len(buffer)
	def load_window_list(self, file_name):
		self.filtered_list = OrderedDict()
		with open(file_name) as f:
			for line in f:
				line_data = line.strip('\n').split()
				if (not line_data[0]  in self.filtered_list.keys()): self.filtered_list[line_data[0]] = OrderedDict()
				self.filtered_list[line_data[0]][line_data[1]] = WindowInfo(*line_data[:4])
	#the following methos aims to fill in all windows information for each of chromosomes
	def fill_sequence_info(self, chromosome_info):
		if self.current_chromosome != chromosome_info[0]:
			self.current_chromosome = chromosome_info[0]
			if (not self.current_chromosome in self.filtered_list.keys()): return
			self.chromosome_windows = iter(sorted(self.filtered_list[self.current_chromosome].items(),key = lambda w: long(w[0])))
			self.current_window = next(self.chromosome_windows)
		
		#break
		offset = long(chromosome_info[2])
		offset_end = long(offset + len(chromosome_info[1]) - 1)
		if long(offset_end) < long(self.current_window[1].window_start):
			return
		if long(self.current_window[1].window_start) < long(offset) and long(self.current_window[1].window_end) > long(offset_end):
			self.current_window[1].sequence += chromosome_info[1]
			return
		if long(offset) <= long(self.current_window[1].window_end) and long(self.current_window[1].window_end) <= long(offset_end) and long(self.current_window[1].window_start < offset):
			self.current_window[1].sequence += chromosome_info[1][ : self.current_window[1].window_end - offset + 1]
			output_file.write(self.current_window[1].print_to_fasta())
			try:
				#while long(current_window[1].window_end) < long(offset):
				self.current_window = next(self.chromosome_windows)
			except StopIteration:
				return#
		try:
			while long(current_window[1].window_start) >= long(offset) and long(current_window[1].window_end) <= long(offset_end):
				self.current_window[1].sequence = chromosome[1][self.current_window[1].window_start - offset : self.current_window[1].window_end - offset + 1]
				output_file.write(self.current_window[1].print_to_fasta())
				#while long(current_window[1].window_end) < long(offset):
				self.current_window = next(self.chromosome_windows)
		except StopIteration:
			return#
		if long(offset) <= long(self.current_window[1].window_start) and long(self.current_window[1].window_start) <= long(offset_end) and long(self.current_window[1].window_end) > long(offset_end):
			self.current_window[1].sequence = chromosome_info[1][self.current_window[1].window_start - offset : ]
			return
		if long(offset) > long(self.current_window[1].window_end):
			try:
				while long(self.current_window[1].window_end) < long(offset):
					self.current_window = next(self.chromosome_windows)
			except:
				return
	#the following method processes entire fasta file
	def load_info_from_fasta(self, fasta_file_name):
		for chromosome_info in self.load_fasta(fasta_file_name):
			self.fill_sequence_info(chromosome_info)
	def save_sequences_to_fasta(self, output_file_name):
		output_file = open(output_file_name, 'w')
		# for each chromosome in list print sequence to output file
		for chromosome_windows in self.filtered_list.values():
			for window in chromosome_windows.values():
				output_file.write(window.print_to_fasta())
		output_file.close()
	#

def main():
	fasta_file_name = 'cheetah_genome.fasta'#'test_input.fasta'
	output_file_name = 'cheetah_exons.fasta' #'test.fasta'#
	bed_file_name = 'top100_variable_exons.bed'#'top_exons.bed'#
	parser = OptionParser()
	
	parser.add_option("-i", "--input", dest="input_file", help="input FILE in .bed format", metavar="FILE")
	parser.add_option("-o", "--output", dest="output_file", help="output file in fasta format", metavar="FILE")
	parser.add_option("-f", "--fasta", dest="reference_file", help = "input FILE in .fasta format", metavar = "FILE")
	
	(options, args) = parser.parse_args()
	if options.input_file: bed_file_name = options.input_file
	if options.output_file: output_file_name = options.output_file
	if options.reference_file: fasta_file_name = options.reference_file
	
	converter = WindowToFastaConverter(100)
	converter.load_window_list(bed_file_name)
	converter.load_info_from_fasta(fasta_file_name)
	converter.save_sequences_to_fasta(output_file_name)
	
if __name__ == "__main__":
	main()
