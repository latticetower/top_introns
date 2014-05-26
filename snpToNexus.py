#!/usr/bin/python
import fileinput
import os.path
from collections import OrderedDict
from window_info2 import WindowInfo
from optparse import OptionParser



class NexusSaver(object):
	def __init(self):
		self.windows_container = OrderedDict()
		self.current_chrom_windows = OrderedDict()
		self.chromosome_name = ""
		self.current_window = list()
	
	def load_windows(self, windows_file_name):
		with open(windows_file_name) as f:
			for st in f:
				window_line = st.split()
				if (not window_line[0] in self.windows_container.keys()):
					self.windows_container[window_line[0]] = OrderedDict()
				self.windows_container[window_line[0]][window_line[3]] = WindowInfo(*window_line[:4])
	#
	def load_snp_info(self, snp_info):
		if self.chromosome_name != snp_info[0]:
			self.chromosome_name = snp_info[0]
			if not self.chromosome_name in self.windows_container.keys():
				print "\n\n--not found {0}\n".format(self.chromosome_name)
				return
			self.current_chrom_windows = iter(sorted(self.windows_container[snp_info[0]].items(), key = lambda w: w[1].window_start))
			self.current_window = next(self.current_chrom_windows)
		window_start = self.current_window[1].window_start
		window_end = self.current_window[1].window_end

		if long(snp_info[1]) - 1 <= long(window_end):
			self.current_window[1].process_line(snp_info)
		else:
			try:
				while long(snp_info[1]) - 1 > long(window_end):
					self.current_window = next(self.current_chrom_windows)
					window_start = self.current_window[1].window_start
					window_end = self.current_window[1].window_end
				current_window[1].process_line(snp_info)
			except:
				return
	def load_from_vcf(self, vcf_file_name):
		with open(vcf_file_name) as f:
			for line in f:
				if line.startswith('#'):
					continue
				snp_info = line.strip('\n').split()
				load_snp_info(snp_info)
	#
	#
	def load_fasta(self, fasta_file_name):
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
	def process_and_save(self, fasta_file_name, output_folder_name):
		for chromosome_info in load_fasta(fasta_file_name):
			#print "chrom info: {0}".format(chromosome)
			self.windows_container[chromosome_info[0]][chromosome_info[1]].sequence = chromosome_info[2]
			self.windows_container[chromosome_info[0]][chromosome_info[1]].print_to_nexus(output_folder_name)



def main():
	parser = OptionParser()
	
	parser.add_option("", "--vcf", dest="vcf_file", help="input FILE in .vcf format", metavar="FILE")
	parser.add_option("-i", "--input", dest="input_file", help="input FILE in .bed format", metavar="FILE")
	parser.add_option("-o", "--output", dest="output_folder", help="output folder for .nexus files", metavar="FOLDER")
	parser.add_option("-f", "--fasta", dest="referencefile", help = "input FILE in .fasta format", metavar = "FILE")
	
	parser.add_option("-v", "--verbose",
										action="store_true", dest="verbose", default=False,
										help="don't print status messages to stdout")
	
	(options, args) = parser.parse_args()
	
	vcf_file_name = 'cheetah_SNP_filtered.vcf'
	exons_file_name = 'top100_variable_exons.bed'
	output_folder_name = 'results/'
	fasta_file_name = 'cheetah_exons.fasta'
	
	if options.vcf_file: vcf_file_name = options.vcf_file
	if options.input_file: exons_file_name = options.input_file
	if options.output_folder: output_folder_name = options.output_folder
	if options.referencefile: fasta_file_name = options.referencefile
	
	nexus_saver = NexusSaver()
	nexus_saver.load_windows(exons_file_name)
	nexus_saver.load_from_vcf(vcf_file_name)
	nexus_saver.process_and_save(fasta_file_name, output_folder_name)
#vcf_file = open(vcf_file_name, 'r')

if __name__ == "__main__":
	main()
