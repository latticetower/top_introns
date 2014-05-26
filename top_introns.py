#!/usr/bin/python
#the main script. should run all subscripts
import subprocess
from data_processor import DataProcessor
from optparse import OptionParser

parser = OptionParser()

#input file options:
parser.add_option("-v", "--vcf", dest="vcf_file_name",
  help="input FILE in .vcf format", metavar="FILE")
parser.add_option("-f", "--fasta", dest="fasta_file_name",
  help="input FILE with reference in .fasta format", metavar="FILE")
parser.add_option("-w", "--windowfile", dest="window_file_name",
  help="input FILE in .bed format", metavar="FILE")

#output file:
parser.add_option("-o", "--output", dest="output_folder",
  help = "output FOLDER for saving the results", metavar = "FOLDER")
#other parameters:
#TODO: decide what to parametrize!
# 1. number of top scored values
parser.add_option("-n", "", dest="top_n", type="int", default=100)
# 2. use normalization or not. by default use
parser.add_option("", "--normalize",
                  action="store_true", dest="normalize", default=False,
                  help="normalize SNP count by window width")
# 3. threshhold - minimum window size
parser.add_option("-s", "--min-size", dest="min_size", type="int", default=1)

(options, args) = parser.parse_args()

vcf_file_name = options.vcf_file_name
fasta_file_name = options.fasta_file_name
window_file_name = options.window_file_name
output_folder_name = options.output_folder

#try:
	#
DataProcessor(window_file_name,
		vcf_file_name,
		fasta_file_name,
		output_folder_name,
		options.top_n, 
		options.min_size,
		options.normalize).process()
print("--------\ndata processed and can be found in folder specified")
#except:
	#print "Got some error while processing input data. Look at README.md and check files you've provided, or call with --help key"
