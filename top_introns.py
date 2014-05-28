#!/usr/bin/env python
#the main script. should run all subscripts
import subprocess
from data_processor import DataProcessor
from optparse import OptionParser
import sys


def main(vcf_file_name, fasta_file_name, window_file_name, output_folder,
          top_n, normalize, min_size):
  result = 0
  if not vcf_file_name:
    print >> sys.stderr, "vcf file was not provided"
    result = 1
  if not fasta_file_name:
    print >> sys.stderr, "file with reference genome was not provided"
    result = 1
  if not window_file_name:
    print >> sys.stderr, "file with window info was not provided"
    result = 1
  if result:
    return result
  try:
    DataProcessor(window_file_name,
      vcf_file_name,
      fasta_file_name,
      output_folder,
      top_n,
      min_size,
      normalize).process()
    print("--------\ndata processed and can be found in folder specified")
  except Exception as e:
    print "Got error: ", e
    result = 1
  return result


if __name__ == "__main__":
  parser = OptionParser()

  #input file options:
  parser.add_option("-v", "--vcf", dest = "vcf_file_name",
    help = "input FILE in .vcf format", metavar = "FILE")
  parser.add_option("-f", "--fasta", dest = "fasta_file_name",
    help = "input FILE with reference in .fasta format", metavar = "FILE")
  parser.add_option("-w", "--windowfile", dest = "window_file_name",
    help = "input FILE in .bed format", metavar = "FILE")

  #output file:
  parser.add_option("-o", "--output", dest = "output_folder",
    help = "output FOLDER for saving the results", metavar = "FOLDER")
  #other parameters:
  # 1. number of top scored values
  parser.add_option("-n", "", dest = "top_n", type = "int", default = 100)
  # 2. use normalization or not. by default use
  parser.add_option("", "--normalize", action = "store_true", dest = "normalize", default = False,
    help = "normalize SNP count by window width")
  # 3. threshhold - minimum window size
  parser.add_option("-s", "--min-size", dest = "min_size", type = "int", default = 1)
  (options, args) = parser.parse_args()

  result = main(options.vcf_file_name,
    options.fasta_file_name,
    options.window_file_name,
    options.output_folder,
    options.top_n,
    options.normalize,
    options.min_size)

  if result:
    parser.print_help()
