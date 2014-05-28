#!/usr/bin/python
import fileinput
import os.path
import sys
from collections import OrderedDict
from window_info import WindowInfo
from optparse import OptionParser

class NexusSaver(object):
  def __init__(self):
    self.windows_container = OrderedDict()
    self.current_chrom_windows = OrderedDict()
    self.chromosome_name = ""
    self.current_window = list()
    self.species_amount = 0
    self.species_ids = list()
    self.window_start = 0
    self.window_end = 0


  def load_windows(self, windows_file_name):
    with open(windows_file_name) as f:
      for st in f:
        window_line = st.split()
        if (not window_line[0] in self.windows_container.keys()):
          self.windows_container[window_line[0]] = OrderedDict()
        self.windows_container[window_line[0]][window_line[3]] = WindowInfo(*window_line[ : 4])


  def load_snp_info(self, snp_info):
    if self.chromosome_name != snp_info[0]:
      self.chromosome_name = snp_info[0]
      if not self.chromosome_name in self.windows_container.keys():
        return
      self.current_chrom_windows = iter(sorted(self.windows_container[snp_info[0]].items(), key = lambda w: w[1].window_start))
      self.current_window = next(self.current_chrom_windows)
      self.window_start = self.current_window[1].window_start
      self.window_end = self.current_window[1].window_end

    if long(snp_info[1]) - 1 <= long(self.window_end):
      self.current_window[1].process_line(snp_info)
    else:
      try:
        while long(snp_info[1]) - 1 > long(self.window_end):
          self.current_window = next(self.current_chrom_windows)
          self.window_start = self.current_window[1].window_start
          self.window_end = self.current_window[1].window_end
        current_window[1].process_line(snp_info)
      except:
        return


  def load_from_vcf(self, vcf_file_name):
    with open(vcf_file_name) as f:
      for line in f:
        if line.find("#CHROM") != -1:
          table_header_line = line.strip('\r\n').split()
          try:
            self.species_amount = len(table_header_line) - table_header_line.index('FORMAT') - 1
            self.species_ids = table_header_line[- self.species_amount : ]
          except ValueError:
            self.species_amount = 0
          continue
        if line.startswith('#'):
          continue
        snp_info = line.strip('\r\n').split()
        self.load_snp_info(snp_info)


  def load_fasta(self, fasta_file_name):
    chromosome_name = ''
    window_name = ''
    with open(fasta_file_name) as f:
      for line in f:
        if line.startswith('>'):
          [chromosome_name, window_name] = line[1:].strip('\r\n').split()[:2]
          print "...Loading new chromosome {0} {1}".format(chromosome_name, window_name)
        else:
          buffer = line.strip('\r\n')
          yield(chromosome_name, window_name, buffer)


  def process_and_save(self, fasta_file_name, output_folder_name):
    for chromosome_info in self.load_fasta(fasta_file_name):
      self.windows_container[chromosome_info[0]][chromosome_info[1]].sequence = chromosome_info[2]
      self.windows_container[chromosome_info[0]][chromosome_info[1]].print_to_nexus(output_folder_name, self.species_amount, self.species_ids)


def main(vcf_file_name, windows_file_name, fasta_file_name, output_folder_name):
  result = 0
  if not vcf_file_name:
    result = 1
  if not windows_file_name:
    result = 1
  if not fasta_file_name:
    result = 1
  if result:
    return result
  nexus_saver = NexusSaver()
  nexus_saver.load_windows(windows_file_name)
  nexus_saver.load_from_vcf(vcf_file_name)
  nexus_saver.process_and_save(fasta_file_name, output_folder_name)
  return result

if __name__ == "__main__":
  parser = OptionParser()

  parser.add_option("", "--vcf", dest = "vcf_file", help = "input FILE in .vcf format", metavar = "FILE")
  parser.add_option("-i", "--input", dest = "input_file", help = "input FILE in .bed format", metavar = "FILE")
  parser.add_option("-o", "--output", dest = "output_folder", help = "output folder for .nexus files", metavar = "FOLDER")
  parser.add_option("-f", "--fasta", dest = "referencefile", help = "input FILE in .fasta format", metavar = "FILE")

  (options, args) = parser.parse_args()

  try:
    result = main(options.vcf_file, options.input_file, options.referencefile, options.output_folder)
    if result:
      parser.print_help()
      exit()
  except Exception as e:
    print >> sys.stderr, e
    exit(1)
