import fileinput
import sys
from collections import OrderedDict
from window_info import WindowInfo
from optparse import OptionParser

class WindowToFastaConverter(object):
  def __init__(self):
    self.current_chromosome = ""
    self.chromosome_windows = list()
    self.current_window = 0


  def load_fasta(self, fasta_file_name):
    chromosome_name = ''
    buffer = ''
    offset = 0 # zero-based offset for current string
    with open(fasta_file_name) as f:
      for line in f:
        if line.startswith('>'):
          buffer = ''
          offset = 0
          chromosome_name = line.strip('\r\n')[1 : ]
        else:
          buffer = line.strip('\r\n')
          yield(chromosome_name, buffer, offset)
          offset += len(buffer)


  def load_window_list(self, file_name):
    self.filtered_list = OrderedDict()
    with open(file_name) as f:
      for line in f:
        line_data = line.strip('\r\n').split()
        if (not line_data[0]  in self.filtered_list.keys()): self.filtered_list[line_data[0]] = OrderedDict()
        self.filtered_list[line_data[0]][line_data[1]] = WindowInfo(*line_data[ : 4])


  #the following methos aims to fill in all windows information for each of chromosomes
  def fill_sequence_info(self, chromosome_info):
    if self.current_chromosome != chromosome_info[0]:
      self.current_chromosome = chromosome_info[0]
      if (not self.current_chromosome in self.filtered_list.keys()):
        self.current_chromosome = ''
        return
      self.chromosome_windows = iter(sorted(self.filtered_list[self.current_chromosome].items(), key = lambda w: long(w[0])))
      self.current_window = next(self.chromosome_windows)

    offset = long(chromosome_info[2])
    offset_end = long(offset + len(chromosome_info[1]) - 1)
    if long(offset_end) < long(self.current_window[1].window_start):
      return
    if long(self.current_window[1].window_start) < long(offset) and long(self.current_window[1].window_end) > long(offset_end):
      self.current_window[1].sequence += chromosome_info[1]
      return
    if long(offset) <= long(self.current_window[1].window_end) and long(self.current_window[1].window_end) <= long(offset_end) and long(self.current_window[1].window_start < offset):
      self.current_window[1].sequence += chromosome_info[1][ : self.current_window[1].window_end - offset + 1]
      try:
        self.current_window = next(self.chromosome_windows)
      except StopIteration:
        return
    try:
      while long(self.current_window[1].window_start) >= long(offset) and long(self.current_window[1].window_end) <= long(offset_end):
        self.current_window[1].sequence = chromosome_info[1][self.current_window[1].window_start - offset : self.current_window[1].window_end - offset + 1]
        self.current_window = next(self.chromosome_windows)
    except StopIteration:
      return
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


def main(bed_file_name, output_file_name, fasta_file_name):
  result = 0
  if not bed_file_name:
    print ".bed file was not provided"
    result = 1
  if not fasta_file_name:
    print ".fasta file name was not provided"
    result = 1
  if result:
    return result
  converter = WindowToFastaConverter()
  converter.load_window_list(bed_file_name)
  converter.load_info_from_fasta(fasta_file_name)
  converter.save_sequences_to_fasta(output_file_name)
  return result


if __name__ == "__main__":
  parser = OptionParser()

  parser.add_option("-i", "--input", dest = "input_file", help = "input FILE in .bed format", metavar = "FILE")
  parser.add_option("-o", "--output", dest = "output_file", help = "output file in fasta format", metavar = "FILE")
  parser.add_option("-f", "--fasta", dest = "reference_file", help = "input FILE in .fasta format", metavar = "FILE")
  (options, args) = parser.parse_args()

  try:
    result = main(options.input_file, options.output_file, options.reference_file)
    if result:
      parser.print_help()
      exit()
  except Exception as e:
    print >> sys.stderr, e
  #
