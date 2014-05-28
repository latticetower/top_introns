#!/usr/bin/python
import fileinput
from collections import OrderedDict
from optparse import OptionParser
from window_info import WindowInfo


class TopNWindowsSelector(object):
  def __init__(self, output_amount, normalize):
    self.output_amount = output_amount
    self.windows_container = OrderedDict()
    self.chromosome_name = ""
    self.normalize = normalize
    self.species_amount = 0
    self.species_ids = list()


  def load_windows_info(self, windows_file_name):
    with open(windows_file_name) as f:
      for st in f:
        window_line = st.split()
        if (not window_line[0]  in self.windows_container.keys()):
          self.windows_container[window_line[0]] = OrderedDict()
        self.windows_container[window_line[0]][window_line[1]] = WindowInfo(*(window_line[ : 4] + [self.species_amount]))


  def process_snp(self, snp_info):
    if self.chromosome_name != snp_info[0]:
      self.chromosome_name = snp_info[0]
      self.current_chrom_windows = iter(self.windows_container[snp_info[0]].items())
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
        self.current_window[1].process_line(snp_info)
      except:
        return


  def load_snp(self, vcf_file_name):
    with open(vcf_file_name) as f:
      for line in f:
        if line.find("#CHROM") != -1:
          table_header_line = line.strip('\r\n').split()
          try:
            self.species_amount = len(table_header_line) - table_header_line.index('FORMAT') - 1
            self.species_ids = table_header_line[ - self.species_amount : ]
          except ValueError:
            self.species_amount = 0
          continue
        if line.startswith('#'):
          continue
        snp_info = line.strip('\r\n').split()
        self.process_snp(snp_info)


  def select_top_n(self, min_size):
    counter_list = [item for chrom_windows in self.windows_container.values() for item in chrom_windows.values() if long(long(item.window_end) - long(item.window_start) + 1) >= long(min_size)]
    self.filtered_list = sorted(counter_list, key = lambda window: window.normalized_counter() if self.normalize else window.counter, reverse = True)
    self.filtered_list =  self.filtered_list[ : self.output_amount]


  def save_to_bed_file(self, output_file_name):
    g = open(output_file_name, 'w')
    for x in self.filtered_list:
      g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
            x.chromosome,
            x.window_start,
            x.window_end,
            x.window_name,
            x.counter,
            x.normalized_counter() if self.normalize else x.counter
          )
        )
    g.close()


def main(vcf_file_name, windows_file_name, output_file_name, min_size, output_amount, normalize):
  result = 0
  if not vcf_file_name: result = 1
  if not windows_file_name: result = 1
  if result:
    return result
  selector = TopNWindowsSelector(output_amount, normalize)
  selector.load_windows_info(windows_file_name)
  selector.load_snp(vcf_file_name)
  selector.select_top_n(min_size)
  selector.save_to_bed_file(output_file_name)
  return 0


if __name__ == "__main__":
  parser = OptionParser()

  parser.add_option("-v", "--vcf", dest = "vcf_file", help = "input FILE in .vcf format", metavar = "FILE")
  parser.add_option("-i", "--input", dest = "input_file", help = "input FILE in .bed format", metavar = "FILE")
  parser.add_option("-o", "--output", dest = "output_file", help = "output FILE in .bed format", metavar = "FILE")
  parser.add_option("", "--normalize", action = "store_true", dest = "normalize", default = False,
    help = "normalize snp counter by sequence length")
  parser.add_option("-n", "", dest = "top_n", type = "int", default = 100)
  parser.add_option("-s", "--min-size", dest = "min_size", type = "int", default = 1)
  (options, args) = parser.parse_args()
  try:
    result = main(options.vcf_file, options.input_file, options.output_file,
      options.min_size, options.top_n, options.normalize)
    if result:
      parser.print_help()
      exit()
  except Exception as e:
    print e
