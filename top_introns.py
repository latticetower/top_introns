#!/usr/bin/python
#the main script. should run all subscripts
import subprocess

from optparse import OptionParser
parser = OptionParser()

#input file options:
parser.add_option("-v", "--vcf", dest="vcf_file_name",
  help="input FILE in .vcf format", metavar="FILE")
parser.add_option("-f", "--fasta", dest="fasta_file_name",
  help="input FILE with reference in .fasta format", metavar="FILE")
parser.add_option("-w", "--windowfile", dest="windows_file",
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
                  action="store_true", dest="normalize", default=True,
                  help="normalize SNP count by window width")

# 3. threshhold - minimum window size
parser.add_option("-s", "--min-size", dest="min_size", type="int", default=1)

# TODO: (4. if there is no window file, use fixed value as window width) - additional
# 5. verbose output or not (by default output is silent)
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

vcf_file_name = options.vcf_file_name
fasta_file_name = options.fasta_file_name
window_file_name = options.window_file_name
output_folder_name = options.output_folder

#ensure output folder existance and create if nessesary
dir_path = os.path.join(output_folder_name)
if not os.path.isdir(dir_path): os.makedirs(dir_path)
bed_file_name_1 = os.path.join(dir_path, "results_top{0}.bed".format(options.top_n))
fasta_file_name_1 = os.path.join(dir_path, "window_top{0}.fasta".format(options.top_n))

windows_filter = popen("snpToBed.py --vcf={0} --input={1} --output={2} -n {3}".format(
        vcf_file_name,
        window_file_name,
        bed_file1_name,
        options.top_n
      ))
print windows_filter.communicate()[0]

#
fasta_processor = popen("exonsToFasta.py --vcf={0} --input={1} --fasta={2} --output={3}".format(
    vcf_file_name,
    bed_file_name_1,
    fasta_file_name,
    fasta_file_name_1
))
print fasta_processor.communicate()[0]
#
nexus_generator = popen("snpToNexus.py --vcf={0} --input={1} --fasta={2} --output={3}".format(
  vcf_file_name,
  bed_file_name_1,
  fasta_file_name_1,
  output_folder
))
print nexus_generator.communicate()[0]

print("--------\ndata processed and can be found in folder specified")
