import os.path
from snp_to_nexus import NexusSaver
from snp_to_bed import TopNWindowsSelector
from windows_to_fasta import WindowToFastaConverter

class DataProcessor(object):
	def __init__(self, window_file_name,
											vcf_file_name,
											fasta_file_name,
											output_folder,
											top_n, 
											min_size,
											normalize):
		self.window_file_name = window_file_name
		self.vcf_file_name = vcf_file_name
		self.fasta_file_name = fasta_file_name
		self.output_folder = output_folder
		self.top_n = top_n 
		self.min_size = min_size
		self.normalize = normalize
		#ensure output folder existance and create if nessesary
		dir_path = os.path.join(self.output_folder)
		if not os.path.isdir(dir_path): os.makedirs(dir_path)
		self.bed_file_name_1 = os.path.join(dir_path, "results_top{0}.bed".format(top_n))
		self.fasta_file_name_1 = os.path.join(dir_path, "window_top{0}.fasta".format(top_n))

		#processes given data
	def process(self):
		self.select_top_n_windows()
		self.convert_to_fasta()
		self.save_to_nexus()
	def select_top_n_windows(self):
		print "Selecting top n windows from .bed file and saving them..."
		from snp_to_bed import TopNWindowsSelector
		selector = TopNWindowsSelector(self.top_n, self.normalize)
		selector.load_windows_info(self.window_file_name)
		selector.load_snp(self.vcf_file_name)
		selector.select_top_n(self.min_size)
		selector.save_to_bed_file(self.bed_file_name_1)
	def convert_to_fasta(self):
		print "Getting fasta sequences for top n windows..."
		from windows_to_fasta import WindowToFastaConverter
		converter = WindowToFastaConverter()
		converter.load_window_list(self.bed_file_name_1)
		converter.load_info_from_fasta(self.fasta_file_name)
		converter.save_sequences_to_fasta(self.fasta_file_name_1)
	def save_to_nexus(self):
		print "Saving results to .nexus files..."
		from snp_to_nexus import NexusSaver
		nexus_saver = NexusSaver()
		nexus_saver.load_windows(self.bed_file_name_1)
		nexus_saver.load_from_vcf(self.vcf_file_name)
		nexus_saver.process_and_save(self.fasta_file_name_1, self.output_folder)
#
		#