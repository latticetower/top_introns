import fileinput
from collections import OrderedDict
# basic window object definition
class WindowInfo(object):
	def __init__(self, chromosome, window_start, window_end, name, species_no = 7):
		self.chromosome = chromosome
		self.window_start = long(window_start)
		self.window_end = long(window_end)
		self.window_name = name
		self.counter = 0
		self.species_no = species_no
		self.species = [ 0 for x in range(species_no)]
		self.line_infos = OrderedDict()
		self.sequence = ''
	#
	def process_line(self, line):
		if long(self.window_start) > long(line[1]) - 1 or long(self.window_end) < long(line[1]) - 1:
			return
		self.counter += 1 #the same as before
		for x in range(self.species_no):
			if line[- self.species_no + x].startswith('0/1'):
				self.species[x] += 1
		self.line_infos[long(line[1])] = line
	#
	def to_str(self):
		return "{0}\t{1}\t{2}".format(
			self.chromosome, 
			self.window_start, 
			self.window_end)
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
		return {
			'R': ['A', 'G'],
			'Y': ['C', 'T'],
			'M': ['A', 'C'],
			'K': ['T', 'G'],
			'W': ['T', 'A'],
			'S': ['C', 'G'],
			'B': ['C', 'T', 'G'],
			'D': ['A', 'T', 'G'],
			'H': ['A', 'T', 'C'],
			'V': ['A', 'C', 'G'],
			'N': ['A', 'C', 'T', 'G']
			}.get(letter, [letter])
	#
	def sequence_for_species(self, i):
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
		output_file_name = "{1}.nexus".format(output_folder_name, self.window_name)
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
			g.write("Species{0} {1}\n".format(i + 1, self.sequence_for_species(i)))
		g.write(";\nEnd;")
		g.close()#for i in range(7):