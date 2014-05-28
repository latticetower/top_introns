top_introns
===========

top_introns is a collection of Python scripts for finding genome sequences that have a large number of SNPs per window length.

Initially it was made for counting SNPs in cheetah genome.

Requirements
------------
python 2.7

Usage
-----

The main script is ```top_introns.py```. It calls other scripts in the following order, providing all necessary command line arguments:

- snp_to_bed.py,

- windows_to_fasta.py,

- snp_to_nexus.py.

As a result, you may expect to recieve the folder containing .nexus files, which describe top n windows.

Sample usage (Windows):

```
python top_introns.py -v .\test_dataset2\test_snp_list.vcf -f .\test_dataset2\test_input.fasta -w .\test_dataset2\test_windows.bed -o results -n 3 --normalize
```
Sample usage (Linux):
```
./top_introns.py --vcf test_dataset2/test_snp_list.vcf --fasta test_dataset2/test_input.fasta -w test_dataset2/test_windows.bed --normalize -o "res"
```
or you can call ```python top_introns.py --help``` and get information about available options and their meaning.

Input files format
------------------
The main script and worker scripts expect the following input data:

1. .fasta file with genome information. It should contain collection of .fasta strings, one for each chromosome, with chromosome name in description line.

2. .vcf file with collection of tab-separated fields, describing collection of SNPs found, with order of fields like in specification:
    http://www.1000genomes.org/node/101, with few specific features:
   + chromosome name should be equal to string given as description for particular string in .fasta file
   + last `m` tab-separated columns should correspond to `m` species genotype information for given SNP (as stated in .vcf format specification).
   + script assumes that SNP lines are sorted in ascending order by chromosome and position in chromosome.

3. .bed-like file with information about windows being examined (list of exons, introns, etc.). File is tab-separated, only first 4 columns are used: `chromosome name`, `start` and `end` positions, `window name`.

Example input files can be found in folders named `test_dataset1` and `test_dataset2`.

Output files format
-------------------
All output files are placed to selected folder and are named as `window name`.nexus (for each of top `n` windows).
The folder also contains:
  + .fasta file with list of windows sequences (taken from reference file)
  + .bed file with list of top `n` windows

These two files are generated and used by worker scripts to store intermediate results.

Calling worker scripts independently
------------------------------------
If you want to call worker scripts independently, call them with "--help" key, to get information about their parameters.

Project home page
-----------------

The latest version of top_introns can be found at

http://github.com/latticetower/top_introns
