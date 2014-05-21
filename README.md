top_introns
===========
top_introns is a collection of python scripts for finding genome sequences that have a large number of SNPs per window length.

Initially it was made for counting SNPs in cheetah genome.

=Dependencies

All scripts are build on top of python 2.7.

=Usage

the main script is ```top_introns.py```. It calls other scripts in the following order, providing all necessary command line arguments:

- snpToBed.py

- exonsToFasta.py

- snpToNexus.py

as a result, you may expect to recieve the folder containing .nexus files, which describe top n windows.

=Input files format

=Project home page
The latest version of top_introns can be found at 

http://github.com/latticetower/top_introns