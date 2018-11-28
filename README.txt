Barcode QC Software



Requirements:

Unix-based terminal
Python 2.7
 - pip installed with:
   biopython
   bisect
   itertools
   collections
   argparse
   csv
   gzip
   datetime
   sqlite3



Installation:

Open a Unix terminal.
Navigate directories to find location of bcqc.tar.gz.
Extract bcqc.tar.gz in a working directory.

	tar -xvzf bcqc.tar.gz

Files included in bcqc/:
	-nermal.py
	-bc_reader.py
	-README.txt
	-spikedb/
		-*.db references


Use:

Navigate to working bcqc directory.
Execute python script.

	python bc_reader.py --i path/to/in_dir/ --o out_suffix --l seed_length

	--h for more help
	
	--i path/to/in_dir/
	The input directory must contain only .fastq.gz sequencing files. Do not
	include I1/I2 index files. The software will automatically filter the R1
	and/or R2 files with respect to paired-ended-ness. For example, if only 
	R1s are present in the input directory, there will not be any respect of
	pairs when filtering tagged reads. However, if R1 and R2 for a library
	are both present, the reads that are tagged in either will be removed in
	both.

	--o out_suffix
	All output reports will contain the out_suffix specified. Option 
	available for labelling purposes.

	--l seed_length
	Seed length. Recommended value: 12. Searching for tags in chunks of 
	seed length. Seed can be length 2, 3, 4, 6, 9, or 12. Less mismatches
	allowed and speed optomized for larger seeds. Any seed lengths under 9
	are not recommended. 
	
	--m (optional) allows 1 mismatch per seed length 
	Default: False.

	--d (optional) distributes multiple reads to each hit.
	Fractions of reads are distributed to tags where there are too many 
	mismatches to determine original tag. Default: False.



Example:

	cd bcqc
	python bc_reader.py --i testFastqs/ --l 12 --o 20180615 --f



Output:
	
The only file is is "prefix_counts.csv". This .csv file contains counts for 
each barcode within individual libraries, and the purity of the barcode within
a set of libraries. 


Version 1.0.0.2
20181128
Dylan Fox
dylan.fox@perkinelmer.com
