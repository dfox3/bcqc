# Barcode QC Software



## Concept:

This script qualifies **_[NEXTFLEX® Unique Dual Index Barcode](https://perkinelmer-appliedgenomics.com/library-preparation-kits/nextflex-ngs-barcodes/nextflex-unique-dual-index-barcodes/)_** purity by aligning NGS data ([fastq reads](https://en.wikipedia.org/wiki/FASTQ_format)) to a barcoded amplicon reference. This aligner uses [Levenshtein](https://en.wikipedia.org/wiki/Levenshtein_distance) automata to account for fuzzy matching of single position mismatched sequences, inserted/deleted sequences, and frame shifted sequences. 



## Installation:

Open a Unix terminal.

`git clone https://github.com/dfox3/bcqc.git`

Alternatively, you can download the zip from [the current repo](https://github.com/dfox3/bcqc) and unzip the directory.


Files included in bcqc/:
 - nermal.py
 - bc_reader.py
 - README.txt
 - spikedb/
    - \*.db references



## Requirements:

Unix-based terminal
Python 2.7
 - pip installed with:
    - biopython

For more information on how to pip install, visit [this guide](https://packaging.python.org/tutorials/installing-packages/#ensure-you-can-run-pip-from-the-command-line) for more information.

 - other imports from [The Python Standard Library](https://docs.python.org/2/library/) (typically pre-installed):
    - bisect
    - itertools
    - collections
    - argparse
    - csv
    - gzip
    - datetime
    - sqlite3

Script may be compatible with Python 3 version, but never tested.



## Use:

Navigate to working bcqc directory.
Execute python script.

```
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

--m (optional) allows 1 mismatch per seed length 
Default: False.

--d (optional) distributes multiple reads to each hit.
Fractions of reads are distributed to tags where there are too many 
mismatches to determine original tag. Default: False.

--l (optional) Seed length. Default: 12 (recommended). 
Minimum: 2. Maximum: 12. Ignore if reads are > 36bp.

--s (optional) Scan length. Default: 32 (recommended). 
Ignore if reads are > 36bp.
```



### Example:

```
cd bcqc
python bc_reader.py --i testFastqs/ --o test
```



### Output:
	
The only file output is "prefix_counts.csv" (the example output will be in 
"test_counts.csv"). This .csv file contains counts for each barcode within 
individual libraries, and the purity of the barcode within a set of libraries. 



### Info:

**20181213**

**Dylan Fox**

**dylan.fox@perkinelmer.com**