import os
import sys
import re
from os import listdir
from os.path import isfile, join
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import nermal
import bisect
import random
import itertools
from collections import defaultdict
import argparse
import csv
import gzip
import string
from datetime import datetime
import sqlite3

#--------------------------------------------------------------------------------------------------
#Command line input parameters
parser = argparse.ArgumentParser(description='Spike-ins Tag Reader')

parser.add_argument('--i',
                    metavar='-input',
                    type=str,
                    required=True,
                    help="Input directory, contains reads in fastq.gz " + \
                    "format. Do not include I1/I2 index file" + \
                    "s. The software will automatically filter the R1 and/" + \
                    "or R2 files with respect to paired-ended-ness. For ex" + \
                    "ample, if only R1s are present in the input directory" + \
                    ", there will not be any respect of pairs when filterin" + \
                    "g tagged reads. However, if R1 and R2 for a library a" + \
                    "re both present, the reads that are tagged in either " + \
                    "will be removed in both.")
parser.add_argument('--o',
                    metavar='-output',
                    type=str,
                    required=True,
                    help="Name of output report.\nAll output reports will " + \
                    "contain the out_suffix specified. Option available fo" + \
                    "r labelling purposes.")
parser.add_argument('--m',
                    dest='m',
                    action='store_true',
                    help='Turns on one mismatch per seed length.')
parser.add_argument('--d',
                    dest='d',
                    action='store_true',
                    help='Turns on the distribution of reads for multiple alignments. (not recommended)')
parser.add_argument('--l',
                    metavar='-len',
                    type=int,
                    help='Seed length. Default: 12 (recommended). Ignore if reads are > 36bp')
parser.add_argument('--s',
                    metavar='-scanlen',
                    type=int,
                    help='Scan length. Default: 32 (recommended). Ignore if reads are > 36bp')


parser.set_defaults(l=12)
parser.set_defaults(s=32)
parser.set_defaults(m=False)
parser.set_defaults(d=False)
#--------------------------------------------------------------------------------------------------

def main():
    start_time = datetime.now()
    print("Parsing arguments")
    options = parser.parse_args()
    seed_length, scan_length, distribute_multiples, only_files, db_name, input_dir, count_out = initializeVariables(options)
    collapsed_ref_names, filt_ref_pos, match_ref_pos = noahLoad(db_name)
    tag_bin_bin, total_seqs = binning(seed_length, scan_length, distribute_multiples, only_files, input_dir, filt_ref_pos, match_ref_pos)
    filt_final_csv, filt_percent_final_csv = toPrintInfo(only_files, collapsed_ref_names, tag_bin_bin, total_seqs)
    toPrint(filt_final_csv, count_out)
    print("Run time: " + str(datetime.now() - start_time))


def initializeVariables(options):
    seed_length = int(options.l)
    scan_length = int(options.s)
    distribute_multiples = options.d
    SEED_LENGTHS = [2,3,4,5,6,7,8,9,11,12]

    if seed_length not in SEED_LENGTHS:
        seed_length = seedFix(seed_length)
    mismatch = 0
    if options.m:
        mismatch = 1

    only_files = [f for f in listdir(options.i) if isfile(join(options.i, f))]
    only_files.sort()

    db_name = "spikedb/noah_spike" + str(seed_length) + "." + str(mismatch) + ".db"
    
    count_out = str(options.o) + "_counts.csv"

    return seed_length, scan_length, distribute_multiples, only_files, db_name, str(options.i), count_out


def binning(seed_length, scan_length, distribute_multiples, only_files, input_dir, filt_ref_pos, match_ref_pos):
    print("Counting tags")
    filter_reads = {}
    tag_bin_bin = {}
    total_seqs = {}
    for o in only_files:
        print("Parsing " + str(o))

        seqs = []
        titles = []
        qualities = []
        
        handle = str(input_dir) + "/" +str(o)

        name_fields = o.split('_')
        paired_end_reads_half_finished = False
        root_name = ""
        for e in xrange(len(name_fields)-2):
            root_name += str(name_fields[e])
            root_name += "_"
        root_name = root_name[:-1]
        print("Root:\t" + str(root_name))
        if root_name not in filter_reads:
            filter_reads[root_name] = []
        else:
            paired_end_reads_half_finished = True

        for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
            seqs.append(sequence)
            titles.append(title)
            qualities.append(quality)
            if paired_end_reads_half_finished == False:    
                filter_reads[root_name].append(False)

        found_count = 0
        miss_count = 0
        small_reads = 0
        sizeable_reads = 0
        processing = []
        score = float(1.0)
        total_reads = len(seqs)
        total_tags = 0
        tag_bin = {"multiple" : 0, "not_aligned" : 0}

        
        for s in range(len(seqs)):
            if len(seqs[s]) <= seed_length:
                small_reads += 1
                
            else:
                sizeable_reads += 1
                primer_scan = True
                found = False
                z = 0
                screened_bases = []
                amplicon_lengths = []
                score = float(1.0)
                candidates = {}
                new_candidates = {}
                cand_path = []
                meta_path = []
                iters = 1

                while primer_scan:
                    seq = seqs[s][z:z+seed_length]
                    temp = []
                    temp_seqs = seqs[s]
                    new_candidates = {}

                    if z == 0:
                        if seq in filt_ref_pos:
                            for y in filt_ref_pos[seq]:
                                for m in match_ref_pos[y]:
                                    if m[3] == "1F" or m[3] == "1R":
                                        temp.append((z,seq,m[3],m[2],m[0]))
                                    if m[3] not in tag_bin:
                                        tag_bin[m[3]] = 0
                                    if m[3] not in candidates:
                                        candidates[m[3]] = {}
                                    if m[0] > 0: 
                                        candidates[m[3]][m[2]+seed_length] = m[0]
                                    else:
                                        candidates[m[3]][m[2]-seed_length] = m[0]                                              

                    else:
                        if seq in filt_ref_pos:
                            temp = []
                            temp_seqs = seqs[s]
                            
                            for y in filt_ref_pos[seq]:
                                for m in match_ref_pos[y]:
                                    if m[3] in candidates:
                                        temp.append((z,seq,m[3],m[2],m[0]))
                                        if m[2] in candidates[m[3]]:
                                            if candidates[m[3]][m[2]] == m[0]:
                                                if m[3] not in new_candidates:
                                                    new_candidates[m[3]] = {}
                                                
                                                if m[0] > 0: 
                                                    new_candidates[m[3]][m[2]+seed_length] = m[0]
                                                else:
                                                    new_candidates[m[3]][m[2]-seed_length] = m[0]

                    if z != 0:
                        candidates = new_candidates
                    cand_path.append(candidates)
                    meta_path.append(temp)
                    if z + seed_length >= scan_length or len(candidates) == 0:
                        primer_scan = False
                        if len(candidates) != 0:
                            found = True
                    z = seed_length * iters
                    iters += 1
                        
                if found:
                    total_tags += 1
                    found_count += 1

                    c_keys = candidates.keys()  
                    c_keys.sort()

                    if len(c_keys) == 1:
                        tag_bin[c_keys[0]] += 1.0
                    elif len(c_keys) > 1:
                        if distribute_multiples:
                            d_mult = 1.0 / float(len(c_keys))
                            for c in c_keys:
                                tag_bin[c] += d_mult
                        tag_bin["multiple"] += 1.0
                    filter_reads[root_name][s] = True
                else:
                    miss_count += 1
                    tag_bin["not_aligned"] += 1.0
                    f_true = False
                    r_true = False
                    hey_look_at_this = False

        tag_keys = tag_bin.keys()
        tag_keys.sort()
        for t in tag_keys:
            if tag_bin[t] != 0:
                print("\t* " + str(t) + ":\t" +str(tag_bin[t]) + "\t% total: " + str(format(float(tag_bin[t])/float(len(seqs))*100.0,'.2f')))
            
        print("\t* total reads: " + str(len(seqs)))
        tag_bin_bin[o] = tag_bin
        total_seqs[o] = len(seqs)
    return tag_bin_bin, total_seqs


def toPrintInfo(only_files, collapsed_ref_names, tag_bin_bin, total_seqs):
    filt_final_csv = [ [f] for f in collapsed_ref_names ]
    filt_final_csv = [["libs"]] + filt_final_csv + [["multiple"],["not_aligned"],["total_reads"]]
    filt_percent_final_csv = [ [f] for f in collapsed_ref_names ]
    filt_percent_final_csv = [["libs"]] + filt_percent_final_csv + [["multiple"],["not_aligned"],["total_reads"]]
    f_keys = [ f[0] for f in filt_final_csv ]
    for p in only_files:
        most_abundant = 0
        most_abundant_hits = 0
        tag_total = 0
        for q in xrange(len(f_keys)):
            if q == 0:
                filt_final_csv[q].append(p)
                filt_percent_final_csv[q].append(p)
            elif q == len(f_keys)-1:
                filt_final_csv[q].append(total_seqs[p])
                filt_percent_final_csv[q].append(total_seqs[p])
            
            else:
                ff = unicode(str(f_keys[q]) + "F", 'utf-8')
                rr = unicode(str(f_keys[q]) + "R", 'utf-8')
                temp_hits = 0
                if rr in tag_bin_bin[p] and ff in tag_bin_bin[p]:
                    filt_final_csv[q].append(tag_bin_bin[p][ff] + tag_bin_bin[p][rr])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][ff] + tag_bin_bin[p][rr])/float(total_seqs[p])*100.0)
                    temp_hits = tag_bin_bin[p][ff] + tag_bin_bin[p][rr]
                elif rr in tag_bin_bin[p]:
                    filt_final_csv[q].append(tag_bin_bin[p][rr])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][rr])/float(total_seqs[p])*100.0)
                    temp_hits = tag_bin_bin[p][rr]     
                elif ff in tag_bin_bin[p]:
                    filt_final_csv[q].append(tag_bin_bin[p][ff])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][ff])/float(total_seqs[p])*100.0)
                    temp_hits = tag_bin_bin[p][ff]
                elif total_seqs[p]!=0 and (f_keys[q] == "multiple" or f_keys[q] == "not_aligned"):
                    filt_final_csv[q].append(tag_bin_bin[p][f_keys[q]])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][f_keys[q]])/float(total_seqs[p])*100.0)
                else:
                    filt_final_csv[q].append(0.0)
                    filt_percent_final_csv[q].append(0.0)
                tag_total += temp_hits
                if temp_hits > most_abundant_hits:
                    most_abundant = q
                    most_abundant_hits = temp_hits

    for i, j in enumerate(filt_final_csv):
        if i != 0 and i < len(f_keys)-4:
            max_val = 0
            sum_vals = 0
            for k, l in enumerate(j):
                if k != 0:
                    if l > max_val:
                        max_val = l
                    sum_vals += l

            filt_final_csv[i].append(sum_vals)
            if sum_vals != 0:
                print(max_val)
                print(sum_vals)
                filt_final_csv[i].append(str(format((float(max_val)/float(sum_vals)) * 100.0, '.2f')))
            else:
                filt_final_csv[i].append("n/a")
        else:
            if i == 0:
                filt_final_csv[i].append("total_bc_reads")
                filt_final_csv[i].append("bc_purity")

    return filt_final_csv, filt_percent_final_csv


def noahLoad(db_name):
    print("Accessing " + str(db_name))
    
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    filt_ref_pos = {}
    cur.execute("SELECT * FROM filter_ref")
    for row in cur.fetchall():
        if row[0] not in filt_ref_pos:
            filt_ref_pos[row[0]] = []
        filt_ref_pos[row[0]].append(row[2])

    cur.execute("SELECT * FROM match_dict")
    match_ref_pos = {}
    ids = {}
    for row in cur.fetchall():
        if row[0] not in match_ref_pos:
            match_ref_pos[row[0]] = []
        if row[5] not in ids:
            ids[row[5]] = 0
        match_ref_pos[row[0]].append((row[2],row[3],row[4],row[5]))

    cur.execute("SELECT id FROM ids")
    ref_names = cur.fetchall()
    ref_names = [ r[0] for r in ref_names ]
    collapsed_ref_names = list(set([ int(''.join(re.split('[a-z]+', str(r), flags=re.IGNORECASE))) if isInt(''.join(re.split('[a-z]+', str(r), flags=re.IGNORECASE))) else r for r in ref_names ]))
    collapsed_ref_names.sort()

    return [ str(c) for c in collapsed_ref_names ], filt_ref_pos, match_ref_pos


def toPrint(filt_final_csv, count_out):
    ff_csv = [ b for a, b in enumerate(filt_final_csv) if a < 97 ]
    with open(count_out, 'w') as of:
        a = csv.writer(of, delimiter=',')
        a.writerows(ff_csv)
        of.close()


def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


def seedFix(seed_length):
    ret_length = 2
    print(seed_length)
    if seed_length < 2:
        print("Seed length is too small; defaulting to 2. Advising seed_length of 12.")
    if seed_length > 12:
        print("Seed length capped at 12; defaulting to 12.")
        ret_length = 12
    return ret_length


#--------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
