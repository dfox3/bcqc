import os
import sys
import re
from os import listdir
from os.path import isfile, join
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
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
import shutil
import operator

#import nermal
#from index_filter import printOut, dictionaryMerge, averageDips

#--------------------------------------------------------------------------------------------------
#Command line input parameters
parser = argparse.ArgumentParser(description='Spike-ins Tag Reader')

parser.add_argument('-i',
                    metavar='-input',
                    type=str,
                    required=True,
                    help="Input directory, contains reads in fastq.gz " + \
                    "format. Include I1/I2 index files if --x is thrown.")
parser.add_argument('-o',
                    metavar='-output',
                    type=str,
                    required=True,
                    help="Name of output report.\nAll output reports will " + \
                    "contain the out_suffix specified. Option available fo" + \
                    "r labelling purposes.")
parser.add_argument('-a',
                    metavar='-ampmap',
                    type=str,
                    required=True,
                    help="CSV of linked barcode,amplicon numbers")
parser.add_argument('-t',
                    metavar='-tag',
                    type=str,
                    default="references/tags.txt",
                    help="Name of tag txt file.\n Default is blah,")
parser.add_argument('-x',
                    dest='x',
                    action='store_true',
                    help='(Recommended) Index filter libraries before quality assessment.')
parser.add_argument('-d',
                    dest='d',
                    action='store_true',
                    help='Distribute mutliple alignments evenly')
parser.add_argument('-m',
                    dest='m',
                    action='store_true',
                    help='Turns on one mismatch per seed length. Levenshtein distance for mismatches.')
parser.add_argument('-l',
                    metavar='-scanlen',
                    type=int,
                    default=18,
                    help='Scan length. Default: 18 (recommended). Ignore if reads are > 36bp')
parser.add_argument('-w',
                    metavar='-windows',
                    type=int,
                    default=12,
                    help='Window length. Default: 12 (recommended)')
parser.add_argument('-s',
                    metavar='-seed',
                    type=int,
                    default=9,
                    help='Seed length. Default: 9 (recommended)')



parser.set_defaults(x=False)
parser.set_defaults(d=False)
parser.set_defaults(m=False)
#--------------------------------------------------------------------------------------------------

Q_SCORES = {'!':0, '"':1, '#':2, '$':3, '%':4, '&':5, "'":6, '(':7, ')':8,
            '*':9, '+':10,',':11,'-':12,'.':13,'/':14,'0':15,'1':16,'2':17,
            '3':18,'4':19,'5':20,'6':21,'7':22,'8':23,'9':24,':':25,';':26,
            '<':27,'=':28,'>':29,'?':30,'@':31,'A':32,'B':33,'C':34,'D':35,
            'E':36,'F':37,'G':38,'H':39,'I':40,'J':41}

REGEX_AID = {"x": "\\w", "b": "^A", "d": "^C", "h":"^G", "u":"^T"}
REGEX_ENCODE = {"A": "b", "C": "d", "G": "h", "T": "u"}

class Matcher(object):
    def __init__(self, l):
        self.l = l
        self.probes = 0
    def __call__(self, w):
        self.probes += 1
        pos = bisect.bisect_left(self.l, w)
        if pos < len(self.l):
            return self.l[pos]
        else:
            return None
class BKNode(object):
    def __init__(self, term):
        self.term = term
        self.children = {}
  
    def insert(self, other):
        distance = levenshtein(self.term, other)
        if distance in self.children:
            self.children[distance].insert(other)
        else:
            self.children[distance] = BKNode(other)
  
    def search(self, term, k, results=None):
        if results is None:
            results = []
        distance = levenshtein(self.term, term)
        counter = 1
        if distance <= k:
            results.append(self.term)
        for i in range(max(0, distance - k), distance + k + 1):
            child = self.children.get(i)
            if child:
                counter += child.search(term, k, results)
        return counter

def main():
    start_time = datetime.now()
    print("Parsing arguments")
    options = parser.parse_args()
    scan_length, mm, only_files, input_dir, count_out, index_filt, out_prefix, db_file, amp_map, windows, seed_length = initializeVariables(options)
    if index_filt:
        input_dir = indexFilter(input_dir, out_prefix)

    tag_bin_bin, total_seqs, only_files, collapsed_ref_names = binTags(input_dir, db_file, mm, scan_length, index_filt, amp_map, only_files, count_out, out_prefix, windows, seed_length, options.d)
    filt_final_csv, filt_percent_final_csv = toPrintInfo(only_files, collapsed_ref_names, tag_bin_bin, total_seqs, amp_map)
    toPrint(filt_final_csv, count_out)

    '''
    word = "AAAAGGGGTTTTCCCCAT"
    print("Zero mismatches")
    temp = generateLevRE(word, 0)
    print(str(len(temp)) + "\n\n"+ str(temp) + "\n\n")
    print("One mismatch")
    temp = generateLevRE(word, 1)
    print(str(len(temp)) + "\n\n"+ str(temp) + "\n\n")
    print("Two mismatches")
    temp = generateLevRE(word, 2)
    print(str(len(temp)) + "\n\n"+ str(temp) + "\n\n")
    #print("Three mismatches")
    #temp = generateLevRE(word, 3)
    #print(str(len(temp)) + "\n\n"+ str(temp) + "\n\n")
    '''


def initializeVariables(options):
    scan_length = int(100000)
    seed_length = int(options.s)
    SEED_LENGTHS = [2,3,4,5,6,7,8,9,11,12]

    if seed_length not in SEED_LENGTHS:
        seed_length = seedFix(seed_length)

    if options.l != -1:
        scan_length = options.l

    mismatch = 0
    if options.m:
        mismatch = 1

    only_files = [f for f in listdir(options.i) if isfile(join(options.i, f))]
    only_files.sort()

    count_out = str(options.o) + "_counts.csv"

    db_tag = str(options.t)
    while db_tag[-1] == "/":
        db_tag = db_tag[:-1]
    db_tag = db_tag.split("/")[-1]
    db_file = str(options.t) + "/" + str(db_tag) + str(seed_length) + "." + str(mismatch) + ".db"

    amp_map = getAmpMap(options.a)

    return scan_length, mismatch, only_files, str(options.i), count_out, options.x, options.o, db_file, amp_map, options.w, seed_length

def getAmpMap(amp_file):
    with open(amp_file, 'rb') as f:
        reader = csv.reader(f)
        amp_list = list(reader)[1:]
        #print({ int(a[0]):int(a[1]) for a in amp_list })
        return { int(a[0]):int(a[1]) for a in amp_list }

def populateWords(csv_name):
    with open(csv_name, 'rb') as f:
        reader = csv.reader(f)
        ret_list = []
        bins_list = list(reader)
        for b in bins_list:
            ret_list.append(b[0])
    return ret_list

def indexFilter(input_dir, out_prefix):
    startTime = datetime.now()
    quality_threshold = 30
    
    print("Tagging Ill Indices")
    filter_indecies = {}
    index_filtered_dir = str(input_dir) + "i"
    if not os.path.exists(index_filtered_dir):
        os.makedirs(index_filtered_dir)
    else:
        shutil.rmtree(index_filtered_dir)
        os.makedirs(index_filtered_dir)

    only_files = [ f for f in listdir(input_dir) if isfile(join(input_dir, f)) ]
    only_files.sort()
    temp = []
    root_names = {}

    for o in only_files:
        qualities = [] 
        name_fields = o.split('_')
        index_file = False
        root_name = ""
        for n in name_fields:
            if n == "I1" or n == "I2":
                index_file = True
                break
            if n == "R1" or n == "R2":
                break
            root_name = str(root_name) + str(n) + "_"
        root_name = root_name[:-1]
        if root_name not in root_names:
            root_names[o] = root_name

        print("Parsing " + str(o))
        if index_file:
            paired_end_indecies_half_finished = False
            if root_name not in filter_indecies:
                filter_indecies[root_name] = []
            else:
                paired_end_indecies_half_finished = True

            handle = str(input_dir) + "/" +str(o)
        
            for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
                qualities.append(quality)
                if paired_end_indecies_half_finished == False:    
                    filter_indecies[root_name].append(False)
            
            for q in xrange(len(qualities)):
                if averageDips(qualities[q], quality_threshold):
                    filter_indecies[root_name][q] = True
        else:
            if root_name not in filter_indecies:
                filter_indecies[root_name] = []
                handle = str(input_dir) + "/" +str(o)
                print("***\nWARNING: Library " + str(o) + " is without corresponding index file(s).\n***")
                for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
                    filter_indecies[root_name].append(False)

        print("# total reads: " + str(len(filter_indecies[root_name])))
        true_filtereds = [ f for f in xrange(len(filter_indecies[root_name])) if
                           filter_indecies[root_name][f] == True ]
        print(" # filtered reads: " + str(len(true_filtereds)))

        if len(filter_indecies[root_name]) == 0:
            print("  % filtered reads: 0")
            temp.append([root_name, len(filter_indecies[root_name]),
                         len(true_filtereds), 0])
        else:
            print("  % filtered reads: " + str(float(len(true_filtereds))/
                                               float(len(filter_indecies[root_name]))*100))
            temp.append([root_name, len(filter_indecies[root_name]),
                         len(true_filtereds), float(len(true_filtereds))/
                         float(len(filter_indecies[root_name]))*100])

    printOut(temp, str(out_prefix)+"_index_filter_report.csv")
        
    print("Printing Groomed Reads")
    for o in only_files:
        name_field = o.split('/')[-1][:-3]
        root_name = root_names[o]
        out_name = str(index_filtered_dir) + "/i" + str(name_field)
        print(out_name)
        f = 0
        handle = str(input_dir) + "/" + str(o)
        n_fields = o.split('_')
        for n in n_fields:
            if n == "R1" or n == "R2":
                with open(out_name, 'w') as nfq:
                    content = ""
                    for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
                        if filter_indecies[root_name][f] == False:
                            nfq.write("@"+str(title)+"\n")
                            nfq.write(str(sequence)+"\n")
                            nfq.write("+\n")
                            nfq.write(str(quality)+"\n")
                        f += 1
                    nfq.close()

    print("Run time: " + str(datetime.now() - startTime))
    return index_filtered_dir


def binTags(input_dir, db_file, mm, scan_length, index_filt, amp_map, only_files, count_out, out_prefix, windows, seed_length, distribute_multiples):
    # compile tag regex and map


    collapsed_ref_names, filt_ref_pos, match_ref_pos = noahLoad(db_file)

    print("Counting tags")
    filter_reads = {}
    tag_bin_bin = {}
    total_seqs = {}
    new_only_files = []

    if index_filt:
        only_files = [f for f in listdir(input_dir) if isfile(join(input_dir, f))]
        only_files.sort()
    for o in only_files:

        seqs = []
        titles = []
        qualities = []
        
        handle = str(input_dir) + "/" +str(o)

        name_fields = o.split('_')
        #print(name_fields)
        if "R1" in name_fields:
            new_only_files.append(o)
            print("Parsing " + str(o))
            paired_end_reads_half_finished = False
            root_name = ""
            for e in range(len(name_fields)-2):
                root_name += str(name_fields[e])
                root_name += "_"
            root_name = root_name[:-1]
            print("Root:\t" + str(root_name))
            if root_name not in filter_reads:
                filter_reads[root_name] = []
            else:
                paired_end_reads_half_finished = True

            if index_filt:
                for (title, sequence, quality) in FastqGeneralIterator(open(handle)):
                    seqs.append(sequence)
                    titles.append(title)
                    qualities.append(quality)
                    if paired_end_reads_half_finished == False:    
                        filter_reads[root_name].append(False)
            else:
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
                if len(seqs[s]) < scan_length:
                    small_reads += 1
                    
                else:
                    sizeable_reads += 1
                    found = False
                    x = 0
                    z = 0
                    candidates = {}
                    new_candidates = {}
                    cand_path = []
                    meta_path = []
                    iters = 0

                    while x < windows and found == False:
                        z = 0
                        iters = 0
                        #print("reseting z:\t" + str(z))
                        primer_scan = True
                        candidates = {}
                        new_candidates = {}
                        cand_path = []
                        meta_path = []
                        while primer_scan:
                            seq = seqs[s][z+x:z+x+seed_length]
                            temp = []
                            temp_seqs = seqs[s]
                            new_candidates = {}
                            #print(temp_seqs)
                            #print(seq)
                            if z == 0:
                                if seq in filt_ref_pos:
                                    for y in filt_ref_pos[seq]:
                                        #print(match_ref_pos[y])
                                        for m in match_ref_pos[y]:
                                            #print((z,seq,m[3],m[2],m[0]))
                                            temp.append((z,seq,m[3],m[2],m[0]))
                                            if m[3] not in tag_bin:
                                                tag_bin[m[3]] = 0
                                            if m[3] not in candidates:
                                                candidates[m[3]] = {}
                                            if m[0] > 0: 
                                                candidates[m[3]][m[2]+seed_length] = m[0]
                                            else:
                                                candidates[m[3]][m[2]-seed_length] = m[0]
                                #else:
                                #    print("Nothing found")

                            else:
                                if seq in filt_ref_pos:
                                    temp = []
                                    temp_seqs = seqs[s]
                                    for y in filt_ref_pos[seq]:
                                        #print(match_ref_pos[y])
                                        for m in match_ref_pos[y]:
                                            if m[3] in candidates:
                                                #print((z,seq,m[3],m[2],m[0]))
                                                temp.append((z,seq,m[3],m[2],m[0]))
                                                if m[2] in candidates[m[3]]:
                                                    if candidates[m[3]][m[2]] == m[0]:
                                                        if m[3] not in new_candidates:
                                                            new_candidates[m[3]] = {}
                                                        
                                                        if m[0] > 0: 
                                                            new_candidates[m[3]][m[2]+seed_length] = m[0]
                                                        else:
                                                            new_candidates[m[3]][m[2]-seed_length] = m[0]
                                #else:
                                #    print("Broken Chain")

                            if z != 0:
                                candidates = new_candidates
                            cand_path.append(candidates)
                            meta_path.append(temp)
                            if z + seed_length >= scan_length or len(candidates) == 0:
                                primer_scan = False
                                if len(candidates) != 0:
                                    #print("found!")
                                    found = True
                            z = seed_length * iters
                            iters += 1
                            #print("z:\t" + str(z))
                        #print("x:\t" + str(x))
                        x += 1
                    #print("s:\t" + str(s))
                    #print("")
                            
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
    return tag_bin_bin, total_seqs, new_only_files, collapsed_ref_names


def toPrintInfo(only_files, collapsed_ref_names, tag_bin_bin, total_seqs, amp_map):
    filt_final_csv = [ [f] for f in collapsed_ref_names ]
    filt_final_csv = [["libs"]] + filt_final_csv + [["multiple"],["not_aligned"],["total_reads"]]
    filt_percent_final_csv = [ [f] for f in collapsed_ref_names ]
    filt_percent_final_csv = [["libs"]] + filt_percent_final_csv + [["multiple"],["not_aligned"],["total_reads"]]
    f_keys = [ f[0] for f in filt_final_csv ]
    for p in only_files:
        most_abundant = 0
        most_abundant_hits = 0
        tag_total = 0
        for q in range(len(f_keys)):
            if q == 0:
                filt_final_csv[q].append(p)
                filt_percent_final_csv[q].append(p)
            elif q == len(f_keys)-1:
                filt_final_csv[q].append(total_seqs[p])
                filt_percent_final_csv[q].append(total_seqs[p])
            
            else:
                k = unicode(str(f_keys[q]), 'utf-8')
                temp_hits = 0
                if k in tag_bin_bin[p] and total_seqs[p] != 0:

                    filt_final_csv[q].append(tag_bin_bin[p][k])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][k])/float(total_seqs[p])*100.0)
                    temp_hits = tag_bin_bin[p][k]                
                elif total_seqs[p] != 0 and (f_keys[q] == "multiple" or f_keys[q] == "not_aligned"):
                    filt_final_csv[q].append(tag_bin_bin[p][f_keys[q]])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][f_keys[q]])/float(total_seqs[p])*100.0)
                else:
                    filt_final_csv[q].append(0.0)
                    filt_percent_final_csv[q].append(0.0)
                tag_total += temp_hits
                if temp_hits > most_abundant_hits:
                    most_abundant = q
                    most_abundant_hits = temp_hits
    print("filt final csv")
    for f in filt_final_csv:
        print(f)

    not_enough_reads = []
    not_pure = []
    dirty_when_empty = []
    max_idxs = []
    total_tags = []
    purities = []

    for i, j in enumerate(filt_final_csv):
        #print(j)
        if i != 0 and i < len(f_keys)-3:
            max_val = 0
            max_idx = -1
            sum_vals = 0
            for k, l in enumerate(j):
                if k != 0:
                    if l > max_val:
                        max_val = l
                        max_idx = k
                    sum_vals += l


            filt_final_csv[i].append(sum_vals)
            total_tags.append(sum_vals)
            max_idxs.append(max_idx)
            #print(max_idx)
           # print(len(max_idxs))
            if sum_vals != 0:
                filt_final_csv[i].append(str(format((float(max_val)/float(sum_vals)) * 100.0, '.2f')))
                purities.append((float(max_val)/float(sum_vals)) * 100.0)
            else:
                filt_final_csv[i].append("n/a")
                purities.append("n/a")
        else:
            if i == 0:
                filt_final_csv[i].append("total_bc_reads")
                filt_final_csv[i].append("bc_purity")

    found_in_amp_map = []
    amp_map_keys = amp_map.keys()
    amp_map_keys.sort()
    tags = set()
    for a in amp_map_keys:
        tags.add(amp_map[a])
    for i, j in enumerate(max_idxs):
        print(i)
        print(max_idxs[i])
        if j in amp_map:
            print(amp_map[j])
        else:
            print("not in amp_map")
        print(purities[i])
        print(0)
        if j in amp_map and i+1 in tags:
            found_in_amp_map.append(j)
            if (i+1 != amp_map[j] or purities[i] < 99.5):
                not_pure.append((j,amp_map[j]))
            if total_tags[i] < 1000:
                not_enough_reads.append("bc#" + str(j))
        else:
            if total_tags[i] > 100:
                dirty_when_empty.append("bc#" + str(j) + ":tag#" + str(i+1))
    found_in_amp_map = set(found_in_amp_map)
    for a in amp_map_keys:
        if a not in found_in_amp_map:
            not_pure.append((a,amp_map[a]))
    not_pure.sort(key=operator.itemgetter(0))
    not_pure = [ "bc#" + str(n[0]) + ":tag#" + str(n[1]) for n in not_pure ]        


    filt_final_csv = [["Not pure:"] + [ n for n in not_pure ], ["Warning - Low reads:"] + [ n for n in not_enough_reads ], ["Warning - Dirty empty:"] + [ n for n in dirty_when_empty ] ] + filt_final_csv
    filt_final_csv[3] = ["tag"] + [ "BC" + str(x+1).zfill(3) for x in range(len(filt_final_csv[3][1:-2])) ] + filt_final_csv[3][-2:]

    return filt_final_csv, filt_percent_final_csv
    
def generateLevRE(word, mm):
    # recursive levenshtein regex generator
    words = {word}
    if mm != 0:
        # substitution
        for w in range(len(word)):
            if word[w] not in REGEX_AID:
                words = words | generateLevRE(word[:w] + REGEX_ENCODE[word[w]] + word[w+1:], mm-1)
        # insertion
        for w in range(1,len(word)-1):
            words = words | generateLevRE(word[:w] + "x" + word[w:], mm-1)
        # deletion
        for w in range(len(word)):          
            if w == 0:
                words = words | generateLevRE(word[1:], mm-1)
            elif w == len(word)-1:
                words = words | generateLevRE(word[:len(word)-(1)], mm-1)
            else:
                words = words | generateLevRE(word[:w] + word[w+1:], mm-1)
    return words


def loadCSV(file_name):
    with open(file_name) as f:
        reader = csv.reader(f)
        data = list(reader)
    return data

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
    #ff_csv = [ b for a, b in enumerate(filt_final_csv) if a < 97 ]
    with open(count_out, 'w') as of:
        a = csv.writer(of, delimiter=',')
        a.writerows(filt_final_csv)
        of.close()

def seedFix(seed_length):
    ret_length = 2
    print(seed_length)
    if seed_length < 2:
        print("Seed length is too small; defaulting to 2. Advising seed_length of 12.")
    if seed_length > 12:
        print("Seed length capped at 12; defaulting to 12.")
        ret_length = 12
    return ret_length

def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def printOut(barcodes, out_name):
    with open(out_name, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(barcodes)

    
def dictionaryMerge(dict1, dict2):

    ret_dict = defaultdict(list)
    for k, v in itertools.chain(dict1.items(), dict2.items()):
        ret_dict[k].append(v)
    return ret_dict


def averageDips(quality_string, threshold):
    # if the average is below a specified level, then return True
    ret_bool = False
    total_q_score = 0
    for i in xrange(len(quality_string)):
        total_q_score += float(Q_SCORES[quality_string[i]])
    average_q_score = total_q_score / float(len(quality_string))
    if average_q_score < threshold:
        ret_bool = True
    return ret_bool


#--------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
