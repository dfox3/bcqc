import os
import sys
from os import listdir
from os.path import isfile, join
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import random
import itertools
import bisect
import argparse
import csv
import gzip
import string
from datetime import datetime

#--------------------------------------------------------------------------------------------------
#Command line input parameters
parser = argparse.ArgumentParser(description='Index Filtering')

parser.add_argument('--o',
                    metavar='-output',
                    required=True,
                    help='Name of output report')
parser.add_argument('--q',
                    metavar='-qual',
                    required=True,
                    help='Quality threshold')
parser.add_argument('--t',
                    metavar='-type',
                    required=True,
                    help='Type of filtering: 0 for median, 1 for average, 2 for any')
parser.add_argument('--i',
                    metavar='-input',
                    required=True,
                    help='Name of input dir')

#--------------------------------------------------------------------------------------------------

Q_SCORES = {'!':0, '"':1, '#':2, '$':3, '%':4, '&':5, "'":6, '(':7, ')':8,
            '*':9, '+':10,',':11,'-':12,'.':13,'/':14,'0':15,'1':16,'2':17,
            '3':18,'4':19,'5':20,'6':21,'7':22,'8':23,'9':24,':':25,';':26,
            '<':27,'=':28,'>':29,'?':30,'@':31,'A':32,'B':33,'C':34,'D':35,
            'E':36,'F':37,'G':38,'H':39,'I':40,'J':41}

def main():
	
    startTime = datetime.now()
    options = parser.parse_args()
    quality_threshold = int(options.q)
    print("Tagging Ill Indices")

    filter_indecies = {}

    only_files = [ f for f in listdir(options.i) if isfile(join(options.i, f)) ]
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
            


       	    handle = str(options.i) + "/" +str(o)
        
            for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
            	qualities.append(quality)
            	if paired_end_indecies_half_finished == False:    
                    filter_indecies[root_name].append(False)
            
            for q in xrange(len(qualities)):
            	# CHANGE THIS TO anyDips, medianDips, or averageDips
                if int(options.t) == 0:
            	    if medianDips(qualities[q], quality_threshold):
                        filter_indecies[root_name][q] = True
                elif int(options.t) == 1:
                    if averageDips(qualities[q], quality_threshold):
                        filter_indecies[root_name][q] = True
                elif int(options.t) == 2:
                    if anyDips(qualities[q], quality_threshold):
                        filter_indecies[root_name][q] = True
        else:
            if root_name not in filter_indecies:
                filter_indecies[root_name] = []
                handle = str(options.i) + "/" +str(o)
                print("***\nWARNING: Library " + str(o) + " is without index.\n***")
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

    printOut(temp, str(options.o)+"_report.csv")
    print("Run time: " + str(datetime.now() - startTime))
         
        
    print("Printing Groomed Reads")
    
    for o in only_files:
    	name_field = o.split('/')[-1][:-3]
    	root_name = root_names[o]
    	out_name = str(options.o) + '_' + str(name_field)
    	print(out_name)
    	f = 0
    	handle = str(options.i) + "/" +str(o)

    	with open(out_name, 'w') as nfq:
            for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
            	if filter_indecies[root_name][f] == False:
                    nfq.write("@"+str(title)+"\n")
                    nfq.write(str(sequence)+"\n")
                    nfq.write("+\n")
                    nfq.write(str(quality)+"\n")
            	f += 1
            nfq.close()

    print("Run time: " + str(datetime.now() - startTime))

def printOut(barcodes, out_name):
    with open(out_name, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(barcodes)

    
def dictionaryMerge(dict1, dict2):

    ret_dict = defaultdict(list)
    for k, v in itertools.chain(dict1.items(), dict2.items()):
        ret_dict[k].append(v)
    return ret_dict

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return float(sortedLst[index])
    else:
        return (float(sortedLst[index]) + float(sortedLst[index + 1]))/2.0

def anyDips(quality_string, threshold):
    # if the quality string dips a specified level, then return True
    ret_bool = False
    i = 0
    while ret_bool == False and i < len(quality_string):
        if Q_SCORES[quality_string[i]] < threshold:
            ret_bool = True
        i += 1
        
    return ret_bool

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

def medianDips(quality_string, threshold):
    # if the median is below a specified level, then return True
    ret_bool = False
    base_qualities = []
    for i in xrange(len(quality_string)):
        base_qualities.append(float(Q_SCORES[quality_string[i]]))
    median_q_score = median(base_qualities)
    if median_q_score < threshold:
        ret_bool = True
    return ret_bool

#-------------------------------------------------------------------------------
#Run program starting from main()
if __name__ == '__main__':
     main()
