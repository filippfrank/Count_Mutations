#!/usr/bin/python
from Bio.Seq import Seq
import os
from collections import defaultdict

#Read the fastq file
def process(lines):
	ks = ['name', 'sequence', 'optional', 'quality']
	return {k: v for k, v in zip(ks, lines)}

os.chdir('/Users/filipp/Library/CloudStorage/OneDrive-EmoryUniversity/Documents - Ortlund Lab - RADx VTF/RADx VTF/HTS/ScreeningData/220518_N_Lib05')

fastq_files = [f for f in os.listdir('./data') 
               if f.endswith('R2_001.fastq')]

os.chdir('data')

print(fastq_files)

n = 4
reads = {}
read_lengths = []

for fn in fastq_files:
    with open(fn, 'r') as fh:
        lines= []
        for line in fh:
            #if len(reads) < 100: #take the first n reads
                lines.append(line.strip('\n'))
                if len(lines) == n:
                    record = process(lines)
                    #sys.stderr.write("Record: %s\n" % (str(record)))
                    #print(record)
                    lines = []
                    read_lengths.append(len(record["sequence"]))
                    reads[record["name"]] = record


#Find barcodes
found_barcodes = defaultdict(int)

for key, value in reads.items():
	sequence = str(Seq(value["sequence"]).reverse_complement())
	str_start = sequence.find("TGTGGCAGAAGAAGCCACGTTAA")
	primer_start = sequence.find("GCGGCCGCGGAT")
	length_ = len("TGTGGCAGAAGAAGCCACGTTAA")
	if str_start != -1 and primer_start != -1: 
		found_barcodes[sequence[str_start+length_:primer_start]] += 1

#Use the lookup table to determine the mutations
barcode_lookup = {}

with open("LookupTable_N.csv", "r") as f:
	next(f)
	for line in f:
		barcode_lookup[line.split(",")[0]] = line.split(",")[1].strip("\n")

#all barcodes
found_ = 0
not_found_ = 0

with open("found.txt", "w") as found, open("not_found.txt", "w") as notfound:
	notfound.write("barcode\tcount\n")
	found.write("barcode\tmutation\tcount\n")
	for key, value in found_barcodes.items():
		if barcode_lookup.get(key):
			found_ += value
			found.write(key + "\t" + barcode_lookup[key] + "\t" + str(value) + "\n")
		else:
			not_found_ += value
			notfound.write(key + "\t" + str(value) + "\n")

