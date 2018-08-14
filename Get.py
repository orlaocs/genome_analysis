# To Get the Coding Regions from a DNA sequence in FASTA Format:

import sys
import re

def find_all_starts(sequence):
    starts = []
    i = sequence.find("ATG")
    while i >= 0:
        starts.append(i)
        i = sequence.find("ATG",i+1)
    return starts

def find_first_in_register_stop(sequence):
    i = 0
    while i < len(sequence)-2 and sequence[i:i+3] not in ("TAA","TAG","TGA"):
        i += 3
    if i < len(sequence)-2:
        return i+3
    else:
        return -1

def all_orfs(filename):
    orfs = []
    with open(filename, "r") as TEXT:
        content = TEXT.readlines()[1:]
        sequence = "".join(content).replace("\r","").replace("\n","").replace("N","")
        TEXT.close()
        starts_inds = find_all_starts(sequence)
        stop_inds = []
        for start in start_inds:
            relative_stop = find_first_in_register_stop(sequence[start:])
            #if relative_stop != -1:
             #   stop = start + relative_stop
                if stop not in stop_inds:
                    orfs.append((relative_stop,start,stop))
                    stop_inds.append(stop)
        orfs = sorted(orfs,reverse=True)
        for i, orf in enumerate(orfs):
            orfs[i] = (orf[1], orf[2])
    with open("orfs.fasta","w") as OUTPUT:
        for x in orfs:
            OUTPUT.write(str(x))
    return orfs

all_orfs(sys.argv[1])
