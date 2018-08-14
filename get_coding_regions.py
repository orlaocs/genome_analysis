# To Get the Open Reading Frames (ORFs) From a DNA Sequence in a FASTA File:
import sys
import re

## (1) To Find All the Start Codons in the Sequence:

def find_all_start_codons(sequence):
    all_start_codons = []
    #The find() function returns the index for the lowest index of the substring, if found.
    next_start_codon = sequence.find("ATG")
    while next_start_codon >= 0:
        #Once the lowest index value of a start codon is found, the value is appended to "all_start_codons".
        all_start_codons.append(next_start_codon)
        #The value for the the "next_start_codon" is increased to find the next start codon.
        next_start_codon = sequence.find("ATG", next_start_codon + 1)
    #Once all the values are found, all the start codons are returned.
    return all_start_codons

## (2) To Find the Next Stop Codon in the Sequence:

def find_first_stop_codon(sequence):
    i = 0 
    #While the number of "i" is less than the length of the sequence minus two & the values in the sequence do not correlate to any of the stop codons, increase "i" by three.
    while i < len(sequence)-2 and sequence[i:i+3] not in ("TAA","TAG","TGA"):
        i += 3
    #If before the end of the sequence, a stop codon is found, return the the stop codon.
    if i < len(sequence)-2:
        return i+3
    #Else failed to find a stop codon.
    else:
        return -1

## (3) To Generate all the Open Reading Frames from the Sequence:

def all_orfs(filename):
    orfs = []
    with open(filename, "r") as TEXT:
        content = TEXT.readlines()[1:]
        sequence = "".join(content).replace("\r","").replace("\n","").replace("N","")
        TEXT.close()
        #To get all index positions of start and stop codons.
        start_indexes = find_all_start_codons(sequence)
        stop_indexes = []
        #For each start codon, find the next stop codon in the sequence
        for start in start_indexes:
            relative_stop_codon = find_first_stop_codon(sequence[start:])
            if relative_stop_codon != -1:
                #Retrieves the index of the stop codon.
                stop = start + relative_stop_codon
                #If already had a stop codon, a longer ORF contains this one.
                if stop not in stop_indexes:
                    orfs.append((relative_stop_codon, start, stop))
                    stop_indexes.append(stop)
        #To get a sorted list of the ORF length
        orfs = sorted(orfs, reverse = False)
        #To remove the values for the lengths of the of each of the ORFs.
        for i, orf in enumerate(orfs):
            orfs[i] = (orf[1], orf[2])
        orfs = sorted(orfs)
    #The write() function requires that the output information be in the form of a string in the form of str().
    with open("orfs.txt","w") as OUTPUT:
        for x in orfs:
            OUTPUT.write(str(x))
    return orfs

## (4) To Generate the DNA sequence of the ORFs from a fasta file:

def get_orf_seq(filename):
    orfs =[]
    with open(filename, "r") as TEXT:
        content = TEXT.readlines()[1:]
        sequence = "".join(content).replace("\n","").replace("\r","").replace("N","")
        TEXT.close()
        start_indexes = find_all_start_codons(sequence)
        stop_indexes = []
        for start in start_indexes:
            relative_stop_codon = find_first_stop_codon(sequence[start:])
            if relative_stop_codon != -1:
                stop = start + relative_stop_codon
                if stop not in stop_indexes:
                    orfs.append((relative_stop_codon, start, stop))
                    stop_indexes.append(stop)
        orfs = sorted(orfs, reverse = False)
        for i, orf in enumerate(orfs):
            orfs[i] = [orf[1], orf[2]]
        orfs = sorted(orfs)
        coding_sequence = []
        #For each row in the "orfs" variable, extract the values for the start codon in column 1 and the values for the stop codon in column 2 and append to the coding_sequence variable with each coding sequence being separated by a line.
        for row in orfs:
            coding_sequence.append((sequence[row[0]:row[1]]))
    with open("coding_sequence.txt", "w") as OUTPUT:
        for x in coding_sequence:
            OUTPUT.write(x)
    return coding_sequence
        
## (5) To Convert the OUTPUT from all_orfs() to a Basic BED Format File:
# for the amount of lines in the sequence generate the same amount of numbers
def convert_to_bed(filename):
    orfs = []
    with open(filename, "r") as TEXT:
        content = TEXT.readlines()[1:]
        sequence = "".join(content).replace("\r","").replace("\n","").replace("N","")
        TEXT.close()
        start_indexes = find_all_start_codons(sequence)
        stop_indexes = []
        for start in start_indexes:
            relative_stop_codon = find_first_stop_codon(sequence[start:])
            if relative_stop_codon != -1:
                stop = start + relative_stop_codon
                if stop not in stop_indexes:
                    orfs.append((relative_stop_codon, start, stop))
                    stop_indexes.append(stop)
        #To get the start and stop positions in a sorted variable.
        start_pos = sorted(start_indexes)
        stop_pos = sorted(stop_indexes)
        #To Calculate the Number of of ORFs:
        num_orfs = len(stop_pos)
        #To generate the gene names:
        all_nums = list(range(1, num_orfs))
        all_genes = []
        for i in all_nums:
            all_genes += ["Gene%s" % (i)]
        #To generate all the chromosome locations and adding a new line character to begin the format of the bed file:
        all_chrom = ["\nchrom1"] * len(stop_pos)
        #To compile all the values into the one variable:
        info_tuple = zip(all_chrom, all_genes, start_pos, stop_pos)
        #To convert tuple to list of list:
        bed_info = [list(element) for element in info_tuple]
        #To convert list of list to a list:
        bed_info = [item for sublist in bed_info for item in sublist]
        bed_info = list(bed_info)
        #To convert to a string:
        bed_info = ",".join(str(i) for i in bed_info)
        #To replace all the commas with a tab and to remove the empty line at the start of the file:
        bed_info = bed_info.replace(",","\t")[1:]
        #The headings are added to the start of the file:
        bed_info = ("chrom\tname\tchromStart\tchromEnd\n" + bed_info)
    #The modifications are outputed as a .bed file.
    with open("orfs.bed","w") as OUTPUT:
        for x in bed_info:
            OUTPUT.write(x)
    return bed_info

## (8) To Get All 6 Possible Frames From an Open Reading Frame:

def all_6_orfs(filename):
    with open(filename) as TEXT:
        sequence = TEXT.readlines()
        sequence = "".join(sequence).replace("TAG","TAGNNN").replace("TGA","TGANNN").replace("TAA","TAANNN")
        TEXT.close()
        ##To get the first orf:
        ORF_1 = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        ORF_1 = [x for x in ORF_1 if "N" not in x]
        ORF_1 = "".join(ORF_1).replace("TAG,","TAG\n").replace("TGA,", "TGA\n").replace("TAA,","TAA\n")
        ORF_1 = ORF_1.split("\n")
        print ORF_1
        ##To get the second orf:
        ORF_2 = [sequence[i:i+3] for i in range(1, len(sequence)-1, 3)]
        ORF_2 = [x for x in ORF_2 if "N" not in x]
        ORF_2 =  "".join(ORF_2).replace("TAG,","TAG\n").replace("TGA,", "TGA\n").replace("TAA,","TAA\n")
        ORF_2 = ORF_2.split("\n")
        ##To get the third orf:
        ORF_3 = [sequence[i:i+3] for i in range(2, len(sequence)-2, 3)]
        ORF_3 = [x for x in ORF_3 if "N" not in x]
        ORF_3 =  "".join(ORF_3).replace("TAG,","TAG\n").replace("TGA,", "TGA\n").replace("TAA,","TAA\n")
        ORF_3 = ORF_3.split("\n")
    with open("orf_1.fasta", "w") as OUTPUT:
        count = 1
        for line in ORF_1:
            line = line.rstrip("\n")
            OUTPUT.write(">" + "gene" + str(count) + "\n")
            OUTPUT.write(line + "\n")
            count = count + 1
        OUTPUT.close()
    #return ORF_1
    with open("orf_2.fasta","w") as OUTPUT:
        count = 1
        for line in ORF_2:
            line = line.rstrip("\n")
            OUTPUT.write(">" + "gene" +  str(count) + "\n")
            OUTPUT.write(line + "\n")
            count = count + 1
        OUTPUT.close()
    #return ORF_2
    with open("orf_3.fasta","w") as OUTPUT:
        count = 1
        for line in ORF_3:
            line = line.rstrip("\n")
            OUTPUT.write(">" + "gene" + str(count) + "\n")
            OUTPUT.write(line + "\n")
            count = count + 1
        OUTPUT.close()
    #return ORF_3

## (9) To Get the Protein Sequence of the ORFs:

def translate_orf(filename):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T',
        'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N',
        'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R',
        'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H',
        'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R',
        'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V',
        'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G',
        'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S',
        'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L',
        'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}
    with open(filename, "r") as TEXT:
        content = TEXT.readlines()[1::2]
        sequence = "".join(content).replace("\n","")
        protein_sequence = ""
        for i in range(0, len(sequence), 3):
            protein_sequence += codon_table[sequence[i:i+3]]
            continue
        protein_sequence = "".join(protein_sequence).replace("_", "_\n")
        protein_sequence = protein_sequence.split()
        #To remove the empty strings in the list of strings:
        protein_sequence = filter(None, protein_sequence)
    with open("protein_seq_orf3.fasta","w") as OUTPUT:
        count = 1
        for line in protein_sequence:
            line = line.rstrip("\n")
            OUTPUT.write(">" + "protein" + str(count) + "\n")
            OUTPUT.write(line + "\n")
            count = count + 1
        OUTPUT.close()
        #print protein_sequence

## (10) To Get the Length of Each Protein and Return it to a New File:
def len_genes(filename):
    with open(filename, "r") as TEXT:
        sequence = TEXT.readlines()[1::2]
        sequence = "".join(sequence).replace("\n","")
        sequence = sequence.split("_")
        len_seq = [len(i) for i in sequence]
        len_seq = ",".join(str(x) for x in len_seq)
        len_seq = len_seq.split(",")
    with open("len_seq_orf3.fasta","w") as OUTPUT:
        count = 1
        for length in len_seq:
            length = length.rstrip(",")
            #OUTPUT.write(">" + "protein" + str(count) + "\n")
            OUTPUT.write(length + "\n")
            count  = count + 1
        OUTPUT.close()
        
#all_orfs(sys.argv[1])
#get_orf_seq(sys.argv[1])
#convert_to_bed(sys.argv[1])
all_6_orfs("coding_sequence.txt")
translate_orf("orf_3.fasta")
len_genes("protein_seq_orf3.fasta")
