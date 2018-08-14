##Kmercounter Assignment:

# (1) Length Function:
# This function computes the total amount of bases in a file, excluding the
# header line if it is written in FASTA format.

def length(filename):
    with open(filename, 'r') as TEXT:
        sequence = TEXT.read()
        TEXT.close()
        A = []
        T = []
        C = []
        G = []
        total_bases = []
        for line in sequence:
            if line[0] == '>':
                continue
            A = sequence.count("A")
            T = sequence.count("T")
            C = sequence.count("C")
            G = sequence.count("G")
            total_bases = A + T + C + G
    return total_bases

# (2) Complement Function:
# This fumction computes the reverse strand of a fasta file.
# When it opens the file, it immediately gets rid of the header line
# Later getting rid of the \r and \n and N symbols.
# Then, it returns the reverse complement of the remaining fasta file.

def complement(filename):
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    with open(filename, 'r') as text:
            content = text.readlines()[1:]
            sequence = "".join(content).replace("\r","").replace("\n","").replace("N","")
            text.close()
            return "".join([complement[base] for base in sequence[::-1]])

# (3) Translate Function:
# NOTE: Start codon: AUG & Stop codons: UAA, UAG, UGA.

from itertools import takewhile

def translate(filename):
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
    with open(filename, 'r') as text:
        content = text.readlines()[1:]
        sequence = "".join(content).replace("\r","").replace("\n","").replace("N","")
        text.close()
        start_codon = sequence.find("ATG")
        start = sequence[start_codon:]
        stop_codons = ("TAA", "TAG", "TGA")
        codons = [start[i:i+3] for i in range(0, len(start), 3)]
        coding_sequence = takewhile(lambda x: x not in stop_codons and len(x) == 3, codons)
        protein_sequence = "".join([codon_table[codon] for codon in coding_sequence])
        return protein_sequence


##https://stackoverflow.com/questions/19521905/translation-dna-to-protein
        
# (4) Di-nucleotide Function:

#then merge the lines
#define a way of identifying the dinucleotides and count how many there are.

##To start:
def dinulc(filename):
    with open(filename, 'r') as text:
        content = text.readlines()[1:]
        sequence = "".join(content).replace("\r","").replace("\n","").replace("N","")
        text.close()
        for line in sequence:
            seq += line.strip()
            return seq
            

##Not working:
def dinulc(filename):
    ecnt = {}
    seq = ""
    with open(filename, "r") as seq_data:
    # First we want to merge the lines so that the result will include the
    # nucleotides at the start ans end of each of the lines. We do this by
    # merging the values into a string of data.
        for line in seq_data:
            seq += line.strip()
        # By slicing, we can disern the dinucleotides in the sequence. Here,
        # we look at both the current nucleotide and the nucleotide following
        # it.
            for i in range(len(seq)-1):
                dinuc = seq[i:i+2]
            #We can then count the amount of dinucleotides and store the
            # the results to the vector; ecnt.
                if dinuc in ecnt:
                    ecnt[dinuc] += 1
                else:
                    ecnt[dinuc] = 1
                return ecnt

##
from collections import defaultdict

filename = "Sample.txt"

dinucleotide_counts = defaultdict(int)

sequence = ""

with open(filename) as file:
    for line in file:
        sequence += line.strip()

for i in range(len(sequence)-1):
    dinucleotide_counts(sequence[i:i+2] += 1

for key, value in sorted(dinucleotide_counts.items()):
                        print(key, value)



####
print(

                    
