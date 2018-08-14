# To Get the Summary Statistics of a Sequence:
#The following results will be generated using the protein_sequence.fasta file.

import re

## (1) To Get the Number of Genes Per Reading Frame:

def genes_per_frame(filename):
    with open(filename, "r") as TEXT:
        sequence = TEXT.readlines()[-2]
        TEXT.close()
        num_genes = re.findall("%s([0-9]*)%s" % (">protein", "\n"), sequence)
        num_genes = "".join(num_genes)
        print num_genes
#Can also use len() to get the amount of items in a file when in a list.

## (2) To Get the Length Distance of Genes:

#First to define a function that will return the length of a sequence:
def get_length(seq):
    return len(seq)

#Then to create a function that will generate a function: "sorted_sequences" based on the length of each of the protein coding sequences:
def len_dist_genes(filename):
    with open(filename, "r") as TEXT:
        sequence = TEXT.readlines()[1::2]
        sequence = "".join(sequence).replace("\n","")
        sequence = sequence.split("_")
        #sort the seqeuences:
        sorted_sequences = sorted(sequence, key=get_length)
        print sorted_sequences

## (3) To Get the Length of Each Protein and return it to a new file:
def len_genes(filename):
    with open(filename, "r") as TEXT:
        sequence = TEXT.readlines()[1::2]
        sequence = "".join(sequence).replace("\n","")
        sequence = sequence.split("_")
        #To remove the empty strings from the list of strings:
        #sequence = filter(None, sequence)
        #To get the length of each string in the list of strings in "sequence":
        len_sequences = [len(i) for i in sequence]
        len_sequences = ",".join(str(x) for x in len_sequences)
        len_sequences = len_sequences.split(",")
        #print len_sequences
    with open("len_protein_sequence.fasta","w") as OUTPUT:
        count = 1
        for length in len_sequences:
            length = length.rstrip(",")
            OUTPUT.write(">" + "protein" + str(count) + "\n")
            OUTPUT.write(length + "\n")
            count = count + 1
        OUTPUT.close()

#genes_per_frame("protein_sequence.fasta")
#len_dist_genes("protein_sequence.fasta")
len_genes("protein_sequence.fasta")
