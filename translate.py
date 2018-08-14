#To Translate a DNA Sequence from a File in FASTA Format and Also to Generate the Exons from that Protein Sequence:
import sys
import re
from itertools import takewhile

# (1) To Generate a Protein Sequence from a DNA Sequence:

def translate(filename):
    #To Create a Dictionary of all the codons and their corresponding Amino Acid.
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
        'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
        'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',      
}
    with open(filename, 'r') as TEXT:
        content = TEXT.readlines()[1:]
        sequence = "".join(content).replace("\r","").replace("\n","").replace("N","")
        TEXT.close()
        protein_sequence = ""
        coding_protein = ""
        #To find the start codon with the lowest index value in the string (IE the first start codon in the string).
        start_codon = sequence.find("ATG")
        #To modify the sequence wherein the sequence starts with the first start codon found.
        start = sequence[start_codon:]
        #For each codon in the sequence, convert it into its corresponding amino acid and add the resulting protein to a new variable; "protein_sequence".
        for i in range(0, len(start)-2, 3):
            protein_sequence += codon_table[start[i:i+3]]
            continue
    #Outputting the vector as a file in my directory:
    with open("protein_sequence.fasta", "w+") as OUTPUT:
        OUTPUT.write(protein_sequence)
    return protein_sequence

# (2) To Generate all the Exons from a Protein Sequence:

def exons(protein_sequence):
    with open(protein_sequence, 'r') as TEXT:
        sequence = TEXT.read()
        TEXT.close()
        pos = ""
        start = "M"
        stop = "-"
        pos = re.findall("%s[A-Z]*%s" % (start, stop), sequence)
#NOTE: Instead, can also use the following for adding back in the start and the stop codons:
# pos = map((lambda x: x + "-"), pos) 
# pos = [i + "-" for i in pos]
#or an even more concise way is to the following:
# pos = ["M" + i + "-" for i in pos]
#However, when we deleted the brackets in the findall() function, it printed the start and stop codon. The original script is as follows:
# pos = re.findall("%s([A-Z])*%s" % (start, stop), sequence)
        #print pos
    return pos

translate(sys.argv[1])
exons("protein_sequence.fasta")

#bed format for co-ordinates of the genes. gene feature map.
#fasta nucleotide, protein
