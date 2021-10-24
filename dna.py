#load python third party modules
import random

#load your files or data here this code will read it from a text file or other file format:
#dna = "ATGCTAGC"

inputfile = 'data/dna_sequence_second.txt'
#inputfile = input("Enter the name of your file: ")

def read_seq_data(inputfile):
    with open(inputfile, "r") as f:
        dna_seq = f.read()
    dna_seq = dna_seq.replace("\n", "")
    dna_seq = dna_seq.replace("\r", "")
    return dna_seq

dna = read_seq_data(inputfile)


# dna class holds every function for dna analysis
class Dna:
    def __init__(self):
        self.dna = dna
        self.length = len(dna)
        self.count_A = dna.count("A")
        self.count_T = dna.count("T")
        self.count_G = dna.count("G")
        self.count_C = dna.count("C")
        self.sort = sorted(dna)
        self.list = list(dna)
        self.tuple = tuple(dna)
        self.set = set(dna)
        self.transcripted = dna.replace("T", "U")
        self.name = name = 'Hok'

    def Seq(self, dna):
        return self.dna

    def count_nucleotides(self, dna):
        output = f"A: {self.count_A}, T:{self.count_T}, G:{self.count_G}, C:{self.count_C}"
        return output

    def reverse_seq(self, dna):
        return dna[::-1]

    def transcription(self, dna):
        return self.transcripted

    def indexing_seq(self, dna):
        for index, letter in enumerate(dna):
            print("%i %s" % (index, letter))

    def length_of_seq(self, dna):
        output = f"length of seq is {self.length}"
        return output

    def gc_content(self, dna):
        output = f"GC Content is {round(float(self.count_C + self.count_G) / self.length * 100)}%"
        return output

    def at_Content(self, dna):
        output = f"AT Content is {round(float(self.count_A + self.count_C) / self.length * 100)}%"
        return output

    def sort_data(self, dna):
        return self.sort

    def get_list(self, dna):
        return self.list

    def get_tuple(self, dna):
        return self.tuple

    def get_set(self, dna):
        return self.set

    def generate_random_neclueotides(self, neclueotides):
        neclueotides = ["A", "T", "G", "C"]
        for x in range(18):
            random_neclueotides = random.choice(neclueotides)
            lis = random_neclueotides
            print(lis, end='')

    def complement_strand(self, dna):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        bases = list(self.dna)
        bases = [complement[base] for base in bases]
        return ''.join(bases)

    def rev_complement(self, dna):
        return seq.complement_strand(self.dna[::-1])

    # Generate protein sequence
    def translate(self, dna):
        table = {"TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
                 "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
                 "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
                 "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
                 "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
                 "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
                 "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
                 "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
                 "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
                 "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
                 "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
                 "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
                 "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
                 "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
                 "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
                 "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
                 }

        protein = ""
        if len(dna) % 3 == 0:
            for i in range(0, len(dna), 3):
                codon = dna[i:i + 3]
                protein += table[codon]
        return protein

    def coding_template_strand(self, dna, complement_strand):
        print(dna + "\n" + f"{''.join(['|' for c in range(len(dna))])}" "\n" + complement_strand)


    def helix_strand(self, dna):
        for helix in dna:
            print(helix)

    def complete_helix_strand(self, complement_strand):
        for hel in complement_strand:
            print(hel)


    def fasta_format(self, name, dna):
        with open("output/fasta_file.fasta", "w") as fasta:
            fasta.write(">" + name)
            fasta.write("\n" + dna)
            fasta.close()
            print("data have been saved.")


seq = Dna()


'''
print(seq.Seq(dna))
print(seq.count_nucleotides(dna))

print(seq.reverse_seq(dna))
print(seq.transcription(dna))

#print(seq.indexing_seq(dna))
#print(seq.generate_random_neclueotides(dna))

print(seq.length_of_seq(dna))
print(seq.gc_content(dna))

print(seq.at_Content(dna))

print(seq.sort_data(dna))
print(seq.get_list(dna))
print(seq.get_tuple(dna))
print(seq.get_set(dna))

print(dna)
print(seq.complement_strand(dna))

print(seq.rev_complement(dna))
print(seq.translate(dna[0:66]))

print(seq.helix_strand(dna))
'''
# print(seq.complete_helix_strand(dna))
#print(dna + "\n" + f"{''.join(['|' for c in range(len(dna))])}" "\n" + seq.complement_strand(dna))
#

# print(dna.split())

# print(seq.fasta_format('test', "AAAAGTGCAACCGATC"))

print(seq.coding_template_strand("TTTC", "AAAG"))