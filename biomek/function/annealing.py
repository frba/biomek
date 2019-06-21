import os

from Bio import SeqIO
from Bio.Seq import Seq

from misc import genbank

COUNT = 0

class Primer:
    """Class representing a sequencing primer

    Parameters
    ----------
    name
      Name (label) of the primer

    sequence
      An ATGC string

    metadata
      A dictionnary {field, value}
    """

    def __init__(self, name=None, sequence=None, metadata=None):
        """Initialize"""
        self.name = name
        self.sequence = sequence
        self.metadata = {} if metadata is None else metadata


    def find_location(self, sequence):
        """Return the (start, end, strand) of the primer in the sequence.

        The convention is that always (start < end), the strand indicates
        the direction of the homology.
        """
        ind = sequence.find(self.sequence)
        strand = 1
        if ind == -1:
            ind = sequence.find(reverse_complement(self.sequence))
            if ind == -1:
                return None
            else:
                # ind = len(sequence) - ind - len(self.sequence)
                strand = -1
        start, end = ind, ind + len(self.sequence)
        self.metadata = start, end, strand
        return start, end, strand


def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.
    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.
    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


def verify_primer_annealing(gb_file, selected_primers):
    for primer in selected_primers:
        if primer.find_location(gb_file.seq) is None:
            # print(str(primer.name) + ' cant be annealing in ' + str(gb_file.name))
            return False
    return True



def get_product_length(gb_file, selected_primers):
    start_fwd = 0
    end_rvs = 0

    for primer in selected_primers:
        start, end, strand = primer.find_location(gb_file.seq)
        # print(primer.name, primer.metadata, end=', ')
        if primer.name.__contains__('_Fp'):
            start_fwd = start
        else:
            end_rvs = end

    length = end_rvs - start_fwd
    if length > len(primer.sequence):
        print('assembly:' + str(gb_file.name))
        for primer in selected_primers:
            print(primer.name, end=', ')
        print('\nLength:' + str(length) + '\n')


def get_gbfiles(constructs_path):
    records = [
        genbank.load_record(os.path.join(constructs_path, f), linear=False)
        for f in sorted(os.listdir(constructs_path)) if f.endswith('.gb')
    ]
    return records


def get_primers(fasta_file):
    """Return a list of Primer objects, by reading a FASTA file."""
    records = SeqIO.parse(fasta_file, 'fasta')
    return [
        Primer(name=rec.name, sequence=str(rec.seq))
        for rec in records
    ]


def get_selected_primers(fasta_file, selected_primers_names):
    records = SeqIO.parse(fasta_file, 'fasta')
    return [
        Primer(name=rec.name, sequence=str(rec.seq))
        for rec in records
        for primer_name in selected_primers_names if rec.name == primer_name
    ]


def get_fwd_primers(fasta_file):
    records = SeqIO.parse(fasta_file, 'fasta')
    return [
        Primer(name=rec.name, sequence=str(rec.seq))
        for rec in records if rec.name.__contains__('_Fp')
    ]


def get_rvs_primers(fasta_file):
    records = SeqIO.parse(fasta_file, 'fasta')
    return [
        Primer(name=rec.name, sequence=str(rec.seq))
        for rec in records if rec.name.__contains__('_Rp')
    ]


#get_primers
available_primers_path = '../input/annealing/available_primers.fa'
primers = get_primers(available_primers_path)


def get_length_from_seleted_primers(selected_primers_names, gb_file):
    selected_primers = get_selected_primers(available_primers_path, selected_primers_names)
    if verify_primer_annealing(gb_file, selected_primers):
        get_product_length(gb_file, selected_primers)


''' get length of a selected pair of primers in gb files'''
# selected_primers_names = ['ConS_scar_Fp','Con1_scar_Rp']
# gb_files = get_gbfiles('../input/annealing/constructors/')
# for gb_file in gb_files:
#     get_length_from_seleted_primers(selected_primers_names, gb_file)


'''get all pairs of primers in available file and return the length using .gb files'''
fwd_primers = get_fwd_primers(available_primers_path)
rvs_primers = get_rvs_primers(available_primers_path)
gb_files = get_gbfiles('../input/annealing/constructors/')
for gb_file in gb_files:
    for fwd in fwd_primers:
        for rvs in rvs_primers:
            selected_primers_names = [fwd.name, rvs.name]
            get_length_from_seleted_primers(selected_primers_names, gb_file)