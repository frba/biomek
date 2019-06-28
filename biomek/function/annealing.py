import os

from Bio.Seq import Seq
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator
from Bio.SeqFeature import FeatureLocation
from Bio import SeqFeature

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
                ind = len(sequence) - ind - len(self.sequence)
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
            return False
    return True


def get_product_length(gb_file, selected_primers):
    start_fwd = 0
    end_rvs = 0
    length = 0
    for primer in selected_primers:
        start, end, strand = primer.find_location(gb_file.seq)
        # print(primer.name, start, end, strand)
        if primer.name.__contains__('_Fp'):
            start_fwd = start
            # print(start_fwd)
        else:
            end_rvs = len(gb_file.seq) - end + len(primer.sequence)
            # print(end_rvs)

    length = end_rvs - start_fwd

    if length < 0:
        length = length + len(gb_file.seq)

    return length
    # if length > len(primer.sequence):
    #     return length
    # else:
    #     return None


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
        # for rec in records if rec.name.__contains__('_Fp') and not rec.name.__contains__('ConE')
    ]


def get_rvs_primers(fasta_file):
    records = SeqIO.parse(fasta_file, 'fasta')
    return [
        Primer(name=rec.name, sequence=str(rec.seq))
        for rec in records if rec.name.__contains__('_Rp')
        # for rec in records if rec.name.__contains__('_Rp') and not rec.name.__contains__('ConS')
    ]


def draw_gbfile(gb_file):
    graphic_record = genbank.MyCustomTranslator().translate_record(gb_file)
    ax, _ = graphic_record.plot(figure_width=10)
    ax.figure.tight_layout()
    ax.figure.savefig("custom_biopython_translator.png")


def detect_tu(selected_primers):

    if selected_primers[0].name.__contains__('_Fp'):
        fwd_primer = selected_primers[0].name
        rvs_primer = selected_primers[1].name
    else:
        fwd_primer = selected_primers[1].name
        rvs_primer = selected_primers[0].name

    if fwd_primer.__contains__('ConS') and rvs_primer.__contains__('Con1'):
        return 'TU1'
    elif fwd_primer.__contains__('ConS') and rvs_primer.__contains__('Con2'):
        return 'TU1-TU2'
    elif fwd_primer.__contains__('ConS') and rvs_primer.__contains__('ConE'):
        return 'TU1-TU2-TU3'
    elif fwd_primer.__contains__('Con1') and rvs_primer.__contains__('Con2'):
        return 'TU2'
    elif fwd_primer.__contains__('Con1') and rvs_primer.__contains__('ConE'):
        return 'TU2-TU3'
    elif fwd_primer.__contains__('Con2') and rvs_primer.__contains__('ConE'):
        return 'TU3'
    else:
        return None


def print_output(results):
    print('Assembly' + '\t' + 'Primer_FWD' + '\t' + 'Primer_RVS' + '\t' + 'TU' + '\t' + 'Size')
    for result in results:
        gb_file, selected_primers, length = result
        tu_type = detect_tu(selected_primers)
        if tu_type is not None:
            if selected_primers[0].name.__contains__('_Fp'):
                print(str(gb_file.name) + '\t'
                      + str(selected_primers[0].name) + '\t'
                      + str(selected_primers[1].name) + '\t'
                      + str(tu_type) + '\t'
                      + str(length))
            else:
                print(str(gb_file.name) + '\t'
                      + str(selected_primers[1].name) + '\t'
                      + str(selected_primers[0].name) + '\t'
                      + str(tu_type) + '\t'
                      + str(length))


#get_primers
available_primers_path = '../input/annealing/available_primers.fa'
primers = get_primers(available_primers_path)


# ''' get length of a selected pair of primers in gb files'''
# selected_primers_names = ['ConS_scar_Fp','ConE_scar_Rp']
# gb_files = get_gbfiles('../input/annealing/constructors/')
# result = []
# for gb_file in gb_files:
#     selected_primers = get_selected_primers(available_primers_path, selected_primers_names)
#     if verify_primer_annealing(gb_file, selected_primers):
#         length = get_product_length(gb_file, selected_primers)
#         if length is not None:
#             result.append([gb_file, selected_primers, length])
#
# print_output(result)
            # draw_gbfile(gb_file)


'''get all pairs of primers in available file and return the length using .gb files'''
fwd_primers = get_fwd_primers(available_primers_path)
rvs_primers = get_rvs_primers(available_primers_path)

gb_files = get_gbfiles('../input/annealing/constructors/')
result = []
for gb_file in gb_files:
    for fwd in fwd_primers:
        for rvs in rvs_primers:
            selected_primers_names = [fwd.name, rvs.name]
            selected_primers = get_selected_primers(available_primers_path, selected_primers_names)
            # if selected_primers[0].name.__contains__('_Fp'):
            #     print(selected_primers[0].name,selected_primers[1].name)
            # else:
            #     print(selected_primers[1].name, selected_primers[0].name)
            if verify_primer_annealing(gb_file, selected_primers):
                length = get_product_length(gb_file, selected_primers)
                if length is not None:
                    result.append([gb_file, selected_primers, length])
                    # print(gb_file.name, selected_primers[0].name, selected_primers[1].name, length)

print_output(result)
                    # print(gb_file.name, selected_primers[0].name, selected_primers[1].name, length)



                    # if selected_primers[0].name.__contains__('_Fp'):
                    #     start = SeqFeature.AfterPosition(int(selected_primers[0].metadata[0]))
                    #     end = SeqFeature.BeforePosition(int(selected_primers[1].metadata[1]))
                    # else:
                    #     start = SeqFeature.AfterPosition(int(selected_primers[1].metadata[0]))
                    #     end = SeqFeature.BeforePosition(int(selected_primers[0].metadata[1]))
                    # print(start, end)
                    # f1 = FeatureLocation(start, end, strand=+1)

                    # for feat in gb_file.features:
                        # if feat.location in f1:
                        #     print(feat)


