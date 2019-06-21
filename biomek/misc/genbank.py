import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

from misc import file


def generate_from_csv(path):
    filein = file.verify(path)
    csv_file = file.create_reader_csv(filein)
    header = next(csv_file)

    for row in csv_file:
        name = row[0]
        description = row[1]
        sequence = row[3]
        generate_from_sequence(sequence, name, description)


def generate_from_sequence(sequence, name, description):
    # Create a sequence
    sequence_string = sequence
    sequence_object = Seq(sequence_string, IUPAC.unambiguous_dna)

    # Create a record
    record = SeqRecord(sequence_object,
                       id='123456789',  # random accession number
                       name=name,
                       description=description)

    # Add annotation
    # feature = SeqFeature(FeatureLocation(start=3, end=12), type='misc_feature')
    # record.features.append(feature)

    # Save as GenBank file
    output_file = open('biomek/output/' + str(name)+'.gb', 'w')
    SeqIO.write(record, output_file, 'genbank')


def load_record(filename, linear=True, name='auto'):
    if filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(('fa', 'fasta')):
        record = SeqIO.read(filename, "fasta")
    else:
        raise ValueError('Unknown format for file: %s' % filename)
    record.linear = linear
    if name == 'auto':
        name = os.path.splitext(os.path.basename(filename))[0]
    record.id = name
    record.name = name.replace(" ", "_")[:35]
    return record