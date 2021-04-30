import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from biomek.misc import file
from dna_features_viewer import BiopythonTranslator


class MyCustomTranslator(BiopythonTranslator):
    """Custom translator implementing the following theme:

    - Color terminators in green, CDS in blue, all other features in gold.
    - Do not display features that are restriction sites unless they are BamHI
    - Do not display labels for restriction sites
    - For CDS labels just write "CDS here" instead of the name of the gene.

    """

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            return "blue"
        elif feature.type == "terminator":
            return "green"
        else:
            return "gold"

    def compute_feature_label(self, feature):
        if feature.type == 'restriction_site':
            return None
        elif feature.type == "CDS":
            return "CDS here"
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        """Do not display promoters. Just because."""
        return [
            feature for feature in features
            if (feature.type != "restriction_site")
            or ("BamHI" in str(feature.qualifiers.get("label", '')))
        ]


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
    sequence_object = Seq(sequence_string)

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


def add_anootation(genbank_file, annotation):
    # Add annotation
    feature = SeqFeature(FeatureLocation(start=3, end=12), type='misc_feature')
    genbank_file.features.append(feature)

    # Save as GenBank file
    output_file = open('biomek/output/' + str(genbank_file)+'.gb', 'w')
    SeqIO.write(genbank_file, output_file, 'genbank')


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