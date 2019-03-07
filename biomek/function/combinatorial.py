from ..misc import file
import itertools, os


def get_sets_in_filepath(reader):
    lists_combination_parts = []
    ''' For each line in file get the list'''
    for line in reader:
        parts = []
        sets = line.strip('\n').split(';')
        '''List of parts'''
        for set in sets:
            parts.append(set.split(','))
        ''' Create the combinations from the list'''
        lists_combination_parts.append(list(itertools.product(*parts)))
    return lists_combination_parts


def create_combinations(filepath):
    filein = file.verify(filepath)

    """Create combinations"""
    lists_combination_parts = get_sets_in_filepath(filein)

    """Write a output file"""
    fileout = file.create('biomek/output/combination_' + str(os.path.basename(filepath)), 'w')
    file.write_combinations(fileout, lists_combination_parts)