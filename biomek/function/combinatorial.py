from ..misc import file, calc
import itertools, os


def get_sets_in_filepath(reader):
    lists_parts = []
    lists_combination_parts = []
    ''' For each line in file get the list'''
    for line in reader:
        parts = []
        sets = line.strip('\n').split(';')
        '''List of parts'''
        for set in sets:
            parts.append(set.split(','))
        ''' Create the single list of parts'''
        lists_parts.append(list(parts))
        ''' Create the combinations from the list'''
        lists_combination_parts.append(list(itertools.product(*parts)))
    return lists_combination_parts, lists_parts


def create_combinations(filepath):
    filein = file.verify(filepath)

    """Create combinations"""
    lists_combination_parts, lists_parts = get_sets_in_filepath(filein)

    """Calculate the num of parts in input file"""
    list_set_num_parts = calc.num_listsparts(lists_parts)

    """Calculate number of combinations"""
    num_combinations_in_list = calc.num_combinations(lists_combination_parts)

    """Write a output file"""
    fileout = file.create('biomek/output/combination_' + str(os.path.basename(filepath)), 'w')
    file.write_combinations(fileout, lists_combination_parts)