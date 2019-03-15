"""
# File to manager the selection among the function for BIOMEK outputs
"""

import sys
from biomek.misc import file
from ..function import normalization as nb
from ..function import spotting as tb
from ..function import combinatorial as ct
from ..function import moclo as mc


def moclo():
    # filepath = input('Inform the filepath (biomek/input/combinat_part.csv): ')
    # database = input('Inform the filepath (biomek/input/database.csv): ')
    # filepath = 'biomek/input/combination_parts_egf.txt'
    filepath = 'biomek/input/combination_parts.txt'
    # database = 'biomek/input/database_egf.csv'
    database = 'biomek/input/database.csv'
    # dispenser_parameters = input('dispenser_parameters: ')
    machine = 0
    """Minimal volume in nl obtain from the machine from source plate"""
    min_vol = 2.5e-9
    """Minimal volume in nl dropped by the machine"""
    res_vol = 2.5e-9
    """Volume in ul cant be reached by the machine"""
    dead_vol = 3
    dispenser_parameters = [machine, min_vol, res_vol, dead_vol]

    # mix_parameters = input('mixer parameters: ')
    part_fmol = 40
    bb_fmol = 80
    total_vol = 10
    buffer = 10
    rest_enz = 10
    lig_enz = 10
    mix_parameters = [part_fmol, bb_fmol, total_vol, buffer, rest_enz, lig_enz]
    out_num_well = '96'
    pattern = '0'

    # out_num_well = input('Inform the number of wells in destination plate: ')
    # plate_type = input('Inform the plate type: ')
    # pattern = input('Pattern by row -> ' + file.colours.RED + '0 ' + file.colours.ENDC +
    # 'Pattern by column -> ' + file.colours.RED + '1 ' + file.colours.ENDC + ': ')

    mc.create_moclo(filepath, database, dispenser_parameters, mix_parameters, int(out_num_well), int(pattern))


def combinatorial_data():
    filepath = input('Inform the filepath (biomek/input/parts.csv): ')
    ct.create_combinations(filepath)


def data_normalization():
    # filepath = input('Inform the filepath (biomek/input/database.csv): ')
    filepath = 'biomek/input/to_be_normalized.csv'
    in_num_well = input('Inform the number of wells in source plate: ')
    out_num_well = input('Inform the number of wells in destination plate: ')
    bb_fmol = input('Inform the fmol for backbone (80 or 40):  ')
    part_fmol = input('Inform the fmol for part (40 or 20):  ')
    nb.create_biomek_dilution_output(filepath, int(in_num_well), int(out_num_well), bb_fmol, part_fmol)


def template():
    num_source_plates = input('Inform the number of source plates: ')
    num_pattern = input('Inform the pattern [1, 2, 3, 4, 6 or 8]: ')
    pattern = input('Pattern by row -> ' + file.colours.RED + '0 ' + file.colours.ENDC
                    + 'Pattern by column -> ' + file.colours.RED + '1 ' + file.colours.ENDC
                    + 'Pattern by Biomek -> ' + file.colours.RED + '2' + file.colours.ENDC
                    + ': ')
    tb.verify_biomek_constraints(int(num_source_plates), int(num_pattern), int(pattern))


def function(choose):
    """
    Function to select the function in the system
    :param choose: int number
    """
    if choose == 0:
        '''Create CSV templates'''
        template()
    elif choose == 1:
        '''Create Normalization CSV File'''
        data_normalization()
    elif choose == 2:
        '''Create Combinatorial CSV File'''
        combinatorial_data()
    elif choose == 3:
        moclo()
    else:
        print(file.colours.RED + 'Invalid option' + file.colours.ENDC)
        sys.exit()


def autoplay():
    """
    Autoplay script that helps to choose the function desired
    :return: number int
    """
    print(file.colours.BLUE + "Biomek Script" + file.colours.ENDC)
    print('Please choose one option:\n')
    choose = input(file.colours.RED + '0' + file.colours.ENDC + ' -> Create a CSV file template\n'
                   + file.colours.RED + '1' + file.colours.ENDC + ' -> Create a Normalization CSV file\n'
                   + file.colours.RED + '2' + file.colours.ENDC + ' -> Create a CSV file with combinatorial\n'
                   + file.colours.RED + '3' + file.colours.ENDC + ' -> Create a CSV file for MOCLO\n'
                   +'Choose > ')

    return int(choose)