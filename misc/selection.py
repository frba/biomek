"""
# File to manager the selection among the function for BIOMEK outputs
"""

import sys
import pyfiglet
from misc import file
from termcolor import cprint
from function import normalization_biomek as nb
from function import template_biomek as tb


def data_normalization():
    filepath = input('Inform the filepath (input/database.csv): ')
    in_num_well = input('Inform the number of wells in source plate: ')
    out_num_well = input('Inform the number of wells in destination plate: ')
    nb.create_Biomek_dilution_output(filepath, in_num_well, out_num_well)


def template():
    num_source_plates = input('Inform the number of source plates: ')
    num_pattern = input('Inform the pattern [1, 2, 4, 6]: ')
    pattern = input('Pattern by row -> ' + file.colours.RED + '0 ' + file.colours.ENDC
                    + 'Pattern by column -> ' + file.colours.RED + '1' + file.colours.ENDC + ': ')
    tb.verify_biomek_constraints(num_source_plates, num_pattern, pattern)


def function(choose):
    if choose == '0':
        '''Create CSV templates'''
        template()
    elif choose == '1':
        '''Create Normalization CSV File'''
        data_normalization()
    else:
        print(file.colours.RED+'Invalid option'+file.colours.ENDC)
        sys.exit()


def autoplay():

    cprint(pyfiglet.figlet_format("Biomek Script"), "blue")
    print('Please choose one option:\n')
    choose = input(file.colours.RED+'0'+file.colours.ENDC+' -> Create CSV files templates for Biomek\n'
                   +file.colours.RED+'1'+file.colours.ENDC+' -> Create a Normalization CSV file to be used on BioMek\n'
                   +'Choose > ')

    return choose