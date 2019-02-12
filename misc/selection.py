"""
# File to manager the selection among the function for BIOMEK outputs
"""

import sys
import pyfiglet
from misc import file
from termcolor import cprint
from function import normalization_biomek as nb
from function import template_biomek as tb


def verify_entry(type, num):
    try:
        '''Verify if the numbers type is correct'''
        num = type(num)
    except ValueError:
        message = str(num) + ' is not a number'
        print(message)
        sys.exit()
    if num <= 0:
        message = 'the value needs to be greater than ' + str(num)
        print(message)
        sys.exit()
    else:
        return num


def combinatorial_data():
    print('Under development')
    sys.exit()


def data_normalization():
    filepath = input('Inform the filepath (input/database.csv): ')
    in_num_well = input('Inform the number of wells in source plate: ')
    out_num_well = input('Inform the number of wells in destination plate: ')
    nb.create_Biomek_dilution_output(filepath, in_num_well, out_num_well)


def template():
    num_source_plates = input('Inform the number of source plates: ')
    num_pattern = input('Inform the pattern [1 to 11]: ')
    pattern = input('Pattern by row -> ' + file.colours.RED + '0 ' + file.colours.ENDC
                    + 'Pattern by column -> ' + file.colours.RED + '1' + file.colours.ENDC + ': ')
    tb.verify_biomek_constraints(num_source_plates, num_pattern, pattern)


def function(choose):
    """
    Function to select the function in the system
    :param choose: int number
    """
    if choose == '0':
        '''Create CSV templates'''
        template()
    elif choose == '1':
        '''Create Normalization CSV File'''
        data_normalization()
    elif choose == '2':
        '''Create Combinatorial CSV File'''
        combinatorial_data()
    else:
        print(file.colours.RED+'Invalid option'+file.colours.ENDC)
        sys.exit()


def autoplay():
    """
    Autoplay script that helps to choose the function desired
    :return: number int
    """
    cprint(pyfiglet.figlet_format("Biomek Script"), "blue")
    print('Please choose one option:\n')
    choose = input(file.colours.RED+'0'+file.colours.ENDC+' -> Create a CSV file template\n'
                   +file.colours.RED+'1'+file.colours.ENDC+' -> Create a Normalization CSV file\n'
                   +file.colours.RED+'2'+file.colours.ENDC+' -> Create a CSV file with combinatorial\n'
                   +'Choose > ')

    return choose