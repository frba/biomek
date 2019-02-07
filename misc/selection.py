#Autoselect function

import sys
import pyfiglet
from misc import file
from termcolor import cprint
from function import normalization_biomek as nb


def data_normalization():
    filepath = input('Inform the filepath (input/database.csv): ')
    in_num_well = input('Inform the number of wells in source plate: ')
    out_num_well = input('Inform the number of wells in destination plate: ')

    nb.create_Biomek_dilution_output(filepath, in_num_well, out_num_well)


def template():
    # TODO: template function for biomek
    print('Not implemented...')
    sys.exit()
    # num_source_plate = input('Inform the number of source plates: ')


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