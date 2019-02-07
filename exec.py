# Concordia Genome Foundry
# Script to create CSV files to be used in  a normalization cvs file to be used on BioMek
# author: Flavia Araujo

# imported packages
from container import plate
from misc import file, calc, selection
import sys


class Experiment:
    def __init__(self, plate_source, well_source, plate_dest, well_dest, volume):
        self.plate_in = plate_source
        self.plate_out = plate_dest
        self.volume = volume


def main():

    if len(sys.argv) == 1:
        choose = selection.autoplay()
        selection.function(choose)

    else:
        print("Insert the expected number of arguments")
        print("#Usage: python exec.py\n")


if __name__ == '__main__':
    try:
        main()

    except KeyboardInterrupt:
        print("\nKeybord Interrupted")