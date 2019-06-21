"""
# Concordia Genome Foundry
# Script to create CSV files to be used on BioMek
# author: Flavia Araujo
"""

# imported packages
from biomek.misc import selection, genbank
import sys
sys.path.append("..")


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
