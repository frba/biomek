"""
# Concordia Genome Foundry
# Script to create CSV files to be used on BioMek
# author: Flavia Araujo
"""


# imported packages
from biomek.misc import selection
from biomek.db import db
import sys


def main():

    if len(sys.argv) == 1:
        # database = db.connect_db()
        choose = selection.autoplay()
        selection.function(choose)



        # cursor = db.get_sample(database, 'pYTK001')
        # cursor = db.get_sample_in_plate(database, 'pYTK001')
        # cursor = db.get_sample_info(database, 'pYTK001')
        # cursor = db.get_samples_from_project(database, '1')

    else:
        print("Insert the expected number of arguments")
        print("#Usage: python exec.py\n")


if __name__ == '__main__':
    try:
        main()

    except KeyboardInterrupt:
        print("\nKeybord Interrupted")
