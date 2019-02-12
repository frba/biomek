"""
#Library to deal with input and output files
"""


import os
import sys

# DATABASE = "../input/database.csv"

# SEP = '\t'
SEP = ','


class colours:
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    BLUE = '\033[94m'
    RED = '\033[93m'


def create(filename, mode):
    """Create a new file"""
    newfile = open(filename, mode)
    return newfile


def get_extension(path):
    """Return the file extension"""
    filename, file_extension = os.path.splitext(path)
    return file_extension


def get_filename(path):
    """Return filename without extension"""
    filename, file_extension = os.path.splitext(path)
    return filename


def verify(path):
    """Verify if the file exist"""
    try:
        filein = open(path, 'r')
        return filein

    except IOError:
        print('File could not be opened.')
        sys.exit(0)


def get_header(filein):
    """Get the first line in file (header)"""
    temp = filein.readline()
    header = temp.split(SEP)
    return header


def set_header(fileout):
    if SEP == ',':
        fileout.write('Source Plate Name,Source Well,Destination Plate Name,Destination Well,Volume\n')
    else:
        fileout.write('Source Plate Name\tSource Well\tDestination Plate Name\tDestination Well\tVolume\n')


def write_result(fileout, result):

    for i in range(0, len(result)):
        part, source_plate_name, source_well_name, water_plate_name, water_plate_well, dest_plate, \
        dest_well, fmol, vol_sample, vol_water, error = result[i]

        if len(error) > 0:
            print(colours.RED + str(error))
        else:
            fileout.write(str(source_plate_name) + '\t' + str(source_well_name) + '\t' + str(dest_plate) + '\t' + str(
                dest_well) + '\t' + str(vol_sample) + str(error) + '\n')

            fileout.write(str(water_plate_name) + '\t' + str(water_plate_well) + '\t' + str(dest_plate) + '\t' + str(
                dest_well) + '\t' + str(vol_water) + '\n')

        # print(str(source_plate_name) + '\t' + str(source_well_name) + '\t' + str(dest_plate) + '\t' + str(
        #     dest_well) + '\t' + str(vol_sample) + str(error))
        # print(str(water_plate_name) + '\t' + str(water_plate_well) + '\t' + str(dest_plate) + '\t' + str(
        #     dest_well) + '\t' + str(vol_water))






