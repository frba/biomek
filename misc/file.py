'''Library to deal with input and output files'''

import os
import sys
DATABASE = "../input/database.csv"

SEP = ','


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
    part, source_plate_name, source_well_name, water_plate_name, water_plate_well, fmol, vol_sample, vol_water, error = result
    print(result)
    if SEP == ',':
        fileout.write(str(source_plate_name)+','+str(source_well_name)+','+str(source_plate_name)+','+str(source_well_name)+','+str(vol_sample)+str(error)+'\n')
        fileout.write(str(water_plate_name)+','+str(water_plate_well)+','+str(source_plate_name)+','+str(source_well_name)+','+str(vol_water)+'\n')
    else:
        fileout.write('Source Plate Name\tSource Well\tDestination Plate Name\tDestination Well\tVolume\tError\n')







