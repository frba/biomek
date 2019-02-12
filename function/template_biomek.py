"""
Functions to create CSV files to be used in Biomek

Source Plate Name,Source Well,Destination Plate Name,Destination Well,Volume
PlateS1,A1,PlateD1,A1,4
PlateS1,A1,PlateD1,B1,4
PlateS1,A2,PlateD1,C1,4
PlateS1,A2,PlateD1,D1,4
"""

from misc import calc, file, selection

from container import plate
import sys

MAX_PLATES = 12
VOLUME = 4
BY_ROW = '0'
BY_COL = '1'


def verify_biomek_constraints(num_source_plates, num_pattern, pattern):
    """
    Calls a function to create a output file
    The output file has the source plate and the distribution of the samples according to the num_pattern and pattern
    :param num_source_plates: int number
    :param num_pattern: int number
    :param pattern: 0 or 1
    """
    ver_num_source = selection.verify_entry(int, num_source_plates)
    ver_pattern = selection.verify_entry(int, num_pattern)
    total_destination = ver_num_source * ver_pattern
    total_plates = ver_num_source + total_destination
    if total_plates > MAX_PLATES:
        print('The total plates (%d) exceeds the Biomek limit of %d' % (total_plates, MAX_PLATES))

    else:
        print('The total plates in Biomek is %d' % total_plates)
        print('The total destination plate(s) is %d and total source plate(s) is %d' % (total_destination, ver_num_source))
        create_output_file(ver_num_source, total_destination, pattern)


def generate_random_names(name, init, end):
    """
    Returns a vector with a main name + number
    :param name: string
    :param init: int number
    :param end: int number
    :return: a vector
    """
    names = []
    for i in range(init, end):
        names.append(str(name) + str(i))
    return names


def create_plate(num_wells, name):
    """
    Returns a named plate type 96 or 384 according to num_wells
    :param num_wells: int [96, 384]
    :param name: string
    :return: object from class Plate
    """
    rows, cols = calc.rows_columns(int(num_wells))
    new_plate = plate.Plate(rows, cols, name)
    return new_plate


def write_on_file_by_col(source_plate, destination_plates, num_pattern, outfile):
    """
    Create a .csv file to be used in Biomek
    :param source_plate: object from Plate Class
    :param destination_plates: a vector with of Plate Class that will receive the samples from source plate
    :param num_pattern: number of repetitions samples get from source plates
    :param outfile: A CSV file to be used in Biomek with the choosed pattern
    with 1 source plate and num_pattern = 2, pattern = bycols, the output file will be like:
    Source Plate Name,Source Well,Destination Plate Name,Destination Well,Volume
    PlateS1,A1,PlateD1,A1,4
    PlateS1,A1,PlateD1,B1,4
    PlateS1,A2,PlateD1,C1,4
    PlateS1,A2,PlateD1,D1,4
    """
    source_wells = source_plate.iterR(num_pattern)
    for plateD in destination_plates:
        dest_wells = plateD.iterC(1)
        while source_wells and dest_wells:
            try:
                wellD = next(dest_wells)
                wellS = next(source_wells)
                # print(source_plate.name + ',' + wellS.name + ',' + plateD.name + ',' + wellD.name + ',' + str(VOLUME))
                outfile.write(str(source_plate.name) + ',' + str(wellS.name) + ',' + str(plateD.name) + ',' + str(
                    wellD.name) + ',' + str(VOLUME) + '\n')
            except StopIteration:
                break


def write_on_file_by_row(source_plate, destination_plates, num_pattern, outfile):
    """
    Create a .csv file to be used in Biomek
    :param source_plate: source_plates: object from Plate Class
    :param destination_plates: a vector with of Plate Class that will receive the samples from source plate
    :param num_pattern: number of repetitions samples get from source plates
    :param outfile: A CSV file to be used in Biomek with the choosed pattern
    with 1 source plate and num_pattern = 2, pattern = byrows, the output file will be like:
    Source Plate Name,Source Well,Destination Plate Name,Destination Well,Volume
    PlateS1,A1,PlateD1,A1,4
    PlateS1,A1,PlateD1,A2,4
    PlateS1,A2,PlateD1,A3,4
    PlateS1,A2,PlateD1,A4,4
    """
    source_wells = source_plate.iterR(num_pattern)
    for plateD in destination_plates:
        dest_wells = plateD.iterR(1)
        while dest_wells:
            try:
                wellD = next(dest_wells)
                wellS = next(source_wells)
                # print(source_plate.name + ',' + wellS.name + ',' + plateD.name + ',' + wellD.name + ',' + str(VOLUME))
                outfile.write(str(source_plate.name) + ',' + str(wellS.name) + ',' + str(plateD.name) + ',' + str(
                    wellD.name) + ',' + str(VOLUME) + '\n')
            except StopIteration:
                break


def create_output_file(total_source, total_destination, pattern):
    """
    Create a random output file name, and plates names
    :param total_source: integer number
    :param total_destination: integer number
    :param pattern: integer number 0 -> BY_ROW or 1 -> BY_COL
    """
    num_pattern = int(total_destination/total_source)
    '''Add the header'''
    if pattern == BY_ROW:
        outfile = file.create('output/template_'+str(total_destination)+'x'+str(total_source)+'_byrow.csv', 'w')
        file.set_header(outfile)
        ''' Create the source plates'''
        for i in range(0, total_source):
            plateS_num = i+1
            source_name = 'PlateS' + str(plateS_num)
            source_plate = create_plate(96, source_name)
            destination_names = generate_random_names('PlateD', num_pattern*i+1, num_pattern*i+num_pattern+1)
            destination_plates = []
            for j in range(0, len(destination_names)):
                destination_plates.append(create_plate(96, destination_names[j]))
            '''Call Function to write the CSV by rows'''
            write_on_file_by_row(source_plate, destination_plates, num_pattern, outfile)
        print(file.colours.BOLD + 'Output File: ' + outfile.name + file.colours.BOLD)

    elif pattern == BY_COL:
        outfile = file.create('output/template_' + str(total_source) + 'x' + str(total_destination) + '_bycol.csv', 'w')
        file.set_header(outfile)
        ''' Create the source plates'''
        for i in range(0, total_source):
            plateS_num = i + 1
            source_name = 'PlateS' + str(plateS_num)
            source_plate = create_plate(96, source_name)
            destination_names = generate_random_names('PlateD', num_pattern * i + 1, num_pattern * i + num_pattern + 1)
            destination_plates = []
            for j in range(0, len(destination_names)):
                destination_plates.append(create_plate(96, destination_names[j]))
            '''Call Function to write the CSV by rows'''
            write_on_file_by_col(source_plate, destination_plates, num_pattern, outfile)
        print(file.colours.BOLD + 'Output File: ' + outfile.name + file.colours.BOLD)
    else:
        print('Invalid option')
        sys.exit()