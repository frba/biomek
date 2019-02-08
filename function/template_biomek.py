from misc import calc, file
from container import plate
import sys

MAX_PLATES = 12
VOLUME = 4
BY_ROW = '0'
BY_COL = '1'


def verify_entry(type, num):
    try:
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


def verify_biomek_constraints(num_source_plates, num_pattern, pattern):
    ver_num_source = verify_entry(int, num_source_plates)
    ver_pattern = verify_entry(int, num_pattern)
    total_destination = ver_num_source * ver_pattern
    total_plates = ver_num_source + total_destination

    if total_plates > MAX_PLATES:
        print('The total plates (%d) exceeds the Biomek limit of %d' % (total_plates, MAX_PLATES))

    else:
        print('The total plates in Biomek is %d' % total_plates)
        print('The total destination plate(s) is %d and total source plate(s) is %d' % (total_destination, ver_num_source))
        create_output_file(ver_num_source, total_destination, pattern)


def generate_random_names(name, num_var):
    names = []
    for i in range(1, num_var+1):
        names.append(str(name) + str(i))
    return names


def create_plate(num_wells, name):
    rows, cols = calc.rows_columns(int(num_wells))
    new_plate = plate.Plate(rows, cols, name)
    return new_plate


# def write_on_file_by_row(plateS, plateD, num_pattern, outfile):




def create_output_file(total_source, total_destination, pattern):
    source_plates = []
    destination_plates = []
    num_pattern = int(total_destination/total_source)

    '''Add the header'''
    if pattern == BY_ROW:
        outfile = file.create('output/template_'+str(total_destination)+'x'+str(total_source)+'.csv', 'w')
        file.set_header(outfile)
        ''' Create the source plates'''
        source_names = generate_random_names('PlateS', total_source)
        for i in range(0, len(source_names)):
            source_plates.append(create_plate(96, source_names[i]))

        ''' Create the destination plates'''
        destination_names = generate_random_names('PlateD', total_destination)
        for j in range(0, len(destination_names)):
            destination_plates.append(create_plate(96, destination_names[j]))

        for plateS in source_plates:
            for plateD in destination_plates:
                source_wells = plateS.iterL(num_pattern)
                dest_wells = plateD.iterL(1)

                while dest_wells:
                    try:
                        wellS = next(source_wells)
                        wellD = next(dest_wells)
                        print(plateS.name + ',' + wellS.name + ',' + plateD.name + ',' + wellD.name + ',' + str(VOLUME))
                        outfile.write(str(plateS.name) + ',' + str(wellS.name) + ',' + str(plateD.name) + ',' + str(
                            wellD.name) + ',' + str(VOLUME) + '\n')
                    except StopIteration as e:
                        print (e)
                        return

        print(file.colours.BOLD + 'Output File: ' + outfile.name + file.colours.BOLD)

    elif pattern == BY_COL:
        outfile = file.create('output/template_' + str(total_source) + 'x' + str(total_destination) + '.csv', 'w')
        file.set_header(outfile)
        ''' Create the source plates'''
        source_names = generate_random_names('PlateS', total_source)
        for i in range(0, len(source_names)):
            source_plates.append(create_plate(96, source_names[i]))

        ''' Create the destination plates'''
        destination_names = generate_random_names('PlateD', total_destination)
        for j in range(0, len(destination_names)):
            destination_plates.append(create_plate(96, destination_names[j]))

        for plateS in source_plates:
            for plateD in destination_plates:
                write_on_file_by_col(plateS, plateD, num_pattern, outfile)

        print(file.colours.BOLD + 'Output File: ' + outfile.name + file.colours.BOLD)
    else:
        print('Invalid option')
        sys.exit()

