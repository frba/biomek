from ..misc import calc, file
from ..container import plate
import itertools, csv


# def populate_destination_plates(plates_out, lists_combination_parts, database):

    # for line in filein:
    #     samp_name, samp_len, samp_conc, volume, plate_name, plate_well = line
    #     for i in range(0, len(plates_in)):
    #         if plates_in[i].name == plate_name:
    #             row, col = calc.wellname_to_coordinates(plate_well)
    #             plates_in[i].wells[row][col].samples.append(plate.Sample(samp_name, samp_len, samp_conc, int(volume)))
    # return plates_in


def add_on_list(list, item):
    for i in range(0, len(list)):
        if list[i] == item:
            return False
    return True


def get_list_no_repetition(lists_combination_parts):
    unique_list = []
    for combinations in lists_combination_parts:
        for set in combinations:
            for item in set:
                if add_on_list(unique_list, item):
                    unique_list.append(item)

    return unique_list


def get_count_unique_list(unique_list, lists_combination_parts):
    count_unique_list = []
    for item in unique_list:
        count = calc.num_times_combination(item, lists_combination_parts)
        count_unique_list.append([item, count])
    return count_unique_list


def verify_samples_database(lists_combination_parts, file, database):
    ''' Get unique samples - no repetition '''
    unique_list = get_list_no_repetition(lists_combination_parts)
    # print(unique_list)

    ''' Verify how many times it appears '''
    count_unique_list = get_count_unique_list(unique_list, lists_combination_parts)
    # print(count_unique_list)

    vol_unique_list = []
    for item in count_unique_list:
        total_vol_available = 0
        part = item[0]
        count = item[1]

        file.seek(0)
        for line in database:
            part_d, length_d, concent_d, vol_d, plate_d, well_d = line
            if part == part_d:
                total_vol_available += float(vol_d)
                print(part, part_d, vol_d)
        vol_unique_list.append([part, total_vol_available])

    print(vol_unique_list)
    # for line in database:
    #     print(line)


def create_plate(num_wells, name):
    """
    Returns a Plate with number of wells and name
    :param num_wells: int number (96, 384)
    :param name: String of plate name
    :return: Plate
    """
    rows, cols = calc.rows_columns(int(num_wells))
    new_plate = plate.Plate(rows, cols, name)
    return new_plate


def create_destination_plates(lists_combination_parts, out_num_well):
    plates_out = []

    count_num_wells = 0
    for combinations in lists_combination_parts:
        count_num_wells += len(combinations)
    num_dest_plates = calc.num_destination_plates(count_num_wells, out_num_well)

    for i in range(0, num_dest_plates):
        plates_out.append(create_plate(out_num_well, 'Destination_' + str(i+1)))
    return plates_out


def get_sets_in_filepath(reader):
    "TODO: Read csv file correctly"
    lists_combination_parts = []
    count_num_wells = 0
    ''' For each line in file get the list'''

    # for line in reader:
    #     set_parts = []
    #     for col in line:
    #         temp = col.split(';')
    #         for i in range(0, len(temp)):
    #             set_parts.append(temp)
    #     print(set_parts)


    '''List of parts'''
    set_parts1 = list([["TU1-5'", "TU1-3'", "TU1-5' gblock"], ["PSMB gblock"], ["pYTK001", "pYTK002", "pYTK003", "pYTK004"]])
    # set_parts2 = list([['TU1-5', 'TU1-3', 'TU1-5 gblock'], ['PSMB gblock'], ['pYTK001', 'pYTK002', 'pYTK003', 'pYTK004']])

    ''' Create the combinations from the list'''
    lists_combination_parts.append(list(itertools.product(*set_parts1)))
    # lists_combination_parts.append(list(itertools.product(*set_parts2)))

    return lists_combination_parts


def create_combinations(database, filepath, out_num_well):

    filein = file.verify(filepath)
    database = file.verify(database)
    database_csv = file.create_reader_CSV(database)
    reader = file.create_reader_CSV(filein)

    """Create combinations"""
    lists_combination_parts = get_sets_in_filepath(reader)


    """Verify the volume of parts"""
    verify_samples_database(lists_combination_parts, database, database_csv)


    """Create a platemap"""
    plates_out = create_destination_plates(lists_combination_parts, out_num_well)


    """Populate plate"""
    # plates_out_filled = populate_destination_plates(plates_out, lists_combination_parts)