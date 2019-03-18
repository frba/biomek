from ..misc import calc, file, parser
from ..container import plate, machine
import sys


def get_plate_with_empty_well(destination_plates):
    for i in range(0, len(destination_plates)):
        if destination_plates[i].get_empty_well_coord() is not None:
            return i


def get_localization_vol(part_name, list_source_wells):
    for part in list_source_wells:
        name, num_part, vol_part, source_plate, source_well = part
        if name == part_name:
            return source_plate, source_well, num_part, vol_part


def populate_destination_plates(plates_out, list_destination_plate, list_source_wells, robot, mix_parameters):
    part_fmol, bb_fmol, total_vol, per_buffer, per_rest_enz, per_lig_enz = mix_parameters
    list_source_mix = []
    result = []
    for set in list_destination_plate:
        i = get_plate_with_empty_well(plates_out)
        row, col = plates_out[i].get_empty_well_coord()
        total_parts_vol = 0
        for part in set:
            source_plate, source_well, num_part, vol_part = get_localization_vol(part, list_source_wells)
            total_parts_vol += vol_part
            #TODO: Add part fmol concentration in list

            # type_part_fmol = calc.fmol_by_parttype(num_part, bb_fmol, part_fmol)
            plates_out[i].wells[row][col].samples.append(
                plate.Sample(part, None, None, None, vol_part))

            result.append(['', '', part, source_plate, source_well, vol_part, '', '', '', plates_out[i].name,
                          plates_out[i].wells[row][col].name, plates_out[i].id])

        '''Calculate buffer and enzimes'''
        buffer_vol = (per_buffer * total_vol)/100
        rest_enz_vol = (per_rest_enz * total_vol)/100
        lig_enz_vol = (per_lig_enz * total_vol)/100
        total_buffer = buffer_vol + rest_enz_vol + lig_enz_vol
        # print(buffer_vol, rest_enz_vol, lig_enz_vol)

        '''Total water volume in well'''
        vol_water = total_vol - (total_buffer + total_parts_vol)
        # print(vol_water)

        list_source_mix.append([buffer_vol, rest_enz_vol, lig_enz_vol, vol_water])

        result.append(['', '', 'buffer', '', '', buffer_vol, '', '', '', plates_out[i].name, plates_out[i].wells[row][col].name, plates_out[i].id])
        result.append(['', '', 'Rest_Enz', '', '', rest_enz_vol, '', '', '', plates_out[i].name, plates_out[i].wells[row][col].name, plates_out[i].id])
        result.append(['', '', 'Ligase_Enz', '', '', lig_enz_vol, '', '', '', plates_out[i].name, plates_out[i].wells[row][col].name, plates_out[i].id])
        result.append(['', '', 'water', '', '', vol_water, '', '', '', plates_out[i].name, plates_out[i].wells[row][col].name, plates_out[i].id])

        ''' Add receipts in destination plate '''
        plates_out[i].wells[row][col].samples.append(plate.Sample('buffer', None, None, None, buffer_vol))
        plates_out[i].wells[row][col].samples.append(plate.Sample('restriction_enzime', None, None, None, rest_enz_vol))
        plates_out[i].wells[row][col].samples.append(plate.Sample('ligase_enzime', None, None, None, lig_enz_vol))
        plates_out[i].wells[row][col].samples.append(plate.Sample('water', None, None, None, vol_water))
    return plates_out, list_source_mix, result


def create_plate(num_wells, name):
    rows, cols = calc.rows_columns(int(num_wells))
    new_plate = plate.Plate(rows, cols, name)
    return new_plate


def create_destination_plates(list_destination_plate, out_num_well):
    plates_out = []
    num_receipts = len(list_destination_plate)
    num_plates = calc.num_destination_plates(num_receipts, out_num_well)

    for i in range(0, num_plates):
        plates_out.append(create_plate(out_num_well, 'Destination_' + str(i+1)))
    return plates_out


def create_entry_list_for_destination_plate(lists_parts, list_part_low_vol):
    """Calculate the number of destination plates"""
    list_destination_plate = []
    for set in lists_parts:
        low_vol = False
        for part in set:
            for low_vol_part in list_part_low_vol:
                name_lvp = low_vol_part[0]
                vol_lvp = low_vol_part[1]
                if part == name_lvp:
                    print('For Contructor: ' + str(set) + ' There is not enough volume for sample: ' + str(
                        part) + '. Required : ' + str(vol_lvp))
                    low_vol = True

        if low_vol is False:
            '''Add constructor to list for destination plate'''
            list_destination_plate.append(set)
    return list_destination_plate


def get_part_info(found_list, name):
    for part in found_list:
        part_name, part_type, part_length, part_conc, part_vol, source_plate, source_well, num_well = part
        if part_name == name:
            part_conc = float(part_conc)
            return part_type, part_length, part_conc
    return None


def calc_part_volumes_in_plate(count_unique_list, plates_in, mix_parameters, dispenser_parameters):
    part_fmol, bb_fmol, total_vol, buffer, rest_enz, lig_enz = mix_parameters
    machine, min_vol, res_vol, dead_vol = dispenser_parameters
    total_vol_parts = []

    for pair in count_unique_list:
        part_name = pair[0]
        count = pair[1]  # Number of times it appears in experiment

        for plate_in in plates_in:
            list_wells = plate_in.get_samples_well(part_name)
            while list_wells:
                try:
                    wellD = next(list_wells)

                    for sample in wellD.samples:
                        if sample.name == part_name:
                            # print(sample.name, sample.type, sample.length, sample.concentration, sample.volume)
                            '''fmol -> ng  of the part to give 80 or 40 fmol of that part'''
                            fmol, concent_fmol = calc.fmol(sample.type, sample.length, bb_fmol, part_fmol)

                            '''Volume of part to get the selected fmol in ng '''
                            vol_part_add = float(fmol) / float(sample.concentration)
                            # print(vol_part_add)

                            '''Rounding the part volume according to machine resolution'''
                            vol_part_add = calc.round_at(vol_part_add, res_vol)
                            # print(vol_part_add)

                            '''Minimal dispense volume'''
                            vol_part_add = max(vol_part_add, min_vol)
                            # print(vol_part_add)

                            total_vol_parts.append([sample.name, sample.type, sample.length, sample.concentration, sample.volume, count, vol_part_add, plate_in.name, wellD.name])

                except StopIteration:
                    break
    return total_vol_parts


def calc_part_volumes(count_unique_list, found_list, mix_parameters, dispenser_parameters):
    part_fmol, bb_fmol, total_vol, buffer, rest_enz, lig_enz = mix_parameters
    machine, min_vol, res_vol, dead_vol = dispenser_parameters

    total_vol_parts = []

    for pair in count_unique_list:
        part_name = pair[0]
        count = pair[1]  # Number of times it appears in experiment
        part_type, part_length, part_conc = get_part_info(found_list, part_name)
        # part_in_plate = get_part_in_plate(plates_in, part_name)

        '''fmol -> ng  of the part to give 80 or 40 fmol of that part'''
        fmol, concent_fmol = calc.fmol(part_type, part_length, bb_fmol, part_fmol)
        # print(part_name, part_conc, part_type, part_length)

        '''Volume of part to get the selected fmol in ng '''
        vol_part_add = float(fmol)/float(part_conc)
        # print(vol_part_add)

        '''Rounding the part volume according to machine resolution'''
        vol_part_add = calc.round_at(vol_part_add, res_vol)
        # print(vol_part_add)

        '''Minimal dispense volume'''
        vol_part_add = max(vol_part_add, min_vol)
        # print(vol_part_add)

        total_vol_parts.append([part_name, count, vol_part_add])

    return total_vol_parts


def calc_mixer_volumes(list_destination_plate, mix_parameters):
    vol_for_mixer = []
    part_fmol, bb_fmol, total_vol, per_buffer, per_rest_enz, per_lig_enz = mix_parameters
    '''Number of wells in destination plate'''
    for set in list_destination_plate:
        '''Calculate buffer and enzimes'''
        buffer_vol = (per_buffer * total_vol) / 100
        rest_enz_vol = (per_rest_enz * total_vol) / 100
        lig_enz_vol = (per_lig_enz * total_vol) / 100
        total_vol_buffer = buffer_vol + rest_enz_vol + lig_enz_vol
        vol_for_mixer.append([buffer_vol, rest_enz_vol, lig_enz_vol, total_vol_buffer])
    return vol_for_mixer


def add_on_list(lista, item):
    for i in range(0, len(lista)):
        if lista[i] == item:
            return False
    return True


def get_list_no_repetition(lists_parts):
    unique_list = []
    for lista in lists_parts:
        for part in lista:
                if add_on_list(unique_list, part):
                    unique_list.append(part)
    return unique_list


def get_count_unique_list(unique_list, lists_parts):
    count_unique_list = []
    for part in unique_list:
        count = calc.num_times_part(part, lists_parts)
        count_unique_list.append([part, count])
    return count_unique_list


def verify_samples_volume(vol_for_part, found_list, robot):
    '''Volume needed of parts for the experiment'''
    list_source_wells = []
    list_part_low_vol = []
    for part in vol_for_part:
        sample_name, sample_type, sample_length, sample_concentration, sample_volume, count, vol_part_add, plate_in_name, wellD_name = part
        total_vol_part = count * vol_part_add

        '''Volume available of parts in database'''
        found = False
        for part_f in found_list:
            name_part_f, part_type, part_length, part_conc, part_vol, source_plate, source_well, num_well = part_f
            if name_part_f == sample_name:
                """ Verifies if the total amount of volume minus the dead volume is enough"""
                available_vol = float(part_vol) - robot.dead_vol
                # print(name_part_f, available_vol, total_vol_part, float(part_vol), robot.dead_vol)
                if available_vol > total_vol_part:
                    found = True
                    list_source_wells.append(part)
                    break
        if found is False:
            # print('Not enough volume for sample: ' + str(sample_name))
            list_part_low_vol.append([sample_name, total_vol_part])
    return list_source_wells, list_part_low_vol


def find_samples_database(unique_list, database, db_reader):
    ''' Verify parts in database '''
    found_list = []
    missing_list = []
    for part in unique_list:
        found = False
        database.seek(0)
        for line in db_reader:
            part_indb, part_type, part_length, part_conc, part_vol, source_plate, source_well, num_well = line
            if part == part_indb and float(part_vol) > 0:
                found = True
                # print(part_indb, source_plate, source_well)
                found_list.append(line)
        if found is False:
            # print(part + ' is missing in database.')
            missing_list.append(part)
    return found_list, missing_list


def get_sets_in_filepath(reader):
    lists_parts = []
    ''' For each line in file get the list'''
    for line in reader:
        set = []
        parts = line.strip("\n").split(',')
        '''List of parts'''
        for part in parts:
            set.append(part)
        ''' Create the single list of parts'''
        lists_parts.append(list(set))
    return lists_parts


def create_and_populate_sources_plate(db_reader, database):
    database.seek(0)
    plates_in = parser.create_source_plates_from_csv(db_reader)
    database.seek(0)
    plates_in = parser.csv_to_source_plates(db_reader, plates_in)
    return plates_in


def create_moclo(filepath, database, dispenser_parameters, mix_parameters, out_num_well, pattern):

    filein = file.verify(filepath)
    database = file.verify(database)
    db_reader = file.create_reader_CSV(database)
    # reader = file.create_reader_CSV(filein)

    name_machine, min_vol, res_vol, dead_vol = dispenser_parameters
    robot = machine.Machine(name_machine, min_vol, res_vol, dead_vol)

    """Create combinations"""
    lists_parts = get_sets_in_filepath(filein)
    # print(lists_parts)

    ''' Get unique samples - no repetition'''
    unique_list = get_list_no_repetition(lists_parts)
    # print(unique_list)

    ''' Verify how many times it appears'''
    count_unique_list = get_count_unique_list(unique_list, lists_parts)
    # print(count_unique_list)

    """Verify the parts on database"""
    found_list, missing_list = find_samples_database(unique_list, database, db_reader)
    # print(found_list)

    if len(missing_list) > 0:
        print('Alert for the missing parts: ' + str(missing_list))
        sys.exit(0)

    else:
        """Create and Populate Source Plates"""
        plates_in = create_and_populate_sources_plate(db_reader, database)

        """Calculate the part volumes"""
        vol_for_part = calc_part_volumes_in_plate(count_unique_list, plates_in, mix_parameters, dispenser_parameters)

        """Verify parts volume in source plate"""
        # list_source_wells, list_part_low_vol = verify_samples_volume(vol_for_part, found_list, robot)
        list_source_wells, list_part_low_vol = verify_samples_volume(vol_for_part, found_list, robot)

        """Create entry list for destination plates"""
        list_destination_plate = create_entry_list_for_destination_plate(lists_parts, list_part_low_vol)
        # print(list_destination_plate)

        """Calculate the master mix volumes"""
        vol_for_mixer = calc_mixer_volumes(list_source_wells, mix_parameters)
        # print(vol_for_mixer)


        if len(list_destination_plate) > 0:
            """Create a destination plates"""
            plates_out = create_destination_plates(list_destination_plate, out_num_well)

            """Populate plate"""
            plates_out_filled, result = populate_destination_plates(plates_out, list_destination_plate, vol_for_mixer, robot, mix_parameters)
            # file.write_plate_by_col(plates_out_filled)

            # for i in range(0, len(result)):
            #     print(result[i])









        # else:
        #     print('Not available samples')
        #     sys.exit()