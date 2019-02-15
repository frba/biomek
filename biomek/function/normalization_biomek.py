"""
# File to calculate the normalization of samples and create a csv file to be used on BIOMEK
"""

from ..container import plate
from ..misc import file, calc
from biomek.function import biomek
import os


def verify_sample_volume(volume_needed, volume_available):
    """
    Returns the minimal sample volume needed if available
    :param volume_needed:
    :param volume_available:
    :return:
    """
    ''' Verifies if the volume needed is available'''
    if volume_needed > volume_available:
        return None
    else:
        # biomek cant deal with liquids below 4ul
        ''' Verifies if the volume needed is bigger then Biomek minimal volume'''
        if volume_needed < biomek.MIN_VOL:
            #Set the needed volume to minimal
            if volume_available >= biomek.MIN_VOL:
                volume_needed = biomek.MIN_VOL
                return volume_needed
            else:
                # print('Not enough sample')
                return None
        #volume_needed less than available and greater or equal to 4ul
        else:
            return volume_needed


def calc_normalization_from_plate(sample, plate_in, plate_water, plate_out, i, j):
    """result = [part, plate_name, well_name, fmol, vol_sample, vol_water, message]"""
    fmol = calc.fmol(sample.get_length(), sample.get_concentration())
    dilut_factor = calc.dilution_factor(fmol, sample.get_concentration())

    '''Verify the total sample volume available'''
    sample_vol_available = sample.get_volume()

    '''Calculate the sample volume to get at least 35ul in total volume'''
    sample_vol_needed = calc.sample_volume(dilut_factor, plate.Well.MIN_VOL)
    sample_vol_verified = verify_sample_volume(sample_vol_needed, sample_vol_available)

    if sample_vol_verified is None:
        message = 'Sample %s: volume needed %s / available volume %s.' % \
                  (sample.name, sample_vol_needed, sample_vol_available)
        r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                  plate_water.wells[i][j].name, '', '', fmol, '', '', message]
        return r_norm

    else:
        ''' Calculate total volume'''
        total_volume = calc.total_volume(sample_vol_verified, dilut_factor)
        if total_volume > plate.Well.MAX_VOL:
            message = 'The total volume %d exceeds the well limit.' % total_volume
            r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                      plate_water.wells[i][j].name, '', '', fmol, '', '', message]
            return r_norm

        else:
            '''Calculate water volume'''
            water_vol_needed = calc.water_volume(total_volume, sample_vol_needed)

            if plate.Sample.MIN_VOL > water_vol_needed > 0:
                ''' the water volume is too low'''
                message = 'The total volume of water %dul to dilute the sample %s is too low.' % (water_vol_needed, sample.name)
                r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                          plate_water.wells[i][j].name, '', '', fmol, '', '', message]
                return r_norm
            else:
                if water_vol_needed < 0:
                    ''' the sample concentration is too low'''
                    message = 'The sample %s concentration %sng/ul is too low to be used.' \
                              % (sample.name, sample.get_concentration())
                    r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                              plate_water.wells[i][j].name, '', '', fmol, '', '', message]
                    return r_norm
                else:
                    r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                              plate_water.wells[i][j].name, '', '', fmol, sample_vol_verified, water_vol_needed, '']
                    return r_norm


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


def populate_plate_sample(plate_in, filein):
    """
    Returns a populate Plate from file
    :param plate_in: Plate Source
    :param filein: csv file
    :return: Plate
    """
    file.get_header(filein)
    for line in filein:
        # temp = line.split(file.SEP)
        samp_name, samp_len, samp_conc, volume, plate_name, plate_well = line
        if plate_in.get_name() is None:
            plate_in.set_name(plate_name)
        else:
            """create a new plate"""
        row, col = calc.wellname_to_coordinates(plate_well)
        plate_in.wells[row][col].samples.append(plate.Sample(samp_name, samp_len, samp_conc, int(volume)))
    return plate_in


def populate_plate_water(plate_water):
    """
    Returns a populate water Plate
    :param plate_water: Plate
    :return: Plate
    """
    plate_water.set_name('Source water')
    for i in range(0, plate_water.num_rows):
        for j in range(0, plate_water.num_cols):
            plate_water.wells[i][j].samples.append(plate.Sample('water', None, None, plate.Well.MAX_VOL))
    return plate_water


def populate_plate_destination(plate_out, norm_result):
    """
    Returns filled destination plate and result
    :param plate_out: object Plate: destination plate for normalized samples
    :param norm_result: a list with normalized results
    :return: Plate, list
    """
    for k in range(0, len(norm_result)):
        sample_name, plate_in_name, plate_in_wells, plate_water_name, \
        plate_water_wells, dest_plate, dest_well, fmol, sample_vol_verified, water_vol_needed, message = \
        norm_result[k]
        if message == '':
            if plate_out.get_empty_well_coord() is not None:
                i, j = plate_out.get_empty_well_coord()
                plate_out.wells[i][j].samples.append(plate.Sample(sample_name, None, None, sample_vol_verified))
                plate_out.wells[i][j].samples.append(plate.Sample('water', None, None, water_vol_needed))
                norm_result[k][5] = plate_out.name
                norm_result[k][6] = plate_out.wells[i][j].name
    return plate_out, norm_result


def create_biomek_dilution_output(path, in_well, out_well):
    """
    Creates a CSV file to be used in Biomek with the volume of sample and water to normalize the sample
    :param path: input file with parts information
    :param in_well: number of well of input plate
    :param out_well: number of well of output plate
    """
    filein = file.verify(path)
    outfile = file.create('output/dilution_' + str(os.path.basename(path)), 'w')
    out_csv = file.createCSV(outfile)
    file.set_header(out_csv)

    """Create a platemap"""
    plate_in = create_plate(in_well, 'Source')
    plate_water = create_plate(in_well, 'Water')
    plate_out = create_plate(out_well, 'Destination')

    """Populate the plate"""
    populate_plate_sample(plate_in, filein)
    populate_plate_water(plate_water)
    norm_result = []
    for i in range(0, plate_in.num_rows):
        for j in range(0, plate_in.num_cols):
            """Calculate the normalization parameters"""
            for sample in plate_in.wells[i][j].samples:
                norm_result.append(calc_normalization_from_plate(sample, plate_in, plate_water, plate_out, i, j))

    plate_out_filled, new_result = populate_plate_destination(plate_out, norm_result)
    """Write the result in a CSV file"""
    file.write_result(out_csv, new_result)
    print(file.colours.BOLD + 'Output File: ' + outfile.name + file.colours.BOLD)
