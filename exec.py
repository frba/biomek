# Concordia Genome Foundry
# Script to output a normalization cvs file to be used on BioMek
# author: Flavia Araujo
# use: python exec.py <input file path> <input num wells> <output num wells>

# imported packages
from container import plate
from misc import file, calc
import sys


class Experiment:
    def __init__(self, plate_source, well_source, plate_dest, well_dest, volume):
        self.plate_in = plate_source
        self.plate_out = plate_dest
        self.volume = volume


def verify_sample_volume(volume_needed, volume_available):
    if volume_needed > volume_available:
        # print('Not enough sample')
        return None
    else:
        # Biomek cant deal with liquids below 4ul
        if volume_needed < plate.Sample.MIN_VOL:
            #Set the needed volume to minimal
            if volume_available >= plate.Sample.MIN_VOL:
                volume_needed = plate.Sample.MIN_VOL
                return volume_needed
            else:
                # print('Not enough sample')
                return None
        #volume_needed less than available and greater or equal to 4ul
        else:
            return volume_needed


def calc_normalization_from_plate(sample, plate_in, plate_water, i, j):
    """result = [part, plate_name, well_name, fmol, vol_sample, vol_water, message]"""
    fmol = calc.fmol(sample.get_length(), sample.get_concentration())
    dilut_factor = calc.dilution_factor(fmol, sample.get_concentration())

    '''Verify the total sample volume available'''
    sample_vol_available = sample.get_volume()

    '''Calculate the sample volume to get at least 35ul in total volume'''
    sample_vol_needed = calc.sample_volume(dilut_factor, plate.Well.MIN_VOL)
    sample_vol_verified = verify_sample_volume(sample_vol_needed, sample_vol_available)

    if sample_vol_verified is None:
        message = 'Sample volume needed %s / available volume %s.' % \
                  (sample_vol_needed, sample_vol_available)
        r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                  plate_water.wells[i][j].name, fmol, '', '', message]
        return r_norm

    else:
        ''' Calculate total volume'''
        total_volume = calc.total_volume(sample_vol_verified, dilut_factor)
        if total_volume > plate.Well.MAX_VOL:
            message = 'The total volume %d exceeds the well limit.' % total_volume
            r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                      plate_water.wells[i][j].name, fmol, '', '', message]
            return r_norm

        else:
            '''Calculate water volume'''
            water_vol_needed = calc.water_volume(total_volume, sample_vol_needed)

            if plate.Sample.MIN_VOL > water_vol_needed:
                ''' the water volume is too low'''
                message = 'The total volume of water %d to dilute the sample is too low.' % water_vol_needed
                r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                          plate_water.wells[i][j].name, fmol, '', '', message]
                return r_norm
            else:
                r_norm = [sample.name, plate_in.name, plate_in.wells[i][j].name, plate_water.name,
                          plate_water.wells[i][j].name, fmol, sample_vol_needed, water_vol_needed, '']
                return r_norm


def create_water_plate(num_wells):
    rows, cols = calc.rows_columns(int(num_wells))
    new_plate = plate.Plate(rows, cols, 'water_plate')
    return new_plate


def create_plate(num_wells):
    rows, cols = calc.rows_columns(int(num_wells))
    new_plate = plate.Plate(rows, cols, None)
    return new_plate


#TODO: Create a better input file structure
def populate_plate_sample(plate_in, filein):
    file.get_header(filein)
    for line in filein:
        temp = line.split(file.SEP)
        samp_name, samp_len, samp_conc, samp_fmol, samp_uL, vol50, plate_name, plate_well = temp
        if plate_in.get_name() is None:
            plate_in.set_name(plate_name)
        else:
            """create a new plate"""
        row, col = calc.wellname_to_coordinates(plate_well)
        plate_in.wells[row][col].samples.append(plate.Sample(samp_name, samp_len, samp_conc, 5))

    return plate_in


def populate_plate_water(plate_water):
    plate_water.set_name('Source water')
    for i in range(0, plate_water.num_rows):
        for j in range(0, plate_water.num_cols):
            plate_water.wells[i][j].name = 'A1'
            plate_water.wells[i][j].samples.append(plate.Sample('water', None, None, plate.Well.MAX_VOL))
    return plate_water


def create_Biomek_dilution_output(path, in_well, out_well):
    filein = file.verify(path)
    outfile = file.create('output/dilution_biomek.csv', 'w')
    file.set_header(outfile)
    if file.get_extension(path) == '.csv':
        # TODO: verify the separator

        """Create a platemap"""
        plate_in = create_plate(in_well)
        plate_water = create_plate(in_well)
        plate_out = create_plate(out_well)

        """Populate the plate"""
        populate_plate_sample(plate_in, filein)
        populate_plate_water(plate_water)

        for i in range(0, plate_in.num_rows):
            for j in range(0, plate_in.num_cols):
                """Calculate the normalization parameters"""
                for sample in plate_in.wells[i][j].samples:
                    norm_result = calc_normalization_from_plate(sample, plate_in, plate_water, i, j)
                    file.write_result(outfile, norm_result)

        filein.close()
        outfile.close()
    else:
        # TODO: support to xlsx files
        print('File could not be opened.')


def main():
    if len(sys.argv) > 3:
        filepath = sys.argv[1]
        input_num_wells = sys.argv[2]
        output_num_wells = sys.argv[3]
        create_Biomek_dilution_output(filepath, input_num_wells, output_num_wells)

    else:
        print("Insert the expected number of arguments")
        print("#Usage: python exec.py input_file #numwell_input #numwell_output\n")


if __name__ == '__main__':
    try:
        main()

    except KeyboardInterrupt:
        print("Interrupted")