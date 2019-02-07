# Library to create a class for different plates

from misc import calc
from container import well


class Plate:

    """
    Parameters
    ----------
    name
      The well's name, for instance "A1"

    rows
      The well's row (a number, starting from 0)

    cols
      The well's column (a number, starting from 0)

    num_wells
      The quantity of wells in the plate (96, 384, 1536)

    wells
      A matrix of wells according the number of cols and rows
    """

    def __init__(self, rows, cols, name):
        self.name = name
        self.num_rows = rows
        self.num_cols = cols
        self.num_wells = rows * cols
        self.wells = [[Well(calc.coordinates_to_wellname(coords=[j, i]), i, j, [None, None, None])
                       for i in range(0, self.num_cols)] for j in range(0, self.num_rows)]

    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name


class Well:

    """
    Parameters
    ----------
    name
      The well's name, for instance "A1"

    row
      The well's row (a number, starting from 0)

    column
      The well's column (a number, starting from 0)

    contents
      A vector of well_contents like Name, Lenght, Concentration.
    """

    def __init__(self, name, row, column):
        self.name = name
        self.row = row
        self.column = column
        self.samples = [Sample(None, None, None, None)]

    def set_name(self, name):
        self.name = name

    def set_row(self, row):
        self.row = row

    def set_column(self, column):
        self.column = column

    def get_row(self):
        return self.row

    def get_column(self):
        return self.column


class Sample:

    """
    Parameters
    ----------
    contents
      Sample: Name, Lenght, Concentration.
    """

    def __init__(self, sample_name, length, concentration, volume):
        self.sample_name = sample_name
        self.length = length
        self.concentration = concentration
        self.volume = volume

    def set_sample_name(self, sample_name):
        self.sample_name = sample_name

    def set_length(self, length):
        self.length = length

    def set_concentration(self,concentration):
        self.concentration = concentration

    def set_volume(self, volume):
        self.volume = volume

    def get_sample_name(self):
        return self.sample_name

    def get_length(self):
        return self.length

    def get_concentration(self):
        return self.concentration

    def get_volume(self):
        return self.volume

    def set_data(self, data):
        self.sample_name = data[0]
        self.length = data[1]
        self.concentration = data[2]
        self.volume = data[3]