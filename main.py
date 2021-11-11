import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILE_PATH = 'IGF1-AAV primer combinations.xls'
WELLS_TO_PLOT = ['A4', 'A5', 'B1', 'B2']
well_letters_dictionary = {'0': 'A', '1': 'B', '2': 'C', '3': 'D', '4': 'E', '5': 'F', '6': 'G', '7': 'H'}


# Get Amplification Data and Clean it Up
amplification_sheet = pd.read_excel(FILE_PATH, sheet_name='Amplification Data')
data_cleaned = amplification_sheet.drop(amplification_sheet.index[0:43])
data_cleaned.rename(columns={'Block Type': 'Well', '96-Well 0.2-mL Block': 'Well Letter',
                             'Unnamed: 2': 'Cycle', 'Unnamed: 3': 'Target Name', 'Unnamed: 4': 'Rn',
                             'Unnamed: 5': 'Delta Rn'}, inplace=True)
data_cleaned = data_cleaned.drop('Well', 1)
data_cleaned.reset_index(inplace=True)


# Get Sample Data and Clean it Up
identity_sheet = pd.read_excel(FILE_PATH, sheet_name='Sample Setup')
data_cleaned2 = identity_sheet.drop(identity_sheet.index[0:43])
data_cleaned2.rename(columns={'96-Well 0.2-mL Block': 'Well', 'Unnamed: 2': 'Sample Name',
                              'Unnamed: 6': 'Target Name'}, inplace=True)
data_cleaned2.reset_index(inplace=True)


class Well:
    def __init__(self, well):
        self.well = well
        self.list_of_delta_rn = data_cleaned[data_cleaned["Well Letter"] == well]['Delta Rn']
        self.delta_rn_array = np.array(self.list_of_delta_rn)
        self.cycle_array = np.arange(start=1, stop=41, step=1)


def plot(well_obj):
    plt.scatter(x=well_obj.cycle_array, y=well_obj.list_of_delta_rn, label=str(well_obj.well), s=6)
    plt.plot(well_obj.cycle_array, well_obj.list_of_delta_rn)


fig = plt.figure(figsize=(12, 7))
for row in range(8):
    for col in range(1, 13):
        well_name = f'{well_letters_dictionary[str(row)]}{col}'
        w = Well(well_name)
        if w.well in WELLS_TO_PLOT:
            plot(Well(well_name))

plt.legend()
plt.xlabel('Cycle', fontsize=16)
plt.ylabel('Delta Rn', fontsize=16)
plt.show()


# Find location of well letter in the second excel sheet and then find sample id and target primer. Add to legend