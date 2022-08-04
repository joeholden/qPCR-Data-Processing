import pandas as pd
from mod import Reaction, fold_change, find_dCT
import statistics
import numpy as np

FILE_PATH = 'sGC DOX and no DOX treatment cell lines alpha 1 and 2.xls'
HOUSEKEEPING_GENE = 'beta actin'
CONTROL = 'cal parent'

# Get Setup Data
identity_sheet = pd.read_excel(FILE_PATH, sheet_name='Sample Setup')
data_cleaned = identity_sheet.drop(identity_sheet.index[0:43])
data_cleaned.rename(columns={'96-Well 0.2-mL Block': 'Well', 'Unnamed: 2': 'Sample Name',
                              'Unnamed: 6': 'Target Name'}, inplace=True)
data_cleaned.reset_index(inplace=True)
wells = data_cleaned['Well']
samples = data_cleaned['Sample Name']
targets = data_cleaned['Target Name']

# extracts pertinent information and adds it to rxn list
rxns = []
for row in range(data_cleaned.shape[0]):
    rxn = (wells[row], samples[row], targets[row])
    rxns.append(rxn)

# Gets all reaction information in a list [[sample, target, [well 1, well], ...]
consolidated = []
for i in rxns:
    for j in rxns:
        if i != j and i[1] == j[1] and i[2] == j[2]:
            consolidated.append([i[1], i[2], [i[0], j[0]]])
for i in consolidated:
    for j in consolidated:
        if i != j and i[0]==j[0] and i[1] == j[1] and set(i[2]) == set(j[2]):
            consolidated.remove(j)

# pops everything into a class object
reactions = []
for i in consolidated:
    reactions.append(Reaction(i[0], i[1], i[2]))


# Get CT Results and Clean it Up
results_sheet = pd.read_excel(FILE_PATH, sheet_name='Results')
data_cleaned_2 = results_sheet.drop(results_sheet.index[0:43])
data_cleaned_2.rename(columns={'96-Well 0.2-mL Block': 'Well', 'Unnamed: 14': 'CT'}, inplace=True)
# data_cleaned_2 = data_cleaned_2.drop([127, 128, 129, 130, 131])
data_cleaned_2.reset_index(inplace=True)

wells_ = data_cleaned_2['Well']
ct_ = data_cleaned_2['CT']

for row in range(data_cleaned_2.shape[0]):
    w = wells_[row]
    ct = ct_[row]
    for r in reactions:
        if w in r.wells:
            r.CT.append(ct)

for r in reactions:
    try:
        r.CT = statistics.mean(r.CT)
    except Exception as e:
        print('Sample Likely Did Not Amplify', r.wells, r.CT, e)


# Dynamic Typed >>>
target_genes = [i for i in set(targets) if type(i) == str]
samples = [i for i in set(samples) if type(i) == str]

# {sample: [{target 1: CT1, target2: CT2,...}, {target 1: dCT1, target2: dCT2,...},
#               {target 1: FC1, target2: FC2,...}]
reactions_sorted_by_sample = {}
for s in samples:
    reactions_sorted_by_sample[s] = [{}, {}, {}]
for r in reactions:
    reactions_sorted_by_sample[r.sample][0][r.target] = r.CT

for key in reactions_sorted_by_sample:
    for rx in reactions_sorted_by_sample[key][0]:
        try:
            dCT = reactions_sorted_by_sample[key][0][rx] - reactions_sorted_by_sample[key][0][HOUSEKEEPING_GENE]
            reactions_sorted_by_sample[key][1][rx] = dCT
        except Exception as e:
            print('Sample Likely Did Not Amplify', e)

for key in reactions_sorted_by_sample:
    for rx in reactions_sorted_by_sample[key][0]:
        try:
            ddCT = reactions_sorted_by_sample[key][1][rx] - reactions_sorted_by_sample[CONTROL][1][rx]
            reactions_sorted_by_sample[key][2][rx] = round(2**(-ddCT), 2)

        except Exception as e:
            print('Sample Likely Did Not Amplify. Will not yield FC for this target', e)



for key in reactions_sorted_by_sample:
    print(key, 'FC', reactions_sorted_by_sample[key][2])


