class Well:
    def __init__(self, well):
        self.well = well
        self.list_of_delta_rn = data_cleaned[data_cleaned["Well Letter"] == well]['Delta Rn']
        self.delta_rn_array = np.array(self.list_of_delta_rn)
        self.cycle_array = np.arange(start=1, stop=41, step=1)
class Reaction:
    def __init__(self, sample, target, wells: list):
        """sample name, PCR target, and duplicate wells to be averaged"""
        self.sample = sample
        self.target = target
        self.wells = wells
        self.CT = []

def plot(well_obj):
    plt.scatter(x=well_obj.cycle_array, y=well_obj.list_of_delta_rn, label=str(well_obj.well), s=6)
    plt.plot(well_obj.cycle_array, well_obj.list_of_delta_rn)


def fold_change(dCT_target, dCT_reference):
    ddCT = dCT_target - dCT_reference
    FC = 2^(-ddCT)
    return FC


def find_dCT(rxn_target, rxn_housekeeping):
    """rxn1/2 are objects of the Reaction class"""
    CT_target = rxn_target.CT
    CT_housekeeping = rxn_housekeeping.CT
    dCT = CT_target - CT_housekeeping
    return dCT
