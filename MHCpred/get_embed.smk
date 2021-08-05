import os, glob, sys
import numpy as np
import pandas as pd
import joblib
from jax_unirep import get_reps
from dataset import DataBundle, MHCDataset

WKDIR = "/data/bases/fangzq/ImmunoRep/"

PEPTIDE_DATA = os.path.join(WKDIR,"IEDB/data/threshold")
MHC_ALLELES = os.path.join(WKDIR,"IEDB/data/allelenames")
MHC_PSEUDO = os.path.join(WKDIR,"IEDB/data/MHC_pseudo.dat")
peptide_ba_el_dir = os.path.join(WKDIR,"IEDB/data/threshold/")


TRAIN = os.path.join(WKDIR, "MHCDataset.train.pkl"),
VALID = os.path.join(WKDIR, "MHCDataset.valid.pkl"),
TEST = os.path.join(WKDIR, "MHCDataset.test.pkl"),

BA = sorted(glob.glob(os.path.join(PEPTIDE_DATA,  "*.ba")))
EL = sorted(glob.glob(os.path.join(PEPTIDE_DATA, "*.el")))

BA_EMBED = [f"{b}.unirep.npy" for b in BA]
EL_EMBED = [f"{b}.unirep.npy" for b in EL]


###################################################################
workdir: WKDIR


rule target:
    input: BA_EMBED, EL_EMBED, TRAIN, VALID, TEST


rule get_ba_embed:
    input: os.path.join(DATA, "{seq}.thr.ba")
    output: os.path.join(DATA, "{seq}.thr.ba.unirep.npy")
    run:
        p = input[0]
        ba = pd.read_table(input[0], header=None, sep=" ", names= ['unknown','epitope','affinity'])
        h_avg, h_final, c_final= get_reps(ba['epitope'].to_list())
        np.save(f"{p}.unirep.npy", h_avg)

rule get_el_embed:
    input: os.path.join(DATA, "{seq}.thr.el")
    output: os.path.join(DATA, "{seq}.thr.el.unirep.npy")
    run:
        p = input[0]
        ba = pd.read_table(input[0], header=None, sep=" ", names= ['unknown','epitope','affinity'])
        h_avg, h_final, c_final= get_reps(ba['epitope'].to_list())
        np.save(f"{p}.unirep.npy", h_avg)


rule to_mhc_dataset:
    input: 
        embeds = BA_EMBED,
        mhc_pseudo = MHC_PSEUDO,
        mhc_allels = MHC_ALLELES,
    output:
        databundle = os.path.join(WKDIR, "databundle.pkl"),
        train = TRAIN,
        valid = VALID,
        test = TEST,
    params:
        data_dir = PEPTIDE_DATA,
        seed = 1234
    run:
        data = DataBundle(alleles= input.mhc_allels, 
                        mhc_psudo= input.mhc_pseudo, 
                        ba_el_dir= params.data_dir )
        # preprocess data
        data.parse_ba()
        data.parse_el()
        data.concat()
        data.train_val_test_split(seed=params.seed)
        joblib.dump(data, output.databundle)
        # drop unecessary attributes, save memory
        data.el_path =None
        data.ba_path = None
        data.pseudo2ba = None
        data.pseudo2el =None
        data.el_embed = None
        data.el_affi = None
        data = joblib.load(config.DATA_BUNDLE)
        train_data = MHCDataset(data, data.train)
        joblib.dump(train_data, output.train)
        valid_data = MHCDataset(data, data.val)
        joblib.dump(valid_data, output.valid)
        test_data = MHCDataset(data, data.test)
        joblib.dump(test_data, output.test)



# if __name__ == "__main__":

#     ba_path = sorted(glob.glob(os.path.join(sys.argv[1], "*.ba")))
#     for p in ba_path:
#         ba = pd.read_table(p, header=None, sep=" ", names= ['unknown','epitope','affinity'])
#         h_avg, h_final, c_final= get_reps(ba['epitope'].to_list())
#         np.save(f"{p}.unirep.npy", h_avg)
#         print(f"finished {p}")

#     el_path = sorted(glob.glob(os.path.join(sys.argv[1], "*.el")))
#     for p in el_path:
#         el = pd.read_table(p, header=None, sep=" ", names= ['unknown','epitope','affinity'])
#         h_avg, h_final, c_final= get_reps(el['epitope'].to_list())
#         np.save(f"{p}.unirep.npy", h_avg)
#         print(f"finished {p}")