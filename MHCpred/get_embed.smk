import os, glob, sys
import numpy as np
import pandas as pd
import glob
from jax_unirep import get_reps


DATA = "/data/bases/fangzq/ImmunoRep/IEDB/data/threshold"
BA = sorted(glob.glob(os.path.join(DATA, "*.ba")))
EL = sorted(glob.glob(os.path.join(DATA, "*.el")))

BA_EMBED = [f"{b}.unirep.npy" for b in BA]
EL_EMBED = [f"{b}.unirep.npy" for b in EL]



rule target:
    input: BA_EMBED, EL_EMBED


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