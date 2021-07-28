import os, glob, sys
import numpy as np
import pandas as pd
import glob
from jax_unirep import get_reps



if __name__ == "__main__":

    ba_path = sorted(glob.glob(os.path.join(sys.argv[1], "*.ba")))
    for p in ba_path:
        ba = pd.read_table(p, header=None, sep=" ", names= ['unknown','epitope','affinity'])
        h_avg, h_final, c_final= get_reps(ba['epitope'].to_list())
        np.save(f"{p}.unirep.npy", h_avg)
        print(f"finished {p}")

    el_path = sorted(glob.glob(os.path.join(sys.argv[1], "*.el")))
    for p in el_path:
        el = pd.read_table(p, header=None, sep=" ", names= ['unknown','epitope','affinity'])
        h_avg, h_final, c_final= get_reps(el['epitope'].to_list())
        np.save(f"{p}.unirep.npy", h_avg)
        print(f"finished {p}")