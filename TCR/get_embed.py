import numpy as np
import pandas as pd
from jax_unirep import get_reps



## TODO
## read TRUST4 output results: *_report.tsv 
# header = ["count", "frequency","CDR3nt", "CDR3aa", "V", "D","J","C", "cid", "cid_full_length"]


## select TRB

## select IGH


## Get unirep embeddings for training
X = trb['CDR3aa'].to_list()
y = trb['label'].to_list()
sequences = X
h_avg, h_final, c_final= get_reps(sequences)
np.save("trb.unirep.npy", h_avg)
# each of the arrays will be of shape (len(sequences), 1900),
# with the correct order of sequences preserved
X = igh['CDR3aa'].to_list()
y = igh['label'].to_list()
sequences = X
h_avg, h_final, c_final= get_reps(sequences)
np.save("igh.unirep.npy", h_avg)