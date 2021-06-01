
import numpy as np
import pandas as pd
from jax_unirep import get_reps


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