
import os, re, sys
import numpy as np
import pandas as pd 
from typing import List, Dict, Tuple


def assign(df: pd.DataFrame, ratio: List[float] = [0.8, 0.1, 0.1]):
    output_train = []
    output_test = []
    output_val = []

    counter_train = 0
    counter_test = 0
    counter_valid = 0
    counter_unique = 0

    peptide_train_set = set()
    peptide_test_set = set()
    peptide_val_set = set()
    
    for i, row in df.iterrows():
        # check if the peptide already exists in any of the three lists
        # if yes, this new feature will be assigned to that list
        peptide = row['seq']

        if (peptide in peptide_train_set):
            output_train.append(i)
            counter_train += 1
        elif (peptide in peptide_val_set):
            output_val.append(i)
            counter_valid += 1
        elif (peptide in peptide_test_set):
            output_test.append(i)
            counter_test += 1
        # if not, this new peptide and its spectrum will be randomly assigned
        else:
            counter_unique += 1
            set_num = np.random.choice(a=3, size=1, p=ratio)
            if set_num == 0:
                output_train.append(i)
                peptide_train_set.add(peptide)
                counter_train += 1
            elif set_num == 1:
                output_val.append(i)
                peptide_val_set.add(peptide)
                counter_valid += 1
            else:
                output_test.append(i)
                peptide_test_set.add(peptide)
                counter_test += 1

    return output_train, output_val, output_test

def main(feature_csv, outdir, ratio = [0.8,0.1,0.1], random_state=42):
    
    filename = os.path.basename(feature_csv)
    df = pd.read_csv(feature_csv)
    # df = df[df['seq'] != ''].reset_index(drop=True)
    df = df.dropna().reset_index(drop=True)
    df1 = df[df['seq'].isna()].reset_index(drop=True)
    # shuffle first 
    df2 = df.sample(frac=1, random_state=random_state)
    # MS/MS data have a special property: 
    # the same peptide could appear multiple times in a dataset with different spectra.

    # It is essential to make sure that the training,
    # validation, and testing sets do not share common peptides.
    train_idx, val_idx, test_idx = assign(df2, ratio)
    train = df2.iloc[train_idx]
    val = df2.iloc[val_idx]
    test = df2.iloc[test_idx]

    # # shuffle data first and split using np.split
    # train, val, test = np.split(df2, [int(ratio[0]*len(df2)), int((ratio[0]+ratio[1])*len(df2))])
    train.to_csv( os.path.join(outdir, filename + ".train.nodup"), index=False)
    val.to_csv( os.path.join(outdir, filename + ".valid.nodup"), index=False)
    test.to_csv( os.path.join(outdir, filename +  ".test.nodup"), index=False)
    df1.to_csv( os.path.join(outdir, filename + ".denovo.nodup"), index=False)

if __name__ == '__main__':
    inp = sys.argv[1]
    outdir = sys.argv[2]
    if len(sys.argv) < 3:
        print("Usage: train_test_val.py input.csv")
        sys.exit(0)
    main(inp, outdir)
