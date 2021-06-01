

import os, sys, glob
import joblib
import pandas as pd

smkpath = "/data/bases/fangzq/ImmunoRep/PointNovo"
DATA = "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
workdir:  DATA

mzML = glob.glob(os.path.join(DATA, "peaks/*gz"))
SAMPLES = [os.path.basename(m).replace(".mzML.gz", "") for m in mzML ]
PERCOLATOR = ["percolator/{s}.percolator.target.peptides.txt"  for s in SAMPLES ]
MGF = expand("mgf/{sample}.mgf", sample=SAMPLES)
FEATURES = expand("mgf/{sample}.features.csv", sample=SAMPLES)
LOCATION = "data/spectrum.mgf.location.pytorch.pkl"

TRAIN_SAMPLES = [s for s in SAMPLES if s.startswith('train_')]
TEST_SAMPLES = [s for s in SAMPLES if s.startswith('test_')]
DATASETS = expand("data/features.csv.{dat}.nodup", sample=TRAIN_SAMPLES, dat=['train','test', 'valid','denovo']) # , 'denovo'])


#################################################################


rule target:
    input: MGF, FEATURES, LOCATION, DATASETS

rule mzML2mgf:
    input: 
        mzml = "peaks/{sample}.mzML.gz",
        perlco = "percolator/{sample}.percolator.target.peptides.txt",
    output:
        mgf = "mgf/{sample}.mgf",
        features = "mgf/{sample}.features.csv"
    params:
        path = smkpath
    shell:
        ## better to use hash, not jobid
        "python {params.path}/preprocess/mzml2mgf.py "
        "{input.mzml} {input.perlco} {output.mgf} {output.features} {jobid}"

rule cat_mgf:
    input: MGF
    output: "data/spectrums.mgf"
    shell:
        "cat {input.mgf} > {output.mgf}"

rule cat_features:
    input:
        features = FEATURES,
    output:
        features = "data/features.csv",
    run:
        feats = []
        for feat in input.features:
            f = pd.read_csv(feat, dtype=str)
            feats.append(f)
        feats = pd.concat(feats)
        feats.to_csv(output.features, index=False)

rule mgf2location:
    """
    This step is really slow, and it will be very slow train, test.
    """
    input: 
        mgf = "data/spectrums.mgf"
    output:
        loc = "data/spectrums.mgf.location.pytorch.pkl"
    run:
        spectrum_location_dict = {}
        line = True
        with open(input.mgf, 'r') as f:
            while line:
                # The tell() method returns the current file position in a file stream.
                current_location = f.tell() 
                line = f.readline()
                if "BEGIN IONS" in line:
                    spectrum_location = current_location
                elif "SCANS=" in line:
                    scan = re.split('[=\r\n]', line)[1]
                    spectrum_location_dict[scan] = spectrum_location

        joblib.dump(spectrum_location_dict, output.loc)

rule train_val_test:
    input:  "data/features.csv"
    output:
        train =  "data/features.csv.train.nodup",
        valid =  "data/features.csv.valid.nodup",
        test =  "data/features.csv.test.nodup",
        denovo = "data/features.csv.denovo.nodup",
    params:
        ratio = [0.8, 0.1, 0.1],
        path = smkpath
    shell:
        "python {params.path}/preprocess/train_val_test.py {input} data"


