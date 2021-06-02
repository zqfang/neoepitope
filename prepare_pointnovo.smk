

import os, sys, glob
import joblib
import pandas as pd

smkpath = "/data/bases/fangzq/ImmunoRep/PointNovo"
DATA = "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
workdir:  DATA # working directory

mzML = glob.glob(os.path.join(DATA, "peaks/*gz"))
SAMPLES = [os.path.basename(m).replace(".mzML.gz", "") for m in mzML ]

PERCOLATOR = ["percolator/{s}.percolator.target.peptides.txt"  for s in SAMPLES ]
MGF = expand("mgf/{sample}.mgf", sample=SAMPLES)
FEATURES = expand("mgf/{sample}.features.csv", sample=SAMPLES)
LOCATION = expand("mgf/{sample}.mgf.location.pytorch.pkl", sample=SAMPLES)

TRAIN_SAMPLES = [s for s in SAMPLES if s.startswith('train_')]
TEST_SAMPLES = [s for s in SAMPLES if s.startswith('test_')]
DATASETS = expand("features/{sample}.features.csv.{dat}.nodup", sample=TRAIN_SAMPLES, dat=['train','test', 'valid','denovo']) # , 'denovo'])

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


rule mgf2location:
    """
    This step is really slow, and it will be very slow when train or test.
    """
    input: 
        mgf = "mgf/{sample}.mgf",
    output:
        loc = "mgf/{sample}.mgf.location.pytorch.pkl"
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



rule correct_mass_shift:
    input:
        features="features.csv",
        label="features.csv.labeled",
    output:
        features="features.csv.mass_corrected",
        label ="features.csv.labeled.mass_corrected" 
    run:
        ppm = calculate_mass_shift_ppm(input.label)
        correct_mass_shift_ppm(input.label, ppm)
        correct_mass_shift_ppm(input.feautures, ppm)


rule train_val_test:
    input:   "mgf/{sample}.features.csv"
    output:
        train =  "features/{sample}.features.csv.train.nodup",
        valid =  "features/{sample}.features.csv.valid.nodup",
        test =  "features/{sample}.features.csv.test.nodup",
        denovo = "features/{sample}.features.csv.denovo.nodup",
    params:
        ratio = [0.8, 0.1, 0.1],
        path = smkpath,
        outdir = "features"
    shell:
        "python {params.path}/preprocess/train_val_test.py {input} {params.outdir}"


