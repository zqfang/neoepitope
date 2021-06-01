

import os, sys, glob
import joblib

DATA = "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
workdir:  DATA

mzML = glob.glob(os.path.join(DATA, "peaks/*gz"))
SAMPLES = [os.path.basename(m).replace(".mzML.gz", "") for m in mzML ]
PERCOLATOR = ["percolator/{s}.percolator.target.peptides.txt"  for s in SAMPLES ]
MGF = expand("mgf/{sample}.mgf", sample=SAMPLES)
FEATURES = expand("mgf/{sample}.features.csv", sample=SAMPLES)
LOCATION = expand("mgf/{sample}.location.dict.pkl", sample=SAMPLES)

TRAIN_SAMPLES = [s for s in SAMPLES if s.startswith('train_')]
TEST_SAMPLES = [s for s in SAMPLES if s.startswith('test_')]
DATASETS = expand("datasets/{sample}.features.csv.{dat}.nodup", sample=TRAIN_SAMPLES, dat=['train','test', 'valid', 'denovo'])


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
    shell:
        "python mzml2mgf.py "
        "{input.mzml} {input.perlco} {output.mgf} {output.features} {jobid}"

rule mgf2location:
    input: 
        mgf = "mgf/{sample}.mgf"
    output:
        loc = "mgf/{sample}.location.dict.pkl"
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
    input:  "mgf/{sample}.features.csv"
    output:
        train =  "datasets/{sample}.features.csv.train.nodup",
        valid =  "datasets/{sample}.features.csv.valid.nodup",
        test =  "datasets/{sample}.features.csv.test.nodup",
        denovo = "datasets/{sample}.features.csv.denovo",
    params:
        ratio = [0.8, 0.1, 0.1]
    shell:
        "python train_val_test.py {input} datasets"