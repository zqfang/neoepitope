import os, glob, sys 
import joblib 
import pandas as pd

# scripts path
smkpath = "/home/fangzq/github/neoepitope"
# working directory
WKDIR = "/data/bases/fangzq/ImmunoRep/data/MSV000082648"

##### INPUTS #####################
MZML = sorted(glob.glob(os.path.join(WKDIR, "peaks/*.mzML.gz")))
SAMPLES = [os.path.basename(mz).replace(".mzML.gz", "") for mz in MZML]
PERCOLATOR = ["percolator/{s}.percolator.target.peptides.txt"  for s in SAMPLES ]

TRAIN_SAMPLES = [s for s in SAMPLES if s.startswith('train_')]
TEST_SAMPLES = [s for s in SAMPLES if s.startswith('test_')]

#### OUTPUTS #####################
# output files write to the working directory 
MGF = expand("mgf/{sample}.mgf", sample=SAMPLES)
FEATURES = expand("mgf/{sample}.features.csv", sample=SAMPLES)
LOCATION = expand("mgf/{sample}.mgf.location.pytorch.pkl", sample=SAMPLES)
LOCDICT = "spectrums.location.pytorch.pkl"

WEIGHTS = ["train/forward_deepnovo.pth", "train/backward_deepnovo.pth"]
DATASETS = expand("features/{sample}.features.csv.{dat}.nodup", sample=TRAIN_SAMPLES, dat=['train','test', 'valid','denovo']) 
FEATURES = expand("features.csv.labeled.mass_corrected.{dat}.nodup", dat=['train','test', 'valid','denovo'])

# ================================================================================
# Step 2: Preprocess
# ================================================================================
workdir: WKDIR
include: "rules/aa_preprocess.py"

rule target:
    input: FEATURES, LOCDICT, WEIGHTS

# ================================================================================
# Step 2.1: Prepare the training data.
# ================================================================================

# Run merge_mgf_file() and merge_feature_file()
# We will get two output files in the same folder: "spectrum.mgf" and "feature.csv".
rule mzML2mgf:
    input: 
        mzml = "peaks/{sample}.mzML.gz",
        perlco = "percolator/{sample}.percolator.target.peptides.txt",
    output:
        mgf = "mgf/{sample}.mgf",
        features = "mgf/{sample}.features.csv"
    params:
        path = smkpath,
        sample_idx = lambda wildcards: SAMPLES.index(wildcards.sample), 
    shell:
        ## better to use hash, not jobid
        "python {params.path}/rules/mzml2mgf.py "
        "{input.mzml} {input.perlco} {output.mgf} {output.features} {params.sample_idx} " # {jobid}


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
        "python {params.path}/rules/train_val_test.py {input} {params.outdir}"


rule spectrum_merge:
    input: expand("mgf/{sample}.mgf", sample=TRAIN_SAMPLES)
    output: "spectrums.mgf"
    shell:
        "cat {input} > {output}"

rule spectrum2loc:
    """
    This step is really slow, and it will be very slow when train or test.
    """
    input: 
        mgf = "spectrums.mgf",
    output:
        loc = "spectrums.location.pytorch.pkl"
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


# It will split the "feature.csv" into 2 files: "feature.csv.labeled" and "feature.csv.unlabeled".
rule feautres_merge_split:
    input: expand("mgf/{sample}.features.csv", sample=TRAIN_SAMPLES)
    output: 
        features="features.csv",
        label="features.csv.labeled",
        unlabel="feautures.csv.unlabeled",
    run:
        df = []
        for p in input:
            df.append(pd.read_csv(p, dtype=str))
        df = pd.concat(df)
        df.to_csv(output.features, index=False)
        # seq = '' => unlablel, else labeled
        df_unlabel = df[df['seq'].isna()]
        df_label = df[df['seq'].notnull()]
        df_unlabel.to_csv(output.unlabel, index=False)
        df_label.to_csv(output.label, index=False)

# Run calculate_mass_shift_ppm() and correct_mass_shift_ppm()
# ================================================================================
# The mass shift is calculated from "feature.csv.labeled".
# The mass shift ppm (part per million) is reported as: "mean_precursor_ppm = 7.07514819678".
# Then mass is corrected for 2 files: "feature.csv.labeled.mass_corrected" and "feature.csv.mass_corrected".
rule correct_mass_shift:
    input:
        features="features.csv",
        label="features.csv.labeled",
    output:
        features="features.csv.mass_corrected",
        label ="features.csv.labeled.mass_corrected" 
    run:
        # from aa_preprocess 
        ppm = calculate_mass_shift_ppm(input.label)
        correct_mass_shift_ppm(input.label, ppm)
        correct_mass_shift_ppm(input.features, ppm)

# Run split_feature_training_noshare()
# ================================================================================
# It will split "feature.csv.labeled.mass_corrected" into train/valid/test sets with "proportion = [0.9, 0.05, 0.05]".
# Those 3 sets do not share common peptides.
rule split_feature_training_noshare:
    input:   "features.csv.labeled.mass_corrected"
    output:
        train =  "features.csv.labeled.mass_corrected.train.nodup",
        valid =  "features.csv.labeled.mass_corrected.valid.nodup",
        test =  "features.csv.labeled.mass_corrected.test.nodup",
        denovo = "features.csv.labeled.mass_corrected.denovo.nodup",
    params:
        ratio = [0.90, 0.05, 0.05],
        path = smkpath,
        outdir = "."
    shell:
        "python {params.path}/rules/train_val_test.py {input} {params.outdir}"

