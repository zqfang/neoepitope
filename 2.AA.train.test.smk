import os, glob, sys 
import joblib 
import pandas as pd

WKDIR = "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
workdir: WKDIR

MZML = sorted(glob.glob(os.path.join(WKDIR, "peaks/*.mzML.gz")))
SAMPLES = [os.path.basename(mz).replace(".mzML.gz", "") for mz in MZML]
WEIGHTS = ["train/forward_deepnovo.pth", "train/backward_deepnovo.pth"]

smkpath = "/home/fangzq/github/neoepitope"


PERCOLATOR = ["percolator/{s}.percolator.target.peptides.txt"  for s in SAMPLES ]
MGF = expand("mgf/{sample}.mgf", sample=SAMPLES)
FEATURES = expand("mgf/{sample}.features.csv", sample=SAMPLES)
LOCATION = expand("mgf/{sample}.mgf.location.pytorch.pkl", sample=SAMPLES)
LOCDICT = "spectrums.location.pytorch.pkl"
TRAIN_SAMPLES = [s for s in SAMPLES if s.startswith('train_')]
TEST_SAMPLES = [s for s in SAMPLES if s.startswith('test_')]
DATASETS = expand("features/{sample}.features.csv.{dat}.nodup", sample=TRAIN_SAMPLES, dat=['train','test', 'valid','denovo']) # , 'denovo'])

FEATURES = expand("features.csv.labeled.mass_corrected.{dat}.nodup", dat=['train','test', 'valid','denovo'])

# ================================================================================
# Step 2: Train personalized PointNovo model.
# ================================================================================
include: "rules/aa_preprocess.py"

rule target:
    input: FEATURES, LOCDICT, WEIGHTS

# This step 2 took about 12 hours on a server with GPU Titan X, 32 GB memory

# Note that you will need to specify the paths to your own data and model folders when you run the Python scripts. 
# The following scripts just show examples of my data and model folders.

# ================================================================================
# Step 2.1: Prepare the training data.
# ================================================================================

# Run merge_mgf_file() and merge_feature_file()
# ======================= UNCOMMENT and RUN ======================================
# ~ folder_path = data_training_dir
# ~ fraction_list = range(0, num_fractions)
# ~ merge_mgf_file(
    # ~ input_file_list=[folder_path + "export_" + str(i) + ".mgf" for i in fraction_list],
    # ~ fraction_list=fraction_list,
    # ~ output_file=folder_path + "spectrum.mgf")
# ~ merge_feature_file(
    # ~ input_file_list=[folder_path + "export_" + str(i) + ".csv" for i in fraction_list],
    # ~ fraction_list=fraction_list,
    # ~ output_file=folder_path + "feature.csv")
# ================================================================================
# We will get two output files in the same folder: "spectrum.mgf" and "feature.csv".
# Both functions also report the number of entries that have been processed: "counter = 694565".
# That number should be the same as the total number of MS/MS spectra from the raw files.


# Run split_feature_unlabel()
# ======================= UNCOMMENT and RUN ======================================
# ~ input_feature_file = data_training_dir + "feature.csv"
# ~ split_feature_unlabel(input_feature_file)
# ================================================================================
# It will split the "feature.csv" into 2 files: "feature.csv.labeled" and "feature.csv.unlabeled".
# It also reports the number of labeled and unlabel features: "num_labeled = 207332" and "num_unlabeled = 487233".
# Note that "207332" is also the number of PSMs reported at FDR 1.0% in Step 1.

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
# ======================= UNCOMMENT and RUN ======================================
# ~ input_feature_file = data_training_dir + "feature.csv.labeled.mass_corrected"
# ~ proportion = [0.90, 0.05, 0.05]
# ~ split_feature_training_noshare(input_feature_file, proportion)
# ================================================================================
# It will split "feature.csv.labeled.mass_corrected" into train/valid/test sets with "proportion = [0.9, 0.05, 0.05]".
# Those 3 sets do not share common peptides.
# Their sizes are also reported.
#   "num_total = 207332"
#   "num_unique = 26656"
#   "num_train = 185823"
#   "num_valid = 10900"
#   "num_test = 10609"

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

# ================================================================================
# Step 2.2: Training DeepNovo model.
# ================================================================================
rule train:
    input:
        spectrums = "spectrums.mgf",
        train = "features.csv.labeled.mass_corrected.train.nodup",
        valid = "features.csv.labeled.mass_corrected.valid.nodup",
        locdict = "spectrums.location.pytorch.pkl",
    output:
        "train/forward_deepnovo.pth",
        "train/backward_deepnovo.pth",
    resources:
        partition='gpu',
        gpus=1,
        cpus=8,
        cpu_mem='8g',
        gpu_mem='GPU_MEM:16GB',
        time_min='47:58:00', # less than 2 days
    params:
        batch_size = 16,
        epoch = 20,
        learning_rate = 0.001,
        model = smkpath,
    shell:
        "python {params.model}/PointNovo/main.py --train "
        "--spectrum {input.spectrums} --location_dict {input.locdict} "
        "--train_feature {input.train} "
        "--valid_feature {input.valid} "        

rule test:
    input:
        spectrums = "spectrums.mgf",
        test = "feature.csv.labeled.mass_corrected.test.nodup",
        fweight = "train/forward_deepnovo.pth",
        bweight = "train/backwad_deepnovo.path",
    output:
        "feature.csv.labeled.mass_corrected.test.nodup.deepnovo_denovo"
    params:
        batch_size = 16,
        epoch = 50,
        learning_rate = 0.001,
        modelpath = smkpath,
    resources:
        partition='gpu',
        gpus=1,
        cpus=8,
        cpu_mem='8g',
        gpu_mem='GPU_MEM:16GB',
        time_min='47:58:00', # less than 2 days
    run:
        shell("python {modelpath}/PoinNovo/main.py --search_denovo "
              "--spectrum {input.spectrums} "
              "--test_feature {input.test} ")

        shell("python {modelpath}/main.py --test "
              "--spectrum {input.spectrums} "
              "--test_feature {input.test} ")
        
# Run DeepNovo testing
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test_true_feeding"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--test_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--test_feature", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --search_denovo"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--denovo_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--denovo_feature", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test"]
# ~ command += ["--target_file", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command += ["--predicted_file", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================

# The testing accuracy at the amino acid (AA) and peptide levels will be reported as following:
#   "precision_AA_mass_db  = 0.8425"
#   "precision_peptide_mass_db  = 0.6430"
