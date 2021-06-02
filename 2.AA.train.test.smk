import os

from deepnovo_preprocess import *
from deepnovo_postprocess import *
import aa_workflow_step_4_2
import aa_workflow_step_5


data_fasta_dir = "data.fasta/"
patient_id = "Mel16"
data_training_dir = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/"
num_fractions = 11
model_dir = "train.mel_16.class_1" # before training, create this empty folder at the same level as Python scripts.


# ================================================================================
# Step 2: Train personalized DeepNovo model.
# ================================================================================

# This step 2 took about 12 hours on a server with GPU Titan X, 32 GB memory

# Note that you will need to specify the paths to your own data and model folders when you run the Python scripts. The following scripts just show examples of my data and model folders.

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

rule spectrum_merge:
    input: expand("{samples}.mgf", sample=SAMPLES)
    output: "spectrums.mgf"
    shell:
        "cat {input} > {output}"

rule feautres_merge_split:
    input: expand("{samples}.features.csv", sample=SAMPLES)
    output: 
        features="features.csv",
        label="features.csv.labeled",
        unlabel="feautures.csv.unlabeled",
    run:
        df = []
        for p in input:
            df.append(pd.read_csv(p))
        df = pd.concat(df)
        df.to_csv(output.features, index=False)
        # seq = '' => unlablel, else labeled
        df_unlabel = df[df['seq'].isna()]
        df_label = df[df['seq'].notnull()]
        df_unlabel.to_csv(output.unlabel, index=False)
        df_label.to_csv(output.label, index=False)
        
# Run calculate_mass_shift_ppm() and correct_mass_shift_ppm()
# ======================= UNCOMMENT and RUN ======================================
# ~ labeled_feature_file = data_training_dir + "feature.csv.labeled"
# ~ ppm = calculate_mass_shift_ppm(labeled_feature_file)
# ~ input_feature_file = data_training_dir + "feature.csv.labeled"
# ~ correct_mass_shift_ppm(input_feature_file, ppm)
# ~ input_feature_file = data_training_dir + "feature.csv"
# ~ correct_mass_shift_ppm(input_feature_file, ppm)
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
        ppm = calculate_mass_shift_ppm(input.label)
        correct_mass_shift_ppm(input.label, ppm)
        correct_mass_shift_ppm(input.feautures, ppm)

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
    input:   "feature.csv.labeled.mass_corrected"
    output:
        train =  "features/{sample}.features.csv.train.nodup",
        valid =  "features/{sample}.features.csv.valid.nodup",
        test =  "features/{sample}.features.csv.test.nodup",
        denovo = "features/{sample}.features.csv.denovo.nodup",
    params:
        ratio = [0.90, 0.05, 0.05],
        path = smkpath,
        outdir = "features"
    shell:
        "python {params.path}/preprocess/train_val_test.py {input} {params.outdir}"

# ================================================================================
# Step 2.2: Training DeepNovo model.
# ================================================================================
rule train:
    input:
        spectrums = "spectrums.mgf",
        train = "feature.csv.labeled.mass_corrected.train.nodup",
        valid = "feature.csv.labeled.mass_corrected.valid.nodup",
    output:
        "train/forward_deepnovo.pth",
        "train/backwad_deepnovo.path",
    params:
        batch_size = 16,
        epoch = 50,
        learning_rate = 0.001,
        modelpath = "PoinNovo",
    shell:
        "python {modelpath}/main.py --train "
        "--spectrum {input.spectrums} "
        "--train_feature {input.train} "
        "--valid_feature {input.valid} "

# Run DeepNovo training
# The training will stop after 10 epoch. The model with best performance on the valid set, "ckpt-16200" is saved in the model folder "train.mel_16.class_1".
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --train"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--train_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--train_feature", data_training_dir + "feature.csv.labeled.mass_corrected.train.noshare"]
# ~ command += ["--valid_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--valid_feature", data_training_dir + "feature.csv.labeled.mass_corrected.valid.noshare"]
# ~ command += ["--reset_step"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================
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
        modelpath = "PoinNovo",
    shell:
        shell("python {modelpath}/main.py --search_denovo "
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


# ================================================================================
# Step 3: Perform personalized de novo sequencing with DeepNovo.
# ================================================================================

# This step 3 took about 5 hours on a server with GPU Titan X, 32 GB memory

# Run DeepNovo de novo sequencing on all features (label and unlabeled)
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --search_denovo"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--denovo_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--denovo_feature", data_training_dir + "feature.csv.mass_corrected"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================
# The de novo results will be written to the file "feature.csv.mass_corrected.deepnovo_denovo".
# The tool will also report the number of features that have been processed:
#   "Total spectra: 694565"
#     "read: 690354"
#     "skipped: 4211"
#       "by mass: 4211"


rule search_denovo:
    input:
        spectrums = "spectrums.mgf",
        denovo = "feature.csv.mass_corrected",
        fweight = "train/forward_deepnovo.pth",
        bweight = "train/backwad_deepnovo.path",
    output:
        "feature.csv.mass_corrected.deepnovo_denovo",
    params:
        batch_size = 16,
        epoch = 50,
        learning_rate = 0.001,
        modelpath = "PoinNovo",
    shell:
        shell("python {modelpath}/main.py --search_denovo "
        "--spectrum {input.spectrums} "
        "--denovo_feature {input.denovo} ")

