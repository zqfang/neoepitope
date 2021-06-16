
import os, glob, sys 


# configfile: "config.yaml"
workdir: config['WORKDIR']

# scripts path
SMKPATH = config['SMKPATH']
# working directory
WKDIR = config['WORKDIR']

##### INPUTS #####################
KNAPSACK = "knapsack.npy"
SPECTRUM = "spectrums.mgf"
LOCDICT = "spectrums.location.pytorch.pkl"
FEATURES = expand("features.csv.labeled.mass_corrected.{dat}.nodup", dat=['train','test', 'valid','denovo'])

##### OUTPUTS ##########################
WEIGHTS = ["checkpoints/forward_deepnovo.pth", "checkpoints/backward_deepnovo.pth"]
DENOVO_PREDICT_ALL = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.deepnovo_denovo.accuracy"
TEST = expand("features.csv.labeled.mass_corrected.test.nodup.deepnovo_denovo.{ext}", ext= ["denovo_only", "accuracy"])

# ================================================================================
# Step 2: Train personalized PointNovo model.
# ================================================================================
include: "rules/aa_postprocess.py"


rule target:
    input: WEIGHTS, DENOVO_PREDICT_ALL

# ================================================================================
# Step 2.2: Training DeepNovo model.
# ================================================================================
rule train:
    input:
        train = "features.csv.labeled.mass_corrected.train.nodup",
        valid = "features.csv.labeled.mass_corrected.valid.nodup",
        locdict = LOCDICT,
        spectrums = SPECTRUM,
        knapsack = KNAPSACK,
    output:
        "checkpoints/forward_deepnovo.pth",
        "checkpoints/backward_deepnovo.pth",
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
        model = SMKPATH,
    shell:
        "python {params.model}/PointNovo/main.py --train --train_dir checkpoints "
        "--spectrum {input.spectrums} --location_dict {input.locdict} "
        "--train_feature {input.train} "
        "--valid_feature {input.valid} "
        "--knapsack {input.knapsack} "        

rule test:
    input:
        test = "features.csv.labeled.mass_corrected.test.nodup",
        locdict = LOCDICT,
        spectrums = SPECTRUM,
        knapsack = KNAPSACK,
        fweight = "checkpoints/forward_deepnovo.pth",
        bweight = "checkpoints/backward_deepnovo.pth",
    output: 
        predict = "features.csv.labeled.mass_corrected.test.nodup.deepnovo_denovo",
        #pred1 = "features.csv.labeled.mass_corrected.test.nodup.deepnovo_denovo.denovo_only",
        accurary = "features.csv.labeled.mass_corrected.test.nodup.deepnovo_denovo.accuracy",
    params:
        # batch_size = 16,
        modelpath = SMKPATH,
    resources:
        partition='gpu',
        gpus=1,
        cpus=8,
        cpu_mem='8g',
        gpu_mem='GPU_MEM:16GB',
        time_min='47:58:00', # less than 2 days
    run:
        # output a "feature.csv.labeled.mass_corrected.test.nodup.deepnovo_denovo"
        shell("python {modelpath}/PoinNovo/main.py  --search_denovo --train_dir checkpoints "
              "--denovo_feature {input.test} "
              "--spectrum {input.spectrums} "
              "--location_dict {input.locdict} "
              "--knapsack {input.knapsack} ") 
        ## ...test.nodup.deepnovo_denovo is a predicted file for test
        shell("python {modelpath}/main.py --test --train_dir checkpoints "
              "--test_feature {input.test} "
              "--predict_feature {output.predict}"
              "--spectrum {input.spectrums} "
              "--location_dict {input.locdict} "
              "--knapsack {input.knapsack} ")
        

# The testing accuracy at the amino acid (AA) and peptide levels will be reported as following:
#   "precision_AA_mass_db  = 0.8425"
#   "precision_peptide_mass_db  = 0.6430"


# ================================================================================
# Step 3: Perform personalized de novo sequencing with DeepNovo.
# ================================================================================

# This step 3 took about 5 hours on a server with GPU Titan X, 32 GB memory

# Run DeepNovo de novo sequencing on all features (label and unlabeled)
# The de novo results will be written to the file "feature.csv.mass_corrected.deepnovo_denovo".

rule search_denovo:
    input:
        locdict = LOCDICT,
        spectrums = SPECTRUM,
        knapsack = KNAPSACK,
        denovo = "features.csv.mass_corrected",
        fweight = "checkpoints/forward_deepnovo.pth",
        bweight = "checkpoints/backward_deepnovo.pth",
    output: 
        "features.csv.mass_corrected.deepnovo_denovo",
    params:
        # batch_size = 16,
        # epoch = 50,
        # learning_rate = 0.001,
        modelpath = SMKPATH,
    shell:
        "python {params.modelpath}/PointNovo/main.py --search_denovo --train_dir checkpoints "
        "--denovo_feature {input.denovo} "
        "--spectrum {input.spectrums} "
        "--location_dict {input.locdict} "
        "--knapsack {input.knapsack} "

# ================================================================================
# Step 4: Quality control.
# Step 4.1: Post-process de novo results to improve their accuracy. 
# ================================================================================

rule select_top_score:
    input:
        test_acc = "features.csv.labeled.mass_corrected.test.nodup.deepnovo_denovo.accuracy",
        denovo = "features.csv.mass_corrected.deepnovo_denovo",
    output:
        top95 = "features.csv.mass_corrected.deepnovo_denovo.top95",
    params:
        accuracy_cutoff = 0.95
    run:
        # This script selects a threshold of de novo confidence scores and uses it to filter de novo results.
        # The score threshold is calculated based on a 95% cutoff of the testing accuracy obtained at rule test above.
        score_cutoff = find_score_cutoff(input.test_acc, params.accuracy_cutoff)
        select_top_score(input.denvo, output.top95, params.score_cutoff)

# After this step we'll get the file "feature.csv.mass_corrected.deepnovo_denovo.top95".
# The score cutoff and the number of selected features will also be reported:
#   "score_cutoff =  -0.5"
#   "total_feature =  690354"
#   "select_feature =  233589"

rule correct_and_filter:
    input:
        top95 = "features.csv.mass_corrected.deepnovo_denovo.top95",
    output:
        i2l = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L",
        consensus = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus",
        minlen = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5",
    params:
        minlen = 5 # filters out sequences of length less than 5 amino acids
    run:
        # converts I (Isoleucine) to L (Leucine) in all de novo peptides, because de novo sequencing is not able to distinguish them 
        convert_I_to_L(input.top95, output.i2l)
        # corrects de novo sequencing errors by grouping predicted sequences of the same mass together and voting the consensus sequence.
        correct_by_consensus(output.i2l, output.consensus)
        # filters out sequences of length less than 5 amino acids.
        filter_by_minlen(output.consensus, output.minlen, params.minlen)

# The numbers of features will be reported as:
#   "total_feature =  233589"
#   "minlen_feature =  223507"
#   "removed_feature =  10082"


# Up to this step, we get the following file: 
#   "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"

rule test_all_labels:
    input:
        spectrums = SPECTRUM,
        test = "features.csv.labeled.mass_corrected.test.nodup",
        features = "features.csv.labeled.mass_corrected",
        predit = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5",
        fweight = "checkpoints/forward_deepnovo.pth",
        bweight = "checkpoints/backward_deepnovo.pth",
    output:
        "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only",
        "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.accuracy",
    params:
        batch_size = 16,
        epoch = 50,
        learning_rate = 0.001,
        modelpath = SMKPATH,
    run:
        shell("python {modelpath}/PointNovo/main.py --test --train_dir checkpoints "
              "--test_feature {input.test} "
              "--predict_feature {output.predict}"
              "--spectrum {input.spectrums} "
              "--location_dict {input.locdict} "
              "--knapsack {input.knapsack} ")
        # We get these results:
        #   "precision_AA_mass_db  = 0.9530"
        #   "precision_peptide_mass_db  = 0.8441"
        shell("python {modelpath}/PointNovo/main.py --test --train_dir checkpoints "
              "--test_feature {input.features} "
              "--predict_feature {output.predict}"
              "--spectrum {input.spectrums} "
              "--location_dict {input.locdict} "
              "--knapsack {input.knapsack} ")
        # We get these results:
        #   "precision_AA_mass_db  = 0.9797"
        #   "precision_peptide_mass_db  = 0.9371"
        # Note that these accuracy results look better than those against the test set because the test set was not used for training the model.
        # The number of de novo only features is also reported as
        #   "predicted_only: 68721"
        # and they are written to the file 
        #   "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"