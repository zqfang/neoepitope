

include: "rules/aa_postprocess.py"


rule target:
    output: "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"



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
        bweight = "train/backwad_deepnovo.pth",
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

# ================================================================================
# Step 4: Quality control.
# Step 4.1: Post-process de novo results to improve their accuracy. 
# ================================================================================

rule select_top_score:
    input:
        test_acc = "feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo.accuracy",
        denovo = "feature.csv.mass_corrected.deepnovo_denovo"
    output:
        top95 = "feature.csv.mass_corrected.deepnovo_denovo.top95"
    params:
        accuracy_cutoff = 0.95
    run:
        # This script selects a threshold of de novo confidence scores and uses it to filter de novo results.
        # The score threshold is calculated based on a 95% cutoff of the testing accuracy obtained at the end of Step 2 above.
        score_cutoff = find_score_cutoff(input.test_acc, params.accuracy_cutoff)
        select_top_score(input.denvo, output.top95, params.score_cutoff)

# After this step we'll get the file "feature.csv.mass_corrected.deepnovo_denovo.top95".
# The score cutoff and the number of selected features will also be reported:
#   "score_cutoff =  -0.5"
#   "total_feature =  690354"
#   "select_feature =  233589"


rule correct_and_filter:
    input:
        top95 = "feature.csv.mass_corrected.deepnovo_denovo.top95"
    output:
        i2l = "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L"
        consensus = "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus"
        minlen = "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen"
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


# We test its accuracy against the test set:
# Run DeepNovo testing
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test"]
# ~ command += ["--target_file", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command += ["--predicted_file", data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================


# Repeat the same testing but now against all labeled features:
# Run DeepNovo testing
# ====================== UNCOMMENT and RUN =======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test"]
# ~ command += ["--target_file", data_training_dir + "feature.csv.labeled.mass_corrected"]
# ~ command += ["--predicted_file", data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================

rule test:
    input:
        spectrums = "spectrums.mgf",
        test = "feature.csv.labeled.mass_corrected.test.nodup",
        features = "feature.csv.labeled.mass_corrected",
        predit = "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen"
        fweight = "train/forward_deepnovo.pth",
        bweight = "train/backwad_deepnovo.path",
    output:
        "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"
    params:
        batch_size = 16,
        epoch = 50,
        learning_rate = 0.001,
        modelpath = "PoinNovo",
    shell:
        shell("python {modelpath}/main.py --test "
        "--spectrum {input.spectrums} "
        "--test_feature {input.test} "
        "--predicted_file {input.predict}")
        # We get these results:
        #   "precision_AA_mass_db  = 0.9530"
        #   "precision_peptide_mass_db  = 0.8441"
        shell("python {modelpath}/main.py --test "
        "--spectrum {input.spectrums} "
        "--test_feature {input.features} "
        "--predicted_file {input.predict}")
        # We get these results:
        #   "precision_AA_mass_db  = 0.9797"
        #   "precision_peptide_mass_db  = 0.9371"
        # Note that these accuracy results look better than those against the test set because the test set was not used for training the model.
        # The number of de novo only features is also reported as
        #   "predicted_only: 68721"
        # and they are written to the file 
        #   "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"