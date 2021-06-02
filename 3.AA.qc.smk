# ================================================================================
# Step 4: Quality control.
# ================================================================================

# ================================================================================
# Step 4.1: Post-process de novo results to improve their accuracy. 
# ================================================================================

# Run select_top_score()
# This script selects a threshold of de novo confidence scores and uses it to filter de novo results.
# The score threshold is calculated based on a 95% cutoff of the testing accuracy obtained at the end of Step 2 above.
# ======================= UNCOMMENT and RUN ======================================
# ~ accuracy_cutoff = 0.95
# ~ accuracy_file = data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo.accuracy"
# ~ score_cutoff = find_score_cutoff(accuracy_file, accuracy_cutoff)
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo"
# ~ output_file = input_file + ".top95"
# ~ select_top_score(input_file, output_file, score_cutoff)
# ================================================================================
# After this step we'll get the file "feature.csv.mass_corrected.deepnovo_denovo.top95".
# The score cutoff and the number of selected features will also be reported:
#   "score_cutoff =  -0.5"
#   "total_feature =  690354"
#   "select_feature =  233589"

rule select_top_score:
    input:
        test_acc = "feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo.accuracy",
        denovo = "feature.csv.mass_corrected.deepnovo_denovo"
    output:
        top95 = "feature.csv.mass_corrected.deepnovo_denovo.top95"
    params:
        accuracy_cutoff = 0.95
    run:
        score_cutoff = find_score_cutoff(input.test_acc, params.accuracy_cutoff)
        select_top_score(input.denvo, output.top95, params.score_cutoff)


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
        convert_I_to_L(input.top95, output.i2l)
        correct_by_consensus(output.i2l, output.consensus)
        filter_by_minlen(output.consensus, output.minlen, params.minlen)

# Run convert_I_to_L()
# This script converts I (Isoleucine) to L (Leucine) in all de novo peptides, because de novo sequencing is not able to distinguish them.
# ======================= UNCOMMENT and RUN ======================================
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95"
# ~ output_file = input_file + ".I_to_L"
# ~ convert_I_to_L(input_file, output_file)
# ================================================================================


# Run correct_by_consensus()
# This script corrects de novo sequencing errors by grouping predicted sequences of the same mass together and voting the consensus sequence.
# ======================= UNCOMMENT and RUN ======================================
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L"
# ~ output_file = input_file + ".consensus"
# ~ correct_by_consensus(input_file, output_file)
# ================================================================================

        
# Run filter_by_minlen()
# This script filters out sequences of length less than 5 amino acids.
# ======================= UNCOMMENT and RUN ======================================
# ~ minlen = 5
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus"
# ~ output_file = input_file + ".minlen" + str(minlen)
# ~ filter_by_minlen(input_file, output_file, minlen)
# ================================================================================
# The numbers of features will be reported as:
#   "total_feature =  233589"
#   "minlen_feature =  223507"
#   "removed_feature =  10082"


# Up to this step, we get the following file: 
#   "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"



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
        shell("python {modelpath}/main.py --test "
        "--spectrum {input.spectrums} "
        "--test_feature {input.features} "
        "--predicted_file {input.predict}")


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
# We get these results:
#   "precision_AA_mass_db  = 0.9530"
#   "precision_peptide_mass_db  = 0.8441"


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
# We get these results:
#   "precision_AA_mass_db  = 0.9797"
#   "precision_peptide_mass_db  = 0.9371"
# Note that these accuracy results look better than those against the test set because the test set was not used for training the model.
# The number of de novo only features is also reported as
#   "predicted_only: 68721"
# and they are written to the file 
#   "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"

