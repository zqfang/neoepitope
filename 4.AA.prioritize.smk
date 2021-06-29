import os, glob, sys 
import pandas as pd


# configfile: "config.yaml"
workdir: config['WORKDIR']

# scripts path
SMKPATH = config['SMKPATH']
# working directory
WKDIR = config['WORKDIR']

##### INPUTS #####################
KNAPSACK =  config['KNAPSACK']
SPECTRUM = "spectrums.mgf"
LOCDICT = "spectrums.location.pytorch.pkl"
PROTEOME = config["PROTEOME_PLUS_CONTAMINATION"]
DENOVO_ONLY = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"
ACCURACY_ALL_LABELED = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.accuracy"
##### OUTPUTS ##########################
DENOVO_PSMS = DENOVO_ONLY+"_2nd.pin.psms"

# ================================================================================
# Step 4.2: Run second round of PEAKS X DB search against the list of database and de novo peptides. 
# ================================================================================

include: "rules/aa_workflow_step_4.py"

# include: "rules/aa_workflow_step_5.py"


num_fractions = 11
patient_id = "Mel16"


##################### ##############################
rule target:
    input:  DENOVO_PSMS, # "aa_workflow.step_5.output_neoantigen_criteria.csv"

# Before running PEAKS, we need to combine database and de novo peptides into a list.
# This script will select unique de novo peptides, filter out those that belong to the human Swiss-Prot protein database, 
# and combine the remaining de novo peptides and the database peptides identified from Step 1 into a fasta file.
rule generate_denovo_peptides:
    input:
        denovo = DENOVO_ONLY,
        fasta = PROTEOME,
        labeled = "features.csv.labeled.mass_corrected",
    output:
        fasta = "combined_peptide_list.fasta"
    run:
        # from aa_workflow_step_4 
        ## FIXME: the output fasta contains modifications !!!
        ## HOW to handle this !!!
        preprocess( denovo_file=input.denovo,
                    db_fasta_file=input.fasta,
                    labeled_feature_file=input.labeled,
                    peptide_list_fasta=output.fasta)


# The numbers of de novo and database peptides are reported as following:
#   "Number of top-scoring denovo peptides: 17318"
#   "num_db_peptides = 25274"
#   "num_denovo_peptides = 6444" (not in database)

# Run PEAKS X DB search with as following:
#   Select the DENOVO node result from Step 1.1, and select PEAKS DB search;
#   Select option "No digestion" for "Digest mode";
#   Select the fasta file "aa_workflow.step_4.peptide_list.fasta" as the only database, no contaminant;
#   Leave other settings the same as in Step 1.1.
# Set FDR 1.0% and export the "DB search psm.csv" file, rename it to "aa_workflow.step_4.psm.csv".


# second run db search
rule get_denovo_only_features:
    input: 
        denovo = DENOVO_ONLY,
        features = "features.csv.mass_corrected",
    output:  DENOVO_ONLY+"_2nd"
    run:
        # denovo_ids = []
        # with open(input.denovo, 'r') as denovo:
        #     for line in denovo:
        #         if line.startswith("feature_id"): continue
        #         idx = line.strip().split()[0]
        #         denovo_ids.append(idx)
        denovo = pd.read_table(input.denovo)['feature_id'].dropna().tolist()
        features = pd.read_csv(input.features, index_col='spec_group_id')
        out = features.loc[denovo,:]
        out.to_csv(output[0])

rule search_db:
    input:
        spectrums = SPECTRUM,
        locdict = LOCDICT,
        features = DENOVO_ONLY+"_2nd",
        fweight = "checkpoints/forward_deepnovo.pth",
        bweight = "checkpoints/backward_deepnovo.pth",
        knapsack = KNAPSACK,
        fasta =  "combined_peptide_list.fasta", # need this database from above output to do the search
    output:
        pin = DENOVO_ONLY+"_2nd.pin",
    params:
        batch_size = 16,
        modelpath = SMKPATH,
    shell:
        # FIXME: if fasta input contains modifications, how to handle this !?
        # run on the test set with min5
        "python {params.modelpath}/PointNovo/main.py --search_db "
        "--weight_dir checkpoints "
        "--denovo_feature {input.features} "
        "--spectrum {input.spectrums} "
        "--location_dict {input.locdict} "
        "--knapsack {input.knapsack} "
        "--fasta {input.fasta} "

rule percolator:
    input: DENOVO_ONLY+"_2nd.pin",
    output:
        psms = DENOVO_ONLY+"_2nd.pin.psms",
        xml = "pout.xml" #temp("pout.xml"),
    params:
        pbin = "/home/fangzq/miniconda/bin",
        tmp = "",
    threads: 1
    shell:
        "{params.pbin}/percolator -X {output.xml} "
        "--protein-enzyme no_enzyme "
        "{input} > {output.psms} "

# Extract de novo peptides from the PSMs of PEAKS X DB search round 2.
# ======================= UNCOMMENT and RUN ======================================
# ~ aa_workflow_step_4_2.postprocess(
    # ~ psm_file = data_training_dir + "aa_workflow.step_4.psm.csv",
    # ~ output_denovo_peptide_file = data_training_dir + "aa_workflow.step_4.output_peptide_list")
# ================================================================================
# The number of de novo peptides is reported as following:
#   "num_denovo_peptides = 1259"

rule extract_denovo_peptides:
    input: DENOVO_ONLY+"_2nd.pin.psms",
    output: "denovo_peptide_final_list"
    run:
        # from aa_workflow_step_4
        postprocess(psm_file = input[0], 
                    output_denovo_peptide_file = output[0])


# ================================================================================
# Step 5: Neoantigen selection. 
# ================================================================================
# "aa_workflow.step_4.output_peptide_list" +
# "aa_workflow.step_4.peptide_list.fasta" => netmhcpan, immunogencity(IEDB)
# need to get snp files
rule neoantigen_selection:
    input:
        psm = "aa_workflow.step_4.psm.csv",
        netmhcpan = "aa_workflow.step_5.netmhcpan.csv", 
        immunogenicity = "aa_workflow.step_5.immunogenicity.csv",
        target_decoy =  PROTEOME,
        labeled = "feature.csv.labeled",
        snp = "aa_workflow.step_5.supp_data5_snp.csv",
        snp_enst = "aa_workflow.step_5.supp_data5_snp_enst.fasta",
    output:
        neoantigen = "aa_workflow.step_5.output_neoantigen_criteria.csv",
        protein_mutation = "aa_workflow.step_5.protein_mutation.csv",
    params:
        snp_sample_id = patient_id, 
    run:
        # from aa_workflow_step_5
        step_5( psm_file=input.psm,
                netmhc_file=input.netmhcpan,
                immunogenicity_file=input.immunogenicity,
                db_fasta_file=input.target_decoy,
                labeled_feature_file=input.labeled,
                snp_file=input.snp,
                snp_enst_fasta=input.snp_enst,
                snp_sample_id=params.snp_sample_id,
                output_neoantigen_criteria=output.neoantigen,
                output_protein_mutation=output.protein_mutation)
