import os, glob, sys 


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

##### OUTPUTS ##########################


# ================================================================================
# Step 4.2: Run second round of PEAKS X DB search against the list of database and de novo peptides. 
# ================================================================================

# Before running PEAKS, we need to combine database and de novo peptides into a list.
# This script will select unique de novo peptides, filter out those that belong to the human Swiss-Prot protein database, 
# and combine the remaining de novo peptides and the database peptides identified from Step 1 into a fasta file.
# ======================= UNCOMMENT and RUN ======================================
# ~ aa_workflow_step_4_2.preprocess(
    # ~ denovo_file=data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only",
    # ~ db_fasta_file=data_fasta_dir + "uniprot_sprot.human.plus_contaminants.fasta",
    # ~ labeled_feature_file=data_training_dir + "feature.csv.labeled.mass_corrected",
    # ~ peptide_list_fasta=data_training_dir + "aa_workflow.step_4.peptide_list.fasta")
# ================================================================================
include: "rules/aa_workflow_step_4.py"
include: "rules/aa_workflow_step_5.py"


num_fractions = 11
patient_id = "Mel16"


##################### ##############################
rule target:
    ouput: "aa_workflow.step_5.output_neoantigen_criteria.csv"

rule generate_denovo_peptides:
    input:
        denovo = "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only",
        fasta = "uniprot_sprot.human.plus_contaminants.fasta",
        labeled = "feature.csv.labeled.mass_corrected",
    output:
        fasta = "aa_workflow.step_4.peptide_list.fasta"
    run:
        # from aa_workflow_step_4 
        preproces( denovo_file=input.denovo,
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
# rule PeaksDBSearch_2:
#     input: 
#         fasta = "aa_workflow.step_4.peptide_list.fasta",
#     output:  "aa_workflow.step_4.psm.csv"
#     shell:
#         "touch {output}"


rule search_db:
    input:
        spectrums = SPECTRUM,
        locdict = LOCDICT,
        test = "features.csv.labeled.mass_corrected.test.nodup",
        features = "features.csv.labeled.mass_corrected",
        predict = "features.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5",
        fweight = "checkpoints/forward_deepnovo.pth",
        bweight = "checkpoints/backward_deepnovo.pth",
        knapsack = KNAPSACK,
    output:
        psm = "aa_workflow.step_4.psm.csv",
        percolator_out = "pout.tsv"
    params:
        batch_size = 16,
        modelpath = SMKPATH,
    shell:
        # run on the test set with min5
        "python {params.modelpath}/PointNovo/main.py --search_db "
        "--train_dir checkpoints "
        "--search_feature {input.features} "
        "--spectrum {input.spectrums} "
        "--location_dict {input.locdict} "
        "--knapsack {input.knapsack} "
        "--db_output {output.percolator_out} "


# Extract de novo peptides from the PSMs of PEAKS X DB search round 2.
# ======================= UNCOMMENT and RUN ======================================
# ~ aa_workflow_step_4_2.postprocess(
    # ~ psm_file = data_training_dir + "aa_workflow.step_4.psm.csv",
    # ~ output_denovo_peptide_file = data_training_dir + "aa_workflow.step_4.output_peptide_list")
# ================================================================================
# The number of de novo peptides is reported as following:
#   "num_denovo_peptides = 1259"

rule extract_denovo_peptides:
    input: "aa_workflow.step_4.psm.csv"
    output: "aa_workflow.step_4.output_peptide_list"
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
        target_decoy =  "uniprot_sprot.human.plus_contaminants.fasta",
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
