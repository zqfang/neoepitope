
# ================================================================================
# Step 4.2: Run second round of PEAKS X DB search against the list of database and de novo peptides. 
# ================================================================================

# Before running PEAKS, we need to combine database and de novo peptides into a list.
# This script will select unique de novo peptides, filter out those that belong to the human Swiss-Prot protein database, and combine the remaining de novo peptides and the database peptides identified from Step 1 into a fasta file.
# ======================= UNCOMMENT and RUN ======================================
# ~ aa_workflow_step_4_2.preprocess(
    # ~ denovo_file=data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only",
    # ~ db_fasta_file=data_fasta_dir + "uniprot_sprot.human.plus_contaminants.fasta",
    # ~ labeled_feature_file=data_training_dir + "feature.csv.labeled.mass_corrected",
    # ~ peptide_list_fasta=data_training_dir + "aa_workflow.step_4.peptide_list.fasta")
# ================================================================================
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

# Extract de novo peptides from the PSMs of PEAKS X DB search round 2.
# ======================= UNCOMMENT and RUN ======================================
# ~ aa_workflow_step_4_2.postprocess(
    # ~ psm_file = data_training_dir + "aa_workflow.step_4.psm.csv",
    # ~ output_denovo_peptide_file = data_training_dir + "aa_workflow.step_4.output_peptide_list")
# ================================================================================
# The number of de novo peptides is reported as following:
#   "num_denovo_peptides = 1259"




# ================================================================================
# Step 5: Neoantigen selection. 
# ================================================================================
# ~ aa_workflow_step_5.step_5(
    # ~ psm_file=data_training_dir + "aa_workflow.step_4.psm.csv",
    # ~ netmhc_file=data_training_dir + "aa_workflow.step_5.netmhcpan.csv",
    # ~ immunogenicity_file=data_training_dir + "aa_workflow.step_5.immunogenicity.csv",
    # ~ db_fasta_file=data_fasta_dir + "uniprot_sprot.human.plus_contaminants.fasta",
    # ~ labeled_feature_file=data_training_dir + "feature.csv.labeled",
    # ~ snp_file=data_training_dir + "aa_workflow.step_5.supp_data5_snp.csv",
    # ~ snp_enst_fasta=data_training_dir + "aa_workflow.step_5.supp_data5_snp_enst.fasta",
    # ~ snp_sample_id=patient_id,
    # ~ output_neoantigen_criteria=data_training_dir + "aa_workflow.step_5.output_neoantigen_criteria.csv",
    # ~ output_protein_mutation=data_training_dir + "aa_workflow.step_5.protein_mutation.csv")
