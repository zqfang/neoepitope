# ================================================================================
# Workflow of neoantigen discovery by personalized de novo sequencing.
# ================================================================================

# Step-by-step instructions based on the following example dataset:

#   Patient Mel-16 (Bassani-Sternberg et al., Nature Communication, 2016)
#   HLA class 1: 12 raw files, 1 failed to run PEAKS

#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_1_A
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_1_B
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_2_A
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_2_B
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_3_A
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_3_B
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_1_A_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_1_B_1, failed
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_2_A_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_2_B_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_3_A_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_3_B_1


# ================================================================================
# Step 1: Build the immunopeptidome of the patient.
# ================================================================================

# This step 1 took about ?? hours on a laptop with 4 CPU cores i7, 16 GB memory

# ================================================================================
# Step 1.1: Run PEAKS X DB search on the raw files with the following parameters:
# ================================================================================

#     Enzyme: None
#     Instrument: Orbi-Orbi
#     Fragment: HCD
#     Acquisition: DDA

#     Parent Mass Error Tolerance: 15.0 ppm
#     Fragment Mass Error Tolerance: 0.05 Da
#     Precursor Mass Search Type: monoisotopic
#     Enzyme: None
#     Digest Mode: Unspecific
#     Max Missed Cleavages: 100
#     Variable Modifications:
#       Oxidation (M): 15.99
#       Deamidation (NQ): 0.98
#     Max Variable PTM Per Peptide: 3
#     Database: uniprot_sprot.human
#     Taxon: All
#     Contaminant Database: contaminants_maxquant
#     Searched Entry: 20488
#     FDR Estimation: Enabled
#     Merge Options: no merge
#     Precursor Options: corrected
#     Charge Options: no correction
#     Filter Options: no filter
#     Process: true
#     Associate chimera: no




# ================================================================================
# Step 1.2: Set FDR 1.0%.
# ================================================================================

# The number of MS/MS spectra is "694565", the number of peptide-spectrum matches (PSMs) is "207332", the number of peptide sequences is "26594".

# ================================================================================
# Step 1.3: Right-click on the DB search node "??", select "Deep Denovo Export".
# ================================================================================

# We will get the following 11 pairs of csv and mgf files in the PEAKS project folder:

#       export_0.csv, export_0.mgf
#       export_1.csv, export_1.mgf
#       export_2.csv, export_2.mgf
#       export_3.csv, export_3.mgf
#       export_4.csv, export_4.mgf
#       export_5.csv, export_5.mgf
#       export_6.csv, export_6.mgf
#       export_7.csv, export_7.mgf
#       export_8.csv, export_8.mgf
#       export_9.csv, export_9.mgf
#       export_10.csv, export_10.mgf

#####################################################################################################################
#
# workflow using comet + percolator for neoantigen discovery
# 
# Refer to MHCquant to get more features and results
#####################################################################################################################
## use openMS
## comet + percolator

## conda install -c conda-forge -c bioconda openms openms-thirdparty pyopenms comet percolator
import os, glob 


# configfile: "config.yaml"
workdir: config['WORKDIR']

# scripts path
SMKPATH = config['SMKPATH']
# working directory
WKDIR = config['WORKDIR']

# MZML = glob.glob("MSV000082648/peaks/test_*.mzML.gz")
# SAMPLES = [os.path.basename(mz).replace(".mzML.gz", "") for mz in MZML]
SAMPLES = ["train_sample_0_ms_run_0"] 

PROT_DB = config['dbsearch']['protein_database']
CONTAMINANT = config['dbsearch']['contaminant']  # download containinants from MSV000084172
TARGET_DECOY = "target.decoy.fasta"
PERCOLATOR = expand("{sample}_perc.mzTab", sample=SAMPLES)

####### OpenMS MHC workflow #########################################################
## for Percolator: Search Engine (comet/MSGFPlus..) -> PeptideIndexer > FalseDiscoveryRate 
##                   -> PSMFeatureExtractor -> PercolatorAdapter -> IDFilter.
#####################################################################################


rule target:
    input: PERCOLATOR

rule DecoyDatabase:
    """
    generate target decopy database
    """
    input: 
        proteome = PROT_DB,
        contams = CONTAMINANT,
    output: 
        TARGET_DECOY
    params:
        enzyme = 'no cleavage', #'unspecific cleavage' is not working,
        decoy_string = "DECOY_",
        decoy_string_position = "prefix",
    threads: 6
    shell:
        "DecoyDatabase -in {input.proteome} {input.contams} "
        "-out {output} -threads {threads} "
        "-decoy_string {params.decoy_string} "
        "-decoy_string_position {params.decoy_string_position} "
        #"-enzyme '{params.enzyme}' " # <- won't work


rule PeakPicker:
    """optional, if your input is high resolution data
    """
    input: "{sample}.mzML", 
    output: temp("{sample}.picked.mzML"),
    params:
        pick_ms_levels="2"
    shell:
        "PeakPickerHiRes -in {input} "
        "-out {output}.picked.mzML "
        "-algorithm:ms_levels {params.pick_ms_levels}"


rule CometAdaptor:
    """
    database search using comet
    """
    input:
        mzML = "{sample}.mzML",
        fasta = TARGET_DECOY, 
    output:
        idXML = protected("{sample}_comet.idXML"),
        #percolator_in = "{sample}.comet.pin.tsv",
    params: 
       enzyme = 'unspecific cleavage'
    threads: 16
    output:
    shell:
        "CometAdapter -in {input.mzML} -database {input.fasta} "
        "-out {output.idXML} " # -pin_out {output.percolator_in}
        "-threads {threads} -enzyme '{params.enzyme}' "
        "-precursor_charge 2:3 "
        "-precursor_mass_tolerance 5 "
        "-precursor_error_units ppm "
        #"-fragment_bin_tolerance 0.02 "
        "-fragment_mass_tolerance 0.02 "
        "-fragment_bin_offset 0 "
        "-num_hits 1 "
        "-max_variable_mods_in_peptide 3 "
        "-activation_method ALL " # CID, ECD, ETD, PDD, HCD, IRMPD
        "-missed_cleavages 0 " # no effect of unspecific cleavage
        #"-allowed_missed_cleavages 0 "
        "-fixed_modifications  "
        "-variable_modifications 'Oxidation (M)' "
        "-digest_mass_range 800:2500 "
        "-spectrum_batch_size 500 "


# IDPosteriorErrorProbability
rule calc_peptide_posterior_error:
    input:
        idxml = "{sample}_comet.idXML"
    output:
        idxml = temp("{sample}_idpep.idXML")
    threads: 1
    log:  'log/idpep_{sample}.log'
    shell:
        "IDPosteriorErrorProbability "
        "-in {input.idxml} -out {output.idxml} -debug 10 "
        "-threads {threads}  2>&1 | tee {log} "

rule PeptideIndexer:
    input: 
        idxml = "{sample}_idpep.idXML",
        fasta = TARGET_DECOY,
    output: 
        "{sample}_pi.idXML"
    params:
        decoy_string = "DECOY_",
        decoy_position = "prefix", 
    log:  'log/pi_{sample}.log'
    threads: 1
    shell:
        "PeptideIndexer -in {input.idxml} -out {output} -fasta {input.fasta} "
        "-IL_equivalent " # note here
        "-decoy_string {params.decoy_string} "
        "-missing_decoy_action warn "
        "-decoy_string_position {params.decoy_position} "
        "-threads {threads} "
        "-enzyme:specificity none "
        #"-enzyme:name 'unspecific cleavage' "
        "2>&1 | tee {log} "

rule FalseDiscoveryRate:
    input:  "{sample}_pi.idXML"
    output: "{sample}_fdr.idXML"
    threads: 1
    log:  'log/fdr_{sample}.log'
    shell:
        "FalseDiscoveryRate -in {input} -out {output} "
        "-protein 'false' " # this argument is important
        "-threads {threads} 2>&1 | tee {log} "
        
        
# filter 
rule filter_fdr_for_idalignment:
    input: "{sample}_fdr.idXML"
    output: "{sample}_fdr_filt.idXML"
    params:
        fdr = 0.0, # fdr threshold
        pep_min = 8,
        pep_max = 12,
    threads: 1
    shell:
        "IDFilter -in {input} -out {output} "
        "-score:pep {params.fdr} -threads {threads} "
        "-delete_unreferenced_peptide_hits "
        "-precursor:length '{params.pep_min}:{params.pep_max}' "
        "-remove_decoys "

# #  compute alignment rt transformation
# # This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them
# # Corrects retention time distortions between maps, using information from peptides identified in different maps
rule align_idx_files:
    input: expand("{sample}_fdr_filt.idXML", sample=SAMPLES),  ## Note select the same sample but with multiple run ?
    output: expand("{sample}_fdr_filt.trafoXML", sample=SAMPLES),  ###
    params:
        max_rt_alignment_shift=300,
    shell:
        "MapAlignerIdentification -in {input} "
        "-trafo_out {output} "
        "-model:type linear "
        "-algorithm:max_rt_shift {params.max_rt_alignment_shift} "

# # 
rule align_mzml_files: #using trafoXMLs
    input:
        mzml= "{sample}.mzML",
        trafo= "{sample}_fdr_filt.trafoXML",
    output:
        "{sample}.aligned.mzML"
    threads: 1
    shell:
        "MapRTTransformer -in {input.mzml} "
        "-trafo_in {input.trafo} "
        "-out {output} " # ${Sample}_${Condition}_${id}_aligned.mzML "
        "-threads {threads} "

rule align_peptideindex_files:
    """this rule is used when quntification is skipped
    """
    input:
        idxml="{sample}_pi.idXML",
        trafo= "{sample}_fdr_filt.trafoXML",
    output:
        "{sample}_pi_aligned.idXML"
    threads: 1
    shell:
        "MapRTTransformer -in {input.idxml} "
        "-trafo_in {input.trafo} "
        "-out {output} "   #" ${Sample}_${Condition}_${id}_idx_aligned.idXML "
        "-threads {threads} "

rule merge_aligned_peptideindex_files:
    input: "{sample}_pi_aligned.idXML"
    output: "{sample}_all_ids_merged.idXML"
    shell:
        "IDMerger -in {input} "
        "-out {output} "
        "-threads {threads} "
        "-annotate_file_origin "
        "-merge_proteins_add_PSMs "


# input from merge_align_idxml_files

rule PSMFeatureExtractor:
    input:  "{sample}_all_ids_merged.idXML", 
    output: "{sample}_all_ids_merged_psm.idXML"
    threads: 12
    shell:
        "PSMFeatureExtractor -in {input} -out {output} -threads {threads}"

rule PercolatorAdapter:
    """
    protein identification from shotgun proteomics datasets
    """
    input: "{sample}_all_ids_merged_psm.idXML",
    output: "{sample}_all_ids_merged_psm_perc.idXML",
    params:
        enzyme = 'no_enzyme',
        decoy_prefix = "DECOY_",
        description_correct_features = 0,
        subet_max_train = 0,
        extra = "-klammer " # or false
    threads: 6
    log:  'log/percolator_{sample}.log'
    shell:
        ## set OMP_NUM_THREADS=${task.cpus} 
        "PercolatorAdapter -in {input} -out {output} "
        "-debug 10 -trainFDR 0.05 -testFDR 0.05 "
        "-decoy_pattern {params.decoy_prefix} "
        "-enzyme {params.enzyme} -threads {threads} "
        "-subset_max_train {params.subet_max_train} "
        "-peptide_level_fdrs "
        "-doc {params.description_correct_features} "
        "-seed 4711 "
        "2>&1 | tee {log}"


# filter by percolator q-value
rule filter_peptides:
    input:  "{sample}_all_ids_merged_psm_perc.idXML",
    output: "{sample}_all_ids_merged_psm_perc_filt.idXML",
    params:
        fdr = 0.0, # fdr threshold
        pep_min = 8,
        pep_max = 12,
    threads: 1
    shell:
        "IDFilter -in {input} -out {output} "
        "-score:pep {params.fdr} -threads {threads} "
        "-delete_unreferenced_peptide_hits "
        "-precursor:length '{params.pep_min}:{params.pep_max}' "
        "-remove_decoys "

rule MzTabExporter:
    """Exports various XML formats to an mzTab file"""
    input: "{sample}_all_ids_merged_psm_perc_filt.idXML",
    output: "{sample}_perc.mzTab",
    threads: 1
    priority: 53
    shell:
        "MzTabExporter -in {input} -out {output} "


#  quantify identifications using targeted feature extraction
rule quantify_identifications_targeted:
    input:
        aligned = "{sample}.aligned.mzML",
        idxml = "{sample}_all_ids_merged_psm_perc_filt.idXML",
    output:
        "{sample}.featureXML",
    params: 
        quantification_fdr = False, #  Assess and assign ids matched between runs with an additional quantification FDR
        quantification_min_prob = 0, # Specify a minimum probability cut off for quantification
    threads: 1
    run:
        if not params.quantification_fdr:
            shell("FeatureFinderIdentification -in {input.aligned} "
                  "-id {input.idxml} -out{output} -threads {threads}")
        else:
            shell("FeatureFinderIdentification -in {input.aligned} "
                  "-id {input.idxml} -id_ext ${id_file_quant} "
                  "-svm:min_prob ${params.quantification_min_prob} "
                  "-out {output} -threads {threads} ")



# # # link extracted features
# rule link_extracted_features:
#     input:  "{sample}.featureXML",
#     output: "{sample}_all_features_merged.consensusXML",
#     threads: 1
#     shell:
#        "FeatureLinkerUnlabeledKD -in {input} -out {output} -threads {threads}"


# # - resolve conflicting ids matching to the same feature
# rule resolve_conflicts:
 
#     input:
#      set val(Sample), file(consensus) from consensus_file

#     output:
#      set val(Sample), file("${Sample}_resolved.consensusXML") into (consensus_file_resolved, consensus_file_resolved_2)

#     when:
#      !params.skip_quantification

#     shell:
#      """
#      IDConflictResolver -in ${consensus} \
#                         -out ${Sample}_resolved.consensusXML \
#                         -threads ${task.cpus}
#      """

# rule export_text:
#     publishDir "${params.outdir}/"
 
#     input:
#      set val(Sample), file(consensus_resolved) from consensus_file_resolved

#     output:
#      set val(Sample), file("${Sample}.csv") into consensus_text

#     shell:
#      """
#      TextExporter -in ${consensus_resolved} \
#                   -out ${Sample}.csv \
#                   -threads ${task.cpus} \
#                   -id:add_hit_metavalues 0 \
#                   -id:add_metavalues 0 \
#                   -id:peptides_only
#      """