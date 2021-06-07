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
#####################################################################################################################
## use openMS
## comet + percolator

## conda install -c conda-forge -c bioconda openms openms-thirdparty pyopenms comet percolator
import os, glob 


workdir: "/data/bases/fangzq/ImmunoRep/data"
# MZML = glob.glob("MSV000082648/peaks/test_*.mzML.gz")
# SAMPLES = [os.path.basename(mz).replace(".mzML.gz", "") for mz in MZML]
SAMPLES = ["test_sample_0_ms_run_0"] 

PROT_DB = "uniprot_sprot.fasta"
CONTAMINANT = "UniProtContams_259_20170206.fasta" # download containinants from MSV000084172
TARGET_DECOY = "target.decoy.fasta"
PERCOLATOR = expand("{sample}_perc.mzTab", sample=SAMPLES)

####### OpenMS MHC workflow #########################################################
## for Fido: Search Engine -> PeptideIndexer -> IDPosteriorProbability -> Fido adapter.
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
        enzyme = 'no cleavage', #'unspecific cleavage',
        decoy_string = "DECOY_",
        decoy_string_position = "prefix",
    threads: 6
    shell:
        "DecoyDatabase -in {input.proteome} {input.contams} "
        "-out {output} -threads {threads} "
        "-decoy_string {params.decoy_string} "
        "-decoy_string_position {params.decoy_string_position} "
        "-enzyme '{params.enzyme}' " # <- won't work


rule PeakPicker:
    """optional"""
    input: "{sample}.mzML", 
    output: temp("{sample}.picked.mzML"),
    params:
        pick_ms_levels="2"
    shell:
        "PeakPickerHiRes -in {input} "
        "-out {output}.picked.mzML "
        "-algorithm:ms_levels ${params.pick_ms_levels}"


rule CometAdaptor:
    """
    database search using comet
    """
    input:
        mzML = "{sample}.picked.mzML",
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
        "-precursor_mass_tolerance 5 "
        "-precursor_error_units ppm "
        "-fragment_bin_tolerance 0.02 "
        "-fragment_bin_offset 0 "
        "-num_hits 1 "
        "-max_variable_mods_in_peptide 3"
        "-activation_method 'ALL' "
        "-missed_cleavages 0 "
        #"-fixed_modifications 'Carbamidomethyl (C)' "
        "-variable_modifications 'Oxidation (M)' "


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
        temp("{sample}_pi.idXML")
    params:
        decoy_string = "DECOY_",
        decoy_position = "prefix", 
    log:  'log/pi_{sample}.log'
    threads: 1
    shell:
        "PeptideIndexer -in {input.idxml} -out {output} -fasta {input.fasta} "
        "-IL_equivalent "
        "-decoy_string {params.decoy_string} "
        "-missing_decoy_action warn "
        "-decoy_string_position {params.decoy_position} "
        "-threads {threads} "
        "-enzyme:specificity none "
        "-enzyme:name 'unspecific cleavage' "
        "2>&1 | tee {log} "

rule FalseDiscoveryRate:
    input:  "{sample}_pi.idXML"
    output: temp("{sample}_fdr.idXML")
    threads: 1
    log:  'log/fdr_{sample}.log'
    shell:
        "FalseDiscoveryRate -in {input} -out {output} "
        "-algorithm:add_decoy_peptides "
        "-protein 'false' " # this argument is important
        "-threads {threads} 2>&1 | tee {log} "
        
        
rule PSMFeatureExtractor:
    input: "{sample}_fdr.idXML"
    output: temp("{sample}_fdr_psm.idXML")
    threads: 12
    shell:
        "PSMFeatureExtractor -in {input} -out {output} -threads {threads}"

rule PercolatorAdapter:
    """
    protein identification from shotgun proteomics datasets
    """
    input: "{sample}_fdr_psm.idXML",
    output: protected("{sample}_perc.idXML")
    params:
        enzyme = 'no_enzyme',
        decoy_prefix = "DECOY_",
        description_correct_features = 0,
        subet_max_train = 0
    threads: 6
    log:  'log/percolator_{sample}.log'
    shell:
        "PercolatorAdapter -in {input} -out {output} "
        "-debug 10 -trainFDR 0.05 -testFDR 0.05 "
        "-decoy_pattern {params.decoy_prefix} "
        "-enzyme {params.enzyme} -threads {threads} "
        "-subset_max_train {params.subet_max_train} "
        "-doc {params.description_correct_features} "
        "2>&1 | tee {log}"


# filter by percolator q-value
rule filter_peptides:
    input: "{sample}_perc.idXML"
    output: "{sample}_perc_filt.idXML"
    params:
        fdr = 0.01, # fdr threshold
    threads: 1
    shell:
        "IDFilter -in {input} -out {output} "
        "-score:pep {params.fdr} -threads {threads} "
        "-delete_unreferenced_peptide_hits "
        "-remove_decoys "

rule MzTabExporter:
    """Exports various XML formats to an mzTab file"""
    input: "{sample}_perc_filt.idXML"
    output: "{sample}_perc.mzTab"
    threads: 1
    priority: 53
    shell:
        "MzTabExporter -in {input} -out {output} "