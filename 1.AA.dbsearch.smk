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


## use openMS 
## conda install -c conda-forge -c bioconda openms openms-thirdparty pyopenms
PROT_DB = "SequenceDatabase_MSV000084172/PA_ucsc_proteome_259contams_viruses_26tumorShared.fasta"
PROT_DB_DECOY = "SequenceDatabase_MSV000084172/PA_ucsc_proteome_259contams_viruses_26tumorShared.decoy.fasta"

PERCOLATOR = ""

rule target:
    input: PERCOLATOR


rule cometAdaptor:
    """
    database search using comet
    """
    input:
        mzML = "{sample}.mzML.gz",
        prot_db_fasta = PROT_DB,
    output:
        idXML = "{sample}.pin.idXML",
        percolator_in = "{sample}.pin.tsv",
    params: 
       enzyme = 'unspecific cleavage'
    threads: 16
    output:
    shell:
        "CometAdapter -in {input.mzML} -database {input.prot_db_fasta} "
        "-pin_out {output.percolator_in} -out {output.idXML} "
        "-threads {threads} -enzyme '{params.enzyme}' "

rule decoyBuild:
    """
    generate decopy database for percolator
    """
    input: PROT_DB
    output: PROT_DB_DECOY
    params:
        enzyme = 'unspecific cleavage'
    threads: 6
    shell:
        "DecoyDatabase -in {input} -out {output} " # -only_decoy
        # "-enzyme '{params.enzyme} ' # <- won't work
        "-threads {threads} "


rule percolatorAdapter:
    """
    protein identification from shotgun proteomics datasets
    """
    input: 
        pin = "{sample}.pin.idXML",
        decoy = PROT_DB_DECOY,
    output:
        "{sample}.percolator.out.idXML""
    params:
        enzyme = 'no_enzyme'
    threads: 6
    shell:
        "PercolatorAdapter -in {input.pin} -in_decoy {input.decoy} "
        "-out {output} "
        "-enzyme {params.enzyme} -threads {threads}"


rule percolatorAdapter:
    """
    protein identification from shotgun proteomics datasets
    """
    input: 
        pin = "{sample}.pin.tsv",
        decoy = PROT_DB_DECOY,
    output:
        "{sample}.percolator.out.xml"
    params:
        enzyme = 'no_enzyme'
    threads: 6
    shell:
        "percolator --xmloutput {ouput} " # or -X
        "--protein-enzyme {params.enzyme} "
        "--num-threads {threads} "
        "{input.pin} " # -in_decoy {input.decoy}