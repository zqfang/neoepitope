


import os, glob 


# configfile: "config.yaml"
workdir: config['WORKDIR']

# scripts path
SMKPATH = config['SMKPATH']

######## INPUTS ##########
MZML = glob.glob(os.path.join(config['WORKDIR'], "peaks/*.mzML.gz"))
SAMPLES = [mz.split("/")[-1].replace(".mzML.gz", "") for mz in MZML]
#SAMPLES = ["train_sample_0_ms_run_0"] 

PROT_DB = config['dbsearch']['prot_db']
FULL_DB = config['dbsearch']['full_db']

CONTAMINANT = config['dbsearch']['contaminant']  # download containinants from MSV000084172
CRUX_TSV = expand("{sample}/bullseye.{sample}.mzML.pid.mgf", sample=SAMPLES)
CRUX_MGF = expand("{sample}/percolator.target.peptides.txt", sample=SAMPLES)
######################### CRUX pipeline for MHC peptide ##########################################
rule target:
    input: CRUX_TSV, CRUX_MGF

# rule PeakPicker:
#     """optional, if your input is high resolution data
#        Profile data ->  centroided MS spectra 
#     """
#     input: "peaks/{sample}.mzML.gz" 
#     output: "commet_percolator/{sample}.picked.mzML",
#     shell:
#         "PeakPickerHiRes -in {input} "
#         "-out {output}.picked.mzML "

rule decompress:
    input:  "peaks/{sample}.mzML.gz"
    output: "commet_percolator/{sample}.mzML"
    shell:
        "zcat {input} > {output}"


rule tide_index:
    """protein database search, fasta input"""
    input: FULL_DB
    output:
        "target.decoy.index/auxlocs",
        "target.decoy.index/pepix",
        "target.decoy.index/protix",
    params:
        index = "target.decoy.index"
        decoy_prefix = "decoy_",
        enzyme = "no-enzyme"
    shell:
        "crux tide-index --decoy-prefix {params.decoy_prefix} "
        "--enzyme {params.enzyme} "
        "--output-dir crux_decoy "
        "--overwrite true "
        "{input} {params.index}"

rule crux:
    input:
        mzml = "commet_percolator/{sample}.mzML",
        aux = "target.decoy.index/auxlocs",
        pepix = "target.decoy.index/pepix",
        protix = "target.decoy.index/protix",
    output: 
        tsv = "{sample}/percolator.target.peptides.txt",
        # the fragmentation spectra for which accurate masses were successfully inferred.
        # pid.mgf will will gave be the input of commet
        hits = "{sample}/bullseye.{sample}.mzML.pid.mgf",
        # the fragmentation spectra for which accurate masses were not inferred.
        no_hits = "{sample}/bullseye.{sample}.mzML.bullseye.nopid.mgf" 
    params:
        index = "target.decoy.index"
        enzyme = 'no_enzyme',
        decoy_prefix = "decoy_",
        decoy_string_position = "prefix",
    threads: 12
    output:
    shell:
        "crux pipeline --bullseye true "
        "--search-engine comet --post-processor percolator "
        "--protein-enzyme {params.enzyme} --num-threads {threads} "
        "--allowed_missed_cleavage 0 "
        "--decoy-prefix {params.decoy_prefix} "
        "--decoy_search 2 --overwrite true "
        "--spectrum-format mgf "
        "--mzid-output true "
        "--pepxml-output true "
        "--output-dir {wildcards.sample}  "
        # "--digest_mass_range '800.0 2500.0 " # for 8-12 MHC I peptide
        "--fragment-tolerance 0.02 "
        # --fragment_bin_tol 0.02 ???
        "--peptide_mass_tolerance 3 "
        #"--fragment_bin_offset 0 "
        "{input.mzml} {params.index} "