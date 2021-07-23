import os, glob 


# configfile: "config.yaml"
workdir: config['WORKDIR']
# working directory
WKDIR = config['WORKDIR']
# scripts path
SMKPATH = config['SMKPATH']

######## INPUTS ##########
MZML = sorted(glob.glob(os.path.join(config['WORKDIR'], "peaks/*.mzML.gz")))
SAMPLES = [mz.split("/")[-1].replace(".mzML.gz", "") for mz in MZML]
#SAMPLES = ["test_sample_0_ms_run_2"] 

PROT_DB = config['dbsearch']['prot_db']
FULL_DB = config['dbsearch']['full_db']
TARGET_DECOY = "target.decoy.fasta"
CONTAMINANT = config['dbsearch']['contaminant']  # download containinants from MSV000084172
CRUX_PERCO = expand("crux/{sample}/bullseye.{sample}.mzML.pid.mgf", sample=SAMPLES)
CRUX_MGF = expand("crux/{sample}/percolator.target.peptides.txt", sample=SAMPLES)
FEATURES = expand("mgf/{sample}.features.csv", sample=SAMPLES)

# SAMPLE indexes
with open(os.path.join(WKDIR,"sample_indexes.txt"), 'w') as out:
    for i, s in enumerate(SAMPLES):
        out.write(f"{s}\tF{i}\n")
######################### CRUX pipeline for MHC peptide ##########################################
rule target:
    input: CRUX_PERCO, CRUX_MGF, FEATURES

# # this is for tide-search
# rule tide_index:
#     """protein database search, fasta input
#        output for tide-search, not for commet
#     """
#     input: FULL_DB
#     output:
#         "target.decoy.index/auxlocs",
#         "target.decoy.index/pepix",
#         "target.decoy.index/protix",
#     params:
#         index = "target.decoy.index"
#         decoy_prefix = "decoy_",
#         enzyme = "no-enzyme"
#     shell:
#         "crux tide-index --decoy-prefix {params.decoy_prefix} "
#         "--enzyme {params.enzyme} "
#         "--output-dir crux_decoy "
#         "--overwrite true "
#         "{input} {params.index}"

# rule Profile2Peak:
#     """optional, if your input is Profile-mode data
#        then run this ->  centroided-mode MS spectra 
#        DecoyDatabase (OpenMS) 
#     """
#     input: "peaks/{sample}.mzML.gz" 
#     output: "commet_percolator/{sample}.picked.mzML",
#     shell:
#         "PeakPickerHiRes -in {input} "
#         "-out {output}.picked.mzML "

rule decompress:
    input:  "peaks/{sample}.mzML.gz"
    output: temp("temp/{sample}.mzML")
    shell:
        "zcat {input} > {output}"

rule Crux:
    """ NOTE::
        1) Need to find best parameters if your input is high/low resolution ms.

        2) The scanNumber and TITLE clarification:
        September 16, 2020 
        the TITLE line of mgf file can consist of two formats. 
        The first format is TITLE=fileName.scanNumber.scanNumber.charge.dta. 
        The second format is 
        TITLE=fileName.scanNumber.scanNumber.charge File:"fileNameWithExtension", NativeID:"controllerType=0 controllerNumber=1 scan=scanNumber"

        BUT, scan column of the `percolator.targets.peptide.txt` is just the index of spectrums in the mgf file !!! 
        
    """
    input:
        mzml = "temp/{sample}.mzML",
        fasta = FULL_DB, #TARGET_DECOY,
    output: 
        tsv = "crux/{sample}/percolator.target.peptides.txt",
        # the fragmentation spectra for which accurate masses were successfully inferred.
        # pid.mgf will be the input of commet
        hits = "crux/{sample}/bullseye.{sample}.mzML.pid.mgf",
        # the fragmentation spectra for which accurate masses were not inferred.
        no_hits = "crux/{sample}/bullseye.{sample}.mzML.nopid.mgf" 
    params:
        enzyme = 'no_enzyme', # immunopeptidomes, but it's for protein-level FDR
        decoy_prefix = "decoy_",
        bin = "/home/fangzq/program/crux-4.0.Linux.x86_64/bin"
    threads: 12
    log:  'log/crux_{sample}.log'
    shell:
        ## crux version must be >= 4.0
        "{params.bin}/crux pipeline "
        "--bullseye true "
        "--search-engine comet "
        "--post-processor percolator "
        "--protein-enzyme {params.enzyme} "
        "--num-threads {threads} "
        #"--allowed_missed_cleavage 0 "
        "--decoy-prefix {params.decoy_prefix} "
        "--decoy_search 1 --overwrite true "
        "--spectrum-format mgf "
        "--output-dir crux/{wildcards.sample}  "
        "--exact-p-value true "
        # "--digest_mass_range '800.0 2500.0' " # for 8-12 MHC I peptide
        # "--fragment-tolerance 0.02 "
        ## high resolution setting for comet, if low res, use default value: 1.0005079
        #"--fragment_bin_tol 0.02 " 
        #"--fragment_bin_offset 0 "
        ## automatically set fragment-torelrance, invoke param-medic to figure out what bin size to use
        "--auto_fragment_bin_tol warn " 
        "--peptide_mass_tolerance 5 "
        "--spectrum_batch_size 500 "
        "{input.mzml} {input.fasta} "
        "2>&1 | tee {log} "


rule crux2features:
    input:
        mgf = "crux/{sample}/bullseye.{sample}.mzML.pid.mgf",
        perc = "crux/{sample}/percolator.target.peptides.txt",
    output:
        mgf = "mgf/{sample}.mgf",
        features = "mgf/{sample}.features.csv"
    params:
        path = SMKPATH,
        sample_idx = lambda wildcards: SAMPLES.index(wildcards.sample), 
    shell:
        ## better to use hash, not jobid
        "python {params.path}/rules/crux2mgf.py "
        "{input.mgf} {input.perc} {output.mgf} {output.features} {params.sample_idx} "

