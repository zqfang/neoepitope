# IDFilter
rule filter_peptides_sc:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/fdr_{datafile}.idXML"
    output:
        idxml = temp("work/{dbsearchdir}/{datafile}/sc_filt_{datafile}.idXML")
    singularity:
        config['singularity']['default']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    threads: 1
    benchmark:
        "work/{dbsearchdir}/{datafile}/pi_filt_ur_{datafile}.benchmark.txt"
    params:
        pepfdr = '-score:pep {0}'.format(config["peptide"]["fdr"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/{datafile}/pi_filt_ur_{datafile}.log'
    shell:
        "IDFilter "
        "-in {input.idxml} "
        "-out {output.idxml} "
        "{params.pepfdr} "
        "-score:prot 0 "
        "-delete_unreferenced_peptide_hits "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "


# ProteinQuantifierSeparate
rule count_spectra_perfile:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/sc_filt_{datafile}.idXML",
        fido = "work/{dbsearchdir}/proteinid/fido_fdr_filt.idXML"
    output:
        csv = temp("csv/{datafile}_{dbsearchdir}_proteinCounts.csv")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "csv/{datafile}_{dbsearchdir}_proteinCounts.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = "csv/{datafile}_{dbsearchdir}_proteinCounts.log"
    singularity:
        config['singularity']['default']
    shell:
        "ProteinQuantifier "
        "-in {input.idxml} "
        "-protein_groups {input.fido} "
        "-out {output.csv} "
        "-top 0 "
        "-average sum "
        "-include_all "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# ProteinQuantifierCombiner
rule combine_spectral_counts:
    input:
        csvs = expand(
            "csv/{sample}_{{dbsearchdir}}_proteinCounts.csv", sample=SAMPLES)
    output:
        csv = "csv/sc_{dbsearchdir}_proteinCounts.csv"
    params:
        names = expand("{sample}", sample=SAMPLES),
        log = 'false'
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "csv/sc_{dbsearchdir}_proteinCounts.benchmark.txt"
    run:
        import pandas as pd
        import numpy as np
        import glob
        import os

        files = input.csvs
        df = pd.read_csv(files[0], sep="\t", skiprows=2)
        colnames = df.columns

        newnames = [colnames[0]]
        for i in range(1, len(colnames)):
            name = colnames[i]
            newnames.append(params.names[0] + name)
        df.columns = newnames

        for i in range(1, len(files)):
            dft = pd.read_csv(files[i], sep="\t", skiprows=2)

            colnames = dft.columns
            newnames = [colnames[0]]
            for j in range(1, len(colnames)):
                name = colnames[j]
                newnames.append(params.names[i] + "_" + name)
            dft.columns = newnames

            df = pd.merge(df, dft, how="outer", on='protein')

        df.to_csv(output.csv)
