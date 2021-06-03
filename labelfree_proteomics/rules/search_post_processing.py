# IDPosteriorErrorProbability
rule calc_peptide_posterior_error:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/dbsearch_{datafile}.idXML"
    output:
        idxml = temp("work/{dbsearchdir}/{datafile}/idpep_{datafile}.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/{datafile}/idpep_{datafile}.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/{datafile}/idpep_{datafile}.log'
    shell:
        "IDPosteriorErrorProbability "
        "-in {input.idxml} "
        "-out {output.idxml} "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# PeptideIndexer
rule index_peptides:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/idpep_{datafile}.idXML",
        fasta = 'work/database/target_decoy_database.fasta'
    output:
        idxml = temp("work/{dbsearchdir}/{datafile}/pi_{datafile}.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/{datafile}/pi_{datafile}.benchmark.txt"
    params:
        decoy_string = "-decoy_string {0}".format(config["database"]["decoy_string"]),
        decoy_string_position = "-decoy_string_position {0}".format(config["database"]["decoy_string_position"]),
        missing_decoy_action = "-missing_decoy_action {0}".format(config["database"]["missing_decoy_action"]),
        enzyme = "-enzyme:name {0}".format(config["digestion"]["enzyme"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/{datafile}/pi_{datafile}.log'
    shell:
        "PeptideIndexer "
        "-in {input.idxml} "
        "-fasta {input.fasta} "
        "-out {output.idxml} "
        "-allow_unmatched "
        "-IL_equivalent "
        "-enzyme:specificity none "
        "-threads {threads} "
        "{params.decoy_string} "
        "{params.decoy_string_position} "
        "{params.missing_decoy_action} "
        "{params.enzyme} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# FalseDiscoveryRate:
rule peptide_fdr:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/pi_{datafile}.idXML"
    output:
        idxml = temp("work/{dbsearchdir}/{datafile}/fdr_{datafile}.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/{datafile}/fdr_{datafile}.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/{datafile}/fdr_{datafile}.log'
    shell:
        "FalseDiscoveryRate "
        "-in {input.idxml} "
        "-out {output.idxml} "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "
