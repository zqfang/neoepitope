# IDMerge
rule merge_peptides_fido:
    input:
        idxmls = expand(
            "work/{{dbsearch}}/{sample}/pi_{sample}.idXML", sample=SAMPLES)
    output:
        idxml = temp("work/{dbsearch}/proteinid/idmerge.idXML"),
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearch}/proteinid/idmerge.benchmark.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = "work/{dbsearch}/proteinid/idmerge.log"
    shell:
        "IDMerger "
        "-in {input.idxmls} "
        "-out {output.idxml} "
        "-annotate_file_origin "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# IDFilter
rule filter_peptides_fido:
    input:
        idxml = "work/{dbsearch}/proteinid/idmerge.idXML",
    output:
        idxml = temp("work/{dbsearch}/proteinid/idfilter.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearch}/proteinid/idfilter.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = "work/{dbsearch}/proteinid/idfilter.log"
    shell:
        "IDFilter "
        "-in {input.idxml} "
        "-out {output.idxml} "
        "-score:pep 1.0 "
        "-score:prot 0 "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# FidoAdapter
rule fido:
    input:
        idxml = "work/{dbsearch}/proteinid/idfilter.idXML"
    output:
        idxml = temp("work/{dbsearch}/proteinid/fido.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearch}/proteinid/fido.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = "work/{dbsearch}/proteinid/fido.log"
    shell:
        "FidoAdapter "
        "-in {input.idxml} "
        "-out {output.idxml} "
        "-fido_executable Fido "
        "-fidocp_executable FidoChooseParameters "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# FalseDiscoveryRate
rule fido_fdr:
    input:
        idxml = "work/{dbsearch}/proteinid/fido.idXML"
    output:
        idxml = temp("work/{dbsearch}/proteinid/fido_fdr.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearch}/proteinid/fido_fdr.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = "work/{dbsearch}/proteinid/fido_fdr.log"
    shell:
        "FalseDiscoveryRate "
        "-in {input.idxml} "
        "-out {output.idxml} "
        "-PSM false "
        "-protein true "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# IDFilter
rule fido_fdr_filter:
    input:
        idxml = "work/{dbsearch}/proteinid/fido_fdr.idXML"
    output:
        idxml = "work/{dbsearch}/proteinid/fido_fdr_filt.idXML"
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearch}/proteinid/fido_fdr_filt.benchmark.txt"
    params:
        profdr = '-score:pep {0}'.format(config["protein"]["fdr"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = "work/{dbsearch}/proteinid/fido_fdr_filt.log"
    shell:
        "IDFilter "
        "-in {input.idxml} "
        "-out {output.idxml} "
        "-score:pep 0.0 "
        "{params.profdr} "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "
