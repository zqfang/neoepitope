# FeatureFinderCentroided
rule find_features_centroid:
    input:
        mzml = "mzml/{datafile}.mzML"
    output:
        featurexml = "work/featurefindercentroided/{datafile}/ffc_{datafile}.featureXML"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/featurefindercentroided/{datafile}/ffc_{datafile}.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/featurefindercentroided/{datafile}/ffc_{datafile}.log'
    singularity:
        config['singularity']['default']
    shell:
         "FeatureFinderCentroided "
         "-in {input.mzml} "
         "-out {output.featurexml} "
         "-threads {threads} "
         "{params.debug} "
         "2>&1 | tee {params.log} "

# IDFilter
rule filter_peptides_ffc:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/fdr_{datafile}.idXML"
    output:
        idxml = temp("work/{dbsearchdir}/{datafile}/ffc_filt_{datafile}.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        'work/{dbsearchdir}/{datafile}/pi_filt_ur_{datafile}.benchmark.txt'
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

#IDMapper
rule map_ffc_features:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/ffc_filt_{datafile}.idXML",
        featurexml = "work/featurefindercentroided/{datafile}/ffc_{datafile}.featureXML"
    output:
        featurexml = temp("work/{dbsearchdir}/{datafile}/ffc_filt_idmap_{datafile}.featureXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        'work/{dbsearchdir}/{datafile}/ffc_filt_idmap_{datafile}.benchmark.txt'
    params:
        centroid_rt = "-feature:use_centroid_rt ",
        centroid_mz = "-feature:use_centroid_mz ",
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/{datafile}/ffc_filt_idmap_{datafile}.log'
    shell:
        "IDMapper "
        "-id {input.idxml} "
        "-in {input.featurexml} "
        "-out {output.featurexml} "
        "{params.centroid_rt} "
        "{params.centroid_mz} "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

#MapAlignerPoseClustering:
rule align_ffc_maps:
    input:
        featurexmls = expand("work/{{dbsearchdir}}/{sample}/ffc_filt_idmap_{sample}.featureXML",sample=SAMPLES)
    output:
        featurexmls =
        temp(expand("work/{{dbsearchdir}}/{sample}/ffc_filt_idmap_align_{sample}.featureXML",sample=SAMPLES))
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        'work/{dbsearchdir}/ffc_filt_idmap_align.benchmark.txt'
    params:
        mz_max_difference = "-algorithm:pairfinder:distance_MZ:max_difference {0}".format(config["precursor"]["tolerance"]),
        mz_unit = "-algorithm:pairfinder:distance_MZ:unit {0}".format(config["precursor"]["units"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffc_filt_idmap_align.log'
    shell:
        "MapAlignerPoseClustering "
        "-in {input.featurexmls} "
        "-out {output.featurexmls} "
        "-threads {threads} "
        "{params.mz_max_difference} "
        "{params.mz_unit} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

#FeatureLinkerUnlabeledQT:
rule link_ffc_maps:
    input:
        featurexmls = expand("work/{{dbsearchdir}}/{sample}/ffc_filt_idmap_align_{sample}.featureXML",sample=SAMPLES)
    output:
        featurexml = temp("work/{dbsearchdir}/ffc_flq.consensusXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        'work/{dbsearchdir}/ffc_flq.benchmark.txt'
    params:
        mz_max_difference = "-algorithm:distance_MZ:max_difference 20",
        mz_unit = "-algorithm:distance_MZ:unit ppm",
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffc_flq.log'
    shell:
        "FeatureLinkerUnlabeledQT "
        "-in {input.featurexmls} "
        "-out {output.featurexml} "
        "-threads {threads} "
        "{params.mz_max_difference} "
        "{params.mz_unit} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

#IDConflictResolverRule
rule resolve_ffc_map_conflicts:
    input:
        featurexml = "work/{dbsearchdir}/ffc_flq.consensusXML"
    output:
        featurexml = temp("work/{dbsearchdir}/ffc_flq_idcr.consensusXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        'work/{dbsearchdir}/ffc_flq_idcr.benchmark.txt'
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffc_flq_idcr.log'
    shell:
        "IDConflictResolver "
        "-in {input.featurexml} "
        "-out {output.featurexml} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

#ConsensusMapNormalizer:
rule normalize_ffc_maps:
    input:
        featurexml = "work/{dbsearchdir}/ffc_flq_idcr.consensusXML"
    output:
        featurexml = temp("work/{dbsearchdir}/ffc_flq_idcr_cmn.consensusXML")
    singularity:
        config['singularity']['default']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    threads:
        2
    benchmark:
        'work/{dbsearchdir}/ffc_flq_idcr_cmn.benchmark.txt'
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffc_flq_idcr_cmn.log'
    shell:
        "ConsensusMapNormalizer "
        "-in {input.featurexml} "
        "-out {output.featurexml} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# ProteinQuantifier
rule quantify_proteins_ffc_maps:
    input:
        featurexml = "work/{dbsearchdir}/ffc_flq_idcr_cmn.consensusXML",
        fido = "work/{dbsearchdir}/proteinid/fido_fdr_filt.idXML"
    output:
        csv = "csv/ffc_{dbsearchdir}_proteinIntensities.csv"
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        'csv/ffc_{dbsearchdir}_proteinIntensities.txt'
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = "csv/ffc_{dbsearchdir}_proteinIntensities.log"
    shell:
        "ProteinQuantifier "
        "-in {input.featurexml} "
        "-protein_groups {input.fido} "
        "-out {output.csv} "
        "-top 0 "
        "-average sum "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "
