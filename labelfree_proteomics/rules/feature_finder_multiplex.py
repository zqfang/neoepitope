# FeatureFinderMultiplex
rule find_features_multiplex:
    input:
        mzml = "mzml/{datafile}.mzML",
    output:
        featurexml = "work/featurefindermultiplex/{datafile}/multiplex_{datafile}.featureXML",
        consensusxml = "work/featurefindermultiplex/{datafile}/multiplex_{datafile}.consensusxml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
	    "work/featurefindermultiplex/{datafile}/multiplex_{datafile}.benchmark.txt"
    params:
        rt_typical = "-algorithm:rt_typical {0}".format(config["lc"]["peak_width"]),
        mz_tolerance = "-algorithm:mz_tolerance {0}".format(config["precursor"]["tolerance"]),
        mz_unit = "-algorithm:mz_unit {0}".format(config["precursor"]["units"]),
        intensity_cutoff = "-algorithm:intensity_cutoff {0}".format(config["precursor"]["threshold"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/featurefindermultiplex/{datafile}/multiplex_{datafile}.log'
    singularity:
        config['singularity']['default']
    shell:
        "FeatureFinderMultiplex "
        "-in {input.mzml} -out {output.consensusxml} "
        "-out_features {output.featurexml} "
	    "{params.rt_typical} "
	    "{params.mz_tolerance} "
	    "{params.mz_unit} "
	    "{params.intensity_cutoff} "
	    "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "


# IDFilter
rule filter_peptides_ffm:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/fdr_{datafile}.idXML"
    output:
        idxml = temp("work/{dbsearchdir}/{datafile}/ffm_filt_{datafile}.idXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
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

# IDMapper
rule map_ffm_features:
    input:
        idxml = "work/{dbsearchdir}/{datafile}/ffm_filt_{datafile}.idXML",
        featurexml = "work/featurefindermultiplex/{datafile}/multiplex_{datafile}.featureXML"
    output:
        featurexml = temp("work/{dbsearchdir}/{datafile}/ffm_filt_idmap_{datafile}.featureXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/{datafile}/ffm_filt_idmap_{datafile}.benchmark.txt"
    params:
        centroid_rt = "-feature:use_centroid_rt ",
        centroid_mz = "-feature:use_centroid_mz ",
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/{datafile}/ffm_filt_idmap_{datafile}.log'
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

# MapAlignerPoseClustering:
rule align_ffm_maps:
    input:
        featurexmls = expand("work/{{dbsearchdir}}/{sample}/ffm_filt_idmap_{sample}.featureXML",sample=SAMPLES)
    output:
        featurexmls =
        temp(expand("work/{{dbsearchdir}}/{sample}/ffm_filt_idmap_align_{sample}.featureXML",sample=SAMPLES))
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/ffm_filt_idmap_align.benchmark.txt"
    params:
        mz_max_difference = "-algorithm:pairfinder:distance_MZ:max_difference 20",
        mz_unit = "-algorithm:pairfinder:distance_MZ:unit ppm",
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffm_filt_idmap_align.log'
    shell:
        "MapAlignerPoseClustering "
        "-in {input.featurexmls} "
        "-out {output.featurexmls} "
        "-threads {threads} "
        "{params.mz_max_difference} "
        "{params.mz_unit} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# FeatureLinkerUnlabeledQT:
rule link_ffm_maps:
    input:
        featurexmls = expand("work/{{dbsearchdir}}/{sample}/ffm_filt_idmap_align_{sample}.featureXML",sample=SAMPLES)
    output:
        featurexml = temp("work/{dbsearchdir}/ffm_flq.consensusXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/ffm_flq.benchmark.txt"
    params:
        mz_max_difference = "-algorithm:distance_MZ:max_difference 20",
        mz_unit = "-algorithm:distance_MZ:unit ppm",
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffm_flq.log'
    shell:
        "FeatureLinkerUnlabeledQT "
        "-in {input.featurexmls} "
        "-out {output.featurexml} "
        "-threads {threads} "
        "{params.mz_max_difference} "
        "{params.mz_unit} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# IDConflictResolverRule
rule resolve_ffm_map_conflicts:
    input:
        featurexml = "work/{dbsearchdir}/ffm_flq.consensusXML"
    output:
        featurexml = temp("work/{dbsearchdir}/ffm_flq_idcr.consensusXML")
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/ffm_flq_idcr.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffm_flq_idcr.log'
    shell:
        "IDConflictResolver "
        "-in {input.featurexml} "
        "-out {output.featurexml} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# ConsensusMapNormalizer:
rule normalize_ffm_maps:
    input:
        featurexml = "work/{dbsearchdir}/ffm_flq_idcr.consensusXML"
    output:
        featurexml = temp("work/{dbsearchdir}/ffm_flq_idcr_cmn.consensusXML")
    singularity:
        config['singularity']['default']
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/{dbsearchdir}/ffm_flq_idcr_cmn.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/{dbsearchdir}/ffm_flq_idcr_cmn.log'
    shell:
        "ConsensusMapNormalizer "
        "-in {input.featurexml} "
        "-out {output.featurexml} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# ProteinQuantifier
rule quantify_proteins_ffm_maps:
    input:
        featurexml = "work/{dbsearchdir}/ffm_flq_idcr_cmn.consensusXML",
        fido = "work/{dbsearchdir}/proteinid/fido_fdr_filt.idXML"
    output:
        csv = "csv/ffm_{dbsearchdir}_proteinIntensities.csv"
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "csv/ffm_{dbsearchdir}_proteinIntensities.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = "csv/ffm_{dbsearchdir}_proteinIntensities.log"
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
