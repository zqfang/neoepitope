# MSFragger
rule msfragger_search:
    input:
        mzml = "mzml/{datafile}.mzML",
        fasta = 'work/database/target_decoy_database.fasta'
    output:
        idxml = "work/msfragger/{datafile}/dbsearch_{datafile}.idXML"
    singularity:
        config['singularity']['default']
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8192
    benchmark:
        "work/msfragger/{datafile}/dbsearch_{datafile}.benchmark.txt"
    params:
        pmt = "-tolerance:precursor_mass_tolerance {0}".format(config["precursor"]["tolerance"]),
        pmu = "-tolerance:precursor_mass_unit {0}".format(config["precursor"]["units"]),
        ptt = "-tolerance:precursor_true_tolerance {0}".format(config["precursor"]["tolerance"]),
        ptu = "-tolerance:precursor_true_unit {0}".format(config["precursor"]["units"]),
        fmt = "-tolerance:fragment_mass_tolerance  {0}".format(config["fragment"]["tolerance"]),
        fmu = "-tolerance:fragment_mass_unit {0}".format(config["fragment"]["units"]),
        e = "-digest:search_enzyme_name {0}".format(config["digestion"]["enzyme"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/msfragger/{datafile}/dbsearch_{datafile}.log'
    shell:
        "MSFraggerAdapter "
        "-in {input.mzml} -out {output.idxml} -database {input.fasta} "
        "-executable "
        "/usr/local/openms_thirdparty/All/MSFragger/MSFragger.jar "
        "-java_heapmemory 8192 "
        "{params.pmt} "
        "{params.pmu} "
        "{params.ptt} "
        "{params.ptu} "
        "{params.fmt} "
        "{params.fmu} "
        "{params.e} "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "
