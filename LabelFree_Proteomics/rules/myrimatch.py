# myrimatch:
rule myrimatch_search:
    input:
        mzml = "mzml/{datafile}.mzML",
        fasta = 'work/database/target_decoy_database.fasta',
    output:
        idxml = "work/myrimatch/{datafile}/dbsearch_{datafile}.idXML"
    singularity:
        config['singularity']['default']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 32000
    threads: 4
    benchmark:
        "work/myrimatch/{datafile}/dbsearch_{datafile}.benchmark.txt"
    params:
        pmt = "-precursor_mass_tolerance {0}".format(config["precursor"]["tolerance"]),
        pmu = "-precursor_mass_tolerance_unit {0}".format(config["precursor"]["units"]),
        fmt = "-fragment_mass_tolerance {0}".format(config["fragment"]["tolerance"]),
        fmu = "-fragment_mass_tolerance_unit {0}".format(config["fragment"]["units"]),
        e = "-CleavageRules {0}".format(config["digestion"]["enzyme"]),
        mc = "-MaxMissedCleavages {0}".format(config["digestion"]["missed_cleavages"]),
        fm = "-fixed_modifications {0}".format(config["modifications"]["fixed"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/myrimatch/{datafile}/dbsearch_{datafile}.log'
    shell:
        "MyriMatchAdapter "
        "-in {input.mzml} "
        "-out {output.idxml} "
        "-database {input.fasta} "
        "-myrimatch_executable myrimatch "
        "{params.pmt} "
        "{params.pmu} "
        "{params.fmt} "
        "{params.fmu} "
        "{params.e} "
        "{params.mc} "
        "{params.fm} "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "
