# XTandem:
rule xtandem:
    input:
        mzml = "mzml/{datafile}.mzML",
        fasta = 'work/database/target_decoy_database.fasta',
    output:
        idxml = "work/xtandem/{datafile}/dbsearch_{datafile}.idXML"
    singularity:
        config['singularity']['default']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    threads: 4
    benchmark:
        "work/xtandem/{datafile}/dbsearch_{datafile}.benchmark.txt"
    params:
        pmt = "-precursor_mass_tolerance {0}".format(config["precursor"]["tolerance"]),
        pmu = "-precursor_error_units {0}".format(config["precursor"]["units"]),
        fmt = "-fragment_mass_tolerance {0}".format(config["fragment"]["tolerance"]),
        fmu = "-fragment_error_units {0}".format(config["fragment"]["units"]),
        e = "-enzyme {0}".format(config["digestion"]["enzyme"]),
        mc = "-missed_cleavages {0}".format(config["digestion"]["missed_cleavages"]),
        fm = "-fixed_modifications {0}".format(config["modifications"]["fixed"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/xtandem/{datafile}/dbsearch_{datafile}.log'
    shell:
        "XTandemAdapter "
        "-in {input.mzml} "
        "-out {output.idxml} "
        "-database {input.fasta} "
        "-xtandem_executable tandem.exe "
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
