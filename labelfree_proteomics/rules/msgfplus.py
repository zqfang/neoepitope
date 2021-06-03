# MSGFPlusIndexDB:
rule msgfplus_db_index:
    input:
        fasta = 'work/database/target_decoy_database.fasta'
    output:
        index = 'work/database/target_decoy_database.canno'
    singularity:
        config['singularity']['default']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/database/MSGFPlusIndexDB.benchmark.txt"
    params:
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/database/MSGFPlusIndexDB.log'
    shell:
        "java -Xmx3500M -cp /usr/local/openms_thirdparty/All/MSGFPlus/MSGFPlus.jar "
        "edu.ucsd.msjava.msdbsearch.BuildSA "
        "-d {input.fasta} "
        "-tda 0 "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "

# MSGFPlus
rule msgfplus_search:
    input:
        mzml = "mzml/{datafile}.mzML",
        fasta = 'work/database/target_decoy_database.fasta',
        index = 'work/database/target_decoy_database.canno'
    output:
        idxml = "work/msgfplus/{datafile}/dbsearch_{datafile}.idXML"
    singularity:
        config['singularity']['default']
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    benchmark:
        "work/msgfplus/{datafile}/dbsearch_{datafile}.benchmark.txt"
    params:
        pmt = "-precursor_mass_tolerance {0}".format(config["precursor"]["tolerance"]),
        pmu = "-precursor_error_units {0}".format(config["precursor"]["units"]),
        e = "-enzyme {0}".format(config["digestion"]["enzyme"]),
        fm = "-fixed_modifications {0}".format(config["modifications"]["fixed"]),
        debug = '-debug {0}'.format(config["debug"]),
        log = 'work/msgfplus/{datafile}/dbsearch_{datafile}.log'
    shell:
        "MSGFPlusAdapter "
        "-in {input.mzml} -out {output.idxml} -database {input.fasta} "
        "-executable "
        "/usr/local/openms_thirdparty/All/MSGFPlus/MSGFPlus.jar "
        "-java_memory 8192 "
        "-isotope_error_range 0,0 "
        "-instrument Q_Exactive "
        "-max_precursor_charge 4 "
        "-add_features "
        "-max_mods 3 "
        "{params.pmt} "
        "{params.pmu} "
        "{params.e} "
        "{params.fm} "
        "-threads {threads} "
        "{params.debug} "
        "2>&1 | tee {params.log} "
