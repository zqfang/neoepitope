import errno, glob, os, os.path, sys

shell.prefix('')
configfile: "config.yaml"

SAMPLES = []
MZMLFILES = glob.glob("{0}/*".format(config["mzml"]["location"]))
allowed_exts = [".mzml",".mzML","MZML"]
for rawfile in MZMLFILES:
    rbase = os.path.basename(rawfile)
    rbase,rext = os.path.splitext(rbase)
    if rext in allowed_exts:
        SAMPLES.append(rbase)

DBFILES = glob.glob("{0}/*".format(config["fasta"]["location"]))
DATABASES = []
allowed_exts = [".fasta",".FASTA"]
for dbfile in DBFILES:
    dbase = os.path.basename(dbfile)
    dbase,dext = os.path.splitext(dbase)
    if dext in allowed_exts:
        DATABASES.append(dbase)

# Setup Targets for Pipeline
rule targets:
    input:
        expand("csv/sc_{dbsearch}_proteinCounts.csv", dbsearch=config["search_engines"]),
        expand("csv/ffm_{dbsearch}_proteinIntensities.csv", dbsearch=config["search_engines"]),
        expand("csv/ffc_{dbsearch}_proteinIntensities.csv", dbsearch=config["search_engines"]),
        expand("csv/ffidi_{dbsearch}_proteinIntensities.csv", dbsearch=config["search_engines"]),

# Construct Concatenated Databases With Decoys
include: "rules/decoy_database.py"

# Remove Compression
include: "rules/file_conversion.py"

# Perform Database Search
for search in config["search_engines"]:
    include: "rules/%s.py" % search

# Perfrom Post Processing of Database Searches
include: "rules/search_post_processing.py"

# Perform Protein Inference
include: "rules/fido.py"

# Perfrom Feature Finding
include: "rules/feature_finder_identification_internal.py"

# Perfrom Feature Finding
include: "rules/feature_finder_multiplex.py"

# Perfrom Feature Finding
include: "rules/feature_finder_centroid.py"

# Determine Spectral Counts
include: "rules/spectral_counts.py"

# Finalize script
include: "rules/onsuccess.py"

