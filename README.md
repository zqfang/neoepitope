

# Neoepitope
Deep learning models (Pytorch) for cancer immunology 

1. Neoantigen discovery based on immunopeptidome: 
    - Database Search: Comet + Percolator workflow -> [MHCquant](https://github.com/Leon-Bichmann/MHCquant)
    - De novo search: [PointNovo](https://github.com/volpato30/PointNovo)
2. Antibody drug discovery:
   - tumor classification using tcr/bcr (inspired by [DeepBCR](https://bitbucket.org/liulab/deepbcr) ): TODO
   - other

3. binding affinity model 
   - 5' scRNA-seq + Tetmer: [tcellmatch](https://github.com/theislab/tcellmatch))
   - MHCII binding: [BERTMHC](https://github.com/s6juncheng/BERTMHC)
   - MHCI and II: [MHCAttnNet](https://github.com/gopuvenkat/MHCAttnNet)

## Dependency
1. proteomics
   - openms
   - pyopenms
   - biopython
   - comet
   - percolator
   - snakemake

2. de novo sequencing 
   - python >= 3.7
   - pytorch >= 1.5
   - scikit-learn
   - python-Levenshtein
   - biopython
   - pyteomics
   - pyopenms
   - cython


## Neoantigen discovery workflow


### 0. Prepare input files

```yaml
SMKPATH: "/data/bases/fangzq/ImmunoRep/neoepitope"
WORKDIR: "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
## mzML files in peaks folder:  MSV000082648/peaks/{sample}.mzML.gz

dbsearch:
  ## proteome database, fasta file 
  prot_db: "/data/bases/fangzq/ImmunoRep/data/PA_ucsc_proteome.fasta"
  contaminant: "/data/bases/fangzq/ImmunoRep/data/UniProtContams_259_20170206.fasta"
  # prot_db + contaminant in one file
  prot_plus_contams: "/data/bases/fangzq/ImmunoRep/data/human.proteome.plus.contams.fasta"
  # prot_db + contaminant + tumor mutation in one file
  full_db: "/data/bases/fangzq/ImmunoRep/data/PA_ucsc_proteome_259contams_viruses_26tumorShared.fasta"
```

### 1. Database search from the very begining 
1. comet + percolator:
* MHCquant (recommended) 
* simplify pipeline of MHCquant
   ```shell
   snakemake -s 1.AA.db.search.smk --configfile config.yaml -p -j 12
   ```

2. run `2.AA.preprocess.smk` to generate train, test, valid data
  - makesure `SMKPATH` is correct

```
snakemake -s 2.AA.preprocess.smk --configfile config.yaml -j 32 -p 
```

### 2. Train and predict

1. complie cython files
```shell
cd PointNovo
# build cython extension
python deepnovo_cython_setup.py build_ext --inplace
```


2. train, test, denovo search + postprocessing
```shell
cd ../
## config gpu resource if needed
snakemake -s 3.AA.train.search.smk --configfile config.yaml -p -j 8
```

### 3. Neoantigen selection
perform second round db search (comet + percolator).  

TODO: use step.1 pipeline to do second round db search

```shell
snakemake -s 4.AA.prioritize.smk --configfile config.yaml -p -j 8
```

## Tumor classificatiton

TODO

## Binding affinity

TODO


## QA
## Contacts