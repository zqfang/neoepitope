

# Neoepitope
Deep learning models (Pytorch) for cancer immunology 

1. Neoantigen discovery based on immunopeptidome: 
    - Database Search: Comet + Percolator workflow -> [MHCquant](https://github.com/Leon-Bichmann/MHCquant)
    - De novo search: [PointNovo](https://github.com/volpato30/PointNovo)
2. Antibody drug discovery:
   - tumor classification using tcr/bcr (inspired by [DeepBCR](https://bitbucket.org/liulab/deepbcr) ): TODO
   - other

3. binding affinity model (inspired by [tcellmatch](https://github.com/theislab/tcellmatch)): TODO

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
   - pytorch >= 1.0
   - scikit-learn
   - python-Levenshtein
   - biopython
   - pyteomics
   - pyopenms
   - cython


## Neoantigen discovery workflow


### 1. Database search from the very begining 
1. comet + percolator:
* MHCquant (recommended) 
* on going pepeline: TODO
   ```shell
   snakemake -s 1.AA.db.search.smk --configfile config.yaml -p -j 12
   ```

2. reformat the percolator output 
  - TODO

### 1.1 Start from the published data
e.g. Tmuor HLA peptides (MSV000082648)
1. download peaks and percolator data from `MSV000082648`
```
   - |- MSV000082648
       |- peaks
       |- percolator
```

2. run `2.AA.preprocess.smk` to generate train, test, valid data
  - makesure these two paths are correct

```yaml
SMKPATH: "/data/bases/fangzq/ImmunoRep/neoepitope"
WORK: "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
```
run
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

```shell
snakemake -s 4.AA.prioritize.smk --configfile config.yaml -p -j 8
```

## Tumor classificatiton

TODO

## Binding affinity

TODO


## QA
## Contacts