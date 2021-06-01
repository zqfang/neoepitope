

# Neoepitope

Deep learning models for immunopeptidome
- neoantigen discovery: PointNovo
- epitope classification: TCR
- binding affinity model: TODO

## Dependency
- python >= 3.7
- pytorch >= 1.0
- biopython
- pyteomics
- cython

For database search you also need to install [percolator](http://percolator.ms/).

## Data Preprocessing

1. download peaks and percolator data from `MSV000082648`
```
   - |- MSV000082648
       |- peaks
       |- percolator
```

2. run `prepare_pointnovo.smk` to generate train, test, valid data
  - makesure these two path is correct

```
smkpath = "/data/bases/fangzq/ImmunoRep/neoepitope"
DATA = "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
```
run
```
snakemake -s prepare_pointnovo.smk -j 32 -p 
```

## Neoantigen discovery

1. train, valid, test
```shell
cd PointNovo
# build cython extension
python deepnovo_cython_setup.py build_ext --inplace
# train
python main.py --train
```

refer to PointNovo

## Epitope classificatiton

refer to DeepBCR

## Binding affinity model

TODO


