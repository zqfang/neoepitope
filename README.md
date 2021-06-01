

# Neoepitope
Deep learning models (Pytorch) for cancer immunology 

1. Neoantigen discovery using immunopeptidome: [PointNovo](https://github.com/volpato30/PointNovo)
2. Antibody drug discovery:
   - tumor classification using tcr/bcr (inspired by [DeepBCR](https://bitbucket.org/liulab/deepbcr) ): TODO
   - other

3. binding affinity model (inspired by [tcellmatch](https://github.com/theislab/tcellmatch)): TODO

## Dependency
- python >= 3.7
- pytorch >= 1.0
- scikit-learn
- biopython
- pyteomics
- cython
- snakemake

For protein database search you also need to install [percolator](http://percolator.ms/).

## Neoantigen discovery



### Data Preprocessing
e.g. Tmuor HLA peptides (MSV000082648)
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

### Train and predict

1. train, valid, test
```shell
cd PointNovo
# build cython extension
python deepnovo_cython_setup.py build_ext --inplace

# set correct data path in config.py first, then
# train
python main.py --train
```

## Tumor classificatiton

TODO

## Binding affinity

TODO


## QA
## Contacts