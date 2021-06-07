

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
   - biopython
   - pyteomics
   - pyopenms
   - cython


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
  - makesure these two paths are correct

```
smkpath = "/data/bases/fangzq/ImmunoRep/neoepitope"
DATA = "/data/bases/fangzq/ImmunoRep/data/MSV000082648"
```
run
```
snakemake -s prepare_pointnovo.smk -j 32 -p 
```

### Train and predict

1. train, valid
```shell
cd PointNovo
# build cython extension
python deepnovo_cython_setup.py build_ext --inplace

# set correct data path in config.py first, then
# train
python main.py --train
```

2. test
```shell
# denovo search first, then test
python main.py --search_denovo
python main.py --test
```

## Tumor classificatiton

TODO

## Binding affinity

TODO


## QA
## Contacts