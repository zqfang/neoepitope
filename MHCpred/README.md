# MHCpred

MHC I affinity predition


## training data download

```shell
wget http://www.cbs.dtu.dk/services/NetMHCpan/data.tar.gz
wget http://www.cbs.dtu.dk/services/NetMHCpan/test.tar.gz
```

## get UniRep embedding
What's [UniRep](https://www.nature.com/articles/s41592-019-0598-1)
```
pip install jax-unirep
snakemake -s get_embeds.smk -j 32
```

## training
```shell
python train.py 
```

## predit
```shell
python eval.py
```