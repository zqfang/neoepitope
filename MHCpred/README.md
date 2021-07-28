# MHCpred

MHC I affinity predition


## training data download

```shell
wget http://www.cbs.dtu.dk/services/NetMHCpan/data.tar.gz
wget http://www.cbs.dtu.dk/services/NetMHCpan/test.tar.gz
```

## get uniRep embeding
```
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