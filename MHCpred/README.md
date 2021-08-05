# MHCpred

MHC I affinity predition


## training data download

```shell
wget http://www.cbs.dtu.dk/services/NetMHCpan/data.tar.gz
wget http://www.cbs.dtu.dk/services/NetMHCpan/test.tar.gz
```

## prepare train, valid, test datasets
We use UniRep to get sequence embeddings

What's [UniRep](https://www.nature.com/articles/s41592-019-0598-1)

Edit the working directory and run:

```
pip install jax-unirep
snakemake -s get_embeds.smk -j 32
```

the `${work_dir}` will have 3 dataset output
- MHCDataset.train.pkl
- MHCDataset.valid.pkl
- MHCDataset.test.pkl


## training

Currently, only `BA` data are used.

```shell

python train.py --data_path ${work_dir}  --batch_size 100
```

## predict
```shell
python eval.py --data_path ${work_dir} --batch_size 100
```