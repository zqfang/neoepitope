
import glob, os
import numpy as np
import pandas as pd

import torch 
from torch.utils.data.dataset import Dataset
from jax_unirep import get_reps
from typing import Dict, Sequence, Tuple, Any, List, Union 


class DataBundle:
    def __init__(self, alleles:str , mhc_psudo: str , ba_el_dir: str, *args, **kwargs):
        """
        read training data 
        """
        self.alleles = pd.read_table(alleles, sep=" ", header=None, names=['alias','full'])
        ## mhc sequence
        self.mhc = []
        with open(mhc_psudo, 'r') as allel:
            for line in allel:
                self.mhc.append(line.strip().split())
        self.mhc = pd.DataFrame(self.mhc, columns=['alias','psudo_seq'])
        self.mhc.drop_duplicates(inplace=True)
        
        # get embeding
        self.pseudo2embed = self._pesudo2embed()
        
        # ba, el paths
        self.ba_path = sorted(glob.glob(os.path.join(ba_el_dir, "*.ba")))
        self.el_path = sorted(glob.glob(os.path.join(ba_el_dir, "*.el")))
        ##  

    def _pesudo2embed(self) -> Dict[str, np.ndarray]: 
        sequences = self.mhc['psudo_seq'].to_list()
        self.pesudo2idx = {seq: i for i, seq in enumerate(sequences)}
        h_avg, h_final, c_final= get_reps(sequences)
        ## 
        self.pseudo_embed = h_avg
        return {seq: arr for seq, arr in zip(sequences, h_avg)}


    @property
    def name2pesudo(self) -> Dict[str, str]:
        return {row['alias']: row['psudo_seq'] for i, row in self.mhc.iterrows()}
    
    @property
    def pesudo2name(self) -> Dict[str, str]:
        return { row['psudo_seq']: row['alias'] for i, row in self.mhc.iterrows()}


    @property
    def name2embed(self) -> Dict[str, np.ndarray]:
        return {self.name2pesudo[name]: self.pesudo2embed[seq] for name, seq in self.name2pesudo.items()} 

    def parse_ba(self):
        self.pseudo2ba = {}
        for p in self.ba_path:
            seq = os.path.basename(p).split(".")[0]
            if seq not in self.pesudo2idx: continue
            ba = pd.read_table(p, header=None, sep=" ", names= ['unknown','epitope','affinity'])
            ## # each of the arrays will be of shape (len(sequences), 1900),
            ## require: pip install jax_unirep
            # h_avg, h_final, c_final= jax_unirep.get_reps(ba['epitope'].to_list())
            h_avg = np.load(f"{p}.unirep.npy")
            ba['pseudo_idx'] = self.pesudo2idx[seq]
            self.pseudo2ba[seq] = [ba, h_avg]

    def parse_el(self):
        self.pseudo2el = {}
        for p in self.el_path:
            seq = os.path.basename(p).split(".")[0]
            if seq not in self.pesudo2idx: continue
            el = pd.read_table(p, header=None, sep=" ", names= ['unknown','epitope','affinity'])
            ## # each of the arrays will be of shape (len(sequences), 1900),
            ## require: pip install jax_unirep
            # h_avg, h_final, c_final= get_reps(el['epitope'].to_list())
            h_avg = np.load(f"{p}.unirep.npy")
            el['pseudo_idx'] = self.pesudo2idx[seq]
            self.pseudo2el[seq] = [el, h_avg] #, h_final, c_final] 
    
    def concat(self):
        self.ba_affi = pd.concat([v[0] for k, v in self.pseudo2ba.items()], ignore_index=True)
        self.ba_embed = np.vstack([v[1] for k, v in self.pseudo2ba.items()])
        self.el_affi = pd.concat([v[0] for k, v in self.pseudo2el.items()], ignore_index=True)
        self.el_embed = np.vstack([v[1] for k, v in self.pseudo2el.items()])

    def train_val_test_split(self, ratio = [0.8, 0.1, 0.1], seed=1234):
        rs = np.random.RandomState(seed=seed)
        indices = list(range(len(self.ba_embed)))
        rs.shuffle(indices)
        ## split
        train, val, test = np.split(indices,
                                    [int(ratio[0]*len(indices)), int((ratio[0] + ratio[1])*len(indices))] )
        self.train = train
        self.val = val
        self.test = test
        return train, val, test

    def get_ba_data(self, idx):
        return self.ba_affi.iloc[idx], self.ba_embed[idx]
        
    def get_el_data(self, idx):
        return self.el_affi.iloc[idx], self.el_embed[idx]
      

class MHCDataset(Dataset):
    def __init__(self, data: DataBundle, partition: Sequence ):
        """
        inputs: path to the training dataset
        partition:  indices of train, val, test.
        """
        # ag: antigen
        ba_affi, ba_embed = data.get_ba_data(partition)
        self.ag_embed = torch.from_numpy(ba_embed.astype(np.float32))
        self.affinity = torch.from_numpy(ba_affi['affinity'].values.astype(np.float32))
        self.mhc_embed = np.vstack([data.pseudo_embed[idx] for idx in ba_affi['pseudo_idx'].to_list()])
        self.mhc_embed = torch.from_numpy(self.mhc_embed.astype(np.float32))

    def __getitem__(self, idx):
        ## multi-class classification need longTensor for y
        # self.targets = self.targets.type(torch.LongTensor)
        return {'mhc_embed': self.mhc_embed[idx], 'ag_embed': self.ag_embed[idx], 'target': self.affinity[idx]}
    def __len__(self):
        return len(self.affinity)

    def _pesudo2embed(self) -> Dict[str, np.ndarray]: 
        sequences = self.mhc['psudo_seq'].to_list()
        self.pesudo2idx = {seq: i for i, seq in enumerate(sequences)}
        h_avg, h_final, c_final= get_reps(sequences)



class MHCEvalDataset(Dataset):
    def __init__(self, data: Union[pd.DataFrame, str], mhc2pesudo: Union[pd.DataFrame, str]):
        if not isinstance(mhc2pesudo, pd.DataFrame):
            mhc2pesudo = pd.read_csv(mhc2pesudo)
            # 3 column df
            # hla, hla_alias, pesudo_seq
        self.mhc2pesudo_dict = {row['hla']:row['pseudo_seq'] for _, row in mhc2pesudo.iterrows()}
        
        if not isinstance(data, pd.DataFrame):
            data = pd.read_table(data, sep=" ", header=None, 
                                 names=['species','hla','seq_len','antigen','target'])
        data['pseudo'] = data['hla'].map(self.mhc2pesudo_dict)
        data = data.dropna()
        self.data = data
        # get embeds
        ag_embed = self._get_embed(data['antigen'].to_list())
        self.ag_embed = torch.from_numpy(ag_embed.astype(np.float32))
        mhc_embed = self._get_embed(data['pseudo'].to_list())
        self.mhc_embed = torch.from_numpy(mhc_embed.astype(np.float32))
        self.target = torch.as_tensor(data['target'].values, dtype=torch.long)

    def __getitem__(self, idx):
        # ag: antigen
        return {'mhc_embed': self.mhc_embed[idx], 'ag_embed': self.ag_embed[idx], 'target': self.target[idx]}

    def __len__(self):
        return len(self.target)

    def _get_embed(self, aa: List[str]) -> np.ndarray: 
        h_avg, h_final, c_final= get_reps(aa)
        return h_avg
