
import numpy as np
import pandas as pd


class TCRDataset(Dataset):
    def __init__(self, inputs, targets):
        inputs = inputs.astype(np.float32)
        targets = targets.astype(np.float32)
        self.inputs = torch.from_numpy(inputs)
        self.targets = torch.from_numpy(targets)
    def __getitem__(self, idx):
        ## multi-class classification need longTensor for y
        self.targets = self.targets.type(torch.LongTensor)
        return {'embed': self.inputs[idx], 'target': self.targets[idx]}
    def __len__(self):
        return len(self.targets)




