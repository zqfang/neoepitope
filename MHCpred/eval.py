import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob,sys,os,json
from datetime import datetime

import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data.dataloader import DataLoader
from .datasets import MHCDataset, DataBundle
from .model import MHCModel
import config


os.makedirs("checkpoints", exist_ok=True)
device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

data = DataBundle(alleles= mhc_allel_filename, 
                  mhc_psudo= pesudo_filename, 
                  ba_el_dir= peptide_ba_el_dir )

data.parse_ba()
data.parse_el()
data.concat()
train, val, test = data.train_val_test_split(seed=1234)
print("Prepare DataLoader ")
train_data = MHCDataset(data, test)
test_loader =  DataLoader(valid_data, batch_size=batch_size, num_workers= num_workers)


print("Build Model")
model = MHCModel(input_size, 1)

criterion = torch.nn.MSELoss() 
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

PATH = "checkpoints/model.best.pth"
checkpoint = torch.load(PATH)
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

model.to(device)
## now test loss
model.eval()
with torch.no_grad():
    test_loss = 0
    for i, embeds  in enumerate(valid_loader):
        inp_mhc, inp_ag, targets = embeds['mhc_embed'], embeds['ag_embed'], embeds['target']
        inp_mhc = inp_mhc.to(device)
        inp_ag = inp_ag.to(device)
        targets = targets.to(device)
        outputs = model(inp_mhc, inp_ag)
        test_loss += criterion(outputs, targets) # 
    test_loss /= len(test_loader)
    print('test loss: %.7f' % test_loss)