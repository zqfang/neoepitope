import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob,sys,os,json
from datetime import datetime
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data.dataloader import DataLoader
from torch.utils.tensorboard import SummaryWriter

from datasets import MHCDataset, DataBundle
from model import MHCModel
import joblib
from config import *


os.makedirs("checkpoints", exist_ok=True)
logger = SummaryWriter(log_dir = "checkpoints")

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

data = DataBundle(alleles= mhc_allel_filename, 
                  mhc_psudo= pesudo_filename, 
                  ba_el_dir= peptide_ba_el_dir )

# preprocess data
data.parse_ba()
data.parse_el()
data.concat()
data.train_val_test_split(seed=1234)
# joblib.dump("/data/bases/fangzq/ImmunoRep/databundle.pkl",filename = data)
# data = joblib.load("/data/bases/fangzq/ImmunoRep/databundle.pkl")
train_data = MHCDataset(data, data.train)
valid_data = MHCDataset(data, data.val)
# test_data = MHCDataset(data, test)

print("Prepare DataLoader")
train_loader = DataLoader(train_data, batch_size=batch_size, num_workers= num_workers) #sampler=train_sampler, num_workers=1 )# sampler=SubsetRandomSampler() )
valid_loader =  DataLoader(valid_data, batch_size=batch_size, num_workers= num_workers)


print("Build Model")
model = MHCModel(input_size, 1)
model.to(device)
criterion = torch.nn.MSELoss() 
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
# let learning_rate decrease by 50% at 500, 1000 and 2000-th epoch
scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, [500, 1000, 2000], gamma=0.5)


print("Start training")
os.makedirs("checkpoints", exist_ok=True)
# Training the Model
last_loss = 1000
for epoch in range(num_epochs):
    model.train()
    running_loss = 0.0
    for embeds in tqdm(train_loader, total=len(train_loader)):
        inp_mhc, inp_ag, targets = embeds['mhc_embed'], embeds['ag_embed'], embeds['target']
        inp_mhc = inp_mhc.to(device)
        inp_ag = inp_ag.to(device)
        targets = targets.to(device)

        ## 
        optimizer.zero_grad()
        outputs = model(inp_mhc, inp_ag)
        loss = criterion(outputs, targets) # 
        loss.backward()
        optimizer.step()
        scheduler.step()

        # print statistics
        running_loss += loss.item()   # detach from tensor

    running_loss /= len(train_loader)
    print('epoch [%d] loss: %.7f' % (epoch, running_loss))
    logger.add_scalar('Loss/train',running_loss, epoch) 
    PATH = f'checkpoints/model.best.pth'
    if running_loss < last_loss:
        last_loss = min(last_loss, running_loss)
        torch.save({'model_state_dict': model.state_dict(),         
                    'optimizer_state_dict': optimizer.state_dict()}, PATH)

    # validate
    model.eval()
    with torch.no_grad():
        valid_loss = 0
        for embeds in tqdm(valid_loader, total=len(valid_loader)):
            inp_mhc, inp_ag, targets = embeds['mhc_embed'], embeds['ag_embed'], embeds['target']
            inp_mhc = inp_mhc.to(device)
            inp_ag = inp_ag.to(device)
            targets = targets.to(device)
            outputs = model(inp_mhc, inp_ag)
            valid_loss += criterion(outputs, targets) # 
        valid_loss /= len(valid_loader)
        print('epoch [%d] loss: %.7f' % (epoch, valid_loss))
        logger.add_scalar('Loss/valid',valid_loss, epoch) 

## now test loss
model.eval()
with torch.no_grad():
    test_loss = 0
    for embeds in tqdm(test_loader, total=len(test_loader)):
        inp_mhc, inp_ag, targets = embeds['mhc_embed'], embeds['ag_embed'], embeds['target']
        inp_mhc = inp_mhc.to(device)
        inp_ag = inp_ag.to(device)
        targets = targets.to(device)
        outputs = model(inp_mhc, inp_ag)
        test_loss += criterion(outputs, targets) # 
    test_loss /= len(test_loader)
    print('final test loss: %.7f' % test_loss)

