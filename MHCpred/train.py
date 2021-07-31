import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob,sys,os,json
from datetime import datetime
from pandas.core.indexing import check_bool_indexer
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data.dataloader import DataLoader
from torch.utils.tensorboard import SummaryWriter

from datasets import MHCDataset, DataBundle
from model import MHCModel

import joblib
import config

LOG_DIR = config.args.log_dir
os.makedirs(LOG_DIR, exist_ok=True)
logger = SummaryWriter(log_dir = LOG_DIR)

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

# data = DataBundle(alleles= mhc_allel_filename, 
#                   mhc_psudo= pesudo_filename, 
#                   ba_el_dir= peptide_ba_el_dir )
print("Build Model")
model = MHCModel(config.input_size, 1)
# # preprocess data
# data.parse_ba()
# data.parse_el()
# data.concat()
# data.train_val_test_split(seed=1234)
# joblib.dump(data, "/data/bases/fangzq/ImmunoRep/databundle.pkl")
# data = joblib.load(config.DATA_BUNDLE)
# train_data = MHCDataset(data, data.train)
# valid_data = MHCDataset(data, data.val)
# test_data = MHCDataset(data, data.test)

PATH = config.args.data_path
print("Load data files")
train_data = joblib.load(os.path.join(PATH, "MHCDataset.train.pkl"))
valid_data = joblib.load(os.path.join(PATH, "MHCDataset.valid.pkl"))
test_data = joblib.load(os.path.join(PATH, "MHCDataset.test.pkl"))

print("Prepare DataLoader")
train_loader = DataLoader(train_data, batch_size=config.batch_size, num_workers= config.num_workers) #sampler=train_sampler, num_workers=1 )# sampler=SubsetRandomSampler() )
valid_loader =  DataLoader(valid_data, batch_size=config.batch_size, num_workers= config.num_workers)
test_loader =  DataLoader(valid_data, batch_size=config.batch_size, num_workers= config.num_workers)


print("Logging model")
embeds = next(iter(valid_loader))
logger.add_graph(model, [embeds['mhc_embed'], embeds['ag_embed']])
# print("Build Model")
# model = MHCModel(input_size, 1)
# model.to(device)
criterion = torch.nn.MSELoss() 
optimizer = torch.optim.Adam(model.parameters(), lr= config.learning_rate)
# let learning_rate decrease by 50% at 500, 1000 and 2000-th epoch
#scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, [500, 1000, 2000], gamma=0.5)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, verbose=True,
                                                        threshold=1e-4, cooldown=10, min_lr=1e-5)
## load pretrain model
weight = '{LOG_DIR}/model.best.pth'
if os.path.exists(weight):
    print(f"Load pretrain model: {weight}")
    checkpoint = torch.load(weight, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])


print("Start training")
# Training the Model
last_loss = 1000
total_steps = 0
model.to(device)
for epoch in range(config.num_epochs):
    model.train()
    running_loss = 0.0
    for i, embeds in enumerate(tqdm(train_loader, total=len(train_loader))):
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
        # scheduler.step()
        # print statistics
        train_loss = loss.item() # to cpu
        running_loss += train_loss
        if (i + 1) % config.steps_per_validation == 0:
            # validate
            model.eval()
            valid_loss = 0
            with torch.no_grad():
                for embeds in tqdm(valid_loader, total=len(valid_loader)):
                    inp_mhc, inp_ag, targets = embeds['mhc_embed'], embeds['ag_embed'], embeds['target']
                    inp_mhc = inp_mhc.to(device)
                    inp_ag = inp_ag.to(device)
                    targets = targets.to(device)
                    outputs = model(inp_mhc, inp_ag)
                    valid_loss += criterion(outputs, targets) # 
            valid_loss /= len(valid_loader)
            scheduler.step(valid_loss)
            print('epoch [%d], step [%d], lr [%.7f],  loss: %.7f' % (epoch, i, optimizer.param_groups[0]['lr'], valid_loss))
            total_steps = len(train_loader)*epoch + i
            logger.add_scalar('Loss/valid', valid_loss, total_steps) 
            logger.add_scalar('Loss/train', train_loss, total_steps) 
          
            if valid_loss < last_loss:
                last_loss = min(last_loss, valid_loss)
                PATH = f'{LOG_DIR}/model.best.pth'
                torch.save({'model_state_dict': model.state_dict(),         
                            'optimizer_state_dict': optimizer.state_dict()}, PATH)

    running_loss /= len(train_loader)
    # lr = scheduler.get_last_lr()[-1] # optimizer.param_groups[0]['lr']
    print('epoch [%d], lr: [%.7f], loss: %.7f' % (epoch, optimizer.param_groups[0]['lr'] , running_loss))


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

