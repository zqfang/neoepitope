import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import spearmanr
import seaborn as sns
import glob,sys,os,json
from datetime import datetime
import joblib
import scipy.stats as stat
import torch
import torch.nn as nn
import torch.nn.functional as F

from tqdm import tqdm
from torch.utils.data.dataloader import DataLoader
from datasets import MHCDataset, DataBundle
from model import MHCModel
import config


device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

# data = DataBundle(alleles= mhc_allel_filename, 
#                   mhc_psudo= pesudo_filename, 
#                   ba_el_dir= peptide_ba_el_dir )

# data.parse_ba()
# data.parse_el()
# data.concat()
# train, val, test = data.train_val_test_split(seed=1234)
# print("Prepare DataLoader ")
# test_data = MHCDataset(data, test)

PATH = config.args.data_path
print("Load data files")
test_data = joblib.load(os.path.join(PATH, "MHCDataset.test.pkl"))
test_loader =  DataLoader(test_data, batch_size=config.batch_size, num_workers= config.num_workers)


print("Load Model")
model = MHCModel(config.input_size, 1)

criterion = torch.nn.MSELoss() 
optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)

checkpoint = torch.load(os.path.join(config.args.log_dir, "model.best.pth"), map_location=device)
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

model.to(device)
## now test loss
model.eval()
preds = []
y = []
spearman = []
pearson = []
spearman_pval = []
pearson_pval = []

test_loss = 0
model.eval()
with torch.no_grad():
    for i, embeds  in enumerate(test_loader):
        inp_mhc, inp_ag, targets = embeds['mhc_embed'], embeds['ag_embed'], embeds['target']
        inp_mhc = inp_mhc.to(device)
        inp_ag = inp_ag.to(device)
        y.append(targets.numpy())
        targets = targets.to(device)
        outputs = model(inp_mhc, inp_ag)
        preds.append(outputs.detach().cpu().numpy())
        loss = criterion(outputs, targets).item() # 
        test_loss += loss
        sr, pval_sr = stat.spearmanr(preds[-1], y[-1])
        pr, pval_pr = stat.pearsonr(preds[-1], y[-1])
        spearman.append(sr)
        pearson.append(pr)
        spearman_pval.append(pval_sr)
        pearson_pval.append(pval_pr)
        print(f"batch {i}, test loss: {loss}, speraman's r: {sr}, pearson's r: {pr}")

test_loss /= len(test_loader)

print('averge test loss: %.7f' % test_loss)

preds = np.stack(preds)
y = np.stack(y)
avg_sp, spval = stat.spearmanr(preds, y)
avg_pr, ppval = stat.pearsonr(preds, y)
print("averge Pearson's r: %.7f, pval: %.7f "% (avg_pr, ppval))
print("averge Spearman's r: %.7f, pval: %.7f "% (avg_sp, spval)

# plot the scatter
fig, ax = plt.subplots(figsize=(4,4))
ax.scatter(preds, y)
ax.set_xlabel("Predit")
ax.set_ylabel("True")
ax.set_title("Pred-True scatterplot")
fig.savefig("test.scatter.png", dpi=300)


