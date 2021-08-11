import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
import glob,sys,os,json
from datetime import datetime
import joblib
import scipy.stats as stat
import torch


from tqdm import tqdm
from torch.utils.data.dataloader import DataLoader
from datasets import MHCDataset, MHCEvalDataset
from model import MHCModel
import config


device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')


print("Load data files")
# PATH = config.args.data_path
# test_data = joblib.load(os.path.join(PATH, "MHCDataset.test.pkl"))
PATH = config.eval_filename
test_data = MHCEvalDataset(data=PATH, mhc2pesudo=config.mhc2psedo_filename)
test_loader =  DataLoader(test_data, batch_size=config.batch_size, num_workers= config.num_workers, shuffle=False)


print("Load Model")
model = MHCModel(config.input_size, 1)

criterion = torch.nn.MSELoss() 
optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)

checkpoint = torch.load(os.path.join(config.args.log_dir, "model.best.pth"), map_location=device)
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

model.to(device)
## prediction
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
        loss = criterion(outputs, targets).item() # 
        outputs = outputs.detach().cpu().numpy()
        # convert output to binary value here
        #
        preds.append(outputs)
        
        test_loss += loss
        sr, pval_sr = stat.spearmanr(preds[-1], y[-1])
        pr, pval_pr = stat.pearsonr(preds[-1], y[-1])
        spearman.append(sr)
        pearson.append(pr)
        spearman_pval.append(pval_sr)
        pearson_pval.append(pval_pr)
        print(f"batch {i}, test loss: {loss:.7f}, speraman's r: {sr:.7f}, pearson's r: {pr:.7f}")

test_loss /= len(test_loader)

print('averge test loss: %.7f' % test_loss)

preds = np.concatenate(preds)
y = np.concatenate(y)
avg_sp, spval = stat.spearmanr(preds, y)
avg_pr, ppval = stat.pearsonr(preds, y)
print("averge Pearson's r: %.7f, pval: %d"% (avg_pr, ppval))
print("averge Spearman's r: %.7f, pval: %d"% (avg_sp, spval))

out = np.stack([y, preds])
np.save("preds.npy", out)
# plot the scatter
fig, ax = plt.subplots(figsize=(4,4))
ax.scatter(preds, y)
ax.set_xlabel("Predit")
ax.set_ylabel("True")
ax.set_title("Pred-True scatterplot")
fig.savefig("test.scatter.png", dpi=300)


