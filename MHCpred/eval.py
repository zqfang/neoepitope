import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
import glob,sys,os,json
from datetime import datetime
import joblib
import scipy.stats as stat
import torch
import math
from sklearn.metrics import confusion_matrix
from tqdm import tqdm
from torch.utils.data.dataloader import DataLoader
from datasets import MHCDataset, MHCEvalDataset
from model import MHCModel
import config


device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

def topK_ppv(out,topK=100):
    """
    if there is a very large evaluation dataset, the topK chosen would be a oversampling of extremes,
    so the random part is to make sure the sampling is not too bias
    """
    
    preds=pd.DataFrame(out,columns=['y','preds'])
    sample_size=topK*10
    if topK*10>len(preds):
        sample_size=len(preds)
    preds = preds.sample(sample_size,random_state=2021).sort_values('preds',ascending=True).iloc[:topK,:]
    tn, fp, fn, tp = confusion_matrix(preds['y'].tolist(), preds['preds'].tolist()).ravel()
    ppv = tp / (tp + fp)
    return ppv
    
def ba_score2binary(ba_score,cutoff):
    """
    convert ba_score back into ic50 values,tipical cutoff for binder and no binder is 5000nM, 
    strong vs weak binding may diff in 100, depends on the goal to test
    """
    ba_score = math.e**((1-ba_score)*np.log(50000))
    ba_score[ba_score<cutoff]=1
    ba_score[ba_score>=cutoff]=0
    return ba_score


PATH = config.args.data_path
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

#only be evaluted when both y and pred are continous value 
avg_sp, spval = stat.spearmanr(preds, y)
avg_pr, ppval = stat.pearsonr(preds, y)
print("averge Pearson's r: %.7f, pval: %d"% (avg_pr, ppval))
print("averge Spearman's r: %.7f, pval: %d"% (avg_sp, spval))

preds_binares = []
preds_header = []
for cutoff in range(50, 501, 50):
    pred_b = ba_score2binary(preds, cutoff) #convert to binary with cutoff
    preds_binares.append(pred_b)
    print(f"BA Cutoff: {cutoff}")
    print(confusion_matrix(y, pred_b))
# y = ba_score2binary(y,100) # if y of the data is not binary
out = np.stack([y, preds] + preds_binares, axis=1)
np.save("preds.npy", out)
# plot the scatter
# fig, ax = plt.subplots(figsize=(4,4))
# ax.scatter(preds, y)
# ax.set_xlabel("Predit")
# ax.set_ylabel("True")
# ax.set_title("Pred-True scatterplot")
# fig.savefig("test.scatter.png", dpi=300)
print(topK_ppv(out[:,[0,2]], 100)) # cutoff 100


