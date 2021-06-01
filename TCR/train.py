import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob,sys,os,json
from datetime import datetime

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data.dataset import Dataset
from torch.utils.data.dataloader import DataLoader
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, accuracy_score

from .datasets import TCRDataset
from .model import ImmunoRNN

class_num = 7
input_size = 1900
batch_size = 1000
num_workers = 1
learning_rate = 0.001
num_epochs = 10


device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')


### TODO
X = np.load("igh.unirep.npy")
X_full = pd.read_table("igh.data.txt")
y = X_full['label']
le = LabelEncoder()
y2 = le.fit_transform(y)


train_data = TCRDataset(X_train, y_train)
test_data = TCRDataset(X_test, y_test)
print("Prepare DataLoader ")
train_loader = DataLoader(train_data, batch_size=batch_size, num_workers= num_workers)#sampler=train_sampler, num_workers=1 )# sampler=SubsetRandomSampler() )
test_loader =  DataLoader(test_data, batch_size=batch_size, num_workers= num_workers)




print("Build Model")
model = ImmunoRNN(input_size, class_num)
model.to(device)
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
# let learning_rate decrease by 50% at 500, 1000 and 2000-th epoch
scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, [500, 1000, 2000], gamma=0.5)


print("Start training")
# Training the Model
model.train()
last_loss = 1000
for epoch in range(num_epochs):
    running_loss = 0.0
    #running_loss2 = 0.0
    for i, embeds  in enumerate(train_loader):
        inputs, targets = embeds['embed'], embeds['target']
        inputs = inputs.to(device)
        targets = targets.to(device)
        optimizer.zero_grad()
        outputs = model(inputs)
        # target size must be the same as ouput size
        loss = criterion(outputs, targets) # 
        loss.backward()
        optimizer.step()
        scheduler.step()
        # print statistics
        running_loss += loss.item()
    if (epoch+1) % 10 == 0:     
        print('epoch [%d] loss: %.7f' % (epoch + 1, running_loss /100))
        PATH = f'checkpoints/tcr.model.epoch.{epoch+1}.pth'
        if running_loss < last_loss:
            last_loss = min(last_loss, running_loss)
            torch.save({'model_state_dict': model.state_dict(),         
                        'optimizer_state_dict': optimizer.state_dict()}, PATH)


model.eval()
with torch.no_grad():
    test_loss = 0
    auc = 0
    acc = 0
    y_test = []
    y_pred_prob = []
    for embeds in test_loader:
        inputs, targets = embeds['embed'], embeds['target']
        inputs = inputs.reshape((-1, input_size)).to(device)
        targets = targets.view(-1).to(device)
        #targets = targets.view(-1).type(torch.int)
        outputs = model(inputs)
        test_loss += criterion(outputs, targets) # 
        # y_pred = F.sigmoid(outputs).detach().cpu().numpy()
        y_pred = F.softmax(outputs, dim=-1).detach().cpu().numpy()
        y_pred_prob.append(y_pred)
        targets = targets.detach().cpu().numpy()
        y_test.append(targets)
        #fpr, tpr, thresholds = roc_curve(targets, y_pred, pos_label=1)
        #precision, recall, pr_threshold = precision_recall_curve(targets, y_pred, pos_label=1)
y_test = np.concatenate(y_test)
y_pred_prob = np.concatenate(y_pred_prob)


fpr, tpr, thresholds = roc_curve(y_test, y_pred, pos_label=1)
precision, recall, pr_threshold = precision_recall_curve(y_test, y_pred, pos_label=1)

fig, ax = plt.subplots(1,2, figsize=(8,4))
ax[0].plot(fpr, tpr, color='darkorange',
           lw=2, label="MLP AUC = {:.2f}".format(auc))
ax[0].legend()
ax[0].plot([0,1],[0,1], linestyle='--', color='navy', lw=2,)
ax[0].set_title(f"ROC")
ax[0].set_xlim([0.0, 1.0])
ax[0].set_ylim([0.0, 1.05])
ax[0].set_xlabel('False Positive Rate')
ax[0].set_ylabel('True Positive Rate')
## precision recall
ax[1].step(recall, precision, color='g', alpha=0.2, where='post')
ax[1].fill_between(recall, precision, alpha=0.2, color='g', step='post')
ax[1].set_xlabel('Recall')
ax[1].set_ylabel('Precision')
ax[1].set_ylim([0.0, 1.0])
ax[1].set_xlim([0.0, 1.0])
ax[1].set_title('Precision-Recall curve')
plt.tight_layout()
plt.show()