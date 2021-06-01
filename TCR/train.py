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



class_num = 7
input_size = 1900
batch_size = 1000
num_workers = 1
learning_rate = 0.001
num_epochs = 10
device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

X = np.load("igh.unirep.npy")
X_full = pd.read_table("igh.data.txt")
y = X_full['label']

le = LabelEncoder()
y2 = le.fit_transform(y)


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

    # dim: (batch, sample)

train_data = TCRDataset(X_train, y_train)
test_data = TCRDataset(X_test, y_test)
print("Prepare DataLoader ")
train_loader = DataLoader(train_data, batch_size=batch_size, num_workers= num_workers)#sampler=train_sampler, num_workers=1 )# sampler=SubsetRandomSampler() )
test_loader =  DataLoader(test_data, batch_size=batch_size, num_workers= num_workers)


# model 1    
class MLP(nn.Module):
    def __init__(self, input_size, class_num):
        super(MLP,self).__init__()
        # number of hidden nodes in each layer (512)
        hidden_1 = 1024
        hidden_2 = 512
        # linear layer (input_size -> hidden_1)
        self.fc1 = nn.Linear(input_size, hidden_1)
        # linear layer (n_hidden -> hidden_2)
        self.fc2 = nn.Linear(hidden_1, hidden_2)
        # linear layer (n_hidden -> 1)
        self.fc3 = nn.Linear(hidden_2, class_num)
        # dropout layer (p=0.2)
        # dropout prevents overfitting of data
        self.droput = nn.Dropout(0.2)
        
    def forward(self,x):
        # add hidden layer, with relu activation function
        x = F.relu(self.fc1(x))
        # add dropout layer
        x = self.droput(x)
         # add hidden layer, with relu activation function
        x = F.relu(self.fc2(x))
        # add dropout layer
        x = self.droput(x)
        # add output layer
        x = self.fc3(x)
        return x

# model 2
class ImmunRepRNN(nn.Module):
    def __init__(self, seq_length, hidden_size, output_size, n_layers=1):
        super(ImmunRepRNN, self).__init__()
        self.in_channel = 1 #input_size
        self.hidden_size = hidden_size
        self.output_size = output_size
        self.n_layers = n_layers
        self.seq_len = seq_length

        self.c1 = nn.Conv1d(in_channels=1, out_channels=hidden_size, kernel_size=3)
        self.p1 = nn.AvgPool1d(2)
        self.c2 = nn.Conv1d(hidden_size, hidden_size, 3)
        self.p2 = nn.AvgPool1d(2)
        self.gru = nn.GRU(hidden_size, hidden_size, n_layers, dropout=0.01, bidirectional=True)
        self.fc1 = nn.Linear(hidden_size*2, output_size) # biRNN: concat the bidirectional embeds

    def forward(self, inputs, hidden=None):
        # input size: Batch, seq_len, input_size 
        batch_size = inputs.size(0)
        
        # Turn (batch_size x seq_len x seq_emb_szie) into (batch_size x seq_embed_size x seq_len) for CNN
        ## conv1d input: (batch, channels, W)
        # inputs = inputs.transpose(1, 2)
        inputs = inputs.view(batch_size, 1, self.seq_len)
        # Run through Conv1d and Pool1d layers
        c = self.c1(inputs)
        p = self.p1(c)
        c = self.c2(p)
        p = self.p2(c)
        # Turn (batch_size x hidden_size x seq_len) back into (seq_len x batch_size x hidden_size) for RNN
        p = p.permute(2, 0, 1)
        #breakpoint()
        p = torch.tanh(p)
        output, hidden = self.gru(p, hidden) ## output for seq2seq, last hidden for classification 
        conv_seq_len = output.size(0)
        hidden = hidden.view(self.n_layers, 2, batch_size, self.hidden_size) # 2 for bidirectional
        last_hidden = hidden[-1]
        # last_hidden_fwd = last_hidden[0]
        # last_hidden_bwd = last_hidden[1]
        # convert to (Batch, num_layers, hidden)
        last_hidden = last_hidden.permute(1,0,2).reshape(batch_size, -1).contiguous() ## concat fwd and bwd
        out = self.fc1(last_hidden)
        return out



print("Build Model")
# model = LogisticRegression(input_size)
model = MLP(input_size, class_num)
model.to(device)

# weight = [class_sample_count[0] / class_sample_count[1]]
# weight = [(len(y) - sum(y))/sum(y)]
# criterion = nn.BCEWithLogitsLoss(pos_weight= torch.FloatTensor(weight).to(device))
# optimizer = torch.optim.SGD(model.parameters(), 
#                             lr=learning_rate, 
#                             momentum=0.9, weight_decay=0.0005) 
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(),lr=learning_rate,)
# let learning_rate decrease by 50% at 500, 1000 and 2000-th epoch
scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, [500, 1000, 2000], gamma=0.5)



print("Start training")
# Training the Model
num_epochs = 1000
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