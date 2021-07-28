import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F



class MHCBase(nn.Module):
    def __init__(self, seq_length, output_size, n_layers=1):
        """
        seq_length: dim of unirep output, e.g. 1900
        """
        super(MHCBase, self).__init__()
        self.in_channel = 1 #input_size
        self.hidden_size = 512
        self.output_size = output_size
        self.n_layers = n_layers
        self.seq_len = seq_length

        self.c1 = nn.Conv1d(in_channels=1, out_channels=self.hidden_size, kernel_size=3)
        self.p1 = nn.AvgPool1d(2)
        self.c2 = nn.Conv1d(self.hidden_size, self.hidden_size, 3)
        self.p2 = nn.AvgPool1d(2)
        self.gru = nn.GRU(self.hidden_size, self.hidden_size, n_layers, dropout=0.01, bidirectional=True)
        self.fc1 = nn.Linear(self.hidden_size*2, output_size) # biRNN: concat the bidirectional embeds

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



# model 1    
class MHCModel(nn.Module):
    def __init__(self, input_size, output_size = 1):
        super(MHCModel,self).__init__()
        # number of hidden nodes in each layer (512)
        self.hidden_size = 512
        self.ag = MHCBase(input_size, self.hidden_size, n_layers=1)
        self.mhc = MHCBase(input_size, self.hidden_size, n_layers=1)
        
        self.hidden_size *= 2
        self.ffn = torch.nn.Sequential(torch.nn.Linear(self.hidden_size, 2 * self.hidden_size),
                                       torch.nn.BatchNorm1d(2 * self.hidden_size),
                                       torch.nn.ReLU(),
                                       torch.nn.Dropout(0.2),
                                       torch.nn.Linear(2 * self.hidden_size, output_size)))
        
    def forward(self, mhc, antigen):
        # add hidden layer, with relu activation function
        mhc = self.mhc(mhc)
        ag = self.ag(antigen)
        x = torch.cat([ag, mhc], dim=1)
        x = self.ffn(x)
        return x

# model 