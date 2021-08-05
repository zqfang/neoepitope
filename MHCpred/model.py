import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F




class MHCConv(nn.Module):
    def __init__(self, seq_length, embed_size= 512):
        """
        seq_length: dim of unirep output, e.g. 1900
        """
        super(MHCConv, self).__init__()
        self.in_channel = 1 #input_size
        self.hidden_size = embed_size
        self.seq_len = seq_length

        self.c1 = nn.Conv1d(in_channels=1, out_channels=self.hidden_size, kernel_size=3)
        self.b1 = nn.BatchNorm1d(self.hidden_size)
        self.p1 = nn.AvgPool1d(2)
        self.c2 = nn.Conv1d(self.hidden_size, self.hidden_size, 3)
        self.b2 = nn.BatchNorm1d(self.hidden_size)
        self.p2 = nn.AvgPool1d(2)

    def forward(self, inputs):
        # input size: Batch, seq_len 
        batch_size = inputs.size(0)
        # Turn (batch_size x seq_len x seq_emb_szie) into (batch_size x seq_embed_size x seq_len) for CNN
        ## conv1d input: (batch, channels, W)
        # inputs = inputs.transpose(1, 2)
        # inputs = inputs.view(batch_size, 1, self.seq_len)
        inputs = inputs.unsqueeze(1) # dim(1)
        # Run through Conv1d and Pool1d layers
        c = self.c1(inputs)
        c = F.relu(self.b1(c))
        p = self.p1(c)
        c = self.c2(p)
        c = F.relu(self.b2(c))
        out = torch.tanh(self.p2(c))
        return out

class EpitopeConv(nn.Module):
    def __init__(self, seq_length, embed_size=512):
        """
        seq_length: dim of unirep output, e.g. 1900
        """
        super(EpitopeConv, self).__init__()
        self.in_channel = 1 #input_size
        self.hidden_size = embed_size
        self.seq_len = seq_length

        self.c1 = nn.Conv1d(in_channels=1, out_channels=self.hidden_size, kernel_size=3)
        self.p1 = nn.AvgPool1d(2)
        self.c2 = nn.Conv1d(self.hidden_size, self.hidden_size, 3)
        self.p2 = nn.AvgPool1d(2)
        self.tanh = nn.Tanh()

    def forward(self, inputs):
        # input size: Batch, seq_len, 
        batch_size = inputs.size(0)
        # Turn (batch_size x seq_len x seq_emb_szie) into (batch_size x seq_embed_size x seq_len) for CNN
        ## conv1d input: (batch, channels, W)
        # inputs = inputs.transpose(1, 2)
        # inputs = inputs.view(batch_size, 1, self.seq_len)
        inputs = inputs.unsqueeze(1)
        # Run through Conv1d and Pool1d layers
        c = self.c1(inputs)
        p = self.p1(c)
        c = self.c2(p)
        out = self.p2(c)
        out = self.tanh(out)
        return out


# model   
class MHCModel(nn.Module):
    def __init__(self, input_size, output_size = 1, n_layers=1):
        super(MHCModel,self).__init__()
        # number of hidden nodes in each layer (512)
        self.hidden_size = 512
        self.n_layers = n_layers
        self.ag = EpitopeConv(input_size, self.hidden_size)
        self.mhc = MHCConv(input_size, self.hidden_size)
              
        self.hidden_size *= 2
        self.gru = nn.GRU(self.hidden_size, self.hidden_size, n_layers, dropout=0, bidirectional=True)
        #self.fc = nn.Linear(self.hidden_size*2, self.hidden_size)
        self.attn = nn.MultiheadAttention(self.hidden_size, num_heads=8) # Note: must self.hidden_size // num_heads == 0
        self.ffn = torch.nn.Sequential(torch.nn.Linear(self.hidden_size*2, self.hidden_size),
                                       #torch.nn.BatchNorm1d(self.hidden_size),
                                       torch.nn.ReLU(),
                                       torch.nn.Dropout(0.2),
                                       torch.nn.Linear(self.hidden_size, output_size))
        
    def forward(self, mhc, antigen, hidden=None):

        # both output are  (batch_size x seq_embed_size x seq_len)
        mhc = self.mhc(mhc)
        ag = self.ag(antigen)
        # concat embedding dim
        x = torch.cat([ag, mhc], dim=1)

        # Turn (batch_size x hidden_size x seq_len) back into (seq_len x batch_size x hidden_size) for RNN
        p = x.permute(2, 0, 1)
        seq_len, batch_size, hidden_size = p.size()
        # out, out_weight = self.attn(p, p, p) # k, q, v is same for transformers
        output, hidden = self.gru(p, hidden) ## output for seq2seq, last hidden for classification, regression
        hidden = hidden.view(self.n_layers, 2, batch_size, hidden_size) # 2 for bidirectional
        last_hidden = hidden[-1]
        # last_hidden_fwd = last_hidden[0]
        # last_hidden_bwd = last_hidden[1]

        # convert to (Batch, num_layers, hidden), then flatten
        last_hidden = last_hidden.permute(1, 0, 2).reshape(batch_size, -1).contiguous() ## concat fwd and bwd
        out = self.ffn(last_hidden)
        return out.squeeze(-1)
# model 


class LayerNorm(nn.Module):
    def __init__(self, d_model, eps = 1e-6):
        super().__init__()
        self.size = d_model
        # create two learnable parameters to calibrate normalisation
        self.alpha = nn.Parameter(torch.ones(self.size))
        self.bias = nn.Parameter(torch.zeros(self.size))
        self.eps = eps
    def forward(self, x):
        norm = self.alpha * (x - x.mean(dim=-1, keepdim=True)) \
        / (x.std(dim=-1, keepdim=True) + self.eps) + self.bias
        return norm
    

class Attention(nn.Module):
    def __init__(self, feature_dim, step_dim, context_dim):
        super(Attention, self).__init__()
        self.feature_dim = feature_dim
        self.step_dim = step_dim
        self.context_dim = context_dim
        self.tanh = nn.Tanh()

        weight = torch.zeros(feature_dim, context_dim)
        nn.init.kaiming_uniform_(weight)
        self.weight = nn.Parameter(weight)
        self.b = nn.Parameter(torch.zeros(step_dim, context_dim))

        u = torch.zeros(context_dim, 1)
        nn.init.kaiming_uniform_(u)
        self.context_vector = nn.Parameter(u)

    def forward(self, x):
        eij = torch.matmul(x, self.weight)
        # eij = [batch_size, seq_len, context_dim]
        eij = self.tanh(torch.add(eij, self.b))
        # eij = [batch_size, seq_len, context_dim]
        v = torch.exp(torch.matmul(eij, self.context_vector))  # dot product
        # v = [batch_size, seq_len, 1]
        v = v / (torch.sum(v, dim=1, keepdim=True))
        # v = [batch_size, seq_len, 1]
        weighted_input = x * v
        # weighted_input = [batch_size, seq_len, 2*hidden_dim]             -> 2 : bidirectional
        s = torch.sum(weighted_input, dim=1)
        # s = [batch_size, 2*hidden_dim]                                   -> 2 : bidirectional
        return s