# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 19:42:51 2021

@author: chenwei
"""
import pandas as pd
import numpy as np


filenames = 'E:/SERT_19800/code_mypaper/MUV/'
data = pd.read_csv(filenames + 'IC50_not_Kd.csv')
#data = pd.read_csv(filenames + 'IC50_samples.csv')
smiles = data['SMILES']
targets = data['Target Sequence']
labels = data['Label']

# 正样本占1/2, 1/4, 1/6
inde = len(labels) // 2
sort_label = labels.sort_values()
mdi = sort_label[inde]

smi = []
for ind in range(len(smiles)):
    
    if labels[ind] > mdi:
        tmp = [smiles[ind], targets[ind], '0']
    else:
        tmp = [smiles[ind], targets[ind], '1']
    smi.append(tmp)
    
# 处理后的数据写入txt文本中
with open(filenames+'bindingdb_data1.txt','a') as f_test:
    for line in smi:
        line=' '.join(line)
        f_test.write(line+'\n')