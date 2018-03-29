
# coding: utf-8

# In[3]:


# Preprocessing - Gene Pair Plotting Script
# Rachel Eimen
# 2/23/18
# Description:  Plots the average of the first and second genes for every cluster; colors them based on cancer type in tissue


import pandas as pd
import matplotlib.pyplot as plt
import math

data = pd.read_csv('Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.GEM.txt',sep='\t')
data = data.transpose()

genelist = pd.read_csv('output_new.txt')


# In[ ]:


#xList = {}  # create gene frequency count for x and y as a dicts
#yList = {}
color = []  # create list to store colors (type) of each tissue

for i in list(data.index):
    if i[0] == 'B':
        color.append(1)
    if i[0] == 'O':
        color.append(2)
    if i[0] == 'L':
        color.append(3)
    if i[0] == 'T':
        color.append(4)
    if i[0] == 'G':
        color.append(5)

dataRows, dataCols = data.shape
#print(data.shape)
geneRows, geneCols = genelist.shape
#print ("genelist length = " + str(geneRows))
#print(genelist.shape) # (rows, cols)

#for i in color:
#    print(color[i])
        
import numpy as np
import json


# In[ ]:


# initialize and fill dictionaries which store the frequency of each gene for x and y

xList = dict() # dictionary
yList = dict() # dictionary

for i in range(1, geneRows):
    strTemp = "g" + str(i)
    xList[strTemp] = 0
    yList[strTemp] = 0

#with open('genePair_output.txt', 'w') as file:
with open('output_new.txt', 'r') as genelist:
    for line in genelist:
        gene = line.split(' = ')  # split by " = "
#        print(gene[0])
        gene_separate = gene[0].split('_') # split by "_"
#        print(gene_separate[0] + " " + gene_separate[1])
        xList[gene_separate[0]] = xList[gene_separate[0]] + 1
        yList[gene_separate[1]] = yList[gene_separate[1]] + 1
#        print(str(xList[gene_separate[0]]) + " " + str(yList[gene_separate[1]]))


# In[ ]:


# read through GEM to determine values for each gene; values which will be multiplied by gene count in x and y

gemDict = dict()

for i in range(1, dataCols):  #(0, dataCols):  # genes
#    strTemp = str(data.ix[0,i])  # gene name
#    i+=1
    strTemp = "g" + str(i)
    geneInfo = [None]*dataRows  # temporary list for each gene
    for j in range(1, dataRows): # tissues
        if ((np.isnan(data.ix[j,i])) | (data.ix[j,i] < 0)):
            geneInfo[j] = 0
        else:
            geneInfo[j] = data.ix[j,i]
    gemDict[strTemp] = geneInfo
    
#print(data.ix[1,1])
        
    


# In[ ]:


#print(gemDict)


# In[ ]:


# average gene expression for each tissue; final x and y values for graph
xFinal = [0]*dataRows  # final values; create two lists filled with 0s; entry for every tissue
yFinal = [0]*dataRows

#xTotal = [0]*dataRows  # total count of gene samples; create two arrays filled with 0s
#yTotal = [0]*dataRows


for i in range(1, dataRows):  # tissues
#    strTemp = "g" + str(i)
    xTotal = 0
    yTotal = 0
    for j in range(1, dataCols):  # genes
        strTemp = "g" + str(j)
        if strTemp in xList:
            xFinal[i] = xFinal[i] + (gemDict[strTemp][i] * xList[strTemp])
            xTotal += xList[strTemp]
        if strTemp in yList:
            yFinal[i] = yFinal[i] + (gemDict[strTemp][i] * yList[strTemp])
            yTotal += yList[strTemp]
#        xTotal[i] = xTotal[i] + xList[strTemp]
#        yTotal[i] = yTotal[i] + yTotal[strTemp]
    if xTotal > 0:
        xFinal[i] = xFinal[i] / xTotal
    if yTotal > 0:
        yFinal[i] = yFinal[i] / yTotal


# In[ ]:


with open('genePairOut.txt', 'w') as file:
    file.write("Colors List")
    for i in color:
        file.write(str(i))
        
    file.write("x Values")
    for i in xFinal:
        file.write(str(i))
        
    file.write("y Values")
    for i in yFinal:
        file.write(str(i))

