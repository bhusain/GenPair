# Preprocessing Script
# Rachel Eimen
# 2/16/18
# Description:  Program reads in GEM, deletes points <= 0, and populates
# an array of size 361 with a condensed version of the data

import pandas as pd
import matplotlib.pyplot as plt
import math
data = pd.read_csv('Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.GEM.txt',sep='\t')

#data = data[data.columns.values[:-1]]
data = data.transpose()



global_max = data.max().max()
global_min = data.min().min()

print(global_max)
print(global_min)

print(data.shape)

import numpy as np
import json

d = dict()

with open('output_new.txt', 'w') as file:
        for i in range (0, 100):        #(0,73599):     # g1 options
                if data.ix[:, i].max() > 0:
                        for j in range(i + 1, 100):       # 73599):       # g2 options
                                if data.ix[:, j].max() > 0:
                                        gName = 'g' + str(i + 1) + '_g' + str(j + 1)
                                        array = np.zeros(361, dtype = int)
                                        for k in range(0, 2016):        # run through each sample for every gene
                                                if (data.ix[k, i] != 'nan'):    # verify that entry is not equal to 'nan'
                                                        if (data.ix[k, i] > 0 and data.ix[k, j] > 0):   # if value is greater than 0
                                                                g1_coord = math.ceil(data.ix[k, i])             # round up
                                                                g2_coord = math.ceil(data.ix[k, j])
                                                                array[int(g1_coord - 1 + (18*g2_coord))] += 1   # add 1 to the array at index = g1 - 1 + (18*g2)
                                        d[gName] = array
        # write array to file
        for key in d:
                file.write(key)
                for i in d[key]:
			o = ' '.join(str(i) for i in d[key])
		file.write(o)

        file.close()

