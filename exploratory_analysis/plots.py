
import numpy as np
import pandas as pd
import os

import seaborn as sns

import matplotlib.pyplot as plt

first = pd.DataFrame.from_csv('/users/ryan/downloads/Default Dataset.csv', header=None, index_col=None)
second = pd.DataFrame.from_csv('/users/ryan/downloads/Default Dataset (1).csv', header=None, index_col=None)
third = pd.DataFrame.from_csv('/users/ryan/downloads/Default Dataset (2).csv', header=None, index_col=None)
fourth = pd.DataFrame.from_csv('/users/ryan/downloads/Default Dataset (3).csv', header=None, index_col=None)
sns.set_style("whitegrid")
plt.figure(figsize=(8,6))
plt.plot(first[0]/1000000,first[1], second[0]/1000000,second[1], third[0]/1000000,third[1],fourth[0]/1000000,fourth[1])
plt.axis([0, 61.5, 0, 50000])
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('Reads (millions)', size=16)
plt.ylabel('Transcripts at 10x coverage',size=16)
plt.legend(['sample 1','sample 2','sample 3','sample 4'], loc='best', fontsize=16)

#plt.set_style("whitegrid")
plt.show()