# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
from os import listdir
from os.path import isfile, join
mypath=r"C:\Users\13695\Documents\MATLAB\Macro\hw3_data"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-3:]=="csv"]
createVar = locals()

for file in onlyfiles:
    createVar["df_"+file[:-4]]=pd.read_csv(file)
    
df_incomeshare= pd.merge(df_GDP, df_PINCOME, on='DATE')

df_incomeshare.plot()

