import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from scipy.optimize import root

#python readCSV.py path rowNum
if len(sys.argv)!=3:
    print("wrong number of arguments")

fileName=sys.argv[1]
rowNum=int(sys.argv[2])

dfStr=pd.read_csv(fileName)

oneRow=dfStr.iloc[rowNum,:]

a1=float(oneRow.loc["a1"])
a2=float(oneRow.loc["a2"])
c1=float(oneRow.loc["c1"])
c2=float(oneRow.loc["c2"])




print("a1"+str(a1)+"a2"+str(a2)+"c1"+str(c1)+"c2"+str(c2))




