import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd


runNum=int(sys.argv[1])
pathData="../../version1Data/1d/run"+str(runNum)+".funcquadraticCartesian/row0"

#this script computes thermal expansion coefficient


TVals=[]
TFileNames=[]

for TFile in glob.glob(pathData+"/T*"):

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))




#sort T files

sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]



def compute_alpha(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: one alpha value
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))

    LFilePath=oneTFile+"/jsonData/jsonL/LData.json"
    with open (LFilePath,"r") as fptr:
        LData=json.load(fptr)
    LVec=LData["L"]


    UFilePath=oneTFile+"/jsonData/jsonU/UData.json"
    with open (UFilePath,"r") as fptr:
        UData=json.load(fptr)

    UVec=UData["U"]

    LVec=np.array(LVec)

    UVec=np.array(UVec)

    LUProd=LVec*UVec

    LUMean=np.mean(LUProd)

    LMean=np.mean(LVec)

    UMean=np.mean(UVec)

    alphaVal=1/(TVal**2*LMean)*(LUMean-LMean*UMean)

    return alphaVal



tStart=datetime.now()

alphaAll=[]

for oneTFile in sortedTFiles:
    alphaTmp=compute_alpha(oneTFile)
    alphaAll.append(alphaTmp)
tEnd=datetime.now()

print("alpha time: ",tEnd-tStart)

plt.figure()

plt.plot(sortedTVals,alphaAll,color="black")
plt.scatter(sortedTVals,alphaAll,color="red")
plt.title("Thermal expansion")
plt.xlabel("$T$")
plt.ylabel("$\\alpha$")
plt.ylim((-0.1,0.1))
plt.savefig(pathData+"/alpha.png")
plt.close()






