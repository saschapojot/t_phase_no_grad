
import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime

import pandas as pd


pathData="../../version1Data/1d/funcquadratic/row0"
TVals=[]
TFileNames=[]


for TFile in glob.glob(pathData+"/T*"):
    TFileNames.append(TFile)

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)

    TVals.append(float(matchT.group(1)))

#sort T files

sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]


def readGAA(oneTFile):
    """

    :param oneTFile:
    :return: fA, plot of fA
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    GAAFile=oneTFile+"/GAA.csv"
    GAAData=pd.read_csv(GAAFile,header=None)
    GAAMat=GAAData.to_numpy()
    return GAAMat


def readGBB(oneTFile):
    """

    :param oneTFile:
    :return: fB, plot of fB
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    GBBFile=oneTFile+"/GBB.csv"
    GBBData=pd.read_csv(GBBFile,header=None)
    GBBMat=GBBData.to_numpy()
    return GBBMat

def readGAB(oneTFile):
    """
    :param oneTFile:
    :return: fAB,
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    GABFile=oneTFile+"/GAB.csv"
    GABData=pd.read_csv(GABFile,header=None)
    GABMat=GABData.to_numpy()
    return GABMat

#parameters

whichT=0
rowNum=0
paramFile="../../version1Input/1d/quadratic/quadraticParams.csv"
dfStr=pd.read_csv(paramFile)
oneRow=dfStr.iloc[rowNum,:]

a1=float(oneRow.loc["a1"])
a2=float(oneRow.loc["a2"])
c1=float(oneRow.loc["c1"])
c2=float(oneRow.loc["c2"])


T=sortedTVals[whichT]
file=sortedTFiles[whichT]
def var_xA(j):
    """

    :param j: unit cell index
    :return: variance of xjA
    """
    val=1/2*j*c1**(-1)*T+1/2*j*c2**(-1)*T

    return val



def var_xB(j):
    """

    :param j: unit cell index
    :return: variance of xjB
    """

    val=1/2*(j+1)*c1**(-1)*T+1/2*j*c2**(-1)*T

    return val


GAAMat=readGAA(file)
GBBMat=readGBB(file)
GABMat=readGAB(file)
N,_=GAAMat.shape
unitCells=[j for j in range(0,N)]
#plot var xA
var_xA_mc=GAAMat.diagonal()

var_xA_theory=[var_xA(j) for j in range(0,N)]


plt.figure()
plt.plot(unitCells,var_xA_theory,color="blue",label="theory")
plt.scatter(unitCells,var_xA_mc,color="red",label="mc")

plt.legend(loc="best")
plt.xlabel("unit cell")
plt.title("var xA")
plt.xticks(unitCells)
plt.savefig(file+"/varxA.png")
plt.close()


#plot var xB
var_xB_mc=GBBMat.diagonal()
var_xB_theory=[var_xB(j) for j in range(0,N)]

plt.figure()
plt.plot(unitCells,var_xB_theory,color="green",label="theory")
plt.scatter(unitCells,var_xB_mc,color="magenta",label="mc")
plt.legend(loc="best")
plt.xlabel("unit cell")
plt.title("var xB")
plt.xticks(unitCells)

plt.savefig(file+"/varxB.png")
plt.close()

#plot cov(xjA,xjB)

cov_AB_mc=GABMat.diagonal()

def covAB(j):
    val=1/2*j*c1**(-1)*T+1/2*j*c2**(-1)*T
    return val

cov_AB_theory=[covAB(j) for j in range(0,N)]

plt.figure()
plt.plot(unitCells,cov_AB_theory,color="black",label="theory")
plt.scatter(unitCells,cov_AB_mc,color="red",label="mc")
plt.legend(loc="best")
plt.xlabel("unit cell")
plt.title("cov(xjA,xjB)")
plt.xticks(unitCells)

plt.savefig(file+"/covAB.png")
plt.close()