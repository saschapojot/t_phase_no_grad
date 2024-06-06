import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd

#This script loads json data and plot


pathData="../../version1Data/1d/funcquadraticRing/row0"
inCsv="../../version1Input/1d/ringQuadratic/ringQuadratic.csv"
dfStr=pd.read_csv(inCsv)

oneRow=dfStr.iloc[0,:]

a1=float(oneRow.loc["a1"])
a2=float(oneRow.loc["a2"])
c1=float(oneRow.loc["c1"])
c2=float(oneRow.loc["c2"])
mA=float(oneRow.loc["mA"])
mB=float(oneRow.loc["mB"])


TVals=[]
TFileNames=[]

for TFile in glob.glob(pathData+"/T*"):

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))






#sort T files

sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]


def pltU(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: U plots
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    UFilePath=oneTFile+"/jsonData/jsonU/UData.json"
    with open(UFilePath, 'r') as fptr:
        data = json.load(fptr)
    UVec=np.array(data["U"])
    print("T="+str(TVal)+", data num="+str(len(UVec)))
    meanU=np.mean(UVec)

    varU=np.var(UVec,ddof=1)
    sigmaU=np.sqrt(varU)

    nbins=500
    fig=plt.figure()
    axU=fig.add_subplot()
    (n0,_,_)=axU.hist(UVec,bins=nbins)
    meanU=np.round(meanU,4)
    sigmaU=np.round(sigmaU,4)

    axU.set_title("T="+str(TVal))
    axU.set_xlabel("$U$")
    axU.set_ylabel("#")
    xPosUText=(np.max(UVec)-np.min(UVec))*1/2+np.min(UVec)
    yPosUText=np.max(n0)*2/3
    axU.text(xPosUText,yPosUText,"mean="+str(meanU)+"\nsd="+str(sigmaU))
    plt.axvline(x=meanU,color="red",label="mean")
    axU.text(meanU*1.1,0.5*np.max(n0),str(meanU)+"$\pm$"+str(sigmaU),color="red")
    axU.hlines(y=0,xmin=meanU-sigmaU,xmax=meanU+sigmaU,color="green",linewidth=15)

    plt.legend(loc="best")

    EHistOut="T"+str(TVal)+"UHist.png"
    plt.savefig(oneTFile+"/"+EHistOut)

    plt.close()

    ### test normal distribution for mean U
    #block mean
    USelectedAll=UVec
    def meanPerBlock(length):
        blockNum=int(np.floor(len(USelectedAll)/length))
        UMeanBlock=[]
        for blkNum in range(0,blockNum):
            blkU=USelectedAll[blkNum*length:(blkNum+1)*length]
            UMeanBlock.append(np.mean(blkU))
        return UMeanBlock
    fig=plt.figure(figsize=(20,20))
    fig.tight_layout(pad=5.0)
    lengthVals=[2,5,7,10]
    for i in range(0,len(lengthVals)):
        l=lengthVals[i]
        UMeanBlk=meanPerBlock(l)
        ax=fig.add_subplot(2,2,i+1)
        (n,_,_)=ax.hist(UMeanBlk,bins=100,color="aqua")
        xPosTextBlk=(np.max(UMeanBlk)-np.min(UMeanBlk))*1/7+np.min(UMeanBlk)
        yPosTextBlk=np.max(n)*3/4
        meanTmp=np.mean(UMeanBlk)
        meanTmp=np.round(meanTmp,3)
        sdTmp=np.sqrt(np.var(UMeanBlk))
        sdTmp=np.round(sdTmp,3)
        ax.set_title("L="+str(l))
        ax.text(xPosTextBlk,yPosTextBlk,"mean="+str(meanTmp)+", sd="+str(sdTmp))
    fig.suptitle("T="+str(TVal))
    plt.savefig(oneTFile+"/T"+str(TVal)+"UBlk.png")
    # plt.savefig(EBlkMeanDir+"/T"+str(TTmp)+"EBlk.png")
    plt.close()


def plt_theta(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: plots of positions
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))

    thetaFilePath=oneTFile+"/jsonData/jsontheta/thetaData.json"

    with open (thetaFilePath,"r") as fptr:
        thetaData=json.load(fptr)

    theta0BVec=thetaData["theta0B"]
    theta1AVec=thetaData["theta1A"]
    theta1BVec=thetaData["theta1B"]


    rFilePath=oneTFile+"/jsonData/jsonr/rData.json"
    with open(rFilePath,"r") as fptr:
        rData=json.load(fptr)

    rVec=rData["r"]

    theta0AVec=[]
    for i in range(0,len(theta0BVec)):
        theta0BTmp=theta0BVec[i]
        theta1ATmp=theta1AVec[i]
        theta1BTmp=theta1BVec[i]
        theta0AVec.append(-mB/mA*theta0BTmp-theta1ATmp-mB/mA*theta1BTmp)
    theta0AVec=np.array(theta0AVec)

    d0A0BVec=[]
    d0B1AVec=[]
    d1A1BVec=[]
    d1B0AVec=[]

    for i in range(0,len(rVec)):
        rTmp=rVec[i]
        theta0ATmp=theta0AVec[i]
        theta0BTmp=theta0BVec[i]
        theta1ATmp=theta1AVec[i]
        theta1BTmp=theta1BVec[i]

        d0A0BVec.append(rTmp*(theta0BTmp-theta0ATmp))

        d0B1AVec.append(rTmp*(theta1ATmp-theta0BTmp))

        d1A1BVec.append(rTmp*(theta1BTmp-theta1ATmp))

        d1B0AVec.append(rTmp*(2*np.pi+theta0ATmp-theta1BTmp))


    #summary of distance between neighboring points
    d0A0BMean=np.mean(d0A0BVec)
    d0B1AMean=np.mean(d0B1AVec)
    d1A1BMean=np.mean(d1A1BVec)
    d1B0AMean=np.mean(d1B0AVec)
    smrDistFileName=oneTFile+"/distSummary.txt"

    dist=[d0A0BMean,d0B1AMean,d1A1BMean,d1B0AMean]
    fptrTxt=open(smrDistFileName,"w")
    fptrTxt.write("T="+str(TVal)+"\n")
    fptrTxt.write("position: "+str(dist)+"\n")
    fptrTxt.close()



    d0A0BVar=np.var(d0A0BVec,ddof=1)
    d0B1AVar=np.var(d0B1AVec,ddof=1)
    d1A1BVar=np.var(d1A1BVec,ddof=1)
    d1B0AVar=np.var(d1B0AVec,ddof=1)

    d0A0BSd=np.sqrt(d0A0BVar)
    d0B1ASd=np.sqrt(d0B1AVar)
    d1A1BSd=np.sqrt(d1A1BVar)
    d1B0ASd=np.sqrt(d1B0AVar)

    distAll=[0,d0A0BMean,d0B1AMean,d1A1BMean,d1B0AMean]

    posAll=np.cumsum(distAll)
    pos_A=[posAll[i] for i in range(0,len(posAll),2)]
    pos_B=[posAll[i] for i in range(1,len(posAll),2)]

    sd_A=[d1B0ASd,d0B1ASd,d1B0ASd]

    sd_B=[d0A0BSd,d1A1BSd]
    print(sd_A)
    print(sd_B)

    plt.figure(figsize=(12, 6))
    plt.ylim(-1, 1)
    for i in range(0,len(pos_A)):
        plt.hlines(y=0,xmin=pos_A[i]-sd_A[i],xmax=pos_A[i]+sd_A[i],color="red",linewidth=2,alpha=0.2)
        plt.text(pos_A[i],0.3,str(i)+"A",color="blue", ha='center')
        plt.text(pos_A[i],-0.1,str(np.round(pos_A[i],4)),color="blue", ha='center',fontsize=8)


    plt.scatter(pos_A,[0]*len(pos_A),color="blue",s=8,label="A")

    for i in range(0,len(pos_B)):
        plt.hlines(y=0,xmin=pos_B[i]-sd_B[i],xmax=pos_B[i]+sd_B[i],color="magenta",linewidth=2,alpha=0.2)
        plt.text(pos_B[i],0.1,str(i)+"B",color="green", ha='center')
        plt.text(pos_B[i],-0.3,str(np.round(pos_B[i],4)),color="green", ha='center',fontsize=8)

    plt.scatter(pos_B,[0]*len(pos_B),color="green",s=8,label="B")
    plt.legend(loc="best")
    plt.gca().get_yaxis().set_visible(False)
    plt.axhline(y=0, color='black', linewidth=0.5,alpha=0.3)
    gridOut="T"+str(TVal)+"grid.pdf"
    plt.title("T="+str(TVal))
    plt.savefig(oneTFile+"/"+gridOut)

    plt.close()







tStatsStart=datetime.now()
for oneTFile in sortedTFiles:
    pltU(oneTFile)
    plt_theta(oneTFile)

tStatsEnd=datetime.now()
print("stats total time: ",tStatsEnd-tStatsStart)