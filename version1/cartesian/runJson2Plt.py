import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd

#This script loads json data and plot for a run

runNum=int(sys.argv[1])
pathData="../../version1Data/1d/run"+str(runNum)+".funcquadraticCartesian/row0"



TVals=[]
TFileNames=[]

for TFile in glob.glob(pathData+"/T*"):

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))


sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]

def pltU(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: U plots, U mean, U var
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
    return [meanU,varU]


def plt_x(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: plots of positions, var(x0A), E(x0A), d0A0BMean,d1A1BMean, x0AVar,x1AVar
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))

    xFilePath=oneTFile+"/jsonData/jsonx/xData.json"

    with open (xFilePath,"r") as fptr:
        xData=json.load(fptr)

    x0AVec=xData["x0A"]
    x0BVec=xData["x0B"]
    x1AVec=xData["x1A"]
    x1BVec=xData["x1B"]


    LFilePath=oneTFile+"/jsonData/jsonL/LData.json"
    with open(LFilePath,"r") as fptr:
        LData=json.load(fptr)

    LVec=LData["L"]





    d0A0BVec=[]
    d0B1AVec=[]
    d1A1BVec=[]
    d1B0AVec=[]
    for i in range(0,len(x0BVec)):
        LTmp=LVec[i]
        x0ATmp=x0AVec[i]
        x0BTmp=x0BVec[i]
        x1ATmp=x1AVec[i]
        x1BTmp=x1BVec[i]
        d0A0BVec.append(x0BTmp-x0ATmp)
        d0B1AVec.append(x1ATmp-x0BTmp)
        d1A1BVec.append(x1BTmp-x1ATmp)
        d1B0AVec.append(x0ATmp-x1BTmp+LTmp)

    # #summary of distance between neighboring points
    d0A0BMean=np.mean(d0A0BVec)
    d0B1AMean=np.mean(d0B1AVec)
    d1A1BMean=np.mean(d1A1BVec)
    d1B0AMean=np.mean(d1B0AVec)

    x0AMean=np.mean(x0AVec)
    x0BMean=np.mean(x0BVec)
    x1AMean=np.mean(x1AVec)
    x1BMean=np.mean(x1BVec)

    unitCell0Pos=(x0AMean+x0BMean)/2
    unitCell1Pos=(x1AMean+x1BMean)/2
    unitCellPos=[unitCell0Pos,unitCell1Pos]

    x0AVar=np.var(x0AVec,ddof=1)
    x0BVar=np.var(x0BVec,ddof=1)
    x1AVar=np.var(x1AVec,ddof=1)
    x1BVar=np.var(x1BVec,ddof=1)

    x0ASd=np.sqrt(x0AVar)
    x0BSd=np.sqrt(x0BVar)
    x1ASd=np.sqrt(x1AVar)
    x1BSd=np.sqrt(x1BVar)

    pos_A=[x0AMean,x1AMean]
    sd_A=[x0ASd,x1ASd]
    pos_B=[x0BMean,x1BMean]
    sd_B=[x0BSd,x1BSd]





    smrDistFileName=oneTFile+"/distSummary.txt"
    #
    dist=[d0A0BMean,d0B1AMean,d1A1BMean,d1B0AMean]
    fptrTxt=open(smrDistFileName,"w")
    fptrTxt.write("T="+str(TVal)+"\n")
    fptrTxt.write("distances: "+str(dist)+"\n")
    fptrTxt.write("sd_A: "+str(sd_A)+"\n")
    fptrTxt.write("sd_B: "+str(sd_B)+"\n")
    fptrTxt.close()




    print("sd_A="+str(sd_A))
    print("sd_B="+str(sd_B))

    plt.figure(figsize=(12, 6))
    plt.ylim(-1, 1)
    for i in range(0,len(pos_A)):
        plt.hlines(y=0,xmin=pos_A[i]-sd_A[i],xmax=pos_A[i]+sd_A[i],color="red",linewidth=2,alpha=0.5)
        plt.text(pos_A[i],0.3,str(i)+"A",color="blue", ha='center')
        plt.text(pos_A[i],-0.1,"x"+str(i)+"A="+str(np.round(pos_A[i],4)),color="blue", ha='center',fontsize=8)
        plt.text(pos_A[i],0.1,"sd(x"+str(i)+"A)="+str(np.round(sd_A[i],4)),color="red", ha='center',fontsize=8)


    plt.scatter(pos_A,[0]*len(pos_A),color="blue",s=8,label="A")

    for i in range(0,len(pos_B)):
        plt.hlines(y=0,xmin=pos_B[i]-sd_B[i],xmax=pos_B[i]+sd_B[i],color="magenta",linewidth=2,alpha=0.5)
        plt.text(pos_B[i],0.3,str(i)+"B",color="green", ha='center')
        plt.text(pos_B[i],-0.1,"x"+str(i)+"B="+str(np.round(pos_B[i],4)),color="green", ha='center',fontsize=8)
        plt.text(pos_B[i],0.1,"sd(x"+str(i)+"B)="+str(np.round(sd_B[i],4)),color="magenta", ha='center',fontsize=8)
    plt.scatter(pos_B,[0]*len(pos_B),color="green",s=8,label="B")
    plt.legend(loc="best")
    plt.gca().get_yaxis().set_visible(False)
    plt.axhline(y=0, color='black', linewidth=0.5,alpha=0.3)
    plt.xticks(unitCellPos,[0,1])
    plt.xlabel("Unit cell")
    gridOut="T"+str(TVal)+"grid.pdf"
    plt.title("T="+str(TVal))
    plt.savefig(oneTFile+"/"+gridOut)

    plt.close()
    return [x0AVar,x0AMean,d0A0BMean,d1A1BMean,x0AVar,x1AVar]






xAVarAll=[]
xAMeanAll=[]
UMeanAll=[]
UVarAll=[]
d0A0BAll=[]
d1A1BAll=[]
varx0AAll=[]
varx1AAll=[]
tStatsStart=datetime.now()
for oneTFile in sortedTFiles:
    UMeanTmp,UVarTmp=pltU(oneTFile)
    UMeanAll.append(UMeanTmp)
    UVarAll.append(UVarTmp)
    varTmp,meanTmp,d0A0BTmp,d1A1BTmp,x0AVar,x1AVar=plt_x(oneTFile)
    xAVarAll.append(varTmp)
    xAMeanAll.append(meanTmp)
    d0A0BAll.append(d0A0BTmp)
    d1A1BAll.append(d1A1BTmp)
    varx0AAll.append(x0AVar)
    varx1AAll.append(x1AVar)

tStatsEnd=datetime.now()
print("stats total time: ",tStatsEnd-tStatsStart)



plt.figure()
plt.plot(sortedTVals,xAVarAll,color="black")
plt.scatter(sortedTVals,xAVarAll,color="red")
plt.title("var(x0A)")
plt.ylabel("var(x0A)")
plt.xlabel("$T$")
plt.savefig(pathData+"/varxA.png")
plt.close()

plt.figure()
plt.plot(sortedTVals,xAMeanAll,color="blue")
plt.scatter(sortedTVals,xAMeanAll,color="green")
plt.title("E(x0A)")
plt.ylabel("E(x0A)")
plt.xlabel("$T$")
plt.savefig(pathData+"/ExA.png")
plt.close()


plt.figure()
plt.plot(sortedTVals,UMeanAll,color="green")
plt.scatter(sortedTVals,UMeanAll,color="darkred")
plt.title("E(V)")
plt.xlabel("$T$")
plt.savefig(pathData+"/EV.png")
plt.close()

plt.figure()
plt.plot(sortedTVals,UVarAll,color="navy")
plt.scatter(sortedTVals,UVarAll,color="violet")
plt.title("var(V)")
plt.xlabel("$T$")
plt.savefig(pathData+"/varV.png")
plt.close()

plt.figure()
plt.scatter(sortedTVals,d0A0BAll,color="deeppink",label="d0A0B")
plt.scatter(sortedTVals,d1A1BAll,color="aqua",label="d1A1B",s=1)
plt.title("Intracell distance")
plt.xlabel("$T$")
plt.ylabel("distance")
plt.ylim((0.5,1.5))
plt.legend(loc="best")
plt.xscale("log")
plt.savefig(pathData+"/intracellDistance.png")
plt.close()

plt.figure()
plt.scatter(sortedTVals,varx0AAll,color="blue",label="var(x0A)")
plt.scatter(sortedTVals,varx1AAll,color="red",label="var(x1A)",s=2)
plt.title("var(x0A) and var(x1A)")
plt.xlabel("$T$")
plt.ylabel("var")
plt.legend(loc="best")
plt.xscale("log")
plt.savefig(pathData+"/varx0Ax1A.png")
plt.close()




