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
pathData="../../version1Data/1d/run"+str(runNum)+".funcquadraticDistRelative/row0"
inCsv="../../version1Input/1d/distRelative/distRelative.csv"
dfStr=pd.read_csv(inCsv)
rowNum=0
oneRow=dfStr.iloc[rowNum,:]

a1=float(oneRow.loc["a1"])
a2=float(oneRow.loc["a2"])
c1=float(oneRow.loc["c1"])
c2=float(oneRow.loc["c2"])
# mA=float(oneRow.loc["mA"])
# mB=float(oneRow.loc["mB"])

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


def plt_dist(oneTFile):
    """

    :param oneTFile: corresponds to one temperature

    :return: y1Mean,y1Var,LMean,LVar,d0A0BMean,d1A1BMean
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))

    distFilePath=oneTFile+"/jsonData/jsondist/distData.json"

    with open (distFilePath,"r") as fptr:
        distData=json.load(fptr)

    y0Vec=distData["y0"]
    z0Vec=distData["z0"]
    y1Vec=distData["y1"]



    LFilePath=oneTFile+"/jsonData/jsonL/LData.json"
    with open(LFilePath,"r") as fptr:
        LData=json.load(fptr)

    LVec=LData["L"]

    LMean=np.mean(LVec)
    LVar=np.var(LVec,ddof=1)
    # LSd=np.sqrt(LVar)



    d0A0BVec=[]
    d0B1AVec=[]
    d1A1BVec=[]
    d1B0AVec=[]

    x0BVec=[]
    x1AVec=[]
    x1BVec=[]
    x0AVec=[]

    x0AInit=0
    for i in range(0,len(y0Vec)):
        LTmp=LVec[i]
        y0Tmp=y0Vec[i]
        z0Tmp=z0Vec[i]
        y1Tmp=y1Vec[i]

        d0A0BVec.append(y0Tmp)
        d0B1AVec.append(z0Tmp)
        d1A1BVec.append(y1Tmp)
        d1B0AVec.append(LTmp-y0Tmp-z0Tmp-y1Tmp)
        y1Vec.append(y1Tmp)

        x0BVec.append(x0AInit+y0Tmp)
        x1AVec.append(x0AInit+y0Tmp+z0Tmp)
        x1BVec.append(x0AInit+y0Tmp+z0Tmp+y1Tmp)
        x0AVec.append(x0AInit+LTmp)


    # #summary of distance between neighboring points
    d0A0BMean=np.mean(d0A0BVec)
    d0B1AMean=np.mean(d0B1AVec)
    d1A1BMean=np.mean(d1A1BVec)
    d1B0AMean=np.mean(d1B0AVec)

    d0A0BVar=np.var(d0A0BVec,ddof=1)
    d0B1AVar=np.var(d0B1AVec,ddof=1)
    d1A1BVar=np.var(d1A1BVec,ddof=1)
    d1B0AVar=np.var(d1B0AVec,ddof=1)

    d0A0BSd=np.sqrt(d0A0BVar)
    d0B1ASd=np.sqrt(d0B1AVar)
    d1A1BSd=np.sqrt(d1A1BVar)
    d1B0ASd=np.sqrt(d1B0AVar)












    y1Mean=np.mean(y1Vec)
    y1Var=np.var(y1Vec,ddof=1)





    smrDistFileName=oneTFile+"/distSummary.txt"
    #
    dist=[d0A0BMean,d0B1AMean,d1A1BMean,d1B0AMean]
    varsAll=[d0A0BVar,d0B1AVar,d1A1BVar,d1B0AVar]
    sdAll=[d0A0BSd,d0B1ASd,d1A1BSd,d1B0ASd]
    fptrTxt=open(smrDistFileName,"w")
    fptrTxt.write("T="+str(TVal)+"\n")
    fptrTxt.write("distances: "+str(dist)+"\n")
    fptrTxt.write("var: "+str(varsAll)+"\n")
    fptrTxt.write("sd: " +str(sdAll)+"\n")
    fptrTxt.close()

    rowNames=["0A0B","0B1A","1A1B","1B0A"]
    dfOut=pd.DataFrame({"dist":dist,"var":varsAll},index=rowNames)

    dfOut.to_csv(oneTFile+"/distVar.csv")

    x0AMean=np.mean(x0AVec)
    x0BMean=np.mean(x0BVec)
    x1AMean=np.mean(x1AVec)
    x1BMean=np.mean(x1BVec)

    x0AVar=np.var(x0AVec,ddof=1)
    x0BVar=np.var(x0BVec,ddof=1)
    x1AVar=np.var(x1AVec,ddof=1)
    x1BVar=np.var(x1BVec,ddof=1)

    x0ASd=np.sqrt(x0AVar)
    x0BSd=np.sqrt(x0BVar)
    x1ASd=np.sqrt(x1AVar)
    x1BSd=np.sqrt(x1BVar)


    pos_A=[x0AInit,x1AMean,x0AMean]
    sd_A=[0,x1ASd,x0ASd]

    pos_B=[x0BMean,x1BMean]
    sd_B=[x0BSd,x1BSd]

    # plt.figure(figsize=(12, 6))
    # plt.ylim(-1, 1)
    # for i in range(0,len(pos_A)):
    #     plt.hlines(y=0,xmin=pos_A[i]-sd_A[i],xmax=pos_A[i]+sd_A[i],color="red",linewidth=2,alpha=0.5)
    #     plt.text(pos_A[i],0.3,str(i)+"A",color="blue", ha='center')
    #     plt.text(pos_A[i],-0.1,"x"+str(i)+"A="+str(np.round(pos_A[i],4)),color="blue", ha='center',fontsize=8)
    #     plt.text(pos_A[i],0.1,"sd(x"+str(i)+"A)="+str(np.round(sd_A[i],4)),color="red", ha='center',fontsize=8)
    #
    #
    #
    # plt.scatter(pos_A,[0]*len(pos_A),color="blue",s=8,label="A")
    # for i in range(0,len(pos_B)):
    #     plt.hlines(y=0,xmin=pos_B[i]-sd_B[i],xmax=pos_B[i]+sd_B[i],color="magenta",linewidth=2,alpha=0.5)
    #     plt.text(pos_B[i],0.3,str(i)+"B",color="green", ha='center')
    #     plt.text(pos_B[i],-0.1,"x"+str(i)+"B="+str(np.round(pos_B[i],4)),color="green", ha='center',fontsize=8)
    #     plt.text(pos_B[i],0.1,"sd(x"+str(i)+"B)="+str(np.round(sd_B[i],4)),color="magenta", ha='center',fontsize=8)
    # plt.scatter(pos_B,[0]*len(pos_B),color="green",s=8,label="B")
    # plt.legend(loc="best")
    # plt.gca().get_yaxis().set_visible(False)
    # plt.axhline(y=0, color='black', linewidth=0.5,alpha=0.3)
    # # plt.xticks(unitCellPos,[0,1])
    # plt.xlabel("Unit cell")
    # gridOut="T"+str(TVal)+"grid.pdf"
    # plt.title("T="+str(TVal))
    # plt.savefig(oneTFile+"/"+gridOut)
    #
    # plt.close()












    return [y1Mean,y1Var,LMean,LVar,d0A0BMean,d1A1BMean]




UMeanValsAll=[]
UVarValsAll=[]
y1MeanValsAll=[]
y1VarValsAll=[]
LMeanValsAll=[]
LVarValsAll=[]
d0A0BMeanValsAll=[]
d1A1BMeanValsAll=[]
tStatsStart=datetime.now()
for oneTFile in sortedTFiles:
    UMeanTmp,UVarTmp=pltU(oneTFile)
    UMeanValsAll.append(UMeanTmp)
    UVarValsAll.append(UVarTmp)
    y1Mean,y1Var,LMean,LVar,d0A0BMean,d1A1BMean=plt_dist(oneTFile)
    y1MeanValsAll.append(y1Mean)
    y1VarValsAll.append(y1Var)
    LMeanValsAll.append(LMean)
    LVarValsAll.append(LVar)
    d0A0BMeanValsAll.append(d0A0BMean)
    d1A1BMeanValsAll.append(d1A1BMean)


tStatsEnd=datetime.now()
print("stats total time: ",tStatsEnd-tStatsStart)

interpolatedTVals=np.linspace(np.min(sortedTVals)*0.9,np.max(sortedTVals)*1.1,30)


def EV(T):
    '''

    :param T: temperature
    :return: asymptotic value of E(V)
    '''
    return 2*T
plt.figure()
EVVals=[EV(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,EVVals,color="green",label="theory")
plt.scatter(sortedTVals,UMeanValsAll,color="darkred",label="mc")
plt.title("E(V)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.savefig(pathData+"/EV.png")
plt.close()


def varV(T):
    """

    :param T: temperature
    :return:  asymptotic value of var(V)
    """
    return 2*T**2
plt.figure()
plt.scatter(sortedTVals,UVarValsAll,color="violet",label="mc")
varVVals=[varV(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,varVVals,color="navy",label="theory")
plt.title("var(V)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.savefig(pathData+"/varV.png")
plt.close()

plt.figure()
plt.scatter(sortedTVals,d0A0BMeanValsAll,color="deeppink",label="d0A0B")
plt.scatter(sortedTVals,d1A1BMeanValsAll,color="aqua",label="d1A1B",s=1)
plt.title("Intracell distance")
plt.xlabel("$T$")
plt.ylabel("distance")
plt.ylim((0.5,1.5))
plt.legend(loc="best")
plt.xscale("log")
plt.savefig(pathData+"/intracellDistance.png")
plt.close()











def EL(T):
    """

    :param T: temperature
    :return: asymptotic value of E(L)
    """
    return 2*a1+2*a2

def varL(T):
    """

    :param T: temperature
    :return: asymptotic value of var(L)
    """

    return  c1**(-1)*c2**(-1)*(c1+c2)*T


plt.figure()
plt.scatter(sortedTVals,LMeanValsAll,color="red",label="mc")
ELVals=[EL(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,ELVals,color="blue",label="theory")
plt.title("E(L)")
plt.ylabel("E(L)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(pathData+"/EL.png")
plt.close()

plt.figure()
plt.scatter(sortedTVals,LVarValsAll,color="magenta",label="mc")
varLVals=[varL(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,varLVals,color="green",label="theory")
plt.title("var(L)")
plt.ylabel("var(L)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(pathData+"/varL.png")
plt.close()


def Ey1(T):
    """

    :param T: temperature
    :return: asymptotic value of E(y1)
    """

    return a1



plt.figure()
plt.scatter(sortedTVals,y1MeanValsAll,color="red",label="mc")
Ey1Vals=[Ey1(T) for T in interpolatedTVals]

plt.plot(interpolatedTVals,Ey1Vals,color="blue",label="theory")

plt.title("E(y1)")
plt.ylabel("E(y1)")
plt.xlabel("$T$")
plt.legend(loc="best")

# plt.ylim((0,1.2))

plt.savefig(pathData+"/Ey1.png")
plt.close()


def vary1(T):
    """

    :param T: temperature
    :return: asymptotic value of var(y1)
    """


    val=1/4*c1**(-1)*(c1+c2)**(-1)*(2*c1+c2)*T\
        +1/4*c1**(-1)*c2*(c1+c2)**(-1)*T\


    return val



plt.figure()
plt.scatter(sortedTVals,y1VarValsAll,color="magenta",label="mc")
vary1Vals=[vary1(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,vary1Vals,color="green",label="theory")

plt.title("var(y1)")
plt.ylabel("var(y1)")
plt.xlabel("$T$")
plt.legend(loc="best")

# plt.ylim((0,1.5))
plt.savefig(pathData+"/vary1.png")
plt.close()


