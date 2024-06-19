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
inCsv="../../version1Input/1d/cartesian/cartesianQuadratic.csv"
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


def plt_x(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: plots of positions, var(x0A), E(x0A), d0A0BMean,d1A1BMean, x0AVar,x1AVar, y1Mean, y1Var LMean,LVar
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

    LMean=np.mean(LVec)
    LVar=np.var(LVec,ddof=1)
    # LSd=np.sqrt(LVar)



    d0A0BVec=[]
    d0B1AVec=[]
    d1A1BVec=[]
    d1B0AVec=[]
    y1Vec=[]
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
        y1Vec.append(x1BTmp-x1ATmp)

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

    y1Mean=np.mean(y1Vec)
    y1Var=np.var(y1Vec,ddof=1)





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
    return [x0AVar,x0AMean,d0A0BMean,d1A1BMean,x0AVar,x1AVar,y1Mean,y1Var,LMean,LVar]






xAVarAll=[]
xAMeanAll=[]
UMeanAll=[]
UVarAll=[]
d0A0BAll=[]
d1A1BAll=[]
varx0AAll=[]
varx1AAll=[]
y1MeanAll=[]
y1VarAll=[]

LMeanAll=[]
LVarAll=[]
tStatsStart=datetime.now()
for oneTFile in sortedTFiles:
    UMeanTmp,UVarTmp=pltU(oneTFile)
    UMeanAll.append(UMeanTmp)
    UVarAll.append(UVarTmp)
    varTmp,meanTmp,d0A0BTmp,d1A1BTmp,x0AVar,x1AVar,y1Mean,y1Var,LMean,LVar=plt_x(oneTFile)
    xAVarAll.append(varTmp)
    xAMeanAll.append(meanTmp)
    d0A0BAll.append(d0A0BTmp)
    d1A1BAll.append(d1A1BTmp)
    varx0AAll.append(x0AVar)
    varx1AAll.append(x1AVar)
    y1MeanAll.append(y1Mean)
    y1VarAll.append(y1Var)
    LMeanAll.append(LMean)
    LVarAll.append(LVar)


tStatsEnd=datetime.now()
print("stats total time: ",tStatsEnd-tStatsStart)
interpolatedTVals=np.linspace(0.9*np.min(sortedTVals),1.1*np.max(sortedTVals),30)
def varx0A(T):
    """

    :param T: temperature
    :return: asymptotic value of var(x0A)
    """
    return 1/2*c1**(-1)*c2**(-1)*(c1+c2)+1/3*(a1+a2)**2-1/16*c1**(-2)*c2**(-2)*((c1+c2)/(a1+a2))**2*T**2


plt.figure()

plt.scatter(sortedTVals,xAVarAll,color="red",label="mc")
varx0AVals=[varx0A(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,varx0AVals,color="black",label="theory")
plt.title("var(x0A)")
plt.ylabel("var(x0A)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,0.5))
plt.savefig(pathData+"/varx0A.png")
plt.close()

def Ex0A(T):
    """

    :param T: temperature
    :return: asymptotic value of E(x0A)
    """
    return 1/4*c1**(-1)*c2**(-1)*(c1+c2)/(a1+a2)*T+a1+a2


plt.figure()

plt.scatter(sortedTVals,xAMeanAll,color="green",label="mc")
Ex0AVals=[Ex0A(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,Ex0AVals,color="blue",label="theory")
plt.title("E(x0A)")
yTicks=[0.1*n for n in range(0,11)]
# plt.yticks(yTicks)
plt.ylabel("E(x0A)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,1))
plt.savefig(pathData+"/Ex0A.png")
plt.close()

def EV(T):
    '''

    :param T: temperature
    :return: asymptotic value of E(V)
    '''
    return 2*T
plt.figure()
EVVals=[EV(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,EVVals,color="green",label="theory")
plt.scatter(sortedTVals,UMeanAll,color="darkred",label="mc")
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
plt.scatter(sortedTVals,UVarAll,color="violet",label="mc")
varVVals=[varV(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,varVVals,color="navy",label="theory")
plt.title("var(V)")
plt.xlabel("$T$")
plt.legend(loc="best")
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
plt.ylim((0,0.5))
plt.xscale("log")
plt.savefig(pathData+"/varx0Ax1A.png")
plt.close()









def EL(T):
    """

    :param T: temperature
    :return: asymptotic value of E(L)
    """
    return 1/2*c1**(-1)*c2**(-1)*(c1+c2)/(a1+a2)*T+2*a1+2*a2

def varL(T):
    """

    :param T: temperature
    :return: asymptotic value of var(L)
    """

    return -1/4*c1**(-2)*c2**(-2)*((c1+c2)/(a1+a2))**2*T**2\
            +c1**(-1)*c2**(-1)*(c1+c2)*T


plt.figure()
plt.scatter(sortedTVals,LMeanAll,color="red",label="mc")
ELVals=[EL(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,ELVals,color="blue",label="theory")
plt.title("E(L)")
plt.ylabel("E(L)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.ylim((0,5.5))
plt.savefig(pathData+"/EL.png")
plt.close()

plt.figure()
plt.scatter(sortedTVals,LVarAll,color="magenta",label="mc")
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
    val=1/4*c1**(-1)*(a1+a2)**(-1)*T+a1
    return val


plt.figure()
plt.scatter(sortedTVals,y1MeanAll,color="red",label="mc")
Ey1Vals=[Ey1(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,Ey1Vals,color="black",label="theory")
plt.title("E(y1)")
plt.ylabel("E(y1)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.ylim((0,1.2))
plt.savefig(pathData+"/Ey1.png")
plt.close()


def vary1(T):
    """

    :param T: temperature
    :return: asymptotic value of var(y1)
    """
    val=-1/16*c1**(-2)*(a1+a2)**(-2)*T**2\
        +1/4*c1**(-1)*(c1+c2)**(-1)*(2*c1+c2)*T\
        +1/4*c1**(-1)*c2*(c1+c2)**(-1)*T
    return val


plt.figure()
plt.scatter(sortedTVals,y1VarAll,color="magenta",label="mc")

vary1Vals=[vary1(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,vary1Vals,color="green",label="theory")

plt.title("var(y1)")
plt.ylabel("var(y1)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,1.2))
plt.savefig(pathData+"/vary1.png")
plt.close()
