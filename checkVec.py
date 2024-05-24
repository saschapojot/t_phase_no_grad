import pickle
import numpy as np
# from datetime import datetime
import statsmodels.api as sm
# from scipy import stats
import glob
import sys
import re
from copy import deepcopy
import warnings
from scipy.stats import ks_2samp


#This script checks if  U values reach equilibrium
sigWrongArg="wrong number of arguments"
sigEq="equilibrium"
sigContinue="continue"
sigStop="stop"


if (len(sys.argv)!=2):
    print(sigWrongArg)
    exit()

pklFile=str(sys.argv[1])

with open(pklFile,"rb") as fptr:
    vec=np.array(pickle.load(fptr))


startingLoop=3000*20000

vecTruncated=vec[-startingLoop:]

NLags=int(np.ceil(len(vecTruncated)*5/6))
eps=1e-3
lagVal=0
same=False
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(vecTruncated)
    except Warning as w:
        same=True


if same==True:
    print(sigStop+" same")
    exit()

else:
    if np.min(np.abs(acfOfVec))>eps:
        print("high correlation")
        # print(np.min(np.abs(acfOfVec)))

        print(sigContinue)
        exit()
    else:
        lagVal=np.where(np.abs(acfOfVec)<=eps)[0][0]
        vecValsSelected=vecTruncated[::lagVal]
        lengthTmp=len(vecValsSelected)
        if lengthTmp%2==1:
            lengthTmp-=1
        vecValsToCompute=vecValsSelected[:lengthTmp]
        lenPart=int(len(vecValsToCompute)/2)
        selectedFromPart0=vecValsToCompute[:lenPart]
        selectedFromPart1=vecValsToCompute[lenPart:]
        result = ks_2samp(selectedFromPart0, selectedFromPart1)
        numDataPoints=len(selectedFromPart0)+len(selectedFromPart1)
        print("len(selectedFromPart0)="+str(len(selectedFromPart0)))
        print("len(selectedFromPart1)="+str(len(selectedFromPart1)))
        msg="lag="+str(lagVal)+"\n"+"K-S statistic: "+str(result.statistic)+"\n"+"P-value: "+str(result.pvalue)+"\n"\
            +"numDataPoints="+str(numDataPoints)+"\n"+"nCounterStart="+str(len(vec)-startingLoop)+"\n"
        if result.pvalue>0.1:
            print(sigEq)
            print(msg)
            exit()
        else:
            print(sigContinue)
            exit()



