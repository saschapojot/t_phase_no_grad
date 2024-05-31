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





# N=10
#
# T=1
#
# den=N*dV(x*2)
# print(0.1*T/den)
# rValsAll=np.linspace(0.5,1.2,100)
# dV1ValsAll=dV1(rValsAll)
# dV2ValsAll=dV2(rValsAll)
# plt.figure()
# plt.plot(rValsAll,dV1ValsAll,color="black")
# plt.title("$dV_{1}$")
# plt.savefig("tmpdV1.png")
# plt.close()
#
# plt.figure()
# plt.plot(rValsAll,dV2ValsAll,color="blue")
# plt.title("$dV_{2}$")
# plt.savefig("tmpdV2.png")
# plt.close()

# sol1=root(dV1,1,method="broyden2",tol=1e-9)
#
# sol2=root(dV2,1,method="broyden2",tol=1e-9)
#
# print(sol1)
# print(sol2)

def V1(r):
    return alpha1/r**p1-beta1/r**q1+r**4

def V2(r):
    return alpha2/r**p2-beta2/r**q2+r**4

def V(r):
    return V1(r)+V2(r)
# rValsAll=np.linspace(0.5*x,2*x,50)
# VValsAll=[V(r) for r in rValsAll]
#
# import matplotlib.pyplot as plt
#
# plt.figure()
# plt.plot(rValsAll,VValsAll,color="black")
# plt.savefig("tmp.png")
# plt.close()