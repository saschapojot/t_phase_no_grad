from pathlib import Path

# TValsAll1=[0.001*n for n in range(1,10,2)]
#
# TValsAll2=[0.01*n for n in range(1,10,2)]
#
# TValsAll3=[0.1*n for n in range(1,6,2)]
#
# TVals=TValsAll1+TValsAll2+TValsAll3
TVals=[0.01,0.02,0.05,0.07,0.09,0.1,0.2,0.3,0.4,0.5,1,2,5]
print(TVals)
outDir="./bash/"
Path(outDir).mkdir(parents=True,exist_ok=True)
for T in TVals:
    bashContents = []
    bashContents.append("#!/bin/bash\n")
    bashContents.append("#SBATCH -n 24\n")
    bashContents.append("#SBATCH -N 1\n")
    bashContents.append("#SBATCH -t 0-20:00\n")
    bashContents.append("#SBATCH -p hebhcnormal01\n")
    bashContents.append("#SBATCH --mem=20GB\n")
    bashContents.append("#SBATCH -o outT"+str(T)+"row0.out\n")
    bashContents.append("#SBATCH -e outT"+str(T)+"row0.err\n")
    bashContents.append("cd /public/home/hkust_jwliu_1/liuxi/Document/cppCode/t_phase_no_grad\n")
    bashContents.append("./runV1CartesianQuadratic "+str(T)+" 0\n")
    fileName=outDir+"/runV1CartesianQuadraticT"+str(T)+".sh"
    with open(fileName,"w+") as fptr:
        for oneline in bashContents:
            fptr.write(oneline)
        fptr.close()
