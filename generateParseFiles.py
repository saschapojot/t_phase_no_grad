from pathlib import Path


TIndsAll=[n for n in range(0,11)]
outDir="./parseBash/"
Path(outDir).mkdir(parents=True,exist_ok=True)
for ind in TIndsAll:
    bashContents = []
    bashContents.append("#!/bin/bash\n")
    bashContents.append("#SBATCH -n 60\n")
    bashContents.append("#SBATCH -N 1\n")
    bashContents.append("#SBATCH -t 0-20:00\n")
    bashContents.append("#SBATCH -p hebhcnormal01\n")
    bashContents.append("#SBATCH --mem=100GB\n")
    bashContents.append("#SBATCH -o outParse"+str(ind)+".out\n")
    bashContents.append("#SBATCH -e outParse"+str(ind)+".err\n")
    bashContents.append("cd /public/home/hkust_jwliu_1/liuxi/Document/cppCode/t_phase_no_grad\n")
    bashContents.append("./runV1LJ2AtomParseBinPBC  0 "+str(ind)+"\n")
    fileName=outDir+"/parse"+str(ind)+".sh"
    with open(fileName,"w+") as fptr:
        for oneline in bashContents:
            fptr.write(oneline)
        fptr.close()
