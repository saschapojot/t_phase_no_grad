from pathlib import Path



TVals=[0.001+0.001*n for n in range(0,5)]
print(TVals)
print(len(TVals))
outDir="./bash/"
Path(outDir).mkdir(parents=True,exist_ok=True)
for T in TVals:
    bashContents = []
    bashContents.append("#!/bin/bash\n")
    bashContents.append("#SBATCH -n 60\n")
    bashContents.append("#SBATCH -N 1\n")
    bashContents.append("#SBATCH -t 0-20:00\n")
    bashContents.append("#SBATCH -p hebhcnormal01\n")
    bashContents.append("#SBATCH --mem=180GB\n")
    bashContents.append("#SBATCH -o outT"+str(T)+"row0.out\n")
    bashContents.append("#SBATCH -e outT"+str(T)+"row0.err\n")
    bashContents.append("cd /public/home/hkust_jwliu_1/liuxi/Document/cppCode/t_phase_no_grad\n")
    bashContents.append("./runV1LJ2AtomPBC "+str(T)+" 0\n")
    fileName=outDir+"/runPBCV1LJ2AtomRow0T"+str(T)+".sh"
    with open(fileName,"w+") as fptr:
        for oneline in bashContents:
            fptr.write(oneline)
        fptr.close()
