#!/bin/bash
#SBATCH -n 8   #Number of cores
#SBATCH -t 40    #Runtime in minutes -t D-HH:MM
#SBATCH -p shared #Partition to submit to 
#SBATCH --mem-per-cpu=2000 #Memory per cpu in MB (see also --mem-per-cpu for memory per cpu in MB; --mem=20000 for #Total memory in MB, BF: 4000)
#SBATCH -o 22_05_31_8core2000MB40min.out    # Standard out goes to this file
#SBATCH -e 22_05_31_8core2000MB40min.err    # Standard err goes to this file
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=erikknall@g.harvard.edu  #Email to which notifications will be sent
# load new modules system
source new-modules.sh

module load intel/21.2.0-fasrc01
module load openmpi/4.1.1-fasrc01
module load lumerical-seas/2021_7bf43e7149-fasrc01


# list your modules to make sure they are there:
module list
export PWDD=`pwd`
for fname in `ls -1v ${PWDD}/*.fsp`
do
    export FullFILENAME=$fname
    
    srun -n $SLURM_NTASKS --mpi=pmix fdtd-engine-ompi-lcl -fullinfo ${FullFILENAME}
     
    date
    
    echo "done with ${FullFILENAME}"
done
exit
