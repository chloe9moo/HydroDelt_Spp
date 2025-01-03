#!/bin/bash
#SBATCH -J gfmodel
#SBATCH -o gfmodel_%j.txt
#SBATCH -p cloud72
#SBATCH --qos cloud
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=voorhees@uark.edu

#on hpc ONLY:
# source /etc/profile.d/lmod.sh
# module load gcc/9.3.1 mkl/19.0.5 R/4.2.2
# module list
# 
# cd $SLURM_SUBMIT_DIR #should be working_dir/
# cp -r * /scratch/$SLURM_JOB_ID #copy contents to scratch dir
# cd /scratch/$SLURM_JOB_ID

#variables
RUNFILE="gradient_forests_rINDEX_NUMBER"
BIOSITES="SET_BIO_FILE"
ENVVARS=(SET_ENV_VAR)
FLOWTYPES=(SET_FLOW_TYPE)

#prep code for run
#sed -i 's|PATH <- getwd()|PATH <- "./working_dir"|' $RUNFILE #set working directory (if necessary)

for flow in "${FLOWTYPES[@]}"; do
  for env in "${ENVVARS[@]}"; do
  
	cp gradient_forest_hpc.R "${RUNFILE}_${flow}_${env}.R"
	
	#set which site x bio file to run (same for all iterations)
	sed -i "s|BIO_FILE_HERE|${BIOSITES}|g" "${RUNFILE}_${flow}_${env}.R"
	#set which flow group to run
	sed -i "s/FLOW_TYPE_HERE/${flow}/g" "${RUNFILE}_${flow}_${env}.R"
	#set which env. variable set to run
	sed -i "s/VAR_SELECTION_COL/${env}/g" "${RUNFILE}_${flow}_${env}.R" 

	#run R script for flow x site combo
	Rscript --vanilla ./"${RUNFILE}_${flow}_${env}.R" >> "out_${RUNFILE}_${flow}_${env}.txt" 2>&1
	echo " "
	wait

	# rm "${flow}_${RUNFILE}"
	done
done

# mv *.R *.sh *.txt /scrfs/storage/voorhees/home/Documents/gf_run_outputs/run_scripts/
# mv output/* /scrfs/storage/voorhees/home/Documents/gf_run_outputs/

