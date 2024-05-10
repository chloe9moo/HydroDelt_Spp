#!/bin/bash
##SBATCH -N 1
##SBATCH --ntasks=1
##SBATCH --mem-per-cpu=80G
##SBATCH -t 08:00:00
##SBATCH -J env_filt
##SBATCH -p normal_q
##SBATCH --account=usgs_rcs
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=chloe9mo@vt.edu

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
	# Rscript --vanilla --no-save ./"${flow}_${RUNFILE}"
	# echo " "
	wait

	# rm "${flow}_${RUNFILE}"
	done
done


