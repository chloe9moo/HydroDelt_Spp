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
RUNFILE="gradient_forests_r1.R"
BIOSITES="20_Traits/site_x_trait_fish_abundance_cont2mn.csv"
ENVVARS="include_baseline" #options: include_baseline, include_baseline_arch, include_alt
FLOWTYPES=("Int" "RO" "GW") 

mkdir ../working_dir
mkdir ../working_dir/20_Traits

#copy for run
cp gradient_forest_hpc.R "../working_dir/${RUNFILE}"
cp "../${BIOSITES}" "../working_dir/${BIOSITES}"
cp ../02_EnvDat/*.csv ../working_dir/

cd ../working_dir

#prep code for run
#sed -i 's|PATH <- getwd()|PATH <- "./working_dir"|' $RUNFILE #set working directory

sed -i "s|BIO_FILE_HERE|${BIOSITES}|" $RUNFILE #set bio type to run

sed -i "s/VAR_SELECTION_COL/${ENVVARS}/g" $RUNFILE #set which env variables to include

for ind in "${!FLOWTYPES[@]}"; do
	cp $RUNFILE "${ind}_${RUNFILE}"
	sed -i "s/FLOW_TYPE_HERE/${FLOWTYPES[$ind]}/" "${ind}_${RUNFILE}" #set which flow group to run
	
	#run R script for flow x site combo
	Rscript --vanilla --no-save ./"${ind}_${RUNFILE}"
	echo " "
	wait

	rm "${ind}_${RUNFILE}"
done


