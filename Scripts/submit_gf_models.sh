#!/bin/bash
#### handle multiple simultaneous runs of gf models ####

#set variables
#comment out whatever site x spp files don't want run
BIOSITES_FILES=(
  #traits
  # "20_Traits/site_x_trait_fish_abundance_cont2mn.csv"
  "20_Traits/site_x_trait_fish_abundance_cont2cat.csv"
  "20_Traits/site_x_trait_fish_abundance_clust.csv"
  # "20_Traits/site_x_trait_bug_abundance_cont2mn.csv"
  # "20_Traits/site_x_trait_bug_abundance_cont2cat.csv"
  # "20_Traits/site_x_trait_bug_abundance_clust.csv"
  #taxonomic
  "01_BioDat/occ_fish_15k_wide_20240412.csv"
  "01_BioDat/occ_bug_15k_wide_20240412.csv"
  )
#changing these here so that the sbatch shell script doesn't need to be edited directly
ENVVARS='"include_baseline" "include_baseline_arch" "include_alt" "include_hydro"' #columns from environmental_variable_info.csv
FLOWTYPES='"Int" "RO" "GW"'

#make working directory
mkdir ../working_dir
mkdir ../working_dir/01_BioDat
mkdir ../working_dir/20_Traits
mkdir ../working_dir/output

#move over datasets (all)
cp ../01_BioDat/*.csv ../working_dir/01_BioDat
cp ../20_Traits/*.csv ../working_dir/20_Traits
cp ../02_EnvDat/*.csv ../working_dir/
cp gradient_forest_hpc.R run_gf_models.sh ../working_dir/

cd ../working_dir

#set variable types in template shell script
sed -i "s/SET_ENV_VAR/${ENVVARS}/" run_gf_models.sh
sed -i "s/SET_FLOW_TYPE/${FLOWTYPES}/" run_gf_models.sh

for ind in "${!BIOSITES_FILES[@]}"; do
  
  RUNFILE="gf_run${ind}.sh"
  
  #copy for individual runs
  cp run_gf_models.sh $RUNFILE
  
  #edit site file to run in shell script
  sed -i "s|SET_BIO_FILE|${BIOSITES_FILES[$ind]}|" $RUNFILE
  sed -i "s/INDEX_NUMBER/${ind}/" $RUNFILE

	#submit sbatch job for each site file
	./$RUNFILE

	# echo " "
	wait

	# rm "${ind}_${RUNFILE}"
done


