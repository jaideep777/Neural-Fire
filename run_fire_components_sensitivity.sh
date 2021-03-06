#!/bin/bash

#VARS=( ltn  gpp   gppm1  pr  ts  cld  vp  pop  rdtot   ftmap11 )
#USEV=(   0    0       1   1   1    0   1    0      0         0 )

FOLDER=output_globe_sensitivity_3
MODEL=AUS_mod1464.6.2_gppm1s_gpp_gppl1_ts_cld_vp

#R1=4
#R2=5

# Get fire_dir (directory where the fire code - and this script - resides)
FIRE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## Generate a unique code number for the model from the variables used.
#MODNUM="${USEV[@]}" 			# join USEV array serially (joins with spaces)
#MODNUM=$((2#${MODNUM// /}))		# remove spaces so MODNUM becomes a binary number, then convert to decimal
#MODEL=${MODEL}_mod$MODNUM.5		# Append decimal number to model name. Thus each model gets unique name
#echo MODEL CODE=$MODEL

##################################################################
##	Update code, folder names etc as per specified variables    #
##################################################################

#N=${#VARS[@]}		# get length of variables array
#VARS_LIST=""		# Create a text list of varnames to be used by C++ nc2asc 
#XID="X_ids = ["		# Create a line of tensorflow code that defines training variables
#for ((i=0; i<$N; ++i)); do
##	echo ${VARS[$i]}
#	if ((${USEV[$i]} == 1)); then
#		VARS_LIST=${VARS_LIST}${VARS[$i]}"\n"
#		V=$(sed "s/_//g" <<< ${VARS[$i]})	# remove underscore from variable name (for folder naming)
#		MODEL=${MODEL}_$V					# append variable to model name
#		XID=${XID}"ID_${VARS[$i]}, " 		# append to the X_ids array for tensorflow code
#	fi
#done
echo MODEL = $MODEL
#XID=$(sed '$ s/.$//' <<< $XID)				# remove the trailing comma
#XID=${XID}"] + ID_ft" 						# Add closing bracket + ID_ft + ID_croplands

#XIDLINE=$(grep -rn "X\_ids \= " tensorflow/nn_const_data_fire_v5_pureNN.py | cut -f1 -d:)	# Get the line that defines X_ids in the tensorflow code file

#echo "Replace Line $XIDLINE:" 
#echo -e "\033[0;31m- " $(grep -r "X\_ids \= " tensorflow/nn_const_data_fire_v5_pureNN.py)
#echo -e "\033[0;32m+ " $XID "\033[0m"

#sed -i "${XIDLINE}s/.*/${XID}/" tensorflow/nn_const_data_fire_v5_pureNN.py	# Replace this line with the newly created X_ids 

#mkdir -p $FOLDER/$MODEL

#echo -e $VARS_LIST > $FOLDER/$MODEL/nn_vars.txt		# Output the list of vars to file (to be used by nc2asc)

#echo -e "r1 = $R1 \nr2 = $R2" > $FOLDER/$MODEL/regions.py

make 

##################################################################
##	Run fire components                                         #
##################################################################

## Aggregate training data
#./nc2asc train params_newdata/params_ip_global.r
#Rscript Rscripts/prepare_train_eval_datasets.R

# Train NN
#cd tensorflow
#. runtf 
#cd ..

# Run trained NN on data
# nc2asc eval syntax: ./nc2asc eval <params_file> <model_dir> <weights_file> <vars_file> <regions file>
./nc2asc eval params_newdata/params_ip_global_eval.r $FOLDER/$MODEL weights_ba.txt nn_vars.txt regions.py

# Plot results
cd Rscripts
Rscript plot_aggregate_maps_timeseries.R model_dir="$MODEL" output_dir="$FOLDER" fire_dir="${FIRE_DIR}"

#cd ..


