#!/usr/bin/bash

if [ -z "$1" ]; then
        echo "Enter name of raw output file as a CLA"
        exit 1
fi
echo "Splitting $1 ..."

#Check to make sure this file has job info (otherwise this job would take forever and I dont feel like dealing with it)
jobLines=`tail -12 $1 | grep "Number" | wc -l`
if [ $jobLines -ne 2 ]; then
	echo " Unable to find job info at the end of file"
	echo -e " (file is incomplete or GIVE_DETAILS = 0 in Constants.h)\nDone"
	exit 1
fi

#Get job info into variables
nConfigs_=(`tail -12 $1 | grep "Matrix Dimensions"`)
nConfigs=${nConfigs_[9]}
nAtoms_=(`tail -12 $1 | grep "Number of sites"`)
nAtoms=${nAtoms_[3]}
nEnvs_=(`tail -12 $1 | grep "Matrix Dimensions"`)
nEnvs=${nEnvs_[11]}

#Get correct line numbers for splitting from grep
#nAtoms*nEnvs+16 is the max number of lines we'll have to search
nTailLines=$(( $nAtoms*$nEnvs + 16))
grepOut=(`tail -$nTailLines $1 | grep -n "\-\-" | cut -d ':' -f1`)

#Check to make sure that there are four hyphen lines in the file otherwise this will break
nHyps=${#grepOut[@]}
if [ $nHyps -ne 4 ]; then
	echo " Was excecting 4 seperators (-) but only found $nHyps"
	echo -e " (the file is corrupted or GIVE_ENVS = 0 in Constants.h)\nDone"
	exit 2
fi

#Do the actual file splitting
##Matrix file
echo -n " Matrix ...... "
head -$(( $nConfigs )) $1 >$1".MAT"
echo "written to $1.MAT"
##Enviornment file
echo -n " Chem envs ... "
tail -$nTailLines $1 | tail -n +$(( ${grepOut[0]} + 1)) | head -n $(( ${grepOut[1]} - ${grepOut[0]} - 2)) >$1".ENV"
echo "written to $1.ENV"
##Job info file
echo -n " Job info .... "
tail -$nTailLines $1 | tail -n +$(( ${grepOut[2]} + 1)) | head -n $(( ${grepOut[3]} - ${grepOut[2]} - 1)) >$1".JOB"
echo "written to $1.JOB"

#Finish up
echo "Done"
