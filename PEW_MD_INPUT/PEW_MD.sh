#!/bin/bash

#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3900M
#SBATCH --time=7-00:00
#SBATCH --job-name PEW_MD
#SBATCH --output=PEW_MD.txt
#SBATCH --mail-type=BEGIN,END,FAIL 

###############################################################################
### DESCRTIPTION ##############################################################
###############################################################################
# This set of scripts comprise an addon module for performing born-oppenheimer
# molecular dynamics in the presence of a time-varying electric field at the 
# tddft level of theory by interfacing with the GAMESS-US quantum chemistry
# package. The initial GAMESS-US input file must be colocated with this file
# in the 'INPUT' directory. This can be run either locally or on a SLURM 
# workload manager (e.g. on ComputeCanada HPCs). GAMESS-US must be available
# by the command 'rungms' for proper functionality of these scripts. The 
# parameters immediately in the block below are to be changed as necessary
# for your dynamics simulation. An example GAMESS-US input for use with this
# program is provided colocated with this script.


###############################################################################
### USER DEFINED PARAMETERS FOR MQC_MD ########################################
###############################################################################
export rstart=2	# dynamics restart variable (1: NORMAL RUN; 2: RESTART; last t)
t_step=1.0		# timestep (in femtoseconds)
t_finl=301    	# total simulation time (in femtoseconds)
w_freq=0.56		# frequency of wave (in PHz)
t_fwhm=100     	# full-width half-max of gaussian laser pulse (in femtoseconds)
max_ef=2.5 		# max magnitude of gaussian laser pulse (in V/nm) 
algorm=beeman   # numerical integration algorithm  Beeman: 'beeman'
				#								   Velcity Verlet:'verlet'
export frozat="0,1,8,9,14,15,18,19,28,31,40,47,49"
###############################################################################
### SET UP ENVIRONMENT ########################################################
###############################################################################
echo "SETTING UP ENVIRONMENT"
# PATH TO MQC_MD DIRECTORY
export PEW_MD=$SCRATCH/PEW_MD_BRZ
# RUNNING ON SMALL WORKSTATION/HPC WORKSTATION WITHOUT SLURM
if [ -z "$SLURM_SUBMIT_DIR" ]; then
	export SLURM_SUBMIT_DIR=$PWD
# RUNNING ON HPC WORKSTATION WITH SLURM
else
	# DIRECTORY TO RUN - $SLURM_SUBMIT_DIR is directory job was submitted from
	cd $SLURM_SUBMIT_DIR
	# MODULES REQUIRED FOR GAMESS
	module load nixpkgs/16.09
	module load intel/2018.3	# INTEL COMPILERS   (DEPENDENCY)
	module load intelmpi/2018.3.222	# INTEL MPI LIBRARY (DEPENDENCY)
	module load gamess-us/20180920-R3 			# GAMESS-US
	module load python/3.7.4
	module load scipy-stack
fi
# ENVIRONMENTAL VARIABLES
export MODL=$PEW_MD/MODEL 			# TOP-PARENT DIRECTORY		
export PARM=$MODL/PARM 				# PARAMETERS DIRECTORY 
export LOGS=$MODL/LOGS 				# GAMESS LOG DIRECTORY
export ICON=$MODL/ICON				#
export TCON=$MODL/TCON				#
export  XYZ=$LOGS/XYZ
export  VEL=$LOGS/VEL
export GRAD=$LOGS/GRAD
export PATH=$HOME/gamess:$PATH		# CUSTOM RUNGMS SCRIPT
export SCRT=$SCRATCH/gamess-scratch    # DICTIONARY FILE LOC.
# IF NOT RESTARTING, REMOVE ALL FILES IN DIRECTORY
if [ $rstart -eq 1 ]; then
	rm -r $MODL
fi
# CREATE DIRECTORES FOR FILES
mkdir -p $MODL
mkdir -p $PARM
mkdir -p $LOGS
mkdir -p $ICON
mkdir -p $TCON
mkdir -p $XYZ
mkdir -p $VEL
mkdir -p $GRAD
mkdir -p $LOGS/GMS
echo "DONE"
# VARIABLE TO CONTROL RESTART
if [ $rstart -eq 1 ]; then
	###########################################################################
	### INITIAL CONDITIONS FROM .DAT FILE OF GAMESS GRADIENT RUN ##############
	###########################################################################
	echo "GAMESS-US GRADIENT RUN"
	rm $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp))
	rungms $(find *.inp) 00 ${SLURM_CPUS_PER_TASK} > $(sed 's/\.inp$/.gms/' <<< $(find *.inp))
	mv $(sed 's/\.inp$/.gms/' <<< $(find *.inp)) $LOGS/GMS/$(sed 's/\.inp$/_0.0.gms/' <<< $(find *.inp))
	echo "DONE"
	echo \
	"RETRIEVING NUCLEAR COORDINATES, CHARGES, MASSES & GRAD. FROM GAMESS .DAT FILE"
	# RETRIEVE ATOMS NAMES FROM GAMESS-US .DAT FILE. SAVE TO PARM DIRECTORY
	awk '/GMAX/{ f=1;r=""; next}f && /END/{f=0} f{ r=(r=="")? $0: r RS $0 }\
	END{ print r }' $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp)) | 
	awk '{ print $1}' > $PARM/name_array.mdl
	# RETRIEVE CHARGES FROM GAMESS-US .DAT FILE. SAVE TO PARM DIRECTORY
	awk '/GMAX/{ f=1;r=""; next}f && /END/{f=0} f{ r=(r=="")? $0: r RS $0 }\
	END{ print r }' $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp)) | 
	awk '{ print $2}' > $PARM/chrg_array.mdl
	# RETRIEVE MASSES OF EACH ATOM FROM MASSES.txt. SAVE TO PARM DIRECTORY
	while read atom; do
		grep -w $atom $PEW_MD/SOURCE/MASSES.txt | awk '{ print $2}'\
		>> $PARM/mass_array.mdl
	done < $PARM/name_array.mdl
	# RETRIEVE NUCLEAR COORDINATES FROM GAMESS-US .DAT FILE. SAVE TO PARM DIRECTORY
	awk '/C1/{ f=1;r=""; next}f && /END/{f=0} f{ r=(r=="")? $0: r RS $0 }\
	END{ print r }' $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp)) | 
	awk '{ print $3, $4, $5 }' | grep '^[[:blank:]]*[^[:blank:]#;]'\
	> $ICON/xyz_array.mdl
	# RETRIEVE ENERGY GRADIENTS FROM GAMESS-US .DAT FILE. SAVE TO PARM DIRECTORY
	awk '/GMAX/{ f=1;r=""; next}f && /END/{f=0} f{ r=(r=="")? $0: r RS $0 }\
	END{ print r }' $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp)) | 
	awk '{ print $3, $4, $5}' > $ICON/dV_current_array.mdl
	echo "DONE"
	# SAVE INITIAL CONDITIONS FOR RESTARTING FROM BEGINNING
	cp *.inp $ICON
	# SET RUNTYP KEYWORD OF GAMESSFILE TO GRADIENT
	RUNTYP=$(grep -Po 'RUNTYP=\K[^ ]+' $(find $ICON/*.inp));
	sed "s/$RUNTYP/GRADIENT/" *.inp > $(find $ICON/*.inp)
	# SET GUESS KEYWORD OF GAMESSFILE TO MOREAD
	GUESS=$(grep -Po 'GUESS=\K[^ ]+' $(find $ICON/*.inp));
	sed "s/$GUESS/MOREAD/" *.inp > $(find $ICON/*.inp)
	# START DYNAMICS FROM 0 FEMTOSECONDS
	export t_init=0	
elif [ $rstart -eq 2 ]; then
	# COPY FILES FOR LAST SUCCESSFUL TIMESTEP INTO TCON
	cp $(find $XYZ -name *$(ls $XYZ | sed "s|xyz_array_||" | sed s/.mdl// | sort -n | tail -1)*) $TCON/xyz_array.mdl
	cp $(find $VEL -name *$(ls $VEL | sed "s|vel_array_||" | sed s/.mdl// | sort -n | tail -1)*) $TCON/vel_array.mdl
	cp $(find $GRAD -name *$(ls $GRAD | sed "s|grad_array_||" | sed s/.mdl// | sort -n | tail -1)*) $TCON/dV_current_array.mdl
	cp $(find $GRAD -name *$(ls $GRAD | sed "s|grad_array_||" | sed s/.mdl// | sort -n | tail -2 | head -1)*) $TCON/dV_bckward_array.mdl
	# START DYAMICS FROM LAST SUCCESSFUL TIMESTEP
	export t_init=$(ls $XYZ | sed "s|xyz_array_||" | sed s/.mdl// | sort -n | tail -1)
else 
	echo "OPTION FOR RSTART VALUE IS UNDEFINED"
	exit 1
fi

###############################################################################
### STARTING PEW MD RUN #######################################################
###############################################################################
# START MIXED QUANTUM-CLASSICAL DYNAMICS
echo "STARTING PEW MD RUN FROM $t_init FEMTOSECONDS"
python $PEW_MD/SOURCE/PEW_MD.py $t_init $t_step $t_finl $w_freq $t_fwhm\
				$max_ef $algorm $frozat
echo "DONE"
