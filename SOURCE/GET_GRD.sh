#!/bin/bash

###############################################################################
### DESCRTIPTION ##############################################################
###############################################################################
#
#
#
#

###############################################################################
### CREATE GAMESS INPUT FILE ##################################################
###############################################################################
# DEFINE ARGUMENTS PASSED FROM PYTHON
EFLD=$1 # EFIELD (a.u.)
T=$2 	# TIME (FEMTOSECONDS)
# APPEND EFIELD KEYWORD AND ITS STRENGTH FROM LASER TO ORIGINAL GAMESS-US INPUT FILE
echo -e ' $EFIELD\n'"   EVEC(3)=$EFLD"'\n   $END\n'"$(cat $ICON/*inp)" > $TCON/$(find $ICON/*.inp -type f -printf "%f\n")
# UPDATE COORDINATES IN INPUT
sed '/C1/q' $TCON/*.inp > $(find $TCON/*.inp).1
paste -d' ' $PARM/name_array.mdl $PARM/chrg_array.mdl $TCON/xyz_array.mdl > $(find $TCON/*.inp).2
cat $(find $TCON/*.inp).1 $(find $TCON/*.inp).2 > $(find $TCON/*.inp)
rm $(find $TCON/*.inp).1 $(find $TCON/*.inp).2
echo ' $END' >> $(find $TCON/*.inp)
sed -n '/ $VEC/,/ $END/p' $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp)) >> $(find $TCON/*.inp)

###############################################################################
### RUN GAMESS ################################################################
###############################################################################
cd $TCON
# CLEAN GAMESS-US .DAT FILE FOR CURRENT TIMESTEP
rm $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp))
# SOLVE TIME-DEPENDENT SCHRODINGER EQUATION WITH GAMESS FOR GRADIENT
rungms $(find *.inp) 00 ${SLURM_CPUS_PER_TASK} > $(sed 's/\.inp$/.gms/' <<< $(find *.inp))

###############################################################################
### REFRESH FILES FOR NEXT TIMESTEP ###########################################
###############################################################################
# RETRIEVE ENERGY GRADIENTS FROM GAMESS-US .DAT FILE. SAVE TO PARM DIRECTORY
awk '/GMAX/{ f=1;r=""; next}f && /END/{f=0} f{ r=(r=="")? $0: r RS $0 }\
END{ print r }' $SCRT/$(sed 's/\.inp$/.dat/' <<< $(find *.inp)) | 
awk '{ print $3, $4, $5}' > $TCON/dV_current_array.mdl
# MOVE GAMESS-US LOG TO LOGS DIRECTORY
mv $TCON/*.gms $LOGS/GMS/$(sed "s/\.gms$/_${T}.gms/" <<< $(find *.gms))
