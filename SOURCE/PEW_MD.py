#!/usr/bin/env python3

'''
DESCRIPTION
xyz (fs)
vel (A/fs)
grad(Ha/Bohr) -> amu*A/(fs)^2
'''

#
import numpy as np
import os
import subprocess as s
import sys
from PEW_IO import copy_2log, load_other_parm, load_parm, save_parm
from UPDATE import xyz_update, efield_update, vel_update, froz_grad

# Parameters for dynamics accepted from MQC_MD.sh script 
t_init = float(sys.argv[1])     # init. time of dynamics (femtoseconds)
t_step = float(sys.argv[2])     # timestep (femtoseconds)
t_finl = float(sys.argv[3])     # total simulation time (femtoseconds)
w_freq = float(sys.argv[4])     # frequency of wave (petahertz)
t_fwhm = float(sys.argv[5])     # FWHM of gauss. laser pulse (femtoseconds)
max_ef = float(sys.argv[6])     # max mag. of gaussian laser pulse (V/nm; GV/m)
algorm = str(sys.argv[7])       # ode integration algorithm
try:
    frozat =sys.argv[8].split(",")  # 
except:
    frozat=""
t_shft = t_fwhm + int(t_fwhm/2) # shift laser pulse to +/+ quadrant

# Conversion Factor for Gradient: Hartree/Bohr -> a.m.u.*A/fs^2
Ha2J    = 4.3597447222*10**(-18)/1  # val J/Ha
Bohr2m  = 1/(5.291722109*10**(-11)) # val Bohr/m
m2Ang   = 1/(1.000014*10**(-10))    # val Ang/m          
kg2amu  = 1/(1.6605402*10**(-27))   # val a.m.u./kg
s2fs    = 1*10**(-15)/1             # val s/fs
c_factor= Ha2J*Bohr2m*m2Ang*kg2amu*s2fs**2

# Conversion Factor for Efield: V/m -> a.u.
Vnm2au = 1/(5.1422082*10**2)        # val V/nm -> a.u.
max_ef = max_ef*Vnm2au

# For roudning dt after using np.arange with float input for step argument
num_dec = sys.argv[2][::-1].find('.')

# Import Required Environmental Variables
PARM = os.environ['PARM']               # Parameters    directory
ICON = os.environ['ICON']               # Initial cond. directory
TCON = os.environ['TCON']               # t_step  cond. directory
XYZ  = os.environ['XYZ']                # 3D coord.arr. directory
VEL  = os.environ['VEL']                # Velocity arr. directory
GRAD = os.environ['GRAD']               # Gradient arr. directory
PEWMD= os.environ['PEW_MD'] + '/SOURCE' # PEW_MD parent directory

# Define Filenames
mass_file = '/mass_array.mdl'           # Nuc. Mass   vector
xyz_file  = '/xyz_array.mdl'            # Nuc. XYZ    array
vel_file  = '/vel_array.mdl'            # Nuc. Vel.   array
grad_file = '/dV_current_array.mdl'     # PES Grad.  array

# Load Files
mass = load_parm(PARM, mass_file)                   # Nuc. Mass   vector
xyz = load_parm(TCON, xyz_file)                     # Nuc. XYZ    array
dV_current = load_parm(TCON, grad_file)*c_factor    # PES Grad.   array
vel, dV_bckward = load_other_parm(ICON,
                                  TCON,
                                  XYZ,
                                  VEL,
                                  GRAD,
                                  xyz_file,
                                  vel_file,
                                  grad_file,
                                  c_factor,
                                  xyz)
# Modify gradient so frozen atoms will not move
dV_current = froz_grad(dV_current, frozat)

# Start Pulsed Electromagnetic Wave Born-Oppenheimer Molecular Dynamics
for dt in np.round(np.arange(t_init+t_step, t_finl, t_step), num_dec):
    # Update Nuclear Coordinates (Beeman Algorithm)
    xyz = xyz_update(xyz, vel, dV_current, dV_bckward, mass, t_step, algorm)
    # Save Nuclear Coordinates at Timestep
    save_parm(TCON, xyz_file, xyz)
    # Calculate Electromagnetic Field Strength After dt Femtoseconds
    efield = efield_update(max_ef, t_fwhm, w_freq, t_shft, dt)
    # Pass laser and timestep to bashscript running gamess; execute
    s.call(PEWMD + "/GET_GRD.sh " + str(efield) + " " + str(dt), shell=True)  
    # Update gradient
    dV_forward = load_parm(TCON, grad_file)*c_factor
    # Modify gradient so frozen atoms will not move
    dV_forward = froz_grad(dV_forward, frozat)
    # Update Velocity (Beeman Algorithm)
    vel = vel_update(vel, dV_forward, dV_current, dV_bckward, mass, t_step, algorm)
    # Save vel at timestep
    save_parm(TCON, vel_file, vel)
    # Move initial conditions into directories for respective properties
    xyz_out = "/xyz_array" + "_" + str(dt) + ".mdl"
    vel_out = "/vel_array" + "_" + str(dt) + ".mdl"
    grad_out= "/grad_array" + "_" + str(dt) + ".mdl"
    copy_2log(TCON, xyz_file, XYZ, xyz_out)
    copy_2log(TCON, vel_file, VEL, vel_out)
    copy_2log(TCON, grad_file, GRAD, grad_out)
    # Prepare gradients for next iteration
    dV_bckward = dV_current
    dV_current = dV_forward
