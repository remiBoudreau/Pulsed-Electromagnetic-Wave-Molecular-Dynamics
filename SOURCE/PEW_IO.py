#!/usr/bin/env python3

'''
DESCRIPTION
'''

#
import numpy as np
import subprocess as s

#
def copy_2log(input_directory, input_file, output_directory, output_file):
    # Move initial conditions into directories for respective properties
    s.call("cp " + input_directory + input_file + " " + output_directory + output_file,
                    shell=True)

#
def load_other_parm(ICON, TCON, XYZ, VEL, GRAD, xyz_file, vel_file, grad_file, 
                    c_factor, xyz):
    bck_grad_file = '/dV_bckward_array.mdl'                    # PES Grad. array(dt-1)    
    try:
        # Restart
        # Load Files
        vel = load_parm(TCON, vel_file)                        # Nuc. Vel. array     
        dV_bckward = load_parm(TCON, bck_grad_file)*c_factor   # PES Grad. array (dt)
        s.call("rm " + TCON + bck_grad_file, shell=True)       # Rem.Grad. array(dt-1)
    except:
        # Fresh Run
        # Load Files
        vel = np.zeros_like(xyz)                               # Nuc. Vel. array
        dV_bckward = np.zeros_like(xyz)                        # PES Grad. array
        # Save Files
        save_parm(ICON, vel_file, vel)                         # Nuc. Vel. array
        save_parm(ICON, bck_grad_file, dV_bckward)             # PES Grad. array
        # Copy initial conditions to 'LOGS' Folder
        copy_2log(ICON, xyz_file, XYZ, '/xyz_array_0.0.mdl')   # Nuc. XYZ  array
        copy_2log(ICON, vel_file, VEL, '/vel_array_0.0.mdl')   # Nuc. Vel. array
        copy_2log(ICON, grad_file, GRAD, '/grad_array_0.0.mdl')# PES Grad. array
    return vel, dV_bckward

#
def load_parm(directory, filename):
    parm = np.loadtxt(directory + filename,
                      delimiter=" ",
                      ndmin=2)
    return parm

#
def save_parm(directory, filename, array):
    np.savetxt(directory + filename,
               array,
               delimiter=" ",
               fmt='%2.15f')