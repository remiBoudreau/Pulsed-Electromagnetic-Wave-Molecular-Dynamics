#!/usr/bin/env python3

'''
DESCRIPTION
'''

# 
import numpy as np

# Calculate Electromagnetic Field Strength After dt Femtoseconds
def efield_update(max_ef, t_fwhm, w_freq, t_shft, dt):
    efield = np.abs(max_ef)*np.exp(-2*np.log(2)*(((dt-t_shft)/t_fwhm)**2))\
            *np.cos(w_freq*(dt-t_shft)) # Gauss. Envoloped Electromagnetic Wave
    return efield

def froz_grad(dV, frozat):
    for atom in frozat:
        dV[int(atom)] = [0, 0, 0]
    return dV

# Update Nuclear Coordinates
def xyz_update(xyz, vel, dV_current, dV_bckward, mass, t_step, algorm):
    # Beeman Algorithm
    if algorm == 'beeman':
        xyz = beeman_xyz(xyz, vel, dV_current, dV_bckward, mass, t_step)
    # Velocity-Verlet
    elif algorm == 'verlet':
        xyz = verlet_xyz(xyz, vel, dV_current, mass, t_step)
    return xyz

# Update Nuclear Velocities
def vel_update(vel, dV_forward, dV_current, dV_bckward, mass, t_step, algorm):
    # Beeman Algorithm
    if algorm == 'beeman':
        vel = beeman_vel(vel, dV_forward, dV_current, dV_bckward, mass, t_step)
    # Velocity Verlet
    elif algorm == 'verlet':
        vel = verlet_vel(vel, dV_forward, dV_current, mass, t_step)
    return vel

# Update Nuclear Coordinates Using Beeman Algorithm
def beeman_xyz(xyz, vel, dV_current, dV_bckward, mass, t_step):
    xyz = xyz + vel*t_step + (-1)*(1/6)*(4*dV_current - dV_bckward)*(t_step**2)/mass
    return xyz

# Update Nuclear Velocities Using Beeman Algorithm
def beeman_vel(vel, dV_forward, dV_current, dV_bckward, mass, t_step):
    vel = vel + (-1)*(1/6)*(2*dV_forward+5*dV_current-dV_bckward)*(t_step)/mass
    return vel

# Update Nuclear Coordinates Using Velocity-Verlet Algorithm
def verlet_xyz(xyz, vel, dV_current, mass, t_step):
    xyz = xyz + vel*t_step + (-1)*(1/2)*(dV_current)*(t_step**2)/mass
    return xyz

# Update Nuclear Velocities Using Velocity-Verlet Algorithm
def verlet_vel(vel, dV_forward, dV_current, mass, t_step):
    vel = vel + (-1)*(1/2)(dV_forward+dV_current)*(t_step)/mass
    return vel
