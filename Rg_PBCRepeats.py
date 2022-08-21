# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 13:24:39 2022

@author: avtei
"""

import glob, os, sys, pandas, copy, time
import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.analysis.base import AnalysisFromFunction
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

def Rg_calc(coordinates):
    center_of_mass = coordinates.mean(axis=0)
    masses = np.ones((coordinates.shape[0]))
    total_mass = masses.sum()
    ri_sq = (coordinates-center_of_mass)**2
    sq = np.sum(ri_sq, axis=1)
    sq_x = np.sum(ri_sq[:,[1,2]], axis=1) # sum over y and z
    sq_y = np.sum(ri_sq[:,[0,2]], axis=1) # sum over x and z
    sq_z = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y
    sq_rs = np.array([sq, sq_x, sq_y, sq_z])
    rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
    return(np.sqrt(rog_sq))

def calc_Rg(gro, traj = None, PBCRepeats = 3, selection = "not resname H2O and not resname W and not resname CL and not resname NA and not resname ION"):
    if type(traj) == type(None):
        U = Universe(gro)
    else:
        U = Universe(gro, traj)
    Rgs = np.ndarray((0,4))
    for iframe, frame in enumerate(U.trajectory):
        box = frame.dimensions[:3]
        protein = U.select_atoms(selection)
        positions = np.ndarray((0,3))

        for PBC_X in np.arange(PBCRepeats):
            for PBC_Y in np.arange(PBCRepeats):
                for PBC_Z in np.arange(PBCRepeats):
                    coords = protein.positions.copy()
                    coords[:,0] += (box[0]*PBC_X)
                    coords[:,1] += (box[1]*PBC_Y)
                    coords[:,2] += (box[2]*PBC_Z)
                    positions = np.vstack((positions, coords))
        Rg = Rg_calc(positions)
        Rgs = np.vstack((Rgs, Rg/PBCRepeats))
    return Rgs

if __name__ == "__main__":
    #gro = sys.argv[1]
    #xtc = sys.argv[2]
    #fname = sys.argv[3]

    PBCRepeats = 3
    skip = 1
    
    top = "."
    folders = [x for x in os.listdir(top) if os.path.isdir(f"{top}/{x}") and ".0_" in x and "MD" not in x]
    
    print(folders)

    for folder in folders:
        gro = f"{top}/{folder}/STRUCTURE.gro"
        xtc = f"{top}/{folder}/TRAJECTORY.xtc"
        fname = f"{top}/{folder}/Rg_PBCRepeat.dat"
        
        if os.path.exists(fname):
            print("Found:", fname)
            continue
        
        
        if not os.path.exists(gro):
            gro = f"{top}/{folder}/GMXSASA.gro"
            xtc = f"{top}/{folder}/GMXSASA.xtc"
        #box = f"{top}/{folder}/cphmd_out.restart.xsc"
        #box = np.diag(np.split(np.loadtxt(box)[1:], 6)[:3])
        pH = top.split("pH")[-1].split("/")[0].split("_")[0]
        
        Rgs = np.ndarray((0,4))
    
        U = Universe(gro, xtc)
        for iframe, frame in enumerate(U.trajectory[::skip]):
            box = frame.dimensions[:3]
            selection = "not resname H2O and not resname W and not resname CL and not resname NA and not resname ION"
            protein = U.select_atoms(selection)
            
            #protein.positions
            
            positions = np.ndarray((0,3))
            
    
            for PBC_X in np.arange(PBCRepeats):
                for PBC_Y in np.arange(PBCRepeats):
                    for PBC_Z in np.arange(PBCRepeats):
                        coords = protein.positions.copy()
                        
                        coords[:,0] += (box[0]*PBC_X)
                        coords[:,1] += (box[1]*PBC_Y)
                        coords[:,2] += (box[2]*PBC_Z)
                            
                        positions = np.vstack((positions, coords))
                        
                            
            Rg = Rg_calc(positions)
            print(folder, iframe*skip, "/", len(U.trajectory), "PBCRepeats:", PBCRepeats, "Rg:", Rg, "Rg/PBCRepeats:", Rg/PBCRepeats)
            Rgs = np.vstack((Rgs, Rg/PBCRepeats))
        np.savetxt(fname, Rgs)