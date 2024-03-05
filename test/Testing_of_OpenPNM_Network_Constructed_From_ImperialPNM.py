# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 15:23:45 2021

@author: Priyanka.P.Agrawal
"""

import numpy as np
import openpnm as op
import pandas as pd
import matplotlib.pyplot as plt

def readNetwork (prefix, path):
    
    #Reading  Imperial Network File into OpenPNM
    proj = op.io.Statoil.import_data(path, prefix, network=None)
    
    print(proj)
    pn = proj.network
    pn.name = 'Alumina'
    print(pn)
    
    
    
    pore_radius_OpenPNM = pn['pore.radius']*1000 #Convert from um to nm
    throat_radius_OpenPNM = pn['throat.radius']*1000 #Convert from um to nm
    throat_length_OpenPNM = pn['throat.length']*1e-3

    #Reading of Imperial Network File
    
    ImperialNetwork_ThroatFile1 = path + prefix+'_link1.dat'
    with open(ImperialNetwork_ThroatFile1, mode='r') as f:
        link1 = pd.read_table(filepath_or_buffer=f,header=None, skiprows=1,sep=' ',skipinitialspace=True, index_col=0)
    link1.columns = ['throat.pore1', 'throat.pore2', 'throat.radius','throat.shape_factor', 'throat.total_length']
    
    ImperialNetwork_ThroatFile2 = path + prefix+'_link2.dat'
    with open(ImperialNetwork_ThroatFile2, mode='r') as f:
        link2 = pd.read_table(filepath_or_buffer=f,header=None, skiprows=1,sep=' ',skipinitialspace=True,index_col=0)
    link2.columns = ['throat.pore1', 'throat.pore2','throat.pore1_length', 'throat.pore2_length',
                             'throat.length', 'throat.volume','throat.clay_volume']
    
    ImperialNetwork_PoreFile = path + prefix+'_node2.dat'
   
    with open(ImperialNetwork_PoreFile, mode='r') as f:
        node2 = pd.read_table(filepath_or_buffer=f,header=None, skiprows=1,sep=' ',skipinitialspace=True, index_col=0)
    node2.columns = ['pore.volume', 'pore.radius', 'pore.shape_factor','pore.clay_volume']
    
    Imperial_PoreRadius = node2['pore.radius']*1000 #Convert from um to nm
    Imperial_ThroatRadius = link1['throat.radius']*1000 #Convert from um to nm
    Imperial_ThroatLength = link2['throat.length']*1e-3 #Convert from picometer to nm
    return pore_radius_OpenPNM, throat_radius_OpenPNM, throat_length_OpenPNM, Imperial_PoreRadius, Imperial_ThroatRadius, Imperial_ThroatLength

def plot (OpenPNM, Imperial,title, bins): 
        range1 = (min(OpenPNM), max(OpenPNM))
        [phist1, nbins] = np.histogram(OpenPNM,bins, range1)
        [phist2, nbins] = np.histogram(Imperial,bins, range1)
        fig1, ax = plt.subplots(figsize=(15, 10), dpi=100)
        ax.plot(nbins[1:bins+1],phist1, label = "OpenPNM_"+title, color = '#E69F00', linewidth = 4)
        ax.plot(nbins[1:bins+1],phist2, label = "Imperial_"+title, color = '#56B4E9', linewidth = 4)
        ax.legend(loc=1,fontsize = 20)
        ax.set_title(title, fontsize = 20)
        ax.set_xlabel(title+" (nm)", fontsize = 20)
        ax.tick_params(which='major', labelsize=20)
        plt.show()


#Location for the Imperial PoreNetwork Files
path = r"C:/Users/Priyanka.P.Agrawal/OneDrive - Shell/CatalystProject/Datasets/Constant/Cropped/SensitivityAnalysis/nRSmoothing/nRSmoothing_3/VoxelSize_0.29e-3/"

#Prefix used while extracting pore network
prefix = "Alumina"

#Extract releveant network properties from both networks
pore_radius_OpenPNM, throat_radius_OpenPNM, throat_length_OpenPNM, Imperial_PoreRadius, Imperial_ThroatRadius, Imperial_ThroatLength = readNetwork (prefix, path)

#Plot three releveant network properties from both networks
plot (pore_radius_OpenPNM, Imperial_PoreRadius,'Pore Radius', 10)
plot (throat_radius_OpenPNM, Imperial_ThroatRadius,'Throat Radius', 10)
plot (throat_length_OpenPNM, Imperial_ThroatLength,'Throat Length', 10)


