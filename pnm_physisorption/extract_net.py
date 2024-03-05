import imageio
import numpy as np
import porespy as ps
import openpnm as op
import os
import scipy.ndimage as spim
import matplotlib.pyplot as plt
import pandas as pd
import json
from porespy.filters import find_peaks, trim_saddle_points, trim_nearby_peaks
from porespy.tools import randomize_colors
from skimage.segmentation import watershed
import subprocess
from shutil import move


def extract_net_OpenPNM (inDir,inJson):
    
    #First load the image and complement it
    fPath = inDir + '\\' + inJson['imName']
    print(fPath)
    
    sigma = inJson['OpenPNM']['sigma']
    print('Sigma is ',sigma)
    
    im2 = np.fromfile(fPath,dtype='int8',sep="");
    im = np.where(im2 > 0, 0, 1)
    im = im.reshape(inJson['imSize']['x'], inJson['imSize']['y'], inJson['imSize']['z'])
    #im = im.reshape(imSize)
    plt.imshow(im[10,:,:],cmap = 'gray')
    
    ########################################

    # Gaussian filtering with sigma
    
    print('Filtering')
    dt_orig = spim.distance_transform_edt(input=im)
    dt = spim.gaussian_filter(input=dt_orig,  sigma=sigma)
    peaks = find_peaks(dt=dt)
    print('Filtering done')

    ## Visualize
    fig, ax = plt.subplots(1,2,figsize=(10,10))
    ax[0].imshow(dt_orig[10,:, :], cmap='gray');
    ax[1].imshow(dt[10,:, :], cmap='gray');
    
    ########################################
    
    ## Trim peaks ##
    print('Trimming')
    print('Initial number of peaks: ', spim.label(peaks)[1])
    peaks = trim_saddle_points(peaks=peaks, dt=dt)
    print('Peaks after trimming saddle points: ', spim.label(peaks)[1])
    peaks = trim_nearby_peaks(peaks=peaks, dt=dt)
    peaks, npeaks = spim.label(peaks)
    print('Peaks after trimming nearby peaks: ', npeaks)

    ########################################
    
    ## Do watershed ##
    print('In Watershed')
    regions = watershed(image=-dt, markers=peaks, mask=dt > 0)
    #regions = randomize_colors(regions)
    print('Watershed done')

    ## Plot watershed ##
    fig, ax = plt.subplots(1,2,figsize=(10,10))
    ax[0].imshow((regions*im)[10,:, :], cmap='gray');
    ax[1].imshow((regions)[10,:, :], cmap='gray');
    
    ########################################
    
    ## Extract Network ##
    print('Network extraction')
    net1 = ps.networks.regions_to_network(regions*im, voxel_size=1)
    proj = op.io.PoreSpy.import_data(net1)
    print(proj)
    
    ########################################
    
    return(proj)

    ########################################
  

def extract_net_Imperial (inDir,inJson):
    #First edit the input.mhd file used in the extraction of pore network from the Imperial code
    #read input.mhd file
    fPath = inDir + '\\' + inJson['Imperial']['inputfilename']
    outPath = inJson['Imperial']['outputDir']
    my_file = open(fPath)
    string_list = my_file.readlines()
    my_file.close()
    print(string_list)
    
    ########################################
    
    #Update input.mhd keywords values with the enteries in inJson file
    string_list[4] = 'DimSize =\t\t'+ str(inJson['imSize']['x']) +' '+ str(inJson['imSize']['y'])+' '+ str(inJson['imSize']['z'])+'\n'
    string_list[5] = 'ElementSpacing = '+str(inJson['Imperial']['imResolution'])+' '+str(inJson['Imperial']['imResolution'])+' '+str(inJson['Imperial']['imResolution'])+'\n'
    string_list[8] = 'ElementDataFile = '+ str(inJson['imName'])+'\n'
    string_list[11] = 'medialSurfaceSettings:'+ str(inJson['Imperial']['medialSurfaceSettings']['clipROutx'])+' '+str(inJson['Imperial']['medialSurfaceSettings']['clipROutyz'])+' '+ str(inJson['Imperial']['medialSurfaceSettings']['midRFrac'])+' '+str(inJson['Imperial']['medialSurfaceSettings']['RMedSurfNoise'])+' '+str(inJson['Imperial']['medialSurfaceSettings']['lenNf'])+' '+str(inJson['Imperial']['medialSurfaceSettings']['vmvRadRelNf'])+' '+str(inJson['Imperial']['medialSurfaceSettings']['nRSmoothing'])+' '+str(inJson['Imperial']['medialSurfaceSettings']['RCorsf']) +' ' +str(inJson['Imperial']['medialSurfaceSettings']['RCors'])+'\n'
    string_list[20] = 'title '+str(inJson['imperialPrefix'])+'\n'
    
    #Save the updated input.mhd file
    my_file = open(fPath, "w")
    new_file_contents = "".join(string_list)
    my_file.write(new_file_contents)
    my_file.close()
    
    readable_file = open(fPath)
    read_file = readable_file.read()
    print(read_file)

    ########################################

    #Run Imperial Code for extracting Network    
    path_To_Executable = inDir + '\\' +  '\\'+ 'pnextract.exe'
    process = subprocess.Popen([path_To_Executable, fPath],shell=False,stdout=subprocess.PIPE)

    # Poll process.stdout to show stdout live
    while True:
      output = process.stdout.readline()
      if process.poll() is not None:
        break
      if output:
        print (output.strip())
    rc = process.poll()
    
    Source = os.pardir + '\\' + str('pnm_physisorption') + '\\'
    #Destination = os.pardir + '\\' + str('test') + '\\'
    Destination = outPath
    FilesGeneratedByImperialCode = ['_link1.dat','_link2.dat','_mbHierarchy0.vtu','_medialSurface0.vtu','_node1.dat','_node2.dat','_poreMBs.tif',
     '_pores.vtu', '_poroats.tif','_radius.tif','_throatHierarchy.vtu','_throatMBs.tif','_throats.tif',
     '_throats.vtu','_throatsBalls.vtu','_VElems.mhd','_VElems.tif']
    
    for fle in FilesGeneratedByImperialCode:
        move(Source + inJson['imperialPrefix']+fle, Destination + inJson['imperialPrefix']+fle)
        
    
    return rc

    