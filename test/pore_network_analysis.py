## analyze_pore_intersections is a function to analyze the overlap of pores 
## with surrounding solid. Inputs are:
## proj -> Two element array, first element should be network and second geom object
##         in OpenPNM format
##
## imSize -> Three element array, with X,Y and Z sizes of the volume
## im -> Image
## reso -> Image resolution
## resultPath -> Address of the output file
## resultFile -> Output excel file name. This file is a dataframe dump with three columns:
##               Pore_label : This is just the id of the pore
##               Intersections: It is the number of voxels within the pore, that are solid
##               NonIntersections: Number of voxels which are empty pore space
## In addition to the resultFile, the code also prints out some information, which can also be 
## obtained by analysis of the excel file:
## Number of pores whose conters are solid voxels
## Number of pores with solid voxels in them
## Number of solid voxels inside the pore space
## (Important) Percentage of solid voxels inside pore space. This is defined as:
## Percentage of intersections = 100*(Total solid voxels in the pores in the whole network)/Total pore voxels

import numpy as np
import pandas as pd
    
def pore_network_analysis(proj,im,imSize,reso,resultPath,resultFile):
    
    if (len(proj) == 2):
        networkType = 'openPNM'
    else:
        networkType = 'imperial'

    print('Network Type is ',networkType)

    if (networkType == 'imperial'):
        net = proj.network
    elif (networkType == 'openPNM'):
        net = proj[0]
        geo = proj[1]

    # TESTING
    #net = proj
    #networkType = 'imperial'
    ## TESTING
    print(net)
    
    nPores = len(net['pore.coords'])
    print('Total pores are ', nPores)
    #print(len(geo['pore.diameter']))
    
    # Initialize lists 
    # Number of pores whose center is solid
    nPoreSolid = []
    # Pores with intersections
    pcAlum = []
    # Solid Voxels inside pores
    pSolidVoxel = []
    
    # Main output data structure. First column is pore index
    # Second is number of solid voxels, third is number of pore(empty) voxels
    poreIntersect = np.zeros((nPores,3))
    
    for i in range(nPores):
        
        # Get Radius #
        if (networkType == 'imperial'):
            pc = net['pore.coords'][i]*(1e-6/reso)
            r = net['pore.radius'][i]*(1e-6/reso)
        if (networkType == 'openPNM'):
            pc = net['pore.coords'][i]
            r = 0.5*geo['pore.diameter'][i]
        
        rSq = r*r
        
        ix = min((pc[0]),imSize[0] - 1)
        iy = min((pc[1]),imSize[1] - 1)
        iz = min((pc[2]),imSize[2] - 1)
        
        
        # First check if centre is alumina #
        if im[int(ix),int(iy),int(iz)] == 255:
            nPoreSolid.append(i)
            
        ## Sphere check, limits are per octant##
        ## If the point is inside the octant all other octants are checked ##
        xmin = int(max(0,ix - r))
        xmax = int(min(imSize[0] - 1,ix + r))
        ymin = int(max(0,iy - r ))
        ymax = int(min(imSize[1] - 1,iy + r ))
        zmin = int(max(0,iz - r ))
        zmax = int(min(imSize[2] - 1,iz + r))
        
        print('i: ',i)
        #print('zmin, ymax:',ymin, ymax)
        
        zrange = range(zmin, int(iz + 1))
        yrange = range(ymin, int(iy + 1))
        xrange = range(xmin, int(ix + 1))
    
        isSolid = 0
        isPore = 0
        
        for z in zrange:
            for y in yrange:
                for x in xrange:
                    if ((x - ix)**2 + (y - iy)**2 + (z - iz)**2 > rSq):
                        #print('Here')
                        continue
                    else:
                       
                        #print('Here in big loop')
                        pList = np.zeros((8,4),dtype = np.uint16)
                        #x2 = xmax - x + 1
                        #y2 = ymax - y + 1
                        #z2 = zmax - z + 1
                        
                        x2 = xmax + xmin - x
                        y2 = ymax + ymin - y
                        z2 = zmax + zmin - z
                        
                        # Check for repeat points #
                        
                        if ((x == x2) and (y == y2) and (z == z2)): # centre of pore
                            numPoints = 1
                            pList[0,1:4] = [x,y,z]
                        elif ( (x == x2) and (y == y2) ): # point on z axis
                            numPoints = 2
                            pList[0,1:4] = [x,y,z]
                            pList[1,1:4] = [x,y,z2]
                        elif ( (x == x2) and (z == z2) ): #point on y axis
                            numPoints = 2
                            pList[0,1:4] = [x,y,z]
                            pList[1,1:4] = [x,y2,z]
                        elif ( (z == z2) and (y == y2) ): # point on x axis
                            numPoints = 2
                            pList[0,1:4] = [x,y,z]
                            pList[1,1:4] = [x2,y,z]
                        elif (x == x2): # point on yz plane
                            numPoints = 4
                            pList[0,1:4] = [x,y,z]
                            pList[1,1:4] = [x,y2,z]
                            pList[2,1:4] = [x,y,z2]
                            pList[3,1:4] = [x,y2,z2]
                        elif (y == y2): # point on xz plane
                            numPoints = 4
                            pList[0,1:4] = [x,y,z]
                            pList[1,1:4] = [x2,y,z]
                            pList[2,1:4] = [x,y,z2]
                            pList[3,1:4] = [x2,y,z2]
                        elif (z == z2): # point on xy plane
                            numPoints = 4
                            pList[0,1:4] = [x,y,z]
                            pList[1,1:4] = [x2,y,z]
                            pList[2,1:4] = [x,y2,z]
                            pList[3,1:4] = [x2,y2,z]
                        else:
                            numPoints = 8 # Point situated strictly within octant
                            # Same z plane
                            pList[0,1:4] = [x,y,z]
                            pList[1,1:4] = [x2,y,z]
                            pList[2,1:4] = [x,y2,z]
                            pList[3,1:4] = [x2,y2,z]
                        
                            # Other z plane
                            pList[4,1:4] = [x2,y,z2]
                            pList[5,1:4] = [x,y2,z2]
                            pList[6,1:4] = [x2,y2,z2]
                            pList[7,1:4] = [x,y,z2]
                        
                        for p in range(numPoints):
                            if (im[pList[p,1],pList[p,2],pList[p,3]]) == 255: # solid voxel check
                            #if (im[pList[p,1],pList[p,2],pList[p,3]]) == 0: # solid voxel check
                                isSolid += 1
                                pList[p,0] = i
                                pSolidVoxel.append(pList[p,:])
                            else:
                                isPore += 1
    
        if (isSolid > 0):
            pcAlum.append(i)
        poreIntersect[i,0] = i
        poreIntersect[i,1] = isSolid
        poreIntersect[i,2] = isPore
                     
    
    print('Number of alumpores',len(nPoreSolid))
    #print(net['pore.coords'][nPoreSolid])
    print('Number of pores with intersections ',len(pcAlum))
    print('Number of points with intersections ',len(pSolidVoxel))

    ## Analysis
    
    totalIntersect = np.sum(poreIntersect[:,1])
    totalNonIntersect = np.sum(poreIntersect[:,2])
    totalVoxels = totalIntersect + totalNonIntersect

    print('Total intersect is: ',totalIntersect)
    print('Total inon-ntersect is: ', totalNonIntersect)
    print('Percentage of intersect is: ',(totalIntersect/totalVoxels))
    
    # Write output
    outPath = resultPath + '\\' + resultFile

    df1 = pd.DataFrame(poreIntersect)
    df2 = pd.DataFrame(pSolidVoxel)
    
    df1.columns = ['Pore_label','Solid','Pore']
    df2.columns = ['Pore_label','x','y','z']
    
    df1.to_excel(outPath,sheet_name='pore_intersections')
    #with pd.ExcelWriter(f1_path,engine='openpyxl',mode='a') as writer:
    #    df2.to_excel(writer,sheet_name='points_intersections')

    
    
    