import os
import warnings
import numpy as np
import openpnm as op
import scipy as sp
import scipy.stats as spst
import matplotlib.pyplot as plt
import pandas as pd
import copy


def makeCubicNet(numNodes,dx,inFile,sheet):
    pdfRaw = pd.read_excel(inFile,sheet_name = 'PoreSizeDistribution')
    #pdfRaw.sort_values(by='pDia')
    #print(pdfRaw)
    
    # SI conversion
    pdf = pdfRaw
    pdf['pDia'] *= 1e-9
    pdf['v'] *= 1e-3
    
    pdfArr = pdf.to_numpy() # Numpy 
    #pdfArr[:,1] /= (np.pi/6)*(pdfArr[:,0]**3) # Convert volume to frequencies
    
    # Area normalization
    pdfArea = np.trapz(pdf['v'],pdf['pDia'])
    pdfArr[:,1] /= pdfArea
    newPdfArea = np.trapz(pdf['v'],pdf['pDia'])
    print('New pdf area is: ',newPdfArea)
    
    numVoxels = numNodes**3
    #urand = np.random.rand(numVoxels,1)
    #finalPoreDia = np.zeros([numVoxels,1])
        
    # Sample from the distribution with second method
    # Set up the geometry
    #np.random.seed(10)
    pn = op.network.Cubic(shape=[numNodes, numNodes, numNodes], spacing=dx, connectivity=6)
    geo = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)

    pdf2 = pdfArr
    pdf3 = copy.deepcopy(pdfArr)
    
    pdf2[:,1] /= np.sum(pdf2[:,1])
    print('Sum voxels before ', np.sum(pdf2[:,1]))
    pdf2[:,1] *= numVoxels
    print('Sum voxels after ', np.sum(pdf2[:,1]))
    #print(pdf2)
    plt.xscale('log')
    fig1 = plt.plot(pdf2[:,0],pdf2[:,1])
    plt.xlabel('Pore diameter (log)')
    plt.ylabel('Count total')
    plt.show()
    
    # pdf2 is the histogram to sample from. The total number of pores is
    # numVoxel. Pores of each size have to be distributed in our network randomly.
    # For this, the pore indices from 0 to numVoxels - 1 are generated and shuffled
    # Sizes are then allocated from pdf2
    
    pInds = np.arange(numVoxels)
    np.random.shuffle(pInds)
    #print(pInds)
    s = pdf2.shape[0]
    #print('s is ',s)
    lim1 = 0
    geo['pore.diameter'] = 0.0
    
    for i in range(s):
        lim2 = lim1 + round(pdf2[i,1])
        lim2 = min(lim2,numVoxels)
        if (lim2 > numVoxels):
            break
        z = pInds[lim1:lim2]
        #print(lim1)
        lim1 = lim2
        geo['pore.diameter'][z[:]] = pdf2[i,0]
        
        
    ## TESTING ##########################
    '''geo['pore.diameter'][pInds[0]] = pdf2[s-4,0]
    geo['pore.diameter'][pInds[1]] = pdf2[s-4,0]
    geo['pore.diameter'][pInds[2]] = pdf2[s-4,0]
    geo['pore.diameter'][pInds[3]] = pdf2[s-4,0]
    geo['pore.diameter'][pInds[4]] = pdf2[s-4,0]
    geo['pore.diameter'][pInds[5]] = pdf2[s-3,0]
    geo['pore.diameter'][pInds[6]] = pdf2[s-3,0]
    geo['pore.diameter'][pInds[7]] = pdf2[s-3,0]
    geo['pore.diameter'][pInds[8]] = pdf2[s-3,0]
    geo['pore.diameter'][pInds[9]] = pdf2[s-2,0]
    geo['pore.diameter'][pInds[10]] = pdf2[s-2,0]
    geo['pore.diameter'][pInds[11]] = pdf2[s-2,0]
    geo['pore.diameter'][pInds[12]] = pdf2[s-1,0]
    geo['pore.diameter'][pInds[13]] = pdf2[s-1,0]'''
    
    
    f = op.models.misc.from_neighbor_pores
    geo.add_model(propname='throat.diameter',
                  model=f, prop='pore.diameter', mode='min')
    f = op.models.geometry.throat_length.classic
    geo.add_model(propname='throat.length',
                  model=f,pore_diameter='pore.diameter' )
    return [pdf2,pn,geo]

    ''' Throats to be sampled from a truncated distribution limited by the min diameter of the connecting
    pore 
    ## Throat diameter testing ##
    #pdf2 = pdfArr
    #pdf2[:,1] /= np.sum(pdf2[:,1])
    pdf2[:,1] /= np.trapz(pdf2[:,1],pdf2[:,0])
    #print('Throat Sum voxels before ', np.sum(pdf2[:,1]))
    print('Area under throat distri ', np.trapz(pdf2[:,1],pdf2[:,0]))
    #pdf2[:,1] *= 9450
    #print('Throat Sum voxels after ', np.sum(pdf2[:,1]))
    
    # Sample throat diameter smaller than either connecting pore diameter
    numThroats = len(pn.Ts)
    print('numThroats ',numThroats)
    
    #The CDF of the pore size distribution
    cdf =np.zeros([s,1])
    cdf = sp.integrate.cumtrapz(pdf2[:,1],pdf2[:,0],initial=0)    
    print(np.shape(cdf))
    # Plot
    plt.xscale('log')
    fig3 = plt.plot(pdf2[:,0],cdf)
    plt.xlabel('Pore diameter (log)')
    plt.ylabel('CDF')
    plt.show()
    
    urand = np.random.rand(numThroats,1)
    
    for throat in range(numThroats):
        # Connecting pores of the throat
        pores = geo['throat.conns'][throat]
        # Throat dia should be smaller than either pore dia
        if (geo['pore.diameter'][pores[0]] < geo['pore.diameter'][pores[1]]):
            minDia = geo['pore.diameter'][pores[0]]
            #minPoreId = pores[0] # id of pore
        else:
            minDia = geo['pore.diameter'][pores[1]]
            #minPoreId = pores[1] # id of pore
            
        pIndex = np.where(pdf2[:,0] == minDia)
        rand1 = urand[i]*cdf[pIndex] # Sampling from truncated distribution
        # r1 is the list of all values in the cdf greater than the random number
        #r1 = [cdfVal for cdfVal in cdf if cdfVal > urand[i]]
        r1 = np.where(cdf >= rand1)
        rMin = np.min(r1)
        geo['throat.diameter'][throat] = pdf2[rMin,0]
           
    f = op.models.geometry.throat_length.classic
    geo.add_model(propname='throat.length',
                  model=f,pore_diameter='pore.diameter' )
    
    plt.xscale('log')
    [y,x] = np.histogram(geo['throat.diameter'],bins = pdf2[:,0])
    fig4 = plt.plot(x[0:s-1],y)
    plt.xlabel('Throat diameter (log)')
    plt.ylabel('Pdf')
    plt.show()
        
    return [pdf2,pn,geo]'''
            
    '''tInds = np.arange(9450)                     
    np.random.shuffle(tInds)
    #print(pInds)
    s = pdf2.shape[0]
    #print('s is ',s)
    lim1 = 0
    geo['throat.diameter'] = 0.0
    
    for i in range(s):
        lim2 = lim1 + round(pdf2[i,1])
        lim2 = min(lim2,9450)
        if (lim2 > 9450):
            break
        z = tInds[lim1:lim2]
        #print(lim1)
        lim1 = lim2
        geo['throat.diameter'][z[:]] = pdf2[i,0]'''
    
    '''# Area under pdf should be 1
    pdfArea = np.trapz(pdf['v'],pdf['pDia'])
    
    #print(pdfArea)
    pdf['v'] *= (1.0/pdfArea)
    newPdfArea = np.trapz(pdf['v'],pdf['pDia'])
    print('New pdf area is: ',newPdfArea)
    plt.xscale('log')
    fig2 = plt.plot(pdf['pDia'],pdf['v'],marker='o',linestyle='None')
    plt.xlabel('Pore diameter (log)')
    plt.ylabel('PDF')
    plt.show()
    #plt.plot(pdf['pDia'],pdf['v'])
    
    # Determine the CDF to draw random numbers
    # Convert to numpy array
    pdfArr = pdf.to_numpy()
    s = pdfArr.shape[0]
    print(s)'''
    
    
        
    '''
    for i in range(numVoxels):
        # r1 is the list of all values in the cdf greater than the random number
        #r1 = [cdfVal for cdfVal in cdf if cdfVal > urand[i]]
        r1 = np.where(cdf > urand[i])
        rMin = np.min(r1)
        if (i == 0):
            print('r1 is ',r1)
            print('rmin is ',rMin)
            
        #finalPoreDia[i] = pdf['pDia'].iloc[rMin]
        finalPoreDia[i] = pdfArr[rMin,0]
    plt.xscale('log')
    fig4 = plt.hist(finalPoreDia,bins=28)
    plt.xlabel('Pore diameter (log)')
    plt.ylabel('Count')
    plt.show()
    
    plt.xscale('log')
    spdf, sBins = np.histogram(finalPoreDia, bins=10, normed=True)
    plt.plot(sBins[0:10],spdf,linestyle='None',marker='*')
    plt.xlabel('Pore diameter (log)')
    plt.ylabel('Count')
    plt.show()'''
        

