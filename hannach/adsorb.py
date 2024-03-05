import os
import warnings
import numpy as np
import openpnm as op
import scipy as sp
import scipy.stats as spst
import matplotlib.pyplot as plt
import pandas as pd
import copy

class ptList:
    def __init__(self):
        self.openPores = []
        self.blockedPores = []
        self.openThroats = []
        self.blockedThroats = []

def harkins_jura(x):
    t = np.sqrt(13.99/(0.034 - np.log(x)))
    return t

def r_crit(x,gamSt,vm,R,T):
    rCrit = -(2*gamSt*vm)/(R*T*np.log(x))
    rCrit *= 1e+7
    return rCrit

#def traverse(net,geo,inletPores,numPores,numThroats,t,rCrit):
def traverse(net,geo,inletPores,p1,numPores,numThroats,t,tPrev,rCrit):
    
    # Set data structures for open and blocked pores, throats
    #openPores = []
    #blockedPores = []
    #openThroats = []
    #blockedThroats = []
    
    currVolN2 = 0.0
    
    isPoreVisited = np.zeros((numPores,1))
    isThroatVisited = np.zeros((numThroats,1))
    
    poreStack = copy.deepcopy(inletPores)
    
    while (len(poreStack)):
        currentPore = poreStack.pop()
        #print('Current pore is: ',currentPore)
        #Debug
        #if (currentPore == 115):
        #    print('poreStack: ',currentPore)
        
        if ( (not isPoreVisited[currentPore]) and (currentPore not in p1.blockedPores)): # Traverse if pore is not visited so far
            isPoreVisited[currentPore] = 1
            rPore = 0.5*geo['pore.diameter'][currentPore]
            rFree = max(rPore - t,0)
            if (rFree < rCrit):
                p1.blockedPores.append(currentPore)
                if (currentPore in p1.openPores):
                    p1.openPores.remove(currentPore)
                
                # Add adsorbed N2
                currVolN2 += (4/3)*((rPore - tPrev)**3)
                #Debug
                #if (currentPore == 115):
                #    print('currVolN2: ',currVolN2)
            else:
                if (currentPore not in p1.openPores):
                    p1.openPores.append(currentPore) # Push current open pore on stack
                
                # Add adsorbed N2
                currVolN2 += (4/3)*((rPore - tPrev)**3 - rFree**3)
                #Debug
                #if (currentPore == 115):
                #    print('currVolN2: ',currVolN2)
                
                neighThroats = net.find_neighbor_throats(pores=[currentPore])
                for throat in neighThroats:
                    if ((not isThroatVisited[throat]) and (throat not in p1.blockedThroats) ):
                        isThroatVisited[throat] = 1
                        rThroat = 0.5*geo['throat.diameter'][throat]
                        rFree = max(rThroat - t,0)
                        if (rFree < rCrit):
                            p1.blockedThroats.append(throat)
                            if (throat in p1.openThroats):
                                p1.openThroats.remove(throat)
                            
                            # Add adsorbed N2
                            currVolN2 += ((rThroat - tPrev)**2)*geo['throat.length'][throat]
                        else:
                            if (throat not in p1.openThroats):
                                p1.openThroats.append(throat)
                            # Add adsorbed N2
                            currVolN2 += ((rThroat - tPrev)**2 - rFree**2)*geo['throat.length'][throat]
                            pores = geo['throat.conns'][throat]
                            z = pores[1] if pores[0] == currentPore else pores[0]
                            poreStack.append(z)
                            
    #return [isPoreVisited,isThroatVisited,openPores,blockedPores,blockedThroats]
    return [isPoreVisited,isThroatVisited, currVolN2]
                            
def traverse_desorb(net,geo,inletPores,p1,numPores,numThroats,t,tPrev,rCrit):

    currVolN2 = 0.0
    
    isPoreVisited = np.zeros((numPores,1))
    isThroatVisited = np.zeros((numThroats,1))
    
    poreStack = copy.deepcopy(inletPores)
    
    while (len(poreStack)):
        currentPore = poreStack.pop()
        #print('Current pore is: ',currentPore)
        #Debug
        #if (currentPore == 115):
        #    print('poreStack: ',currentPore)
        justAdded = 0
        
        if ( (not isPoreVisited[currentPore])): # Traverse if pore is not visited so far
            isPoreVisited[currentPore] = 1
            rPore = 0.5*geo['pore.diameter'][currentPore]
            rFree = max(rPore - t,0)
            if ((currentPore in p1.blockedPores) and (rFree > rCrit)): # Check if a blocked pore can be unblocked
                p1.blockedPores.remove(currentPore)
                p1.openPores.append(currentPore)
                justAdded = 1
                
            if ((currentPore in p1.openPores)):
                if (justAdded):
                    currVolN2 += (4/3)*((rPore - t)**3) # Pore was blocked during adsorption but open during desorption
                else:
                    currVolN2 += (4/3)*(rFree**3 - (rPore - tPrev)**3) # Pore was open before, new released volume is a shell
                
                neighThroats = net.find_neighbor_throats(pores=[currentPore])
                for throat in neighThroats:
                    justAddedThroat = 0
                    if ((not isThroatVisited[throat])):
                        isThroatVisited[throat] = 1
                        rThroat = 0.5*geo['throat.diameter'][throat]
                        rFree = max(rThroat - t,0)
                        if ((throat in p1.blockedThroats) and (rFree > rCrit)): # Check if a blocked pore can be unblocked
                            p1.blockedThroats.remove(throat)
                            p1.openThroats.append(throat)
                            justAddedThroat = 1
                            
                        if (throat in p1.openThroats):
                            if (justAddedThroat):
                                currVolN2 += ((rThroat - t)**2)*geo['throat.length'][throat]
                            else:
                                currVolN2 += (rFree**2 - (rThroat - tPrev)**2)*geo['throat.length'][throat]
                            pores = geo['throat.conns'][throat]
                            z = pores[1] if pores[0] == currentPore else pores[0]
                            poreStack.append(z)
                            
    return [isPoreVisited,isThroatVisited, currVolN2]                       
            

def adsorb(net,geo,inJson):
     
    maxDim = max(inJson['imSize']['x'],inJson['imSize']['y'],inJson['imSize']['z'])
    if  inJson['network_extraction_type'] == 'OpenPNM':
        imResolution = inJson['OpenPNM']['imResolution']
    elif inJson['network_extraction_type'] == 'Imperial':
        imResolution = inJson['Imperial']['imResolution']
        imResolution *= 1e-6
    else: 
        imResolution = inJson['cubic']['imResolution']

    print('Resolution: ', imResolution)

    # Gas constant
    R = 8.314
    T = 77
    
    numPres = 100
    
    ### Setup the pressures, thickness and rCrit
    #pres = np.linspace(0.1,1,10)
    pres = np.linspace(1/numPres,0.99,numPres)
    volN2 = np.zeros((numPres,1))
    volN2Desorb = np.zeros((numPres - 1,1))
    
    # Thickness needs special handing. We have to send current and
    # previous thickness to calculate the incremental adsorbed volume 
    # at the current pressure.
    t = np.zeros((numPres+1,1))
    t[0] = 0
    for i in range(1,numPres+1):
        t[i]= np.sqrt(13.99/(0.034 - np.log(pres[i-1])))*1e-10
    #t[1:numPres+1] = np.sqrt(13.99/(0.034 - np.log(pres)))*1e-10
    #print('Thickness Harkins-Jura in nm: ',t*1e+9)
    
    # Surface tension of N2 at 77K
    gamSt = 8.94e-3
    # Molar Volume of N2 at 77K
    molWt = 28e-3
    rho = 808
    vm = molWt/rho
    print('Molar volume of N2: ',vm)
    rCrit = -(2*gamSt*vm)/(R*T*np.log(pres))
    print('rCrit in nm: ',rCrit*1e+9)
    
    ## Testing
    #t = 0.1
    #rCrit = 0.1
    
    #plt.plot(pres,t)
    
    # inletPores has indices of pores which have x coordinate = dx*0.5, i.e. the left face
    inletPores = []
    
    numPores = len(net.Ps)
    numThroats = len(net.Ts)
    
    # Setup lists
    p1 = ptList()
    

    '''for i in range(len(net['pore.coords'])):
        if (net['pore.coords'][i][0] == dx*0.5):
            inletPores.append([i])
        elif (net['pore.coords'][i][1] == dx*0.5):
            inletPores.append([i])
        elif (net['pore.coords'][i][2] == dx*0.5):
            inletPores.append([i])
        elif (net['pore.coords'][i][0] == dx*14.5):
            inletPores.append([i])
        elif (net['pore.coords'][i][1] == dx*14.5):
            inletPores.append([i])
        elif (net['pore.coords'][i][2] == dx*14.5):
            inletPores.append([i]) '''

    for i in range(len(net['pore.coords'])):
        if (net['pore.coords'][i][0] <= imResolution):
            inletPores.append([i])
        elif (net['pore.coords'][i][1] <= imResolution):
            inletPores.append([i])
        elif (net['pore.coords'][i][2] <= imResolution):
            inletPores.append([i])
        #elif (net['pore.coords'][i][0] >= imResolution*(maxDim - 0.5)):
        #    inletPores.append([i])
        #elif (net['pore.coords'][i][1] >= imResolution*(maxDim - 0.5)):
        #    inletPores.append([i])
        #elif (net['pore.coords'][i][2] >= imResolution*(maxDim - 0.5)):
        #    inletPores.append([i])


    print('Number of inlet pores ',len(inletPores))
    print('inletpores are ', inletPores) 
    
    # Adsorption
    for i in range(numPres):
        [isPoreVisited,isThroatVisited,volN2[i]] = traverse(net,geo,inletPores,p1,numPores,numThroats,t[i+1],t[i],0.5*rCrit[i])
        #print('Open pores : ',i, len(p1.openPores))
        #print('Blocked pores : ',i, len(p1.blockedPores))
        #print('Blocked throats : ',i, len(p1.blockedThroats))
        #print('Open pores : ',i, len(p1.openPores))
        #print('Open throats : ',i, len(p1.openThroats))
    
    finalVolN2 = np.cumsum(volN2)*np.pi
    finalVolN2_desorb = np.zeros(numPres);
    finalVolN2_desorb[-1] = finalVolN2[-1]
    
    # Desorption
    numDesorb = numPres - 1
    for i in range(numDesorb):
        [isPoreVisited,isThroatVisited,volN2Desorb[numDesorb - i - 1]] = traverse_desorb(net,geo,inletPores,p1,numPores,numThroats,t[numDesorb - i - 1],t[numDesorb - i],rCrit[numDesorb - i - 1])
        #print('Open pores : ',i, len(p1.openPores))
        print('Blocked pores : ',i, len(p1.blockedPores))
        print('Blocked throats : ',i, len(p1.blockedThroats))
        print('Open pores : ',i, len(p1.openPores))
        print('Open throats : ',i, len(p1.openThroats))
        
        finalVolN2_desorb[numDesorb - i - 1] = finalVolN2_desorb[numDesorb - i] - (volN2Desorb[numDesorb - i - 1]*np.pi)
    
    #finalVolN2_desorb *= np.pi

        
    
    #print('volN2: ',volN2)
    #print('finalVolN2: ',finalVolN2)
    
    #plt.plot(pres,finalVolN2)
    #[isPoreVisited,isThroatVisited,openPores,blockedPores,blockedThroats] = traverse(net,geo,inletPores,numPores,numThroats,t,rCrit)
    #[isPoreVisited,isThroatVisited] = traverse(net,geo,inletPores,p1,numPores,numThroats,t,rCrit)
    
    return [pres,finalVolN2,finalVolN2_desorb,p1]
    
            
    
