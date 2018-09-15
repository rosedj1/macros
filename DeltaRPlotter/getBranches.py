#######################################
##### READ NTUPLE AND MAKE HISTOS #####
#######################################

import ROOT as root
import numpy as np
import root_numpy as rnp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## enter batch mode in root (so python can access displays)
root.gROOT.SetBatch(True)


def openBranch(kinemVar, branch, numEvents):
     
    counter = 0
    print "\n", "#"*20, "\n", "\t", kinemVar, "\n", "#"*20
    #for k in range(len(branch)-1):
    for k in range(numEvents):
        print "Event number:", k, branch[k]
        counter += 1

    print "Total number of events looped over:", counter
    print "Num events =", numEvents

def analyzeBranch(kinemVar, branch, numEvents, numParticles, showInfo=False):

    counter = 0
    ## Create starting array with same shape as branches
    branchArray = np.array([0] * numParticles )

    if numEvents == -1: numEvents = len(branch)
    #for k in range(len(branch)-1):
    for k in range(len(branch[:numEvents])):
        ## Make sure that event is filled with correct num of particles
        if len(branch[k]) != numParticles: continue
        branchArray = np.vstack( (branchArray, branch[k]) ) 
        counter += 1

    ## Remove the original [0,...,0] array 
    ## the 0 in the delete method indicates the 0th array
    branchArray = np.delete(branchArray, 0, axis=0)
    if showInfo == True:
        print "\n", "#"*30, "\n\t", kinemVar, "\n", "#"*30
        print "Total events:", numEvents, "\n", "Events analyzed:", counter
        #print "This is also hopefully the total number of events looped over:", np.shape(branchArray)[0]
        print "Branch array looks like:", "\n", branchArray, "\t" 
        print "Dimensions:", np.shape(branchArray)

    return branchArray

def deltaR(dEta, dPhi):
    dR = np.sqrt(dEta**2 + dPhi**2)
    return dR

## Make sure that: 0 <= deltaPhi <= pi
def deltaPhi(phi1, phi2):
    result = abs(phi1-phi2)
    while result > np.pi:
        result = 2*np.pi - result
    return result
    
def Plotter(x_array, y_array, x1_array=None, y1_array=None):
    plt.scatter(x_array, y_array, c='b')
    plt.scatter(x1_array, y1_array, c='r')
    plt.title('H-->ZdZ-->4l')
    plt.xlabel('eta')
    plt.ylabel('phi')
    #plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    plt.show()

def histoPlotter():
    pass

#openBranch("Z_eta", Z_eta_branch)
#openBranch("Z_phi", Z_phi_branch)
#openBranch("lep_eta", lep_eta_branch, numEvents)
#H_eta_array      = analyzeBranch("H_eta", H_eta_branch, numEvents, 1, showInfo=showInfo)
#H_phi_array     = analyzeBranch("H_phi", H_phi_branch, numEvents, 1, showInfo=showInfo)

## Ultimately want to make a histogram of deltaR values
def makeDeltaRList(inputFile, numEvents):

    ## read data from file into a numpy array
    lep_eta_branch  = rnp.root2array(inputFile, inputTree, "lep_eta")
    lep_phi_branch  = rnp.root2array(inputFile, inputTree, "lep_phi")
    #H_eta_branch    = rnp.root2array(inputFile, inputTree, "H_eta")
    #H_phi_branch    = rnp.root2array(inputFile, inputTree, "H_phi")
    #Z_eta_branch    = rnp.root2array(inputFile, inputTree, "Z_eta")
    #Z_phi_branch    = rnp.root2array(inputFile, inputTree, "Z_phi")

    lep_eta_array  = analyzeBranch("lep_eta", lep_eta_branch, numEvents, 4, showInfo=showInfo)
    lep_phi_array   = analyzeBranch("lep_phi", lep_phi_branch, numEvents, 4, showInfo=showInfo)
    flat_lep_eta_array = list(lep_eta_array.flatten() )
    flat_lep_phi_array = list(lep_phi_array.flatten() )

    final_delta_R_list = []

    for lep in range(len(flat_lep_eta_array)):

        ## Make sure lepton doesn't look at itself for eta or phi values
        delta_eta_list = []
        for eta2 in flat_lep_eta_array[:lep] + flat_lep_eta_array[(lep+1):] :
            if showInfo == True: print "\nlep", lep, "is looking at the following eta:", eta2
            dEta = abs( flat_lep_eta_array[lep] - eta2 )
            if showInfo == True: print "dEta is:", dEta
            delta_eta_list.append(dEta)
        if showInfo == True: print "lep", lep, "delta_eta_list:", delta_eta_list

        delta_phi_list = []
        for phi2 in flat_lep_phi_array[:lep] + flat_lep_phi_array[(lep+1):] :
            if showInfo == True: print "\nlep", lep, "is looking at the following phi:", phi2
            dPhi = deltaPhi( flat_lep_phi_array[lep], phi2 )
            if showInfo == True: print "dPhi is:", dPhi
            delta_phi_list.append(dPhi)
        if showInfo == True: print "lep", lep, "delta_phi_list:", "\n", delta_phi_list

        delta_R_list = []
        if showInfo == True: 
            print "\nlength of delta_eta_list:", len(delta_eta_list)
            print "length of delta_phi_list:", len(delta_phi_list)
        for k in range(len(delta_eta_list)):
            dR = deltaR(delta_eta_list[k], delta_phi_list[k])
            delta_R_list.append(dR)
        if showInfo == True:
            print "delta_R_list for lep", lep, ":", delta_R_list
            print "lep", lep, "was found to have min deltaR of:", min(delta_R_list), "\n"
        final_delta_R_list.append( min(delta_R_list)  )

    if showInfo == True: 
        print "Total number of leptons analyzed:", len(flat_lep_eta_array)
        print "Length of delta_R_list:", len(final_delta_R_list)
        print "Match?", len(flat_lep_eta_array) == len(final_delta_R_list)
        print "final delta_R_list:", final_delta_R_list

    return final_delta_R_list


inputFile = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180706/ZD_UpTo0j_MZD15_Eps1e-2_klo.root"
inputTree = "Ana/passedEvents"
numEvents = 400
showInfo = True
bins = 100

#def analyzeBranch(kinemVar, branch, numEvents, numParticles, showInfo=False):
lep_eta_branch  = rnp.root2array(inputFile, inputTree, "lep_eta")
lep_phi_branch  = rnp.root2array(inputFile, inputTree, "lep_phi")
analyzeBranch("lep_eta", lep_eta_branch, numEvents, 4, showInfo=showInfo)
analyzeBranch("lep_phi", lep_phi_branch, numEvents, 4, showInfo=showInfo)
#dRlist1 = makeDeltaRList(inputFile1, numEvents)
#dRlist2 = makeDeltaRList(inputFile2, numEvents)
#dRlist3 = makeDeltaRList(inputFile3, numEvents)
#dRlist4 = makeDeltaRList(inputFile4, numEvents)
#dRlist5 = makeDeltaRList(inputFile5, numEvents)
#
#sns.distplot(dRlist1, bins, kde=False, color='g')
#sns.distplot(dRlist2, bins, kde=False, color='k')
#sns.distplot(dRlist3, bins, kde=False, color='b')
#sns.distplot(dRlist4, bins, kde=False, color='r')
#sns.distplot(dRlist5, bins, kde=False, color='m')
#plt.show()



#for lep in range(len(flat_lep_eta_array)):
#
#    ## Make sure lepton doesn't look at itself for eta or phi values
#    delta_eta_list = []
#    for eta2 in flat_lep_eta_array[:lep] + flat_lep_eta_array[(lep+1):] :
#        if showInfo == True: print "\nlep", lep, "is looking at the following eta:", eta2
#        dEta = abs( flat_lep_eta_array[lep] - eta2 )
#        if showInfo == True: print "dEta is:", dEta
#        delta_eta_list.append(dEta)
#    if showInfo == True: print "lep", lep, "delta_eta_list:", delta_eta_list
#
#    delta_phi_list = []
#    for phi2 in flat_lep_phi_array[:lep] + flat_lep_phi_array[(lep+1):] :
#        if showInfo == True: print "\nlep", lep, "is looking at the following phi:", phi2
#        dPhi = deltaPhi( flat_lep_phi_array[lep], phi2 )
#        if showInfo == True: print "dPhi is:", dPhi
#        delta_phi_list.append(dPhi)
#    if showInfo == True: print "lep", lep, "delta_phi_list:", "\n", delta_phi_list
#
#    delta_R_list = []
#    if showInfo == True: 
#        print "\nlength of delta_eta_list:", len(delta_eta_list)
#        print "length of delta_phi_list:", len(delta_phi_list)
#    for k in range(len(delta_eta_list)):
#        dR = deltaR(delta_eta_list[k], delta_phi_list[k])
#        delta_R_list.append(dR)
#    if showInfo == True:
#        print "delta_R_list for lep", lep, ":", delta_R_list
#        print "lep", lep, "was found to have min deltaR of:", min(delta_R_list), "\n"
#    final_delta_R_list.append( min(delta_R_list)  )
#
#if showInfo == True: 
#    print "Total number of leptons analyzed:", len(flat_lep_eta_array)
#    print "Length of delta_R_list:", len(final_delta_R_list)
#    print "Match?", len(flat_lep_eta_array) == len(final_delta_R_list)
#    print "final delta_R_list:", final_delta_R_list
####final_delta_R_list = []
####
####for lep in range(len(flat_lep_eta_array)):
####
####    delta_eta_list = []
####    for eta2 in flat_lep_eta_array[(lep+1):]:
####        print "lep", lep, "is looking at the following eta:", eta2
####        dEta = abs( flat_lep_eta_array[lep] - eta2 )
####        delta_eta_list.append(dEta)
####    print "lep", lep, "delta_eta_list:", "\n", delta_eta_list
####
####    delta_phi_list = []
####    for phi2 in flat_lep_phi_array[(lep+1):]:
####        print "lep", lep, "is looking at the following phi:", phi2
####        dPhi = deltaPhi( flat_lep_phi_array[lep], phi2 )
####        print "dPhi is:", dPhi
####        delta_phi_list.append(dPhi)
####    print "lep", lep, "delta_phi_list:", "\n", delta_phi_list
####
####    delta_R_list = []
####    print "length of delta_eta_list:", len(delta_eta_list)
####    print "length of delta_phi_list:", len(delta_phi_list)
####    for k in range(len(delta_eta_list)):
####        dR = deltaR(delta_eta_list[k], delta_phi_list[k])
####        delta_R_list.append(dR)
####    print "delta_R_list:", delta_R_list
####    print "lep", lep, "was found to have min deltaR of:", min(delta_R_list)
####    final_delta_R_list.append( min(delta_R_list)  )
####
####print final_delta_R_list


####delta_phi_list = []
####for phi1 in flat_lep_phi_array:
####    for phi2 in flat_lep_phi_array[ ((list(flat_lep_phi_array).index(phi1))+1):] :
####        dPhi = deltaPhi(phi1, phi2)
####        delta_phi_list.append(dPhi)
####print delta_phi_list
####
####delta_eta_list = []
####for eta1 in flat_lep_eta_array:
####    for eta2 in flat_lep_eta_array[ ((list(flat_lep_eta_array).index(eta1))+1):]:
####        dEta = abs(eta1 - eta2)
####        delta_eta_list.append(dEta)
####print delta_eta_list, "\n"
####
####delta_R_list = []
####print len(delta_eta_list)
####for k in range(len(delta_eta_list)):
####    dR = deltaR(delta_eta_list[k], delta_phi_list[k])
####    delta_R_list.append(dR)
####print delta_R_list
####print min(delta_R_list)

## Make deltaR_array
#deltaR_array = []
#deltaR_array = deltaR_array.append(deltaR(???))
#plt.scatter(H_eta_array, H_phi_array, c='r', s=50)

#Plotter(lep_eta_array, lep_phi_array)
#Plotter(lep_eta_array, lep_phi_array, H_eta_array, H_phi_array)

#fig, axes = plt.subplots()
#axes.plot(H_phi_array, H_eta_array, label='Higgs', color='r')
#axes.plot(lep_phi_array, lep_eta_array, label='4 leptons', color='b')
#axes.legend()
#fig.savefig('firstattempt.pdf')

#analyzeBranch("Z_eta", Z_eta_branch, 1)
#analyzeBranch("Z_phi", Z_phi_branch, 1)







##### GARBAGE ####
def getIndicesPassEtaCut(kinemVar, branch, numEvents):

    counter = 0
    indicesPassedEtaCut = []
    
    #for k in range(len(branch)-1):
    for k in range(numEvents):
        pass
        ## If eta collection and fails eta cut: continue
        ## mapping gives a boolean list
        #print map( lambda x:abs(x)<=2.4 , list(branch[k]) )
        #if ( ('eta' in kinemVar.split('_') ) and ( sum(map( lambda x: (abs(x) >= 2.4), branch[k] ))!=0 ) ): continue
        ## Save event number k
        #else:
        #    indicesPassedEtaCut.append(k)
        #    counter += 1

    print "\n", "Total number of", kinemVar, "events that passed eta cut:", counter 
    print indicesPassedEtaCut
    return indicesPassedEtaCut

