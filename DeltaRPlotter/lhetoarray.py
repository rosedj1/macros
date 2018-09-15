import numpy as np
import root_numpy

##### PURPOSE #####
## Convert any kinematic variable branch from a processed LHE file 
## (cmsgrid_final.lhesyntax2.root) into a numpy array for easy processing. 
## Return the array. 
###################

## inputTree = str tree name from ROOT file
## kinemVar = str name of ROOT branch
## numEvents = -1 for all events in branch
## Possible kinemVar names in LHE file:
    #weight
    #id[3-6]
    #eta[3-6]
    #pt[3-6]
    #pto[1-4]
    #mass4l
    #pT4l
    #eta4l
    #rapidity4l
    #massZ1
    #massZ2
    #costheta1
    #costheta2
    #costhetastar
    #phi
    #phi1

def lheToArray(inputFile,inputTree,kinemVar,numEvents,showInfo=False):
    
    ## Turn ROOT TTree branch into numpy array
    ## Combine arrays such that each row is a single event
    if kinemVar == "id":
        brancharray1  = root_numpy.root2array(inputFile, inputTree, "id3")
        brancharray2  = root_numpy.root2array(inputFile, inputTree, "id4")
        brancharray3  = root_numpy.root2array(inputFile, inputTree, "id5")
        brancharray4  = root_numpy.root2array(inputFile, inputTree, "id6")
        brancharray = np.vstack((brancharray1, brancharray2, brancharray3, brancharray4)).T
    elif kinemVar == "eta":
        brancharray1  = root_numpy.root2array(inputFile, inputTree, "eta3")
        brancharray2  = root_numpy.root2array(inputFile, inputTree, "eta4")
        brancharray3  = root_numpy.root2array(inputFile, inputTree, "eta5")
        brancharray4  = root_numpy.root2array(inputFile, inputTree, "eta6")
        brancharray = np.vstack((brancharray1, brancharray2, brancharray3, brancharray4)).T
    elif kinemVar == "pt":
        brancharray1  = root_numpy.root2array(inputFile, inputTree, "pt3")
        brancharray2  = root_numpy.root2array(inputFile, inputTree, "pt4")
        brancharray3  = root_numpy.root2array(inputFile, inputTree, "pt5")
        brancharray4  = root_numpy.root2array(inputFile, inputTree, "pt6")
        brancharray = np.vstack((brancharray1, brancharray2, brancharray3, brancharray4)).T
    elif kinemVar == "pto":
        brancharray1  = root_numpy.root2array(inputFile, inputTree, "pto1")
        brancharray2  = root_numpy.root2array(inputFile, inputTree, "pto2")
        brancharray3  = root_numpy.root2array(inputFile, inputTree, "pto3")
        brancharray4  = root_numpy.root2array(inputFile, inputTree, "pto4")
        brancharray = np.vstack((brancharray1, brancharray2, brancharray3, brancharray4)).T
    else: 
    ## Then kinemVar has only one kind of branch in ROOT tree
        brancharray = root_numpy.root2array(inputFile, inputTree, kinemVar)
        brancharray = brancharray.reshape(len(brancharray),1)

    ## Print to screen
    if showInfo == True:
        ## Print it out and see what it looks like
        print "\n", 8*"#"+len(kinemVar)*"#"+8*"#"
        print "\t", kinemVar
        print 8*"#"+len(kinemVar)*"#"+8*"#"
        print "The whole branch looks like:\n",brancharray,"\n",brancharray.shape,"\n"
        
        counter = 0
        if numEvents == -1:
            ## Print ALL events
            for k in range(len(brancharray)):
                print "Event number",k,":", brancharray[k]
                counter += 1
        else:
            for k in range(numEvents):
                print "Event number",k,":", brancharray[k]
                counter += 1
        
        print "Total number of events looped over:", counter
        print "Num events called for:", numEvents

    return brancharray
