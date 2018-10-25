import root_numpy

## inputTree = str tree name from ROOT file
## kinemVar = str name of ROOT branch
## numEvents = -1 for all events in branch
def lookAtBranch(inputFile, inputTree, kinemVar, numEvents):
     
    ## Turn ROOT TTree branch into numpy array
    brancharray  = root_numpy.root2array(inputFile, inputTree, kinemVar)

    ## Print it out and see what it looks like
    counter = 0
    print "\n", 8*"#"+len(kinemVar)*"#"+8*"#"
    print "\t", kinemVar
    print 8*"#"+len(kinemVar)*"#"+8*"#"
    print "The whole branch looks like:\n",brancharray,"\n",brancharray.shape,"\n"
    
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
