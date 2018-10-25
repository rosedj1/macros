import numpy as np
import root_numpy
from lookatbranch import lookAtBranch

inputFile = "HIG-RunIIFall17MiniAOD-00079_99.root"
inputTree = "Ana/passedEvents"
kinemVar = "lep_eta"
numEvents = 20
lookAtBranch(inputFile, inputTree, kinemVar, numEvents)
