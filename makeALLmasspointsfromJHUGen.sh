#!/bin/bash
#############################################################################################
## PURPOSE: Run JHUGen on all Zdark mass points in the for loop below. 
##          This creates an LHE file for each mass point. 
##          Then skim each LHE file using ZZD_lhe.C. 
##          This creates a skimmed root file, useful for plotting.
## SYNTAX:  ./<script.sh>
## NOTES:   Check each command line option in the `./JHUGen ...' command
##          Also check the locations of all relevant files, like ZZD_lhe_template.C, 
## AUTHOR:  Jake Rosenzweig
## Date:    2018-10-01
######################################################################################

### PARAMETERS
epsilon=0.01
inputFile="stableHiggsfrompwg_15000events.lhe"
outputFile="new_ggHZZd4l_MZd"

cmsenv

#for zdmass in 1 2 3 4 7 10 15 20 25 30 35; do
for zdmass in 25; do
    ## PChannel=0: g+g
    ## 
    ./JHUGen ReadLHE=$inputFile Collider=1 VegasNc0=20000 PChannel=0 Process=0 Unweighted=1 DecayMode1=0 DecayMode2=0 OffshellX=0 ghz1=0,0 ghzzp1=1,0 ezp_El_left=1,0 ezp_Mu_left=1,0 Interf=0 MZprime=${zdmass} GaZprime=$( bc -l <<< "$zdmass * $epsilon * $epsilon" ) DataFile=${outputFile}${zdmass}.lhe

    #mv ggHZZd4l_MZd${zdmass}.lhe ../../
    #cd ../../

    ## Prepare the LHE skimmer for the current file
    sed -i "s/ZDMASS/${zdmass}/g" ZZD_lhe_template.C
    ## The -q option tells ROOT to quit when done
    root -q -l ZZD_lhe_template.C
    echo "Done processing Zdmass ${zdmass}"
    ## Start from the ZZD_lhe.C template again
    cp ../ZZD_lhe_template.C .
    #cd ./JHUGeneratorpackage/JHUGenerator/
done
