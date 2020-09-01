#!/bin/bash
#####################################################################################################
## PURPOSE: Generate a txt file with LFN of input files to be read by a crab_cfg.py file.
## SYNTAX:  ./<script.sh> <new_file_name> <start_num> <stop_num> <input_root_files>
## NOTES:      <new_file_name>    = name of newly generated txt file to hold input file paths.
##             <start_num>        = root files are usually numbered, so this is the starting number
##             <stop_num>         = stop number
##             <input_root_files> = LFN (file path: `/store/user/...`) of input root files
##                                  Don't include the `.root` at the end!
## AUTHOR:  Jake Rosenzweig
## DATE:    2019-04-01 (miss you, Mom)
## UPDATED: 2019-06-28
#####################################################################################################

#____________________________________________________________________________________________________
# User Parameters
filename="$1"
startnum="$2"
stopnum="$3"
inputfilepath="$4" # don't include .root

#____________________________________________________________________________________________________
# Main
if [[ -e "$filename" ]]; then rm "$filename"; fi

touch $filename

for num in $(seq ${startnum} ${stopnum}); do
    echo "${inputfilepath}${num}.root" >> $filename
done
