# Installation

Framework for running on root ntuples for various CMS physics analysis. 

Make sure that UF-PyNTupleRunner is up-to-date:
```bash
cd <path/to/UF-PyNTupleRunner/>
git pull origin master
cd Core
git pull 
```

Then clone the Statistical FrameWork repo into your UF-PyNTupleRunner directory:
```
cd <path/to/UF-PyNTupleRunner/>
git clone git@github.com:lucien1011/DarkZ-StatFW.git
cd DarkZ-StatFW
git submodule init
git submodule update
```

# Preparation for ZZd from scratch

Before running the statistic framework, one needs to prepare the histograms for stats processing:
* fine mZ2 distributions in TH1 format, this is done by running the configuration file 
```
cd <path/to/UF-PyNTupleRunner/>
source setup.sh /raid/raid7
UFNTuple DarkZ/para_input_Run2016_cfg.py
UFNTuple DarkZ/para_input_Run2017_cfg.py
UFNTuple DarkZ/para_input_Run2018_cfg.py 
```
Example output can be found here:
```
/raid/raid7/lucien/Higgs/DarkZ/ParaInput/DarkPhotonSelection_m4l100To170_Nominal/2019-08-14_m4lSR-m4lSB_HZZd-ppZZd_Run2016/
```
Interpolation of signal: 
* Run DarkZ/Script/interpolate_sig.py in UF-NTupleRunner. 
* This script takes fine mZ2 histograms from the previous step as inputs.
```
python interpolate_sig.py
```
* Interpolation of Z+X background: this is done by running make_ZPlusX_ParaInput_Run201X.py and fit_Shape_Run201X_cfg.py. 
Running fit_Shape_Run201X_cfg.py will make a TF1 for each category (e.g. 4e,4mu,2e2mu,2mu2e). These TF1 will be used 
by make_ZPlusX_ParaInput_Run201X.py to make TH1 histograms in proper format to be read by the statistic framework.
```
UFNTuple plot_Shape_SR_Run2016_cfg.py
python fit_Shape_Run201X_cfg.py
python make_ZPlusX_ParaInput_Run201X.py
```
* Systematic file: Located under Config/. These text files control how much uncertainties you want to put for your analysis.

# Running ZZd
A shell script has been prepared to run all steps. You have to enter the folder paths from the 
preparation step in this shell script.
```bash
source setup.sh
source mk_comb_all.sh
```
There are three steps in this script. First it makes a set of data cards for each signal model (
One folder for one signal model, and inside this folder you will have a set of data cards for each category and maybe each year too). 
Second, all the data cards in each signal model folder will be combined together to make one giant data card.
And from this giant data card it will further make a workspace (A workspace is a RooFit class which contains everything 
you need to perform statistical analysis, e.g. likelihood model, parameters etc). Then the final step will perform 
statistical analysis. You can choose to do each statisitical method yourself. Right now the framework 
supports performing upper limits and significance.
