import sys, os, string, re, pwd, commands, ast, optparse, shlex, time
from array import array
from math import *
from decimal import *
from sample_shortnames import *

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)
    
### Define function for parsing options
def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    # input options
    parser.add_option('-d', '--dir',    dest='SOURCEDIR',  type='string',default='./', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--unfoldModel',dest='UNFOLD',type='string',default='ggH_powheg15_JHUgen_125', help='Name of the unfolding model for central value')
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "inclusive", "pT", "eta", "Njets"')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.0',   help='Mass value for theory prediction')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='Use results from fixed fraction fit, default is False')     
    parser.add_option('',   '--setLog', action='store_true', dest='SETLOG', default=False, help='set plot to log scale y, default is False')     
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')    
    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)
                        
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()
sys.argv = grootargs

if (not os.path.exists("plots")):
    os.system("mkdir plots")

from ROOT import *

from tdrStyle import *
setTDRStyle()
        
datamodel = opt.UNFOLD

sys.path.append('./datacardInputs')
#sys.path.append('./old_lowstathres_unc')
#sys.path.append('./hres_histat')
#sys.path.append('./hres_histat_newbaseline')
#sys.path.append('./hres_histat_newbaseline_bugfix')
#sys.path.append('./hres_histat_newbaseline_bugfix2')
sys.path.append('./hres_histat_newbaseline_bugfix3')
sys.path.append('./lhscan_may18')


def plotXS(obsName, obs_bins):    

    _temp = __import__('inputs_sig_'+obsName, globals(), locals(), ['acc'], -1)
    acc = _temp.acc 
    _temp = __import__('higgs_xsbr', globals(), locals(), ['higgs_xs','higgs4l_br'], -1)
    higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br
    if (opt.FIXFRAC): floatfix = '_fixfrac'
    else: floatfix = ''
    if (obsName == "mass4l"):
        _temp = __import__('resultsXS_'+obsName+'_v3'+floatfix, globals(), locals(), ['modelNames', 'asimovDataModelName', 'resultsXS', 'modelIndUncert'], -1)
        modelNames = _temp.modelNames
        asimovDataModelName = _temp.asimovDataModelName
        resultsXS = _temp.resultsXS
        modelIndUncert = _temp.modelIndUncert
        _temp = __import__('resultsXS_'+obsName+'_v2'+floatfix, globals(), locals(), ['modelNames', 'asimovDataModelName', 'resultsXS', 'modelIndUncert'], -1)
        modelNames_v2 = _temp.modelNames
        asimovDataModelName_v2 = _temp.asimovDataModelName
        resultsXS_v2 = _temp.resultsXS
        modelIndUncert_v2 = _temp.modelIndUncert                                
    else:
        _temp = __import__('resultsXS_'+obsName+'_v3'+floatfix, globals(), locals(), ['modelNames', 'asimovDataModelName', 'resultsXS', 'modelIndUncert'], -1)
        modelNames = _temp.modelNames
        asimovDataModelName = _temp.asimovDataModelName
        resultsXS = _temp.resultsXS
        modelIndUncert = _temp.modelIndUncert
    _temp = __import__('resultsXS_LHScan_'+obsName+'_v3'+floatfix, globals(), locals(), ['resultsXS_LHScan'], -1)
    resultsXS_LHScan = _temp.resultsXS_LHScan
    
    acc_ggH_powheg = {}
    pdfunc_ggH_powheg = {}
    qcdunc_ggH_powheg = {}
    _temp = __import__('ggH_powheg_'+obsName.replace('_reco','_gen')+'_4e_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
    acc_ggH_powheg['4e'] = _temp.acc
    pdfunc_ggH_powheg['4e'] = _temp.pdfUncert
    qcdunc_ggH_powheg['4e'] = _temp.qcdUncert
    _temp = __import__('ggH_powheg_'+obsName.replace('_reco','_gen')+'_4mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
    acc_ggH_powheg['4mu'] = _temp.acc
    pdfunc_ggH_powheg['4mu'] = _temp.pdfUncert
    qcdunc_ggH_powheg['4mu'] = _temp.qcdUncert
    _temp = __import__('ggH_powheg_'+obsName.replace('_reco','_gen')+'_2e2mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
    acc_ggH_powheg['2e2mu'] = _temp.acc
    pdfunc_ggH_powheg['2e2mu'] = _temp.pdfUncert
    qcdunc_ggH_powheg['2e2mu'] = _temp.qcdUncert
                
    acc_ggH_minloHJ = {}
    pdfunc_ggH_minloHJ = {}
    qcdunc_ggH_minloHJ = {}
    _temp = __import__('ggH_minloHJ_weightedSum_'+obsName.replace('_reco','_gen')+'_4e_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
    acc_ggH_minloHJ['4e'] = _temp.acc
    pdfunc_ggH_minloHJ['4e'] = _temp.pdfUncert
    qcdunc_ggH_minloHJ['4e'] = _temp.qcdUncert
    _temp = __import__('ggH_minloHJ_weightedSum_'+obsName.replace('_reco','_gen')+'_4mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
    acc_ggH_minloHJ['4mu'] = _temp.acc
    pdfunc_ggH_minloHJ['4mu'] = _temp.pdfUncert
    qcdunc_ggH_minloHJ['4mu'] = _temp.qcdUncert
    _temp = __import__('ggH_minloHJ_weightedSum_'+obsName.replace('_reco','_gen')+'_2e2mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
    acc_ggH_minloHJ['2e2mu'] = _temp.acc
    pdfunc_ggH_minloHJ['2e2mu'] = _temp.pdfUncert
    qcdunc_ggH_minloHJ['2e2mu'] = _temp.qcdUncert

#    if (obsName=="pT4l" or obsName=="rapidity4l"):
    if (not ("jet" in obsName)):
        useHRes = True
        HRes = 'HRES'
    else:
        useHRes = False
        HRes = 'powheg'

    if obsName == 'mass4l': HRes = 'HRes'
    
    if (useHRes):
        acc_ggH_HRes = {}
        pdfunc_ggH_HRes = {}
        qcdunc_ggH_HRes = {}
        _temp = __import__('ggH_HRES_'+obsName.replace('_reco','_gen')+'_4e_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
        acc_ggH_HRes['4e'] = _temp.acc
        pdfunc_ggH_HRes['4e'] = _temp.pdfUncert
        qcdunc_ggH_HRes['4e'] = _temp.qcdUncert
        _temp = __import__('ggH_HRES_'+obsName.replace('_reco','_gen')+'_4mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
        acc_ggH_HRes['4mu'] = _temp.acc
        pdfunc_ggH_HRes['4mu'] = _temp.pdfUncert
        qcdunc_ggH_HRes['4mu'] = _temp.qcdUncert
        _temp = __import__('ggH_HRES_'+obsName.replace('_reco','_gen')+'_2e2mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
        acc_ggH_HRes['2e2mu'] = _temp.acc
        pdfunc_ggH_HRes['2e2mu'] = _temp.pdfUncert
        qcdunc_ggH_HRes['2e2mu'] = _temp.qcdUncert
    else:
        acc_ggH_HRes = {}
        pdfunc_ggH_HRes = {}
        qcdunc_ggH_HRes = {}
        _temp = __import__('ggH_powheg_'+obsName.replace('_reco','_gen')+'_4e_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
        acc_ggH_HRes['4e'] = _temp.acc
        pdfunc_ggH_HRes['4e'] = _temp.pdfUncert
        qcdunc_ggH_HRes['4e'] = _temp.qcdUncert
        _temp = __import__('ggH_powheg_'+obsName.replace('_reco','_gen')+'_4mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
        acc_ggH_HRes['4mu'] = _temp.acc
        pdfunc_ggH_HRes['4mu'] = _temp.pdfUncert
        qcdunc_ggH_HRes['4mu'] = _temp.qcdUncert
        _temp = __import__('ggH_powheg_'+obsName.replace('_reco','_gen')+'_2e2mu_unc', globals(), locals(), ['acc','pdfUncert','qcdUncert'], -1)
        acc_ggH_HRes['2e2mu'] = _temp.acc
        pdfunc_ggH_HRes['2e2mu'] = _temp.pdfUncert
        qcdunc_ggH_HRes['2e2mu'] = _temp.qcdUncert

    # cross sections
    ggH_powheg15_JHUgen = []
    ggH_minloHJJ = []
    ggH_powheg = []
    ggH_powheg_unc_hi = []
    ggH_powheg_unc_lo = []
    ggH_minloHJ = []
    ggH_minloHJ_unc_hi = []
    ggH_minloHJ_unc_lo = []
    ggH_HRes = []
    ggH_HRes_unc_hi = []
    ggH_HRes_unc_lo = []
    # NNLO theory unc
    ggH_powheg15_JHUgen_NNLOunc_hi = []
    ggH_powheg15_JHUgen_NNLOunc_lo = []
    ggH_powheg_NNLOunc_hi = []
    ggH_powheg_NNLOunc_lo = []
    ggH_minloHJJ_NNLOunc_hi = []
    ggH_minloHJJ_NNLOunc_lo = []
    ggH_minloHJ_NNLOunc_hi = []
    ggH_minloHJ_NNLOunc_lo = []
    ggH_HRes_NNLOunc_hi = []
    ggH_HRes_NNLOunc_lo = []
    # NLO theory unc
    ggH_powheg15_JHUgen_NLOunc_hi = []
    ggH_powheg15_JHUgen_NLOunc_lo = []
    ggH_powheg_NLOunc_hi = []
    ggH_powheg_NLOunc_lo = []
    ggH_minloHJ_NLOunc_hi = []
    ggH_minloHJ_NLOunc_lo = []
    # XH unc
    XH = []
    XH_unc = []
    # Data
    data = []
    data_hi = []
    data_lo = []
    asimovdata = []
    # Systematic unc.
    systematics_hi = []
    systematics_lo = []
    modeldep_hi = []
    modeldep_lo = []
    data_hi_allunc = []
    data_lo_allunc = []
    
    #process ggH qqH WH ZH ttH bkg_qqzz bkg_ggzz bkg_zjets
    #pdf_gg lnN 1.0720 - - - 1.0780 - 1.0710 -
    #pdf_qqbar lnN - 1.0270 1.0350 1.0350 - 1.0342 - -
    #pdf_hzz4l_accept lnN 1.02 1.02 1.02 1.02 1.02 - - -
    #QCDscale_ggH lnN 1.0750 - - - - - - -
    #QCDscale_qqH lnN - 1.0020 - - - - - -
    #QCDscale_VH lnN - - 1.0040 1.0155 - - - -
    #QCDscale_ttH lnN - - - - 1.0655 - - -
    #QCDscale_ggVV lnN - - - - - - 1.2435 -
    #BRhiggs_hzz4l lnN 1.02 1.02 1.02 1.02 1.02 - - -
    unc_theory_ggH_hi = sqrt(0.072**2+0.075**2+0.02**2+0.02**2)
    unc_theory_ggH_lo = sqrt(0.078**2+0.069**2+0.02**2+0.02**2)
    unc_theory_XH_hi  = sqrt(0.027**2+0.02**2+0.002**2+0.02**2)
    unc_theory_XH_lo  = unc_theory_XH_hi
    unc_VBF = sqrt(0.027**2+0.02**2+0.002**2+0.02**2)
    unc_WH = sqrt(0.035**2+0.02**2+0.004**2+0.02**2)
    unc_ZH = sqrt(0.035**2+0.02**2+0.0155**2+0.02**2)
    unc_ttH =  sqrt(0.078**2+0.02**2+0.0655**2+0.02**2)

    unc_acc = 0.02
    unc_br = 0.02

    #unc_pdf_ggH_hi = 0.075
    unc_pdf_ggH_hi = 0.072
    unc_pdf_ggH_lo = 0.069
    unc_pdf_VBF = 0.027
    unc_pdf_WH = 0.023
    unc_pdf_ZH = 0.025
    unc_pdf_ttH = 0.081

    #unc_qcd_ggH_hi = 0.072
    unc_qcd_ggH_hi = 0.075
    unc_qcd_ggH_lo = 0.078
    unc_qcd_VBF = 0.002
    unc_qcd_WH = 0.01
    unc_qcd_ZH = 0.031
    unc_qcd_ttH = 0.0655

    nBins=len(obs_bins)
    for obsBin in range(nBins-1):

        # theory cross sections
        ggH_powheg15_JHUgen.append(0.0)
        ggH_minloHJJ.append(0.0)
        ggH_powheg.append(0.0)
        ggH_powheg_unc_hi.append(0.0)
        ggH_powheg_unc_lo.append(0.0)
        ggH_minloHJ.append(0.0)
        ggH_minloHJ_unc_hi.append(0.0)
        ggH_minloHJ_unc_lo.append(0.0)
        ggH_HRes.append(0.0)
        ggH_HRes_unc_hi.append(0.0)
        ggH_HRes_unc_lo.append(0.0)
        # NNLO theory unc
        ggH_powheg15_JHUgen_NNLOunc_hi.append(0.0)
        ggH_powheg15_JHUgen_NNLOunc_lo.append(0.0)
        ggH_powheg_NNLOunc_hi.append(0.0)
        ggH_powheg_NNLOunc_lo.append(0.0)
        ggH_minloHJJ_NNLOunc_hi.append(0.0)
        ggH_minloHJJ_NNLOunc_lo.append(0.0)
        ggH_minloHJ_NNLOunc_hi.append(0.0)
        ggH_minloHJ_NNLOunc_lo.append(0.0)
        ggH_HRes_NNLOunc_hi.append(0.0)
        ggH_HRes_NNLOunc_lo.append(0.0)
        # NLO theory unc
        ggH_powheg15_JHUgen_NLOunc_hi.append(0.0)
        ggH_powheg15_JHUgen_NLOunc_lo.append(0.0)
        ggH_powheg_NLOunc_hi.append(0.0)
        ggH_powheg_NLOunc_lo.append(0.0)
        ggH_minloHJ_NLOunc_hi.append(0.0)
        ggH_minloHJ_NLOunc_lo.append(0.0)
        # XH 
        XH.append(0.0)
        XH_unc.append(0.0)
        # Data
        data.append(0.0)
        data_hi.append(0.0)
        data_lo.append(0.0)
        asimovdata.append(0.0)
        # Systematic unc
        modeldep_hi.append(0.0)
        modeldep_lo.append(0.0)
        systematics_hi.append(0.0)
        systematics_lo.append(0.0)
        data_hi_allunc.append(0.0)
        data_lo_allunc.append(0.0)

        for channel in ['4e','4mu','2e2mu']:

            XH_fs = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']            

            XH[obsBin]+=XH_fs
            #XH_unc[obsBin]+= unc_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            #XH_unc[obsBin]+= unc_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            #XH_unc[obsBin]+= unc_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            #XH_unc[obsBin]+= unc_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']

            # branching ratio uncertainty
            XH_unc_fs = (unc_br*XH_fs)**2
            # acceptance uncertainty
            XH_unc_fs += (unc_acc*XH_fs)**2

            # qcd scale
            XH_qcdunc_fs = (unc_qcd_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2                
            XH_unc_fs += XH_qcdunc_fs

            # pdf
            XH_qqpdfunc_fs = (unc_pdf_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                              +unc_pdf_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                              +unc_pdf_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2            
            XH_unc_fs += XH_qqpdfunc_fs
            
            # add pdf uncertainty for ttH to total XH uncertainty
            XH_unc_fs += (unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            # total XH uncertainty
            XH_unc[obsBin]+=sqrt(XH_unc_fs)

            # ggH cross sections
            ggH_xsBR = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]

            ggH_powheg15_JHUgen[obsBin]+=ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            ggH_powheg[obsBin]+=ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            ggH_minloHJJ[obsBin]+=ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
            ggH_minloHJ[obsBin]+=ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            #ggH_HRes[obsBin]+=ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            ggH_HRes[obsBin]+=acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]

            # for total uncertainty, correlate br and acc uncertainties across all channels (XH+ggH)
            total_NNLOunc_fs_powheg15_JHUgen_hi = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2
            total_NNLOunc_fs_powheg15_JHUgen_lo = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2
            total_NNLOunc_fs_powheg15_JHUgen_hi += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2
            total_NNLOunc_fs_powheg15_JHUgen_lo += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2

            total_NNLOunc_fs_powheg_hi =  (unc_br*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_powheg_lo =  (unc_br*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_powheg_hi +=  (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_powheg_lo +=  (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2

            total_NNLOunc_fs_minloHJJ_hi = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2
            total_NNLOunc_fs_minloHJJ_lo = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2
            total_NNLOunc_fs_minloHJJ_hi += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2
            total_NNLOunc_fs_minloHJJ_lo += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']))**2

            total_NNLOunc_fs_minloHJ_hi = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_minloHJ_lo = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_minloHJ_hi += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_minloHJ_lo += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2

            #total_NNLOunc_fs_HRes_hi = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            #total_NNLOunc_fs_HRes_lo = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            #total_NNLOunc_fs_HRes_hi += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            #total_NNLOunc_fs_HRes_lo += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_HRes_hi = (unc_br*(XH_fs+acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_HRes_lo = (unc_br*(XH_fs+acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_HRes_hi += (unc_acc*(XH_fs+acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2
            total_NNLOunc_fs_HRes_lo += (unc_acc*(XH_fs+acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]))**2

            # NLO and NNLO are the same at this point
            total_NLOunc_fs_powheg15_JHUgen_hi = total_NNLOunc_fs_powheg15_JHUgen_hi
            total_NLOunc_fs_powheg15_JHUgen_lo = total_NNLOunc_fs_powheg15_JHUgen_lo
            total_NLOunc_fs_powheg_hi = total_NNLOunc_fs_powheg_hi
            total_NLOunc_fs_powheg_lo = total_NNLOunc_fs_powheg_lo
            total_NLOunc_fs_minloHJ_hi = total_NNLOunc_fs_minloHJ_hi
            total_NLOunc_fs_minloHJ_lo = total_NNLOunc_fs_minloHJ_lo
            

            # add ggH qcd uncertainties (uncorrelated with anything else)            
            #NNLO
            total_NNLOunc_fs_powheg15_JHUgen_hi += XH_qcdunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_lo += XH_qcdunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_powheg15_JHUgen_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            
            total_NNLOunc_fs_powheg_hi += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_powheg_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2

            total_NNLOunc_fs_minloHJJ_hi += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJJ_lo += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJJ_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_minloHJJ_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            
            total_NNLOunc_fs_minloHJ_hi += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJ_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_minloHJ_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2

            total_NNLOunc_fs_HRes_hi += XH_qcdunc_fs
            total_NNLOunc_fs_HRes_lo += XH_qcdunc_fs
            #total_NNLOunc_fs_HRes_hi += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
            #                               *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            #total_NNLOunc_fs_HRes_lo += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
            #                               *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_HRes_hi += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                           *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_HRes_lo += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                           *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2

            
            #NLO
            total_NLOunc_fs_powheg15_JHUgen_hi += XH_qcdunc_fs
            total_NLOunc_fs_powheg15_JHUgen_lo += XH_qcdunc_fs
            total_NLOunc_fs_powheg15_JHUgen_hi += (qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                                   *ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NLOunc_fs_powheg15_JHUgen_lo += (qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                                   *ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            
            total_NLOunc_fs_powheg_hi += XH_qcdunc_fs
            total_NLOunc_fs_powheg_lo += XH_qcdunc_fs
            total_NLOunc_fs_powheg_hi += (qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NLOunc_fs_powheg_lo += (qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2

            total_NLOunc_fs_minloHJ_hi += XH_qcdunc_fs
            total_NLOunc_fs_minloHJ_lo += XH_qcdunc_fs
            total_NLOunc_fs_minloHJ_hi += (qcdunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                           *ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NLOunc_fs_minloHJ_lo += (qcdunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                           *ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            
            # add pdf unc, anti correlate ggH and ttH
            #NNLO
            total_NNLOunc_fs_powheg15_JHUgen_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_powheg15_JHUgen_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_powheg_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_powheg_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_minloHJJ_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJJ_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJJ_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_minloHJJ_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_minloHJ_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_minloHJ_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_HRes_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_HRes_lo += XH_qqpdfunc_fs
            #total_NNLOunc_fs_HRes_hi += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
            #                              *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            #                              -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            #total_NNLOunc_fs_HRes_lo += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
            #                              *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            #                              -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_HRes_hi += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_HRes_lo += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            
            #NLO
            total_NLOunc_fs_powheg15_JHUgen_hi += XH_qqpdfunc_fs
            total_NLOunc_fs_powheg15_JHUgen_lo += XH_qqpdfunc_fs
            total_NLOunc_fs_powheg15_JHUgen_hi += (pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                                   *ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                   -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NLOunc_fs_powheg15_JHUgen_lo += (pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                                   *ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                   -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NLOunc_fs_powheg_hi += XH_qqpdfunc_fs
            total_NLOunc_fs_powheg_lo += XH_qqpdfunc_fs
            total_NLOunc_fs_powheg_hi += (pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NLOunc_fs_powheg_lo += (pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NLOunc_fs_minloHJ_hi += XH_qqpdfunc_fs
            total_NLOunc_fs_minloHJ_lo += XH_qqpdfunc_fs
            total_NLOunc_fs_minloHJ_hi += (pdfunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NLOunc_fs_minloHJ_lo += (pdfunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            # finally total uncertainty (different final states are correlated)
            # NNLO
            ggH_powheg15_JHUgen_NNLOunc_hi[obsBin]+=sqrt(total_NNLOunc_fs_powheg15_JHUgen_hi)
            ggH_powheg15_JHUgen_NNLOunc_lo[obsBin]+=sqrt(total_NNLOunc_fs_powheg15_JHUgen_lo)
            ggH_powheg_NNLOunc_hi[obsBin]+=sqrt(total_NNLOunc_fs_powheg_hi)
            ggH_powheg_NNLOunc_lo[obsBin]+=sqrt(total_NNLOunc_fs_powheg_lo) 
            ggH_minloHJJ_NNLOunc_hi[obsBin]+=sqrt(total_NNLOunc_fs_minloHJJ_hi)
            ggH_minloHJJ_NNLOunc_lo[obsBin]+=sqrt(total_NNLOunc_fs_minloHJJ_lo) 
            ggH_minloHJ_NNLOunc_hi[obsBin]+=sqrt(total_NNLOunc_fs_minloHJ_hi)
            ggH_minloHJ_NNLOunc_lo[obsBin]+=sqrt(total_NNLOunc_fs_minloHJ_lo) 
            ggH_HRes_NNLOunc_hi[obsBin]+=sqrt(total_NNLOunc_fs_HRes_hi)
            ggH_HRes_NNLOunc_lo[obsBin]+=sqrt(total_NNLOunc_fs_HRes_lo) 
            # NLO
            ggH_powheg15_JHUgen_NLOunc_hi[obsBin]+=sqrt(total_NLOunc_fs_powheg15_JHUgen_hi)
            ggH_powheg15_JHUgen_NLOunc_lo[obsBin]+=sqrt(total_NLOunc_fs_powheg15_JHUgen_lo)
            ggH_powheg_NLOunc_hi[obsBin]+=sqrt(total_NLOunc_fs_powheg_hi)
            ggH_powheg_NLOunc_lo[obsBin]+=sqrt(total_NLOunc_fs_powheg_lo) 
            ggH_minloHJ_NLOunc_hi[obsBin]+=sqrt(total_NLOunc_fs_minloHJ_hi)
            ggH_minloHJ_NLOunc_lo[obsBin]+=sqrt(total_NLOunc_fs_minloHJ_lo) 
                       
            # OLD WAY
            ggH_powheg_unc_rel_hi = sqrt(pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2
                                         +qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2)
            ggH_powheg_unc_rel_lo = sqrt(pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2
                                         +qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2)
            ggH_powheg_unc_hi[obsBin]+= ggH_powheg_unc_rel_hi*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            ggH_powheg_unc_lo[obsBin]+= ggH_powheg_unc_rel_lo*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            
            ggH_minloHJ_unc_rel_hi = sqrt(pdfunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2
                                          +qcdunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2)
            ggH_minloHJ_unc_rel_lo = sqrt(pdfunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2
                                          +qcdunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2)
            ggH_minloHJ_unc_hi[obsBin]+= ggH_minloHJ_unc_rel_hi*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            ggH_minloHJ_unc_lo[obsBin]+= ggH_minloHJ_unc_rel_lo*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            # END OLD WAY

        ggH_powheg15_JHUgen[obsBin]+=XH[obsBin]
        ggH_minloHJJ[obsBin]+=XH[obsBin]
        ggH_powheg[obsBin]+=XH[obsBin]
        ggH_minloHJ[obsBin]+=XH[obsBin]
        ggH_HRes[obsBin]+=XH[obsBin]

        # OLD WAY
        ggH_powheg_unc_hi[obsBin] += XH_unc[obsBin]
        ggH_powheg_unc_lo[obsBin] += XH_unc[obsBin]
        ggH_minloHJ_unc_hi[obsBin] += XH_unc[obsBin]
        ggH_minloHJ_unc_lo[obsBin] += XH_unc[obsBin]
        # END OLD WAY

        if (opt.UNBLIND):
            data[obsBin] = resultsXS[datamodel+"_"+obsName+"_genbin"+str(obsBin)]["central"]        
            #data_hi[obsBin] = resultsXS[datamodel+"_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            #data_lo[obsBin] = -1.0*resultsXS[datamodel+"_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]        
            data_hi[obsBin] = resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            data_lo[obsBin] = -1.0*resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]
        else:
            data[obsBin] = resultsXS["AsimovData_"+obsName+"_genbin"+str(obsBin)]["central"]        
            #data_hi[obsBin] = resultsXS["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            #data_lo[obsBin] = -1.0*resultsXS["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]        
            data_hi[obsBin] = resultsXS["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            data_lo[obsBin] = -1.0*resultsXS["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]

        if (opt.UNBLIND):
            modeldep_hi[obsBin] = modelIndUncert["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            modeldep_lo[obsBin] = -1.0*modelIndUncert["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"] 
            #systematics_hi[obsBin] = sqrt(max(0.0,resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2))
            #systematics_lo[obsBin] = sqrt(max(0.0,resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2))
            systematics_hi[obsBin] = sqrt(max(0.0,resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2))
            systematics_lo[obsBin] = sqrt(max(0.0,resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2))
        else:
            modeldep_hi[obsBin] = modelIndUncert[asimovDataModelName+"_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            modeldep_lo[obsBin] = -1.0*modelIndUncert[asimovDataModelName+"_"+obsName+"_genbin"+str(obsBin)]["uncerDn"] 
            #modeldep_hi[obsBin] = modelIndUncert["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            #modeldep_lo[obsBin] = -1.0*modelIndUncert["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]             
            systematics_hi[obsBin] = sqrt(max(0.0,resultsXS["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2))
            systematics_lo[obsBin] = sqrt(max(0.0,resultsXS["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2))
            #systematics_hi[obsBin] = sqrt(max(0.0,resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2))
            #systematics_lo[obsBin] = sqrt(max(0.0,resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2))

        data_hi_allunc[obsBin] = sqrt(data_hi[obsBin]**2+modeldep_hi[obsBin]**2)
        data_lo_allunc[obsBin] = sqrt(data_lo[obsBin]**2+modeldep_lo[obsBin]**2)

    if (obsName=="mass4l"):
        for channel in ['2e2mu','4mu','4e']:

            # theory cross sections
            ggH_powheg15_JHUgen.append(0.0)
            ggH_minloHJJ.append(0.0)
            ggH_powheg.append(0.0)
            ggH_powheg_unc_hi.append(0.0)
            ggH_powheg_unc_lo.append(0.0)
            ggH_minloHJ.append(0.0)
            ggH_minloHJ_unc_hi.append(0.0)
            ggH_minloHJ_unc_lo.append(0.0)
            ggH_HRes.append(0.0)
            ggH_HRes_unc_hi.append(0.0)
            ggH_HRes_unc_lo.append(0.0)
            # NNLO theory unc
            ggH_powheg15_JHUgen_NNLOunc_hi.append(0.0)
            ggH_powheg15_JHUgen_NNLOunc_lo.append(0.0)
            ggH_powheg_NNLOunc_hi.append(0.0)
            ggH_powheg_NNLOunc_lo.append(0.0)
            ggH_minloHJJ_NNLOunc_hi.append(0.0)
            ggH_minloHJJ_NNLOunc_lo.append(0.0)
            ggH_minloHJ_NNLOunc_hi.append(0.0)
            ggH_minloHJ_NNLOunc_lo.append(0.0)
            ggH_HRes_NNLOunc_hi.append(0.0)
            ggH_HRes_NNLOunc_lo.append(0.0)
            # NLO theory unc
            ggH_powheg15_JHUgen_NLOunc_hi.append(0.0)
            ggH_powheg15_JHUgen_NLOunc_lo.append(0.0)
            ggH_powheg_NLOunc_hi.append(0.0)
            ggH_powheg_NLOunc_lo.append(0.0)
            ggH_minloHJ_NLOunc_hi.append(0.0)
            ggH_minloHJ_NLOunc_lo.append(0.0)
            # XH 
            XH.append(0.0)
            XH_unc.append(0.0)
            # Data
            data.append(0.0)
            data_hi.append(0.0)
            data_lo.append(0.0)
            # Systematic unc
            modeldep_hi.append(0.0)
            modeldep_lo.append(0.0)
            systematics_hi.append(0.0)
            systematics_lo.append(0.0)
            data_hi_allunc.append(0.0)
            data_lo_allunc.append(0.0)

            if (channel=='2e2mu'): bin = 1
            if (channel=='4mu'): bin = 2
            if (channel=='4e'): bin = 3

            obsBin=0
            
            XH_fs = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0']

            XH[bin]+=XH_fs

            #XH_unc[bin]+= unc_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            #XH_unc[bin]+= unc_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            #XH_unc[bin]+= unc_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            #XH_unc[bin]+= unc_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0']

            # branching ratio uncertainty
            XH_unc_fs = (unc_br*XH_fs)**2
            # acceptance uncertainty
            XH_unc_fs += (unc_acc*XH_fs)**2

            # qcd scale
            XH_qcdunc_fs = (unc_qcd_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2                
            XH_unc_fs += XH_qcdunc_fs

            # pdf
            XH_qqpdfunc_fs = (unc_pdf_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin0_recobin0']
                              +unc_pdf_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0']
                              +unc_pdf_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2            
            XH_unc_fs += XH_qqpdfunc_fs
            
            # add pdf uncertainty for ttH to total XH uncertainty
            XH_unc_fs += (unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2

            # total XH uncertainty
            XH_unc[bin]+=sqrt(XH_unc_fs)

            # ggH cross sections
            ggH_xsBR = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]

            ggH_powheg15_JHUgen[bin]+=ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            ggH_powheg[bin]+=ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            ggH_minloHJJ[bin]+=ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin0_recobin0']            
            ggH_minloHJ[bin]+=ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]            
            ggH_HRes[bin]+=acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]            

             # for total uncertainty, correlate br and acc uncertainties across all channels (XH+ggH)
            total_NNLOunc_fs_powheg15_JHUgen_hi = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_powheg15_JHUgen_lo = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_powheg15_JHUgen_hi += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_powheg15_JHUgen_lo += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2

            total_NNLOunc_fs_powheg_hi =  (unc_br*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_powheg_lo =  (unc_br*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_powheg_hi +=  (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_powheg_lo +=  (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2

            total_NNLOunc_fs_minloHJJ_hi = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_minloHJJ_lo = (unc_br*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_minloHJJ_hi += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_minloHJJ_lo += (unc_acc*(XH_fs+ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2

            total_NNLOunc_fs_minloHJ_hi = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_minloHJ_lo = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_minloHJ_hi += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_minloHJ_lo += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2

            #total_NNLOunc_fs_HRes_hi = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            #total_NNLOunc_fs_HRes_lo = (unc_br*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            #total_NNLOunc_fs_HRes_hi += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            #total_NNLOunc_fs_HRes_lo += (unc_acc*(XH_fs+ggH_xsBR*acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_HRes_hi = (unc_br*(XH_fs+acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_HRes_lo = (unc_br*(XH_fs+acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_HRes_hi += (unc_acc*(XH_fs+acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
            total_NNLOunc_fs_HRes_lo += (unc_acc*(XH_fs+acc_ggH_HRes[channel]['ggH_HRes_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin0']))**2
           
            # add ggH qcd uncertainties (uncorrelated with anything else)            
            #NNLO
            total_NNLOunc_fs_powheg15_JHUgen_hi += XH_qcdunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_lo += XH_qcdunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_powheg15_JHUgen_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            
            total_NNLOunc_fs_powheg_hi += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_powheg_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2

            total_NNLOunc_fs_minloHJJ_hi += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJJ_lo += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJJ_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_minloHJJ_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            
            total_NNLOunc_fs_minloHJ_hi += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJ_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_minloHJ_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2

            total_NNLOunc_fs_HRes_hi += XH_qcdunc_fs
            total_NNLOunc_fs_HRes_lo += XH_qcdunc_fs
            #total_NNLOunc_fs_HRes_hi += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
            #                               *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            #total_NNLOunc_fs_HRes_lo += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
            #                               *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_HRes_hi += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                           *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2
            total_NNLOunc_fs_HRes_lo += (qcdunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                           *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)])**2

            # add pdf unc, anti correlate ggH and ttH
            #NNLO
            total_NNLOunc_fs_powheg15_JHUgen_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg15_JHUgen_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_powheg15_JHUgen_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc['ggH_powheg15_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_powheg_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_powheg_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_minloHJJ_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJJ_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJJ_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_minloHJJ_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc['ggH_minloHJJ_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                                                    -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_minloHJ_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_minloHJ_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            total_NNLOunc_fs_HRes_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_HRes_lo += XH_qqpdfunc_fs
            #total_NNLOunc_fs_HRes_hi += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
            #                              *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            #                              -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            #total_NNLOunc_fs_HRes_lo += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
            #                              *ggH_xsBR*acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            #                              -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_HRes_hi += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2
            total_NNLOunc_fs_HRes_lo += (pdfunc_ggH_HRes[channel]["ggH_"+HRes+"_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *acc_ggH_HRes[channel]['ggH_'+HRes+'_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_pythia_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'])**2

            # finally total uncertainty (different final states are correlated)
            # NNLO
            ggH_powheg15_JHUgen_NNLOunc_hi[bin]+=sqrt(total_NNLOunc_fs_powheg15_JHUgen_hi)
            ggH_powheg15_JHUgen_NNLOunc_lo[bin]+=sqrt(total_NNLOunc_fs_powheg15_JHUgen_lo)
            ggH_powheg_NNLOunc_hi[bin]+=sqrt(total_NNLOunc_fs_powheg_hi)
            ggH_powheg_NNLOunc_lo[bin]+=sqrt(total_NNLOunc_fs_powheg_lo) 
            ggH_minloHJJ_NNLOunc_hi[bin]+=sqrt(total_NNLOunc_fs_minloHJJ_hi)
            ggH_minloHJJ_NNLOunc_lo[bin]+=sqrt(total_NNLOunc_fs_minloHJJ_lo) 
            ggH_minloHJ_NNLOunc_hi[bin]+=sqrt(total_NNLOunc_fs_minloHJ_hi)
            ggH_minloHJ_NNLOunc_lo[bin]+=sqrt(total_NNLOunc_fs_minloHJ_lo) 
            ggH_HRes_NNLOunc_hi[bin]+=sqrt(total_NNLOunc_fs_HRes_hi)
            ggH_HRes_NNLOunc_lo[bin]+=sqrt(total_NNLOunc_fs_HRes_lo) 

            # OLD way
            ggH_powheg_unc_rel_hi = sqrt(pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2
                                         +qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2)
            ggH_powheg_unc_rel_lo = sqrt(pdfunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2
                                         +qcdunc_ggH_powheg[channel]["ggH_powheg_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2)        
            ggH_minloHJ_unc_rel_hi = sqrt(pdfunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2
                                          +qcdunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']**2)
            ggH_minloHJ_unc_rel_lo = sqrt(pdfunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2
                                          +qcdunc_ggH_minloHJ[channel]["ggH_minloHJ_125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']**2)
            ggH_powheg_unc_hi[bin]+= ggH_powheg_unc_rel_hi*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            ggH_powheg_unc_lo[bin]+= ggH_powheg_unc_rel_lo*ggH_xsBR*acc_ggH_powheg[channel]['ggH_powheg_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]            
            ggH_minloHJ_unc_hi[bin]+= ggH_minloHJ_unc_rel_hi*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]
            ggH_minloHJ_unc_lo[bin]+= ggH_minloHJ_unc_rel_lo*ggH_xsBR*acc_ggH_minloHJ[channel]['ggH_minloHJ_125_'+channel+'_'+obsName.replace('_reco','_gen')+'_genbin'+str(obsBin)]

            ggH_powheg_unc_hi[bin] += XH_unc[bin]
            ggH_powheg_unc_lo[bin] += XH_unc[bin]
            ggH_minloHJ_unc_hi[bin] += XH_unc[bin]
            ggH_minloHJ_unc_lo[bin] += XH_unc[bin]
            # End OLD way

            ggH_powheg15_JHUgen[bin]+=XH[bin]
            ggH_minloHJJ[bin]+=XH[bin]
            ggH_powheg[bin]+=XH[bin]
            ggH_minloHJ[bin]+=XH[bin]
            ggH_HRes[bin]+=XH[bin]

            data[bin] = resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["central"]
            data_hi[bin] = resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
            data_lo[bin] = -1.0*resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
            
            if (opt.UNBLIND):
                modeldep_hi[bin] = modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
                modeldep_lo[bin] = -1.0*modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
                systematics_hi[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                systematics_lo[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)
            else:
                #modeldep_hi[bin] = modelIndUncert_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
                #modeldep_lo[bin] = -1.0*modelIndUncert_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
                #systematics_hi[bin] = sqrt(resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                #systematics_lo[bin] = sqrt(resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)
                modeldep_hi[bin] = modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
                modeldep_lo[bin] = -1.0*modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
                systematics_hi[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                systematics_lo[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)

            data_hi_allunc[bin] = sqrt(data_hi[bin]**2+modeldep_hi[bin]**2)
            data_lo_allunc[bin] = sqrt(data_lo[bin]**2+modeldep_lo[bin]**2)
            

        
    print 'data',data
    sumdata = 0.0
    for i in range(len(data)):
        sumdata+=data[i]
    print obsName,'sum data',sumdata
    print 'data_hi',data_hi
    print 'data_lo',data_lo
    print 'ggH HRes + XH',ggH_HRes
    print 'NNLO ggH HRes + XH',ggH_HRes_NNLOunc_hi
    print 'NNLO ggH HRes + XH',ggH_HRes_NNLOunc_lo
    print 'ggH powheg15+JHUgen + XH',ggH_powheg15_JHUgen
    print 'NLO ggH_powheg15_JHUgen_hi',ggH_powheg15_JHUgen_NLOunc_hi
    print 'NLO ggH_powheg15_JHUgen_lo',ggH_powheg15_JHUgen_NLOunc_lo
    print 'NNLO ggH_powheg15_JHUgen_hi',ggH_powheg15_JHUgen_NNLOunc_hi
    print 'NNLO ggH_powheg15_JHUgen_lo',ggH_powheg15_JHUgen_NNLOunc_lo
    print 'ggH_minloHJJ',ggH_minloHJJ
    print 'ggH_powheg',ggH_powheg
    print 'NLO ggH_powheg_hi',ggH_powheg_NLOunc_hi
    print 'NLO ggH_powheg_lo',ggH_powheg_NLOunc_lo
    print 'NNLO ggH_powheg_hi',ggH_powheg_NNLOunc_hi
    print 'NNLO ggH_powheg_lo',ggH_powheg_NNLOunc_lo
    print 'OLD ggH_powheg_hi',ggH_powheg_unc_hi
    print 'OLD ggH_powheg_lo',ggH_powheg_unc_lo
    print 'ggH_minloHJ',ggH_minloHJ
    print 'OLD ggH_minloHJ_hi',ggH_minloHJ_unc_hi
    print 'OLD ggH_minloHJ_lo',ggH_minloHJ_unc_lo
    print 'OLD theory precent. up',unc_theory_ggH_hi,'theory precent. down',unc_theory_ggH_lo
    print 'XH',XH
    print 'XH_unc',XH_unc
    print 'modedlep_hi',modeldep_hi
    print 'modeldep_lo',modeldep_lo
    print 'systematics_hi',systematics_hi
    print 'systematics_lo',systematics_lo
    
    if (obsName=="mass4l"):

        a_observable  = array('d',[0.5+i for i in range(0,4)])
        v_observable    = TVectorD(len(a_observable),a_observable)
        a_dobservable = array('d',[0.5 for i in range(0,4)])
        v_dobservable = TVectorD(len(a_dobservable),a_dobservable)
        
        a_zeros = array('d',[0.0 for i in range(0,4)])
        v_zeros = TVectorD(len(a_zeros),a_zeros)
        a_twos = array('d',[0.2*a_dobservable[i] for i in range(0,4)])
        v_twos = TVectorD(len(a_twos),a_twos)
        
        a_ggH_powheg15_JHUgen = array('d',[ggH_powheg15_JHUgen[0],ggH_powheg15_JHUgen[1],ggH_powheg15_JHUgen[2],ggH_powheg15_JHUgen[3]] )        
        v_ggH_powheg15_JHUgen = TVectorD(len(a_ggH_powheg15_JHUgen),a_ggH_powheg15_JHUgen)
        #a_ggH_powheg15_JHUgen_hi = array('d',[unc_theory_ggH_hi*(a_ggH_powheg15_JHUgen[i]-XH[i])+XH_unc[i] for i in range(len(a_ggH_powheg15_JHUgen))])                             
        a_ggH_powheg15_JHUgen_hi = array('d',[ggH_powheg15_JHUgen_NNLOunc_hi[i] for i in range(len(a_ggH_powheg15_JHUgen))])                             
        v_ggH_powheg15_JHUgen_hi = TVectorD(len(a_ggH_powheg15_JHUgen_hi),a_ggH_powheg15_JHUgen_hi)
        #a_ggH_powheg15_JHUgen_lo = array('d',[unc_theory_ggH_lo*(a_ggH_powheg15_JHUgen[i]-XH[i])+XH_unc[i] for i in range(len(a_ggH_powheg15_JHUgen))])
        a_ggH_powheg15_JHUgen_lo = array('d',[ggH_powheg15_JHUgen_NNLOunc_lo[i] for i in range(len(a_ggH_powheg15_JHUgen))])
        v_ggH_powheg15_JHUgen_lo = TVectorD(len(a_ggH_powheg15_JHUgen_lo),a_ggH_powheg15_JHUgen_lo)

        print 'a_ggH_powheg15JHUgen',a_ggH_powheg15_JHUgen
        print 'a_ggH_powheg15JHUgen_hi',a_ggH_powheg15_JHUgen_hi
        print 'a_ggH_powheg15JHUgen_lo',a_ggH_powheg15_JHUgen_lo

        a_ggH_minloHJJ = array('d',[ggH_minloHJJ[0],ggH_minloHJJ[1],ggH_minloHJJ[2],ggH_minloHJJ[3]] )
        v_ggH_minloHJJ = TVectorD(len(a_ggH_minloHJJ),a_ggH_minloHJJ)
        #a_ggH_minloHJJ_hi = array('d',[unc_theory_ggH_hi*(a_ggH_minloHJJ[i]-XH[i])+XH_unc[i] for i in range(len(a_ggH_minloHJJ))])
        a_ggH_minloHJJ_hi = array('d',[ggH_minloHJJ_NNLOunc_hi[i] for i in range(len(a_ggH_minloHJJ))])
        v_ggH_minloHJJ_hi = TVectorD(len(a_ggH_minloHJJ_hi),a_ggH_minloHJJ_hi)
        #a_ggH_minloHJJ_lo = array('d',[unc_theory_ggH_lo*(a_ggH_minloHJJ[i]-XH[i])+XH_unc[i] for i in range(len(a_ggH_minloHJJ))])
        a_ggH_minloHJJ_lo = array('d',[ggH_minloHJJ_NNLOunc_lo[i] for i in range(len(a_ggH_minloHJJ))])
        v_ggH_minloHJJ_lo = TVectorD(len(a_ggH_minloHJJ_lo),a_ggH_minloHJJ_lo)                                                

        print 'a_ggH_minloHJJ',a_ggH_minloHJJ
        print 'a_ggH_minloHJJ_hi',a_ggH_minloHJJ_hi
        print 'a_ggH_minloHJJ_lo',a_ggH_minloHJJ_lo

        a_ggH_powheg = array('d',[ggH_powheg[i] for i in range(len(ggH_powheg))])
        v_ggH_powheg = TVectorD(len(a_ggH_powheg),a_ggH_powheg)
        #a_ggH_powheg_unc_hi =  array('d',[unc_theory_ggH_hi*(ggH_powheg[i]-XH[i])+XH_unc[i] for i in range(len(ggH_powheg))])
        #a_ggH_powheg_unc_lo =  array('d',[unc_theory_ggH_lo*(ggH_powheg[i]-XH[i])+XH_unc[i] for i in range(len(ggH_powheg))])
        a_ggH_powheg_unc_hi =  array('d',[ggH_powheg_NNLOunc_hi[i] for i in range(len(ggH_powheg))])
        a_ggH_powheg_unc_lo =  array('d',[ggH_powheg_NNLOunc_lo[i] for i in range(len(ggH_powheg))])
        v_ggH_powheg_unc_hi = TVectorD(len(a_ggH_powheg_unc_hi),a_ggH_powheg_unc_hi)
        v_ggH_powheg_unc_lo = TVectorD(len(a_ggH_powheg_unc_lo),a_ggH_powheg_unc_lo)

        print 'a_ggH_powheg',a_ggH_powheg15_JHUgen
        print 'a_ggH_powheg_hi',a_ggH_powheg_unc_hi
        print 'a_ggH_powheg_lo',a_ggH_powheg_unc_lo

        a_ggH_minloHJ = array('d',[ggH_minloHJ[i] for i in range(len(ggH_minloHJ))])
        v_ggH_minloHJ = TVectorD(len(a_ggH_minloHJ),a_ggH_minloHJ)
        #a_ggH_minloHJ_unc_hi =  array('d',[unc_theory_ggH_hi*(ggH_minloHJ[i]-XH[i])+XH_unc[i] for i in range(len(ggH_minloHJ_unc_hi))])
        #a_ggH_minloHJ_unc_lo =  array('d',[unc_theory_ggH_lo*(ggH_minloHJ[i]-XH[i])+XH_unc[i] for i in range(len(ggH_minloHJ_unc_lo))])
        a_ggH_minloHJ_unc_hi =  array('d',[ggH_minloHJ_NNLOunc_hi[i] for i in range(len(ggH_minloHJ_unc_hi))])
        a_ggH_minloHJ_unc_lo =  array('d',[ggH_minloHJ_NNLOunc_lo[i] for i in range(len(ggH_minloHJ_unc_lo))])
        v_ggH_minloHJ_unc_hi = TVectorD(len(a_ggH_minloHJ_unc_hi),a_ggH_minloHJ_unc_hi)
        v_ggH_minloHJ_unc_lo = TVectorD(len(a_ggH_minloHJ_unc_lo),a_ggH_minloHJ_unc_lo)

        print 'a_ggH_minloHJ',a_ggH_minloHJ
        print 'a_ggH_minloHJ_hi',a_ggH_minloHJ_unc_hi
        print 'a_ggH_minloHJ_lo',a_ggH_minloHJ_unc_lo

        a_ggH_HRes = array('d',[ggH_HRes[i] for i in range(len(ggH_HRes))])
        v_ggH_HRes = TVectorD(len(a_ggH_HRes),a_ggH_HRes)
        #a_ggH_HRes_unc_hi =  array('d',[unc_theory_ggH_hi*(ggH_HRes[i]-XH[i])+XH_unc[i] for i in range(len(ggH_HRes_unc_hi))])
        #a_ggH_HRes_unc_lo =  array('d',[unc_theory_ggH_lo*(ggH_HRes[i]-XH[i])+XH_unc[i] for i in range(len(ggH_HRes_unc_lo))])
        a_ggH_HRes_unc_hi =  array('d',[ggH_HRes_NNLOunc_hi[i] for i in range(len(ggH_HRes_unc_hi))])
        a_ggH_HRes_unc_lo =  array('d',[ggH_HRes_NNLOunc_lo[i] for i in range(len(ggH_HRes_unc_lo))])
        v_ggH_HRes_unc_hi = TVectorD(len(a_ggH_HRes_unc_hi),a_ggH_HRes_unc_hi)
        v_ggH_HRes_unc_lo = TVectorD(len(a_ggH_HRes_unc_lo),a_ggH_HRes_unc_lo)

        print 'a_ggH_HRes',a_ggH_HRes
        print 'a_ggH_HRes_hi',a_ggH_HRes_unc_hi
        print 'a_ggH_HRes_lo',a_ggH_HRes_unc_lo

        a_XH = array('d',[XH[0],XH[1],XH[2],XH[3]])
        v_XH = TVectorD(len(a_XH),a_XH)

        a_XH_hi = array('d',[XH_unc[i] for i in range(len(a_XH))])
        v_XH_hi = TVectorD(len(a_XH_hi),a_XH_hi)
        
        a_XH_lo = array('d',[XH_unc[i] for i in range(len(a_XH))])
        v_XH_lo = TVectorD(len(a_XH_lo),a_XH_lo)

        a_data = array('d',[data[0],data[1],data[2],data[3]])
        v_data = TVectorD(len(a_data),a_data)
        a_data_hi = array('d',[data_hi[0],data_hi[1],data_hi[2],data_hi[3]])
        v_data_hi = TVectorD(len(a_data_hi),a_data_hi)
        a_data_lo = array('d',[data_lo[0],data_lo[1],data_lo[2],data_lo[3]])
        v_data_lo = TVectorD(len(a_data_lo),a_data_lo)
        
        a_systematics_hi = array('d',[systematics_hi[0],systematics_hi[1],systematics_hi[2],systematics_hi[3]])
        v_systematics_hi = TVectorD(len(a_systematics_hi),a_systematics_hi)
        a_systematics_lo = array('d',[systematics_lo[0],systematics_lo[1],systematics_lo[2],systematics_lo[3]])
        v_systematics_lo = TVectorD(len(a_systematics_lo),a_systematics_lo) 

        a_modeldep_hi = array('d',[modeldep_hi[0],modeldep_hi[1],modeldep_hi[2],modeldep_hi[3]])
        v_modeldep_hi = TVectorD(len(a_modeldep_hi),a_modeldep_hi)
        a_modeldep_lo = array('d',[modeldep_lo[0],modeldep_lo[1],modeldep_lo[2],modeldep_lo[3]])
        v_modeldep_lo = TVectorD(len(a_modeldep_lo),a_modeldep_lo) 

        v_data_hi_allunc = TVectorD(len(data_hi_allunc), array('d',[data_hi_allunc[i] for i in range(len(data_hi_allunc))]))
        v_data_lo_allunc = TVectorD(len(data_lo_allunc), array('d',[data_lo_allunc[i] for i in range(len(data_lo_allunc))]))
                        
    else:
        
        a_observable  = array('d',[0.5*(float(obs_bins[i])+float(obs_bins[i+1])) for i in range(len(obs_bins)-1)])
        v_observable  = TVectorD(len(a_observable),a_observable)
        a_dobservable = array('d',[0.5*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable = TVectorD(len(a_dobservable),a_dobservable)
        
        a_zeros = array('d',[0.0 for i in range(len(obs_bins)-1)])
        v_zeros = TVectorD(len(a_zeros),a_zeros)
        #a_twos = array('d',[0.2*a_dobservable[i] for i in range(len(obs_bins)-1)])
        a_twos = array('d',[0.015*(float(obs_bins[len(obs_bins)-1])-float(obs_bins[0])) for i in range(len(obs_bins)-1)])
        v_twos = TVectorD(len(a_twos),a_twos)
                                                                                        
        a_ggH_powheg15_JHUgen = array('d',[ggH_powheg15_JHUgen[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg15_JHUgen))])
        v_ggH_powheg15_JHUgen = TVectorD(len(a_ggH_powheg15_JHUgen),a_ggH_powheg15_JHUgen)
        # flat uncertainty NNLO
        #a_ggH_powheg15_JHUgen_hi = array('d',[ggH_powheg15_JHUgen_NNLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg15_JHUgen))])
        #v_ggH_powheg15_JHUgen_hi = TVectorD(len(a_ggH_powheg15_JHUgen_hi),a_ggH_powheg15_JHUgen_hi)
        #a_ggH_powheg15_JHUgen_lo = array('d',[ggH_powheg15_JHUgen_NNLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg15_JHUgen))])
        #v_ggH_powheg15_JHUgen_lo = TVectorD(len(a_ggH_powheg15_JHUgen_lo),a_ggH_powheg15_JHUgen_lo)
        # bin-by-bin uncertainty NLO
        #a_ggH_powheg15_JHUgen_hi =  array('d',[(ggH_powheg_unc_hi[i]*ggH_powheg15_JHUgen[i]/ggH_powheg[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_hi))])
        #a_ggH_powheg15_JHUgen_lo =  array('d',[(ggH_powheg_unc_lo[i]*ggH_powheg15_JHUgen[i]/ggH_powheg[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_lo))])
        a_ggH_powheg15_JHUgen_hi =  array('d',[ggH_powheg15_JHUgen_NLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg15_JHUgen_NLOunc_hi))])
        a_ggH_powheg15_JHUgen_lo =  array('d',[ggH_powheg15_JHUgen_NLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg15_JHUgen_NLOunc_lo))])
        v_ggH_powheg15_JHUgen_hi = TVectorD(len(a_ggH_powheg15_JHUgen_hi),a_ggH_powheg15_JHUgen_hi)
        v_ggH_powheg15_JHUgen_lo = TVectorD(len(a_ggH_powheg15_JHUgen_lo),a_ggH_powheg15_JHUgen_lo)
    
        a_ggH_minloHJJ = array('d',[ggH_minloHJJ[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJJ))])
        v_ggH_minloHJJ = TVectorD(len(a_ggH_minloHJJ),a_ggH_minloHJJ)
        #a_ggH_minloHJJ_hi = array('d',[unc_theory_ggH_hi*ggH_minloHJJ[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJJ))])
        a_ggH_minloHJJ_hi = array('d',[ggH_minloHJJ_NNLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJJ))])
        v_ggH_minloHJJ_hi = TVectorD(len(a_ggH_minloHJJ_hi),a_ggH_minloHJJ_hi)
        #a_ggH_minloHJJ_lo = array('d',[unc_theory_ggH_lo*ggH_minloHJJ[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJJ))])
        a_ggH_minloHJJ_lo = array('d',[ggH_minloHJJ_NNLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJJ))])
        v_ggH_minloHJJ_lo = TVectorD(len(a_ggH_minloHJJ_lo),a_ggH_minloHJJ_lo)
                        
        a_ggH_powheg = array('d',[ggH_powheg[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg))])
        v_ggH_powheg = TVectorD(len(a_ggH_powheg),a_ggH_powheg)
        #a_ggH_powheg_unc_hi =  array('d',[ggH_powheg_unc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_hi))])
        #a_ggH_powheg_unc_lo =  array('d',[ggH_powheg_unc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_lo))])
        a_ggH_powheg_unc_hi =  array('d',[ggH_powheg_NLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_hi))])
        a_ggH_powheg_unc_lo =  array('d',[ggH_powheg_NLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_lo))])
        v_ggH_powheg_unc_hi = TVectorD(len(a_ggH_powheg_unc_hi),a_ggH_powheg_unc_hi)
        v_ggH_powheg_unc_lo = TVectorD(len(a_ggH_powheg_unc_lo),a_ggH_powheg_unc_lo)

        a_ggH_minloHJ = array('d',[ggH_minloHJ[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ))])
        v_ggH_minloHJ = TVectorD(len(a_ggH_minloHJ),a_ggH_minloHJ)
        #a_ggH_minloHJ_unc_hi =  array('d',[ggH_minloHJ_unc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ_unc_hi))])
        #a_ggH_minloHJ_unc_lo =  array('d',[ggH_minloHJ_unc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ_unc_lo))])
        a_ggH_minloHJ_unc_hi =  array('d',[ggH_minloHJ_NLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ_unc_hi))])
        a_ggH_minloHJ_unc_lo =  array('d',[ggH_minloHJ_NLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ_unc_lo))])
        v_ggH_minloHJ_unc_hi = TVectorD(len(a_ggH_minloHJ_unc_hi),a_ggH_minloHJ_unc_hi)
        v_ggH_minloHJ_unc_lo = TVectorD(len(a_ggH_minloHJ_unc_lo),a_ggH_minloHJ_unc_lo)

        a_ggH_HRes = array('d',[ggH_HRes[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes))])
        v_ggH_HRes = TVectorD(len(a_ggH_HRes),a_ggH_HRes)
        #a_ggH_HRes_unc_hi =  array('d',[ggH_HRes_unc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes_unc_hi))])
        #a_ggH_HRes_unc_lo =  array('d',[ggH_HRes_unc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes_unc_lo))])
        a_ggH_HRes_unc_hi =  array('d',[ggH_HRes_NNLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes_unc_hi))])
        a_ggH_HRes_unc_lo =  array('d',[ggH_HRes_NNLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes_unc_lo))])
        v_ggH_HRes_unc_hi = TVectorD(len(a_ggH_HRes_unc_hi),a_ggH_HRes_unc_hi)
        v_ggH_HRes_unc_lo = TVectorD(len(a_ggH_HRes_unc_lo),a_ggH_HRes_unc_lo)


        a_XH = array('d',[XH[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH))])
        v_XH = TVectorD(len(a_XH),a_XH)
        
        #a_XH_hi = array('d',[unc_theory_XH_hi*XH[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH))])
        #a_XH_lo = array('d',[unc_theory_XH_lo*XH[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH))])
        a_XH_hi = array('d',[XH_unc[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH_unc))])
        a_XH_lo = array('d',[XH_unc[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH_unc))])
        v_XH_hi = TVectorD(len(a_XH_hi),a_XH_hi)        
        v_XH_lo = TVectorD(len(a_XH_lo),a_XH_lo)
        
        a_data = array('d',[data[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data))])
        v_data = TVectorD(len(a_data),a_data)
        a_data_hi = array('d',[data_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data_hi))])
        v_data_hi = TVectorD(len(a_data_hi),a_data_hi)
        a_data_lo = array('d',[data_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data_lo))])
        v_data_lo = TVectorD(len(a_data_lo),a_data_lo)
        
        a_systematics_hi = array('d',[(systematics_hi[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(systematics_hi))])
        v_systematics_hi = TVectorD(len(a_systematics_hi),a_systematics_hi)
        a_systematics_lo = array('d',[(systematics_lo[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(systematics_lo))])
        v_systematics_lo = TVectorD(len(a_systematics_lo),a_systematics_lo) 

        a_modeldep_hi = array('d',[(modeldep_hi[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(modeldep_hi))])
        v_modeldep_hi = TVectorD(len(a_modeldep_hi),a_modeldep_hi)
        a_modeldep_lo = array('d',[(modeldep_lo[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(modeldep_lo))])
        v_modeldep_lo = TVectorD(len(a_modeldep_lo),a_modeldep_lo)
        
        v_data_hi_allunc = TVectorD(len(data_hi_allunc), array('d',[data_hi_allunc[i] for i in range(len(data_hi_allunc))]))
        v_data_lo_allunc = TVectorD(len(data_lo_allunc), array('d',[data_lo_allunc[i] for i in range(len(data_lo_allunc))]))                                        
        
    g_ggH_powheg15_JHUgen = TGraphAsymmErrors(v_observable,v_ggH_powheg15_JHUgen,v_dobservable,v_dobservable,v_ggH_powheg15_JHUgen_lo,v_ggH_powheg15_JHUgen_hi)
    g_ggH_powheg15_JHUgen.SetFillStyle(3245);
    g_ggH_powheg15_JHUgen.SetFillColor(ROOT.kAzure)    
    g_ggH_powheg15_JHUgen.SetLineColor(ROOT.kAzure)    
    g_ggH_powheg15_JHUgen.SetLineWidth(2)    
    g_ggH_powheg15_JHUgen.SetMarkerColor(ROOT.kAzure)    

    g_ggH_powheg = TGraphAsymmErrors(v_observable,v_ggH_powheg,v_dobservable,v_dobservable,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
    g_ggH_powheg.SetFillStyle(3254);
    g_ggH_powheg.SetFillColor(ROOT.kMagenta)
    g_ggH_powheg.SetLineColor(ROOT.kMagenta)
    g_ggH_powheg.SetMarkerColor(ROOT.kMagenta)
        
    g_ggH_minloHJ = TGraphAsymmErrors(v_observable,v_ggH_minloHJ,v_dobservable,v_dobservable,v_ggH_minloHJ_unc_lo,v_ggH_minloHJ_unc_hi)
    g_ggH_minloHJ.SetFillStyle(3254);
    g_ggH_minloHJ.SetFillColor(ROOT.kOrange+2)
    g_ggH_minloHJ.SetLineColor(ROOT.kOrange+2)
    g_ggH_minloHJ.SetLineWidth(2)
    g_ggH_minloHJ.SetMarkerColor(ROOT.kOrange+2)
    
    g_ggH_minloHJJ = TGraphAsymmErrors(v_observable,v_ggH_minloHJJ,v_dobservable,v_dobservable,v_ggH_minloHJJ_lo,v_ggH_minloHJJ_hi)
    g_ggH_minloHJJ.SetFillStyle(3245);
    g_ggH_minloHJJ.SetFillColor(ROOT.kOrange+4)
    g_ggH_minloHJJ.SetLineColor(ROOT.kOrange+4)
    g_ggH_minloHJJ.SetMarkerColor(ROOT.kOrange+4)

    g_ggH_HRes = TGraphAsymmErrors(v_observable,v_ggH_HRes,v_dobservable,v_dobservable,v_ggH_HRes_unc_lo,v_ggH_HRes_unc_hi)
    g_ggH_HRes.SetFillStyle(3254);
    g_ggH_HRes.SetFillColor(ROOT.kMagenta)
    g_ggH_HRes.SetLineColor(ROOT.kMagenta)
    g_ggH_HRes.SetLineWidth(2)
    g_ggH_HRes.SetMarkerColor(ROOT.kMagenta)

    g_XH = TGraphAsymmErrors(v_observable,v_XH,v_dobservable,v_dobservable,v_XH_lo,v_XH_hi)
    g_XH.SetFillColor(ROOT.kGreen+3)
    g_XH.SetLineColor(ROOT.kGreen+3)

    g_data = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_data_lo,v_data_hi)
    g_data.SetMarkerColor(ROOT.kBlack)
    g_data.SetLineColor(ROOT.kBlack)
    g_data.SetLineWidth(2)
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(1.4)

    g_data_e0 = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_zeros,v_zeros)
    g_data_e0.SetMarkerColor(ROOT.kBlack)
    g_data_e0.SetLineColor(ROOT.kBlack)
    g_data_e0.SetLineWidth(2)
    g_data_e0.SetMarkerStyle(20)
    g_data_e0.SetMarkerSize(1.4)

    if (useHRes):
        v_ratio_data = TVectorD(len(data), array('d',[data[i]/ggH_HRes[i] for i in range(len(data))]))
        v_ratio_data_hi = TVectorD(len(data), array('d',[data_hi[i]/ggH_HRes[i] for i in range(len(data))]))
        v_ratio_data_lo = TVectorD(len(data), array('d',[data_lo[i]/ggH_HRes[i] for i in range(len(data))]))

        v_ratio_HRes = TVectorD(len(ggH_HRes), array('d',[ggH_HRes[i]/ggH_HRes[i] for i in range(len(ggH_HRes))]))
        v_ratio_HRes_hi = TVectorD(len(ggH_HRes), array('d',[ggH_HRes_NNLOunc_hi[i]/ggH_HRes[i] for i in range(len(ggH_HRes))]))
        v_ratio_HRes_lo = TVectorD(len(ggH_HRes), array('d',[ggH_HRes_NNLOunc_lo[i]/ggH_HRes[i] for i in range(len(ggH_HRes))]))

        v_ratio_powheg15_JHUgen = TVectorD(len(ggH_powheg15_JHUgen), array('d',[ggH_powheg15_JHUgen[i]/ggH_HRes[i] for i in range(len(ggH_powheg15_JHUgen))]))
        v_ratio_powheg15_JHUgen_hi = TVectorD(len(ggH_powheg15_JHUgen), array('d',[ggH_powheg15_JHUgen_NLOunc_hi[i]/ggH_HRes[i] for i in range(len(ggH_powheg15_JHUgen))]))
        v_ratio_powheg15_JHUgen_lo = TVectorD(len(ggH_powheg15_JHUgen), array('d',[ggH_powheg15_JHUgen_NLOunc_lo[i]/ggH_HRes[i] for i in range(len(ggH_powheg15_JHUgen))]))


        v_ratio_minloHJ = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ[i]/ggH_HRes[i] for i in range(len(ggH_minloHJ))]))
        v_ratio_minloHJ_hi = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ_NLOunc_hi[i]/ggH_HRes[i] for i in range(len(ggH_minloHJ))]))
        v_ratio_minloHJ_lo = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ_NLOunc_lo[i]/ggH_HRes[i] for i in range(len(ggH_minloHJ))]))

    else:

        v_ratio_data = TVectorD(len(data), array('d',[data[i]/ggH_minloHJ[i] for i in range(len(data))]))
        v_ratio_data_hi = TVectorD(len(data), array('d',[data_hi[i]/ggH_minloHJ[i] for i in range(len(data))]))
        v_ratio_data_lo = TVectorD(len(data), array('d',[data_lo[i]/ggH_minloHJ[i] for i in range(len(data))]))

        v_ratio_HRes = TVectorD(len(ggH_HRes), array('d',[ggH_HRes[i]/ggH_minloHJ[i] for i in range(len(ggH_HRes))]))
        v_ratio_HRes_hi = TVectorD(len(ggH_HRes), array('d',[ggH_HRes_NNLOunc_hi[i]/ggH_minloHJ[i] for i in range(len(ggH_HRes))]))
        v_ratio_HRes_lo = TVectorD(len(ggH_HRes), array('d',[ggH_HRes_NNLOunc_lo[i]/ggH_minloHJ[i] for i in range(len(ggH_HRes))]))

        v_ratio_powheg15_JHUgen = TVectorD(len(ggH_powheg15_JHUgen), array('d',[ggH_powheg15_JHUgen[i]/ggH_minloHJ[i] for i in range(len(ggH_powheg15_JHUgen))]))
        v_ratio_powheg15_JHUgen_hi = TVectorD(len(ggH_powheg15_JHUgen), array('d',[ggH_powheg15_JHUgen_NLOunc_hi[i]/ggH_minloHJ[i] for i in range(len(ggH_powheg15_JHUgen))]))
        v_ratio_powheg15_JHUgen_lo = TVectorD(len(ggH_powheg15_JHUgen), array('d',[ggH_powheg15_JHUgen_NLOunc_lo[i]/ggH_minloHJ[i] for i in range(len(ggH_powheg15_JHUgen))]))


        v_ratio_minloHJ = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ[i]/ggH_minloHJ[i] for i in range(len(ggH_minloHJ))]))
        v_ratio_minloHJ_hi = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ_NLOunc_hi[i]/ggH_minloHJ[i] for i in range(len(ggH_minloHJ))]))
        v_ratio_minloHJ_lo = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ_NLOunc_lo[i]/ggH_minloHJ[i] for i in range(len(ggH_minloHJ))]))
        
    g_ratio_data = TGraphAsymmErrors(v_observable,v_ratio_data,v_zeros,v_zeros,v_ratio_data_lo,v_ratio_data_hi)
    g_ratio_data.SetMarkerColor(ROOT.kBlack)
    g_ratio_data.SetLineColor(ROOT.kBlack)
    g_ratio_data.SetLineWidth(2)
    g_ratio_data.SetMarkerStyle(20)
    g_ratio_data.SetMarkerSize(1.4)

    g_ratio_HRes = TGraphAsymmErrors(v_observable,v_ratio_HRes,v_dobservable,v_dobservable,v_ratio_HRes_lo,v_ratio_HRes_hi)
    g_ratio_HRes.SetFillStyle(3254);
    g_ratio_HRes.SetFillColor(ROOT.kMagenta)
    g_ratio_HRes.SetLineColor(ROOT.kMagenta)
    g_ratio_HRes.SetLineWidth(2)
    g_ratio_HRes.SetMarkerColor(ROOT.kMagenta)

    g_ratio_powheg15_JHUgen = TGraphAsymmErrors(v_observable,v_ratio_powheg15_JHUgen,v_dobservable,v_dobservable,v_ratio_powheg15_JHUgen_lo,v_ratio_powheg15_JHUgen_hi)
    g_ratio_powheg15_JHUgen.SetFillStyle(3245);
    g_ratio_powheg15_JHUgen.SetFillColor(ROOT.kAzure)    
    g_ratio_powheg15_JHUgen.SetLineColor(ROOT.kAzure)    
    g_ratio_powheg15_JHUgen.SetLineWidth(2)    
    g_ratio_powheg15_JHUgen.SetMarkerColor(ROOT.kAzure) 

    g_ratio_minloHJ = TGraphAsymmErrors(v_observable,v_ratio_minloHJ,v_dobservable,v_dobservable,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
    g_ratio_minloHJ.SetFillStyle(3254);
    g_ratio_minloHJ.SetFillColor(ROOT.kOrange+2)
    g_ratio_minloHJ.SetLineColor(ROOT.kOrange+2)
    g_ratio_minloHJ.SetLineWidth(2)
    g_ratio_minloHJ.SetMarkerColor(ROOT.kOrange+2)


    g_modeldep = TGraphAsymmErrors(v_observable,v_data,v_twos,v_twos,v_modeldep_lo,v_modeldep_hi)
    #g_modeldep = TGraphAsymmErrors(v_observable,v_data,v_twos,v_twos,v_modeldep_lo,v_modeldep_hi)
    g_systematics = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_systematics_lo,v_systematics_hi)

    #g_modeldep.SetLineWidth(5)
    #g_modeldep.SetMarkerColor(ROOT.kGray)
    g_modeldep.SetFillColor(ROOT.kGray)
    g_modeldep.SetLineColor(ROOT.kGray)

    g_systematics.SetLineWidth(5)
    g_systematics.SetMarkerColor(ROOT.kRed)
    g_systematics.SetLineColor(ROOT.kRed)
    g_systematics.SetFillColor(ROOT.kRed)

    g_data_allunc = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_data_lo_allunc,v_data_hi_allunc)
    g_data_allunc.SetMarkerColor(ROOT.kBlack)
    g_data_allunc.SetLineColor(ROOT.kBlack)
    g_data_allunc.SetLineWidth(1)
    g_data_allunc.SetMarkerStyle(20)
    g_data_allunc.SetMarkerSize(1.2)
    
    if (obsName=="pT4l"):
        label="p_{T}^{H}"
        unit="GeV"
    elif (obsName=="massZ2"):
        label = "m(Z_{2})"
        unit = "GeV"
    elif (obsName=="massZ1"):
        label = "m(Z_{1})"
        unit = "GeV"
    elif (obsName=="nJets" or obsName=="njets_reco_pt30_eta4p7"):
        label = "N(jets) |#eta|<4.7"
        unit = ""
    elif (obsName=="njets_reco_pt30_eta2p5"):
        label = "N(jets) |#eta|<2.5"
        unit = ""
    elif (obsName=="pt_leadingjet_reco_pt30_eta4p7"):
        label = "p_{T}(jet)"
        unit = "GeV"
    elif (obsName=="pt_leadingjet_reco_pt30_eta2p5"):
        label = "p_{T}(jet) |#eta|<2.5"
        unit = "GeV"
    elif (obsName=="absrapidity_leadingjet_reco_pt30_eta4p7"):
        label = "|y(jet)|"
        unit = ""
    elif (obsName=="absrapidity_leadingjet_reco_pt30_eta2p5"):
        label = "|y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=="absdeltarapidity_hleadingjet_reco_pt30_eta4p7"):
        label = "|y(H)-y(jet)|"
        unit = ""
    elif (obsName=="absdeltarapidity_hleadingjet_reco_pt30_eta2p5"):
        label = "|y(H)-y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=="rapidity4l"):
        label = "|y^{H}|"
        unit = ""
    elif (obsName=="cosThetaStar"):
        label = "|cos#theta*|"
        unit = ""
    elif (obsName=="cosTheta1"):
        label = "|cos#theta_{1}|"
        unit = ""
    elif (obsName=="cosTheta2"):
        label = "|cos#theta_{2}|"
        unit = ""
    elif (obsName=="Phi"):
        label = "|#Phi|"
        unit = ""
    elif (obsName=="Phi1"):
        label = "|#Phi_{1}|"
        unit = ""
    elif (obsName=="mass4l"):
        label = "inclusive"
        unit = ""
    else:
        label = obsName
        unit = ""
    
    c = TCanvas("c",obsName, 1400, 1400)
    if(opt.SETLOG): c.SetLogy()
    c.SetBottomMargin(0.35)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.18)

    if (obsName=="mass4l"):
        dummy = TH1D("dummy","dummy", 4, 0, 4)
        for i in range(1,5):
            dummy.SetBinContent(i,2.5*max(a_ggH_powheg15_JHUgen))
    else:
        if ("jet" in obsName and (not obsName.startswith("njets"))):
            dummy = TH1D("dummy","dummy", int(float(obs_bins[nBins-1])-float(obs_bins[1])), float(obs_bins[1]), float(obs_bins[nBins-1]))
            for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[1]))):
                dummy.SetBinContent(i,2.5*max(a_ggH_powheg15_JHUgen))
        else:
            dummy = TH1D("dummy","dummy", int(float(obs_bins[nBins-1])-float(obs_bins[0])), float(obs_bins[0]), float(obs_bins[nBins-1]))
            for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[0]))):
                dummy.SetBinContent(i,2.5*max(a_ggH_powheg15_JHUgen))
    if(opt.SETLOG):
        dummy.SetMaximum(55.0*max(max(a_data),(max(a_ggH_powheg15_JHUgen))))
    else:
        if (obsName=="mass4l"): dummy.SetMaximum(1.6*(max(max(a_data),(max(a_ggH_powheg15_JHUgen)))+max(a_data_hi)))
        else: dummy.SetMaximum(1.5*(max(max(a_ggH_powheg15_JHUgen),(max(a_data)+max(a_data_hi)))))
    if (opt.SETLOG): dummy.SetMinimum(0.0501*max(min(a_data),(min(a_ggH_powheg15_JHUgen))))
    else: dummy.SetMinimum(0.0001)
    dummy.SetLineColor(0)
    dummy.SetMarkerColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.GetXaxis().SetLabelSize(0.0)
    dummy.GetYaxis().SetLabelSize(0.04)
    if (opt.SETLOG and (obsName.startswith('njets') or obsName.startswith('pt_leading'))):
        dummy.SetMaximum(200.0*max(max(a_data),(max(a_ggH_powheg15_JHUgen))))
        dummy.GetXaxis().SetTitle("")
    else:
        dummy.GetXaxis().SetTitle("")
    if (obsName.startswith('njets')):
        dummy.GetXaxis().SetTitle('')
        dummy.GetXaxis().SetLabelSize(0.0)
        dummy.GetXaxis().SetBinLabel(1,'')
        dummy.GetXaxis().SetBinLabel(2,'')
        dummy.GetXaxis().SetBinLabel(3,'')
        dummy.GetXaxis().SetBinLabel(4,'')
    dummy.GetXaxis().SetTitleSize(0.0)

    if (label=="inclusive"):
        dummy.GetYaxis().SetTitle("#sigma_{fid.} [fb]")
    elif (unit==""):
        dummy.GetYaxis().SetTitle("d#sigma_{fid.}/d"+label+" [fb]")
    else:
        dummy.GetYaxis().SetTitle("d#sigma_{fid.}/d"+label+" [fb/"+unit+"]")
    if (obsName.startswith('njets')):
        dummy.GetYaxis().SetTitle("d#sigma_{fid.} [fb]")
        
    #dummy.GetYaxis().SetTitleOffset(1.5)
    dummy.GetYaxis().SetTitleOffset(1.4)
    dummy.Draw("hist")
    #g_ggH_minloHJJ.Draw("2same")
    #g_ggH_minloHJJ.Draw("psame")
    #g_ggH_powheg.Draw("2same")
    #g_ggH_powheg.Draw("psame")
    g_ggH_minloHJ.Draw("2same")
    g_ggH_minloHJ.Draw("psame")
    g_ggH_powheg15_JHUgen.Draw("2same")
    g_ggH_powheg15_JHUgen.Draw("psame")
    if (useHRes):
        g_ggH_HRes.Draw("2same")
        g_ggH_HRes.Draw("psame")

    g_XH.Draw("2same")
    g_modeldep.Draw("2same")
    #g_modeldep.Draw("psame")
    g_data.Draw("psameZ")
    g_systematics.Draw("psameZ")
    g_data_e0.Draw("psame")
    #g_data_allunc.Draw("psame")
    #g_data_allunc.Draw("psame[]")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.9, 0.94,"19.7 fb^{-1} (8 TeV)")
    latex2.SetTextSize(0.8*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.19, 0.94, "CMS")
    latex2.SetTextSize(0.7*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    #latex2.DrawLatex(0.28, 0.945, "Unpublished")
    #latex2.DrawLatex(0.28, 0.945, "Preliminary")
    
    legend = TLegend(.2,.72,.93,.90)
    legend . SetNColumns(2)
    #legend . SetTextSize(0.026)
    #legend . SetTextSize(0.029) # less info for XH
    legend . SetTextSize(0.025) # shrink margins
    if (opt.UNBLIND): legend . AddEntry(g_data , "Data (stat.#oplussys. unc.)", "ep")
    else: legend . AddEntry(g_data , "Asimov Data (stat.#oplussys. unc.)", "ep")
    legend . AddEntry(g_ggH_powheg15_JHUgen , "gg#rightarrowH (POWHEG+JHUGen) + XH", "lf") 
    #legend . AddEntry(g_ggH_powheg , "gg#rightarrowH (powheg) + XH", "lf") 
    legend . AddEntry(g_systematics,"Systematic uncertainty","l") 
    #legend . AddEntry(g_ggH_minloHJJ , "gg#rightarrowH (minlo HJJ) + XH", "lf")
    legend . AddEntry(g_ggH_minloHJ , "gg#rightarrowH (MiNLO HJ) + XH", "lf")
    legend . AddEntry(g_modeldep,"Model dependence","f") 
    if (useHRes):
        legend . AddEntry(g_ggH_HRes , "gg#rightarrowH (HRes) + XH", "lf")        
        legend . AddEntry(g_XH , "", "")
        legend . AddEntry(g_XH , "XH = VBF + VH + ttH", "l")
    else:
        legend . AddEntry(g_XH , "XH = VBF + VH + ttH", "l")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw()

    dummy.Draw("axissame")
    
    pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0)
    pad.SetTopMargin(0.65)
    pad.SetRightMargin(0.04)
    pad.SetLeftMargin(0.18)
    pad.SetFillColor(0)
    pad.SetGridy(1)
    pad.SetFillStyle(0)
    pad.Draw()
    pad.cd(0)

    if (obsName=="mass4l"):
        dummy2 = TH1D("dummy2","dummy2", 4, 0, 4)
        for i in range(1,5):
            dummy2.SetBinContent(i,1.02)
        dummy2.GetXaxis().SetLabelSize(0.08)
        dummy2.GetXaxis().SetBinLabel(1,'4l')
        dummy2.GetXaxis().SetBinLabel(2,'2e2#mu')
        dummy2.GetXaxis().SetBinLabel(3,'4#mu')
        dummy2.GetXaxis().SetBinLabel(4,'4e')        
    else:
        if ("jet" in obsName and (not obsName.startswith("njets"))):
            dummy2 = TH1D("dummy2","dummy2", int(float(obs_bins[nBins-1])-float(obs_bins[1])), float(obs_bins[1]), float(obs_bins[nBins-1]))
            for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[1]))):
                dummy2.SetBinContent(i,1.02)
        else:
            dummy2 = TH1D("dummy2","dummy2", int(float(obs_bins[nBins-1])-float(obs_bins[0])), float(obs_bins[0]), float(obs_bins[nBins-1]))
            for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[0]))):
                dummy2.SetBinContent(i,1.02)
        dummy2.GetXaxis().SetLabelSize(0.04)

    dummy2.SetLineColor(0)
    dummy2.SetMarkerColor(0)
    dummy2.SetLineWidth(0)
    dummy2.SetMarkerSize(0)
    if (obsName.startswith('njets')):
        dummy2.GetXaxis().SetTitle(label)
        dummy2.GetXaxis().SetLabelSize(0.08)
        dummy2.GetXaxis().SetBinLabel(1,'0')
        dummy2.GetXaxis().SetBinLabel(2,'1')
        dummy2.GetXaxis().SetBinLabel(3,'2')
        dummy2.GetXaxis().SetBinLabel(4,'#geq 3')
    elif (label=="inclusive"):
        dummy2.GetXaxis().SetTitle("")
    elif (unit==""):
        dummy2.GetXaxis().SetTitle(label)
    else:
        dummy2.GetXaxis().SetTitle(label+" ["+unit+"]")

    dummy2.GetYaxis().SetLabelSize(0.03)
    dummy2.GetYaxis().SetNdivisions(10);
    dummy2.GetXaxis().SetNdivisions(510)
    dummy2.SetMaximum(1.15*(max(v_ratio_data)+max(v_ratio_data_hi)))
    dummy2.SetMinimum(0.0)
    dummy2.Draw("hist")
    dummy2.GetYaxis().CenterTitle()
    dummy2.GetYaxis().SetTitleSize(0.04)
    dummy2.GetYaxis().SetTitleOffset(1.5)
    if (useHRes):
        dummy2.GetYaxis().SetTitle('Ratio to HRes')
    else:
        dummy2.GetYaxis().SetTitleSize(0.03)
        dummy2.GetYaxis().SetTitleOffset(2.0)        
        dummy2.GetYaxis().SetTitle('Ratio to MiNLO HJ')

    if (useHRes):
        g_ratio_HRes.Draw("2same")
        g_ratio_HRes.Draw("psame")
    g_ratio_minloHJ.Draw("2same")
    g_ratio_minloHJ.Draw("psame")
    g_ratio_powheg15_JHUgen.Draw("2same")
    g_ratio_powheg15_JHUgen.Draw("psame")
    g_ratio_data.Draw("psameZ")

    #if (opt.SETLOG): set_log = '_HResTESTlogscale'
    #else: set_log = 'HResTEST'
    
    if (opt.SETLOG): set_log = '_logscale'
    else: set_log = ''

    if (not opt.UNBLIND): set_log = set_log + '_asimov'
 
    c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'.pdf')
    c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'.png')
    c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'.root')
    c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'.C')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.pdf')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.png')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.root')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.C')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.pdf')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.png')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.root')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.C')
   
obs_bins = opt.OBSBINS.split("|")
if (not (obs_bins[0] == '' and obs_bins[len(obs_bins)-1]=='')):
    print 'BINS OPTION MUST START AND END WITH A |'
obs_bins.pop()
obs_bins.pop(0)
#if float(obs_bins[len(obs_bins)-1])>200.0:
#    obs_bins[len(obs_bins)-1]='200.0'
if (opt.OBSNAME=="nJets" or opt.OBSNAME.startswith("njets")):
    obs_bins[len(obs_bins)-1]='4'

plotXS(opt.OBSNAME, obs_bins)  

