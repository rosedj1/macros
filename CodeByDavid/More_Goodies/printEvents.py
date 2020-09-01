from ROOT import *

#infile = TFile('Sync_94X.root','READ')
infile = TFile('Sync_94X_test.root','READ')
#infile = TFile('/raid/raid7/dsperka/Run2/Zprime/SingleDoubleMuon_Run2017-17Nov2017-v1_NoDuplicates.root','READ')
ext=''
#ext='_shortest'
ext='_70to110'
#ext='_150to300'
#ext='_110to150'
#ext='_110to150andabove500'
#ext='_norefit'
#ext='_masses'
#ext='_test'

t = infile.Get('Ana/passedEvents')
#t = infile.Get('passedEvents')
nevents = t.GetEntriesFast()

n4e=0
n4mu=0
n2e2mu=0

with open('eventsBUF_Sync94X'+ext+'.txt','w') as f:
#with open('eventsBUF_Sync94X'+ext+'_data.txt','w') as f:

  for i in xrange(nevents):
    t.GetEntry(i)

    #if (not (t.passedZ4lSelection==1 and t.massZ2>12.0)): continue
    if (not (t.passedFullSelection==1)): continue

    if (ext=='_70to110' and not (t.mass4l>70.0 and t.mass4l<110.0)): continue
    if (ext=='_150to300' and not (t.mass4l>150.0 and t.mass4l<300.0)): continue
    if (ext=='_110to150' and not (t.mass4l>110.0 and t.mass4l<150.0)): continue
    if (ext=='_500' and not (t.mass4l>500.0)): continue
    if (ext=='_110to150andabove500' and not (t.mass4l>500.0 or (t.mass4l>110.0 and t.mass4l<150.0))): continue

    if (t.njets_pt30_eta4p7>0): 
      #ptjet1=t.jet_pt[t.jet_iscleanH4l[0]]
      #qg1=t.jet_QGTagger[t.jet_iscleanH4l[0]]
      ptjet1=t.jet_pt[t.jet1index]
      qg1=t.jet_QGTagger[t.jet1index]
    else: 
      ptjet1=-1.0
      qg1=-1.0

    if (t.njets_pt30_eta4p7>1): 
      ptjet2=t.jet_pt[t.jet2index]
      qg2=t.jet_QGTagger[t.jet2index]
    else: 
      ptjet2=-1.0
      qg2=-1.0

    if (t.finalState==1): 
      n4mu+=1  
      ch="4mu"
    if (t.finalState==2): 
      n4e+=1  
      ch="4e"
    if (t.finalState>=3): 
      n2e2mu+=1
      ch="2e2mu"
    
    #if (not(ch=="4mu")): continue

    if (ext=='_shortest'):
      f.write('{}:{}:{}'.format(t.Run,t.LumiSect,t.Event)+'\n')
    elif (ext=='_short'):
      f.write('{}:{}:{}:{:.2f}:{:.2f}:{:.2f}:{:.2f}:{:.2f}:{:.2f}'.format(t.Run,t.LumiSect,t.Event,t.mass4l,t.mass4lErr,t.massZ1,t.massZ2,t.mass4lREFIT,t.mass4lErrREFIT)+'\n')
    elif (ext=='_norefit'):
      f.write('{}:{}:{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:d}:{:.2f}:{:.2f}:{}'.format(t.Run,t.LumiSect,t.Event,t.mass4l,t.massZ1,t.massZ2,t.D_bkg_kin,t.D_bkg,t.Dgg10_VAMCFM,t.Djet_VAJHU,t.D_g4,t.njets_pt30_eta4p7,ptjet1,ptjet2,t.EventCat)+'\n')
    else:
      f.write('{}:{}:{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:d}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{}:{:.2f}:{:.2f}:{:.3f}'.format(t.Run,t.LumiSect,t.Event,t.mass4l,t.massZ1,t.massZ2,t.D_bkg_kin,t.D_bkg,t.D_VBF,t.D_g4,t.D_VBF1j,t.D_HadWH,t.D_HadZH,t.njets_pt30_eta4p7,ptjet1,ptjet2,qg1,qg2,t.D_VBF_QG,t.D_VBF1j_QG,t.D_HadWH_QG,t.D_HadZH_QG,t.met,t.EventCat,t.mass4lREFIT,t.mass4lErrREFIT,t.genWeight*t.pileupWeight*t.dataMCWeight)+'\n')


  f.close()

print nevents,'events in the file'
print '4e:',n4e,'4mu:',n4mu,'2e2mu:',n2e2mu
