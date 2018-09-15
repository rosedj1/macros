import numpy as np
                                       
## Make sure that: 0 <= deltaPhi <= pi 
def dPhiCalc(phi1, phi2):              
    dPhi = abs(phi1-phi2)            
    while result > np.pi:              
        dPhi = 2*np.pi - dPhi      
    return dPhi

def dRCalc(eta1, phi1, eta2, phi2):                
    dEta = eta1 - eta2
    dPhi = dPhiCalc(phi1, phi2)
    dR = np.sqrt(dEta*dEta + dPhi*dPhi)    
    return dR                          
