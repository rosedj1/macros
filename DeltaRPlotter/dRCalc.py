import numpy as np
                                       
## Make sure that: -Pi <= deltaPhi <= Pi 
def dPhiCalc(phi1, phi2):              
    dPhi = phi1 - phi2            
    while dPhi > np.pi:     dPhi -= 2*np.pi
    while dPhi < -np.pi:    dPhi += 2*np.pi
    return dPhi

def dRCalc(eta1, phi1, eta2, phi2):                
    dEta = eta1 - eta2
    dPhi = dPhiCalc(phi1, phi2)
    dR = np.sqrt(dEta*dEta + dPhi*dPhi)    
    return dR                          
