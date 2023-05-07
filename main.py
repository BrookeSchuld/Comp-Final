import ROOT as RT
from ROOT import TLorentzVector
import numpy as np
import sys
import math
from scipy import optimize
import plotting_functions
import matplotlib.pyplot as plt

#Function that gets called if the particle meets the specifications to be identified. 
def id(evt,vtx,part_id,trk,theory_de_dx,theory_residual):
    de_dx,res_range=plotting_functions.residual(evt,part_id)  #gets the residual range and de/dx of the function
    p=np.linspace(10,600,10000)
    chi_squared=[]
    for i in range(len(theory_de_dx)):
        chi=[]
        for j in range(len(res_range)):
            n=(np.argmin(np.abs(res_range[j]-theory_residual[i])))     #Find the expected residual given the de/dx and compares it to the given residual 
            
            chi.append(np.abs((theory_de_dx[i][n]-de_dx[j])**2/theory_de_dx[i][n]))
        chi_squared.append(np.sum((chi))/len(de_dx)
    return chi_squared                                                 # return  reduced chi squared for each of the three masses


#extract root trees
edep_tree=RT.TChain("EDepSimEvents")
grtk_tree=RT.TChain("DetSimPassThru/gRooTracker")
filelist=[sys.argv[x] for x in range(1,len(sys.argv))]
for file in filelist:
    edep_tree.Add(file)
    grtk_tree.Add(file)
nevt=edep_tree.GetEntries()
tagged=[]                                             #array to keep track of the tagging  
masses=np.array([938.27,135.0,497.61])                #proton,pion,kaon (MeV)
p=np.linspace(10,600,5000)
theory_de_dx=[]
theory_residual=[]
for mass in masses:                                    #calculate the de/dx and residual range for each particle type
    de_dx,res=plotting_functions.res_range(p,mass)
    theory_de_dx.append(np.array(de_dx))
    theory_residual.append(np.array(res))
f=open("tagged.txt","w")
for evt in range(nevt):                              #loop over each event
    print(evt) 					
    edep_tree.GetEntry(evt)
    grtk_tree.GetEntry(evt)
    vtx=edep_tree.Event.Primaries[0]
    traj=edep_tree.Event.Trajectories
    particles=[2212,211,310]                         #PDG Codes for each particle type
    for trk in traj:                                 #loop over particles in the event
        if abs(trk.GetPDGCode())==11:continue        #ignore electrons
        part_id=trk.GetTrackId()
        part_trk=trk
        #Three conditions to see if the particle is a good candidate
        if plotting_functions.is_point_contained(traj[part_id].Points[-1].GetPosition().Vect())==False or plotting_functions.is_point_contained(traj[part_id].Points[0].GetPosition().Vect())==False: continue
        re_co,de_dx=plotting_functions.reco(evt,part_id)       #Is the particle containted?        
        if len(de_dx)<5: continue                              #Is the track long enough to get hits
        if re_co/(trk.GetInitialMomentum().E()-trk.GetInitialMomentum().M())<.9:continue      #is the reconstructed energy enough to rule out another interaction 
       
        likelihood=id(evt,vtx,part_id,trk,theory_de_dx,theory_residual)
        if min(likelihood)>2: continue                         #Is there a good match 
        max_like=np.argmin(likelihood)                         #which match is best
        particle=particles[max_like]
        with open("tagged.txt","a") as f:                      #write the particle id and it's tag
            f.write(str(part_id))
            f.write(",")
            f.write(str(particle))
            f.write("\n")
     
positive_id=np.sum(tagged)/len(tagged)
print(positive_id,len(tagged))
f.close()        
        
