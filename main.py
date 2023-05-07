import ROOT as RT
from ROOT import TLorentzVector
import numpy as np
import sys
#import lar_functions as lar
import math
from scipy import optimize
import plotting_functions
import matplotlib.pyplot as plt
def f(p,mass,de_dx):
    return plotting_functions.bethe_bloch(p,mass)-de_dx
def g(p,mass,res):
    x=np.linspace(0,p,1000)
    de_dx,residual=plotting_functions.res_range(x,mass)
    return np.array(residual)-res
def characteristics(evt,vtx,part_id,trk,theory_de_dx,theory_residual):
    masses=np.array([938.27,135.0,497.61]) #proton,pion,kaon
    de_dx,res_range=plotting_functions.residual(evt,part_id)
    p=np.linspace(10,600,10000)
    percent_error=[]
    for i in range(len(theory_de_dx)):
        expected_de_dx=[]
        for j in range(len(res_range)):
            n=(np.argmin(np.abs(res_range[j]-theory_residual[i])))
            
            expected_de_dx.append(np.abs((theory_de_dx[i][n]-de_dx[j])/de_dx[j]))
        percent_error.append(np.sum(np.abs(expected_de_dx)))
    return percent_error

edep_tree=RT.TChain("EDepSimEvents")
grtk_tree=RT.TChain("DetSimPassThru/gRooTracker")
filelist=[sys.argv[x] for x in range(1,len(sys.argv))]
for file in filelist:
    edep_tree.Add(file)
    grtk_tree.Add(file)
nevt=edep_tree.GetEntries()
tagged=[]
masses=np.array([938.27,135.0,497.61]) #proton,pion,kaon
p=np.linspace(10,600,5000)
theory_de_dx=[]
theory_residual=[]
for mass in masses: 
    de_dx,res=plotting_functions.res_range(p,mass)
    theory_de_dx.append(np.array(de_dx))
    theory_residual.append(np.array(res))
print(theory_residual[0])
f=open("tagged.txt","w")
for evt in range(3275,nevt):
    print(evt) 					
    edep_tree.GetEntry(evt)
    grtk_tree.GetEntry(evt)
    vtx=edep_tree.Event.Primaries[0]
    primary_pdg=[x.GetPDGCode() for x in vtx.Particles]
    traj=edep_tree.Event.Trajectories
    particles=[2212,211,310]
    for trk in traj:
        if abs(trk.GetPDGCode())==11:continue
        part_id=trk.GetTrackId()
        part_trk=trk
        if plotting_functions.is_point_contained(traj[part_id].Points[-1].GetPosition().Vect())==False or plotting_functions.is_point_contained(traj[part_id].Points[0].GetPosition().Vect())==False: continue
        re_co,de_dx=plotting_functions.reco(evt,part_id)
        if len(de_dx)<5: continue
        #print(re_co/(trk.GetInitialMomentum().E()-trk.GetInitialMomentum().M()),trk.GetPDGCode())
        if re_co/(trk.GetInitialMomentum().E()-trk.GetInitialMomentum().M())<.9:continue
       
        likelyhood=characteristics(evt,vtx,part_id,trk,theory_de_dx,theory_residual)
        if min(likelyhood)>3: continue
        #print(trk.GetPDGCode(), likelyhood)
        max_like=np.argmin(likelyhood)
        particle=particles[max_like]
        with open("tagged.txt","a") as f:
            f.write(str(abs(trk.GetPDGCode())))
            f.write(",")
            f.write(str(particle))
            f.write("\n")
        #f=open("tagging.txt","a")
        #f.write(particle, abs(trk.GetPDGCode())) 
        #f.close()
positive_id=np.sum(tagged)/len(tagged)
print(positive_id,len(tagged))
#print(particle_type.count(2212))
#print(particle_type.count(211))
f.close()        
        
