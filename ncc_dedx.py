import ROOT as RT
from ROOT import TLorentzVector
import numpy as np
import sys
import lar_functions as lar
import math
import plotting_functions

kkBlue=RT.TColor(9000,   0/255., 119/255., 187/255.)                ##Add colors from ROOT
kkOrange  = RT.TColor(9003, 238/255., 119/255.,  51/255.)

edep_tree=RT.TChain("EDepSimEvents")
grtk_tree=RT.TChain("DetSimPassThru/gRooTracker")
filelist=[sys.argv[x] for x in range(1,len(sys.argv))]
for file in filelist:
    edep_tree.Add(file)
    grtk_tree.Add(file)

def is_point_contained(pos):                 ##Check to see if the position is within the 2x2 
    if abs(pos[0])>670: return False
    if abs(pos[1]-430)>670:return False
    if abs(pos[2])>670:return False
    return True
    
def statistics(evt,vtx,part_id,trk,mass):                          #Gets kinematics for a specific particle
    energy=[]
    vtx_pos=vtx.GetPosition().Vect()
    position=traj[part_id].Points[-1].GetPosition().Vect()
    distance=(vtx_pos-position).Mag()/10
    for j in traj[part_id].Points:
        energy.append(np.sqrt(j.GetMomentum().Mag()**2+mass**2))
    DE=np.sum(-np.diff(energy))
    KE=(trk.GetInitialMomentum().E()-trk.GetInitialMomentum().M())
    KE_final=(energy[-2]/mass-1)*mass
    return distance,DE,KE_final


def reco(part_id):                                                 #Gets the energy deposited in the detector for a given particle
    reco_energy=0
    de_dx=[]
    for seg in edep_tree.Event.SegmentDetectors:
                nChunks=len(seg[1])
                for n in range(nChunks):
                    key_contrib=seg[1][n].GetContributors()[0]
                    #print(is_2x2_contained(seg[1][n].GetStop()))
                    if part_id==key_contrib:
                        reco_energy+=seg[1][n].GetEnergyDeposit()
                        seg_length=(seg[1][n].GetStart()-seg[1][n].GetStop()).Mag()/10
                        de_dx.append(-seg[1][n].GetEnergyDeposit()/seg_length)
                    #for i in proton_children:
                    #    if i==key_contrib:
                    #if proton_id==key_contrib:
                     #       reco_energy+=seg[1][n].GetEnergyDeposit()
    return reco_energy,de_dx
def residual(part_id):                                                                 #Gets the residual range for a given particle
    for seg in edep_tree.Event.SegmentDetectors:
        nChunks=len(seg[1])
        for n in range(nChunks):
            key_contrib=seg[1][n].GetContributors()[0]
            if part_id==key_contrib:
                seg_length=(seg[1][n].GetStart()-seg[1][n].GetStop()).Mag()/10
                de_dx=-seg[1][n].GetEnergyDeposit()/seg_length
                residual=(seg[1][n].GetStart()-traj[part_id].Points[-1].GetPosition()).Mag()/10
                #if residual>seg_length: res_range.Fill(residual,de_dx)
                res_range.Fill(residual,de_dx)
    return   
############################################################################## Main Code #################################################
proton_mass=938.27 ##MeV/c^2
pion_mass=135.0
nevt=edep_tree.GetEntries()
################ Initialized Variables #################
Proton_Energy_Dep=[]
Pion_Energy_Dep=[]
Proton_KE=[]
Pion_KE=[]
Proton_Reco=[]
Pion_Reco=[]
Proton_Dedx=[]
Pion_Dedx=[]

res_range=RT.TH2D("h","residual range",100,0,30,50,0,80)

for evt in range(nevt):
    edep_tree.GetEntry(evt)
    grtk_tree.GetEntry(evt)
    vtx=edep_tree.Event.Primaries[0]
    primary_pdg=[x.GetPDGCode() for x in vtx.Particles]
    traj=edep_tree.Event.Trajectories
######################## PROTONS ############################
    if primary_pdg==[14,2212]:                            #Checking for ncc events
        print("proton at event",evt)
        for trk in traj: 
            if trk.GetPDGCode()==2212:                    #Pick only the first proton
                proton_id=trk.GetTrackId()
                proton_trk=trk
                break
        #if is_point_contained(traj[proton_id].Points[-1].GetPosition().Vect())==False:
        #        break
        dist,DE,KE=(statistics(evt,vtx,proton_id,trk,proton_mass))
        #residual(proton_id)
        Proton_KE.append(KE)
        Proton_Energy_Dep.append(DE)
        reco_energy,de_dx=reco(proton_id)
        Proton_Reco.append(reco_energy)
        Proton_Dedx.append(de_dx)
########################## PIONS #########################################
    pion_traj=tuple(x for x in vtx.Particles if x.GetPDGCode() in [-211,211])
    for i in primary_pdg: 
        if abs(i)==211:
            for trk in traj: 
                if abs(trk.GetPDGCode())==211:
                    pion_id=trk.GetTrackId()
                    pion_trk=trk
                    break
            if is_point_contained(traj[pion_id].Points[-1].GetPosition().Vect())==False: 
                break
            dist,DE,KE=statistics(evt,vtx,pion_id,pion_trk,pion_mass)
            Pion_Energy_Dep.append(DE)
            residual(pion_id)
            Pion_KE.append(KE)
            reco_energy,de_dx=reco(pion_id)
            Pion_Reco.append(reco_energy)
            Pion_Dedx.append(de_dx)
######## PROTON GRAPHS ##################################          
Proton_Energy=RT.TGraph()
Proton_DE=RT.TGraph()
can1=RT.TCanvas("can","can",1000,800)
can1.cd()
for i in range(len(Proton_KE)):
    Proton_Energy.SetPoint(i,i,Proton_KE[i])
    Proton_DE.SetPoint(i,i,Proton_Energy_Dep[i])
Proton_Energy.SetMarkerColor(9003)
Proton_Energy.GetXaxis().SetTitle("Event Number")
Proton_Energy.GetYaxis().SetTitle("Energy(MeV)")
Proton_DE.Draw("apl")
Proton_Energy.Draw("* same")
RT.gPad.Update()
can1.SaveAs("Proton_Energy_Dep.png")

########### PION GRAPHS ################################
can2=RT.TCanvas("can2","can2",1000,800)
Pion_Energy=RT.TGraph()
Pion_DE=RT.TGraph()
can2.cd()
for i in range(len(Pion_KE)):
    Pion_Energy.SetPoint(i,i,Pion_KE[i])
    Pion_DE.SetPoint(i,i,Pion_Energy_Dep[i])
Pion_Energy.SetLineColor(9003)
Pion_Energy.GetXaxis().SetTitle("Event Number")
Pion_Energy.GetYaxis().SetTitle("Energy(MeV)")
Pion_DE.Draw("apl")
Pion_Energy.Draw("* same")
RT.gPad.Update()
can2.SaveAs("Pion_Energy_Dep.png")

########### BETHE_BLOCH ####################################
pion_mass=139.57 ##MeV
proton_mass= 938.27##MeV
p=np.linspace(20,6000)
proton_de_dx=[]
pion_de_dx=[]
for i in p:
    proton_de_dx.append(plotting_functions.bethe_bloch(i,proton_mass))
    pion_de_dx.append(plotting_functions.bethe_bloch(i,pion_mass))
################### DE_DX #########################################

can3=RT.TCanvas("can3","can3",1000,800)
#Pion_Dep=RT.TGraph()
Pion_Dep=RT.TH2D("h1","Pion de/dx",100,0,6000,50,0,80)
#Proton_Dep=RT.TGraph()
for i in range(len(Pion_KE)):
    Pion_Dep.Fill(Pion_KE[i],np.mean(Pion_Dedx[i][-3:-1]))
Pion_Dep.Draw("colz")
Pion_Curve=RT.TGraph(len(pion_de_dx),(gamma-1)*pion_mass,np.array(pion_de_dx))
Pion_Curve.Draw("same")
Pion_Dep.GetXaxis().SetTitle("Kinetic Energy(MeV)")
Pion_Dep.GetYaxis().SetTitle("de_dx(MeV/cm)")
RT.gPad.Update()
can3.SaveAs("Pion_DEDX_hist.png")

pion_range=plotting_functions.res_range(p,pion_mass,pion_de_dx)
can4=RT.TCanvas("can4","can4",1000,800)
Proton_Dep=RT.TH2D("h2","Proton de/dx",100,0,6000,50,0,80)
for i in range(len(Proton_KE)):
    Proton_Dep.Fill(Proton_KE[i],np.mean(Proton_Dedx[i][-3:-1]))
Proton_Dep.Draw("colz")
Proton_Curve=RT.TGraph(len(proton_de_dx),(gamma-1)*proton_mass,np.array(proton_de_dx))
Proton_Curve.Draw("same")
Proton_Dep.GetXaxis().SetTitle("Kinetic Energy(MeV)")
Proton_Dep.GetYaxis().SetTitle("de_dx(MeV/cm)")
RT.gPad.Update()
can4.SaveAs("Proton_DEDX_hist.png")

############################# Residuals ############################
can5=RT.TCanvas("can5","can5",1000,800)
res_range.Draw("colz")
Pion_Range_Curve=RT.TGraph(len(de_dx),np.array(pion_range),np.array(de_dx))
Pion_Range_Curve.Draw("same")
res_range.GetXaxis().SetTitle("range (cm)")
res_range.GetYaxis().SetTitle("de_dx(MeV/cm")
RT.gPad.Update()
can5.SaveAs("res_range.png")
