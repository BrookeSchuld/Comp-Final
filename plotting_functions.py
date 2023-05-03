import numpy as np
import math
import ROOT as RT
from ROOT import TLorentzVector
import sys

#Function takes the momentum array and mass for a particle and returns the expected bethe bloch de/dx value as an array in MeV/cm 
def bethe_bloch(p,mass):
    electron_mass=.511 #MeV
    c=3*10**8
    LAr_density=1.38 #g/cm^3
    K=.307075 #MeE/mol cm^2
    z=1
    LAr_mass=39.948 #g/mol
    LAr_number=18
    I=188*10**(-6)      #Needs to be in MeV despite what the PDG seems to say
    B=2.0*np.log(10)
    X0=.201
    X1=3.00
    a  = 0.196
    m  = 3.000
    C  = 5.217
    p=np.array(p)    
    beta=p/np.sqrt(p**2+mass**2)
    gamma=1/np.sqrt(1-beta**2)
    de_dx=[]
    for i in range(len(p)):   #the density correction factor is different depending on the value of the momentum
        Wmax=2*electron_mass*beta[i]**2*gamma[i]**2/(1+2*gamma[i]*electron_mass/c**2/mass+(electron_mass/c**2/mass)**2)
        delta=0
        X = np.log10(beta[i] * gamma[i])
        if X0<X and X<X1: delta=(B*X+a*np.power(X1-X,m)-C)
        if X>X1: delta=(B*X-C)
        de_dx.append(K*LAr_density*z**2*LAr_number/LAr_mass/beta[i]**2*(1/2*np.log(2*electron_mass*beta[i]**2*gamma[i]**2*Wmax/I**2)-beta[i]**2-.5*delta))
    return de_dx

#Function takes the momentum and mass of a particle and returns the both the de/dx and the residual range in cm
def res_range(p,mass):
    de_dx=np.array(bethe_bloch(p,mass))   #Requires integrating the bethe bloch equation
    p=np.array(p)
    E=np.sqrt(p**2+mass**2)
    delE=np.diff(E)
    residual=[]
    for i in range(len(p)):
        integrand=np.trapz(1/np.array(de_dx[0:i])*delE[0:i])
        residual.append(integrand)
    return (de_dx,residual)               
       
def is_point_contained(pos):                 ##Check to see if the position is within the 2x2. These numbers are from the specific geometry and are hard coded
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
       
def reco(evt,part_id):                                                 #Gets the energy deposited in the detector for a given particle
    
    edep_tree.GetEntry(evt)
    grtk_tree.GetEntry(evt)
    traj=edep_tree.Event.Trajectories
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
                        if seg_length==0: continue
                        de_dx.append(-seg[1][n].GetEnergyDeposit()/seg_length)
                    #for i in proton_children:
                    #    if i==key_contrib:
                    #if proton_id==key_contrib:
                     #       reco_energy+=seg[1][n].GetEnergyDeposit()
    return reco_energy,de_dx   

def residual(evt,part_id):                                                                 #Gets the residual range for a given particle
    edep_tree.GetEntry(evt)
    grtk_tree.GetEntry(evt)
    traj=edep_tree.Event.Trajectories
    de_dx_array=[]
    residual_array=[]
    for seg in edep_tree.Event.SegmentDetectors:
        nChunks=len(seg[1])
        for n in range(nChunks):
            key_contrib=seg[1][n].GetContributors()[0]
            if part_id==key_contrib:
                seg_length=(seg[1][n].GetStart()-seg[1][n].GetStop()).Mag()/10
                de_dx=-seg[1][n].GetEnergyDeposit()/seg_length
                residual=(seg[1][n].GetStart()-traj[part_id].Points[-1].GetPosition()).Mag()/10
                if 30<abs(residual)<seg_length: continue
                de_dx_array.append(de_dx)
                residual_array.append(np.abs(residual))
    return np.array(de_dx_array), np.array(residual_array)

#Main part of the code that accesses the information in the input file before passing it the called functions
edep_tree=RT.TChain("EDepSimEvents")                                   
grtk_tree=RT.TChain("DetSimPassThru/gRooTracker")
filelist=[sys.argv[x] for x in range(1,len(sys.argv))]
for file in filelist:
    edep_tree.Add(file)
    grtk_tree.Add(file)

