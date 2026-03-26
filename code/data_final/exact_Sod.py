# -*- coding: utf-8 -*-
import numpy as np

######################################################################
def pressure(pini,roz,pz,cz,gamma):
    g1=(gamma-1.)/(2.*gamma)
    g2=(gamma+1.)/(2.*gamma)
    g4=2.0/(gamma-1.0)
    g5=2.0/(gamma+1.)
    g6=(gamma-1.0)/(gamma+1.)
    if pini <= pz:
        # detente
        prat=pini/pz
        f=g4*cz*(prat**g1-1.0)
        df=(1.0/(roz*cz))*prat**(-g2)
    else:
        # ondes de chocs
        az=g5/roz
        bz=g6*pz
        qrt=np.sqrt(az/(bz+pini))
        f=(pini-pz)*qrt
        df=(1.0-0.5*(pini-pz)/(bz+pini))*qrt
    return f,df
            

######################################################################
# Resolution du pbm de Riemann
def solriemann(rol,ul,pl,ror,ur,pr,xl,xr,gamma,cl,cr):
    #calcul des valeurs de la pression et de la vitesse 
    # dans la region intermediaire (etoile) via methode de Newton
    #   utilisation d une methoe de newton
    xeps=1.0E-06
    nitmax=40
    du=ur-ul
    # calcul de la vitesse critique
    g4=2.0/(gamma-1.0)
    ducrit=g4*(cl+cr)-du
    # test pour la naissance du vide
    if (ducrit <= 0.0):
        print("Apparition du vide")
        exit()
    # resolution
    # initialisation
    pini = inivalue(rol,ul,pl,ror,ur,pr,gamma,cl,cr)
    # newton
    p0=pini
    keepon= True
    nit=1
    while (nit <= nitmax) and keepon:
        fl,dfl = pressure(pini,rol,pl,cl,gamma)
        fr,dfr = pressure(pini,ror,pr,cr,gamma)
        pini=pini-(fl+fr+du)/(dfl+dfr)
        dpini=2.0*np.abs(pini-p0)/np.abs(pini+p0)
        # test de convergence         
        if dpini <= xeps:
            keepon = False
        # securite
        if keepon:
            if pini <= 0.0:
                pini=xeps
            p0=pini
        nit=nit+1
    # calcul pression et vitesse
    pet=pini
    uet=0.5*(ul+ur+fr-fl)
  
    return pet,uet
            
######################################################################
#    calcul de la valeur initiale
def inivalue(rol,ul,pl,ror,ur,pr,gamma,cl,cr):
    g1=(gamma-1.)/(2.*gamma)
    g3=2.0*gamma/(gamma-1.0)
    g5=2.0/(gamma+1.)
    g6=(gamma-1.0)/(gamma+1.)
    g7=0.5*(gamma-1.0)
    xeps=1.0E-06
    qmax=2.
    # valeur du solver PVRS de riemann
    pv=0.5*(pl+pr)-0.125*(ur-ul)*(rol+ror)*(cl+cr)
    pmin=np.minimum(pl,pr)
    pmax=np.maximum(pl,pr)
    qrat=pmax/pmin

    zdrap=(pmin <= pv) and (pv <= pmax)
    zdrap=(qrat <= qmax) and zdrap
    if (pmin <= pv) and (pv <= pmax) and (qrat<=qmax):
        pini=max(xeps,pv)
    else:
        if pv < pmin:
            # solution correspondant a deux detentes
            pnu=cl+cr-g7*(ur-ul)
            pde=cl/pl**g1+cr/pr**g1
            pini=(pnu/pde)**g3
        else:
            # solution correspondant a deux chocs
            gel=np.sqrt((g5/rol)/(g6*pl+max(xeps,pv)))
            ger=np.sqrt((g5/ror)/(g6*pr+max(xeps,pv)))
            pini=(gel*pl+ger*pr-(ur-ul))/(gel+ger)
            pini=max(xeps,pini)

    gel=np.sqrt((g5/rol)/(g6*pl+max(xeps,pv)))
    ger=np.sqrt((g5/ror)/(g6*pr+max(xeps,pv)))
    #write(*,*) 'pts ',(gel*pl+ger*pr-(ur-ul))/(gel+ger)
    pnu=cl+cr-g7*(ur-ul)
    pde=cl/pl**g1+cr/pr**g1
    pini=(pnu/pde)**g3
    #write(*,*) 'ptr ',(pnu/pde)**g3

    return pini



######################################################################
# determination de la solution en differents points
def sample(rol,ul,pl,ror,ur,pr,cr,cl,pet,uet,vs,gamma):
    g1=(gamma-1.)/(2.*gamma)
    g2=(gamma+1.)/(2.*gamma)
    g3=2.0*gamma/(gamma-1.0)
    g4=2.0/(gamma-1.0)
    g5=2.0/(gamma+1.)
    g6=(gamma-1.0)/(gamma+1.)
    g7=0.5*(gamma-1.0)
    g8=1.0/gamma
    g9=gamma-1.0
    if vs <= uet:
        # point situe a gauche de la surface de glissement
        if pet <= pl:
            # detente entre l et *
            shl=ul-cl
            # position du point par rapport a la detente
            if vs <= shl:
                ros=rol
                us=ul
                ps=pl
            else:
                cml=cl*(pet/pl)**g1
                stl=uet-cml
                if vs > stl:
                    ros=rol*(pet/pl)**g8
                    us=uet
                    ps=pet
                else:
                    us=g5*(cl+g7*ul+vs)
                    cs=g5*(cl+g7*(ul-vs))
                    ros=rol*(cs/cl)**g4
                    ps=pl*(cs/cl)**g3
        else:
            # choc entre l et *
            pml=pet/pl
            sl=ul-cl*np.sqrt(g2*pml+g1)
            if vs <= sl:
                # etat gauche
                ros=rol
                us=ul
                ps=pl
            else:
                # derriere le choc
                ros=rol*(pml+g6)/(pml*g6+1.0)
                us=uet
                ps=pet               
    else:
        # point a droite de la detente
        if pet > pr:
            # choc entre * et r
            pmr=pet/pr
            sr=ur+cr*np.sqrt(g2*pmr+g1)
            if vs >= sr:
                ros=ror
                us=ur
                ps=pr
            else:
                ros=ror*(pmr+g6)/(pmr*g6+1.0)
                us=uet
                ps=pet
        else:
            shr=ur+cr
            if vs >= shr:
                ros=ror
                us=ur
                ps=pr
            else:
                cmr=cr*(pet/pr)**g1
                sstr=uet+cmr
                if vs <= sstr:
                    ros=ror*(pet/pr)**g8
                    us=uet
                    ps=pet
                else:
                    us=g5*(-cr+g7*ur+vs)
                    cs=g5*(cr-g7*(ur-vs))
                    ros=ror*(cs/cr)**g4
                    ps=pr*(cs/cr)**g3
    return ps,us,ros

######################################################################
# Solution pbm Riemann          
def riesolv(rol,ul,pl,ror,ur,pr,xl,xr,mprof,gamma,timeout,solution):
    # Initialisation
    g1=(gamma-1.)/(2.*gamma)
    g2=(gamma+1.)/(2.*gamma)
    g3=2.0*gamma/(gamma-1.0)
    g4=2.0/(gamma-1.0)
    g5=2.0/(gamma+1.)
    g6=(gamma-1.0)/(gamma+1.)
    g7=0.5*(gamma-1.0)
    g8=1.0/gamma
    g9=gamma-1.0
    # vitesses du son
    cl=np.sqrt(gamma*pl/rol)
    cr=np.sqrt(gamma*pr/ror)
    # Resolution pbm Riemann
    pet,uet = solriemann(rol,ul,pl,ror,ur,pr,xl,xr,gamma,cl,cr)
    # Calcul des profils
    dx=(xr-xl)/mprof
    x_membrane = 0.5
    for i in range(mprof+1):
        xs=xl+i*dx
        vs = (xs - x_membrane) / timeout
        ps,us,ros = sample(rol,ul,pl,ror,ur,pr,cr,cl,pet,uet,vs,gamma)
        cs=np.sqrt(gamma*ps/ros)
        E_massique = (ps / (ros * (gamma - 1.0))) + 0.5 * (us**2)
        xpos = xs
        density = ros
        velocity = us
        pression = ps
        temperature = ps/ros/g9
        mach = us/cs
        entropie = pression/density**gamma
        solution.append((xpos,density,mach, pression,temperature,velocity,entropie))
    return solution

	
# --- CONFIGURATION DU CAS ---
time = 0.20 # par defaut
cas = 7
gamma = 1.4
R = 287.
xdeb = 0.0
xfin = 1.0
npt = 201

# --- SELECTION DES CONDITIONS INITIALES ---
if cas == 1: # Sod subsonique
    rhoL, uL, pL = 1.0, 0.0, 1.0
    rhoR, uR, pR = 0.125, 0.0, 0.1
elif cas == 2: # Sod supersonique
    rhoL, uL, pL = 1.0, 0.75, 1.0
    rhoR, uR, pR = 0.125, 0.0, 0.1
elif cas == 3: # Apparition de vide (double detente)
    rhoL, uL, pL = 1.0, -2.0, 0.4
    rhoR, uR, pR = 1.0, 2.0, 0.4
elif cas == 4: # Ligne de glissement stationnaire
    rhoL, uL, pL = 10.0, 0.0, 10.0
    rhoR, uR, pR = 0.1, 0.0, 10.0
elif cas == 5: # Choc stationnaire 
    rhoL, uL, pL = 1.0, 3.5496478, 1.0
    rhoR, uR, pR = 3.85714285, 0.920279072, 10.33333333
elif cas == 6: # # Blast Wave (Forte pression)
    rhoL, uL, pL = 1.0, 0.0, 1000.0
    rhoR, uR, pR = 1.0, 0.0, 0.01
    time = 0.012
elif cas == 7: # Collision de chocs 
    rhoL, uL, pL = 5.99924, 19.5975, 460.894
    rhoR, uR, pR = 5.99242, -6.19633, 46.0950
    time = 0.035
else:
    print("Cas non defini")
    exit()

print(f"Calcul de la solution exacte - Cas {cas}")
print(f"Etat Gauche : {rhoL}, {uL}, {pL}")
print(f"Etat Droit  : {rhoR}, {uR}, {pR}")
print(f"Temps final : {time}")

# --- CALCUL DE LA SOLUTION ---
cl = np.sqrt(gamma * pL / rhoL)
cr = np.sqrt(gamma * pR / rhoR)

# Resolution du probleme de Riemann
pet, uet = solriemann(rhoL, uL, pL, rhoR, uR, pR, xdeb, xfin, gamma, cl, cr)

# Generation des profils
sol_list = []
riesolv(rhoL, uL, pL, rhoR, uR, pR, xdeb, xfin, npt, gamma, time, sol_list)
result = np.asarray(sol_list)

# --- SAUVEGARDE ---
#header = "x density mach pressure temperature velocity entropy"
np.savetxt("data_final/Exact_Sod.dat", result, fmt='%.8e')

print("Fichier data_final/Exact_Sod.dat enregistre avec succes.")
print("end")