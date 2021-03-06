from sim1 import terre
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('lines', linewidth=4)
plt.rcParams['lines.linewidth'] = 2


#on va calculer T à t=+infini(~20My) pour 4 CI différentes :
# 1.  t0 = 0My  , avec fusion des materiaux
# 2.  t0 = 1My  , avec fusion des materiaux
# 3.  t0 = 0My  , sans fusion des materiaux
# 4.  t0 = 1My  , sans fusion des materiaux

cst = { 'tau_al':  2.26E13, #s
	'H0'    :  1.5E-7, #W kg-1
	'T_neb' :  300   , #K
	'sigma' :  5.67E-8,#W m-2 K-4
	'phi'   :  0.18  , # sans unité
	'kT'    :  11.48 , #W K-1 m-1
	'kT_m'  :  50    , #W K-1 m-1
	'kT_s'  :  3     , #W K-1 m-1
	'Cp'    :  939   , #J K-1 kg-1
	'Cp_m'  :  450   , #J K-1 kg-1
	'Cp_s'  :  1200  , #J K-1 kg-1
	'rho'   :  4028  , #kg m-3 
	'rho_m' :  7800  , #kg m-3  
	'rho_s' :  3200  , #kg m-3  
	'Tfus_m':  1261  , #K
	'Tfus_s':  1408  , #K
	'L_m'   :  2.5E5 , #J kg-1
	'L_s'   :  5.0E5 , #J kg-1
	#'L_m'   :  1 , #J kg-1
	#'L_s'   :  1 , #J kg-1
	}

def graphe1():
    print("Calcule la figure 1 de l'article")
    dt=0.0001
    n = int(1/(dt*0.717))

    t1 = terre(Ri=500000,Ti=300,dt=dt,size=1000,cst=cst)
    print("graphe 1")
    for _ in range(n):
            t1.update_P()
            #t.P[:]=0
            t1.step()
            t1.fusion()
            print("t: {:.3} My     ".format(t1.t*0.717),end="\r")
    print()

    cst['L_m'] =  1  # L_m très petit ~= pas d'influence de la fusion
    cst['L_s'] =  1  #J kg-1

    t3 = terre(Ri=500000,Ti=300,dt=dt,size=1000,cst=cst)
    print("graphe 3")
    for _ in range(n):
            t3.update_P()
            #t.P[:]=0
            t3.step()
            t3.fusion()
            print("t: {:.3} My     ".format(t3.t*0.717),end="\r")
    print()


    plt.figure(figsize=(8,6))
    plt.plot(t3.r[:-1]/t3.r[-2],t3.T[:-1]*t3.T0,'b',label="Sans chaleur latente")
    plt.plot(t1.r[:-1]/t1.r[-2],t1.T[:-1]*t1.T0,'r',label="Avec chaleur latente")
    plt.ylabel(r"T (en K)")
    plt.xlabel(r"r/R")
    #plt.title(r"Temp\'erature apr\`es : {0:.2}My".format(t1.t*0.717))
    plt.legend(loc="best")
    plt.savefig("graph_sim1_fig1.pdf",bbox_inches='tight')

def graphe2():
    print("Calcule un graphe qui compare les 2 deux régimes de diffusion selon le rayon")
    dt=0.0001
    n = int(1/(dt*0.717))

    t1 = terre(Ri=5000,Ti=300,dt=dt,size=1000,cst=cst)
    print("graphe 1")
    for _ in range(n):
            t1.update_P()
            #t.P[:]=0
            t1.step()
            t1.fusion()
            print("t: {:.3} My     ".format(t1.t*0.717),end="\r")
    print()

    print("graphe 2")
    t2 = terre(Ri=30000,Ti=300,dt=dt,size=1000,cst=cst)
    for _ in range(n):
            t2.update_P()
            t2.step()
            t2.fusion()
            print("t: {:.3} My     ".format(t2.t*0.717),end="\r")
    print()

    t3 = terre(Ri=200000,Ti=300,dt=dt,size=1000,cst=cst)
    print("graphe 3")
    for _ in range(n):
            t3.update_P()
            #t.P[:]=0
            t3.step()
            t3.fusion()
            print("t: {:.3} My     ".format(t3.t*0.717),end="\r")
    print()

    plt.figure(figsize=(8,6))
    plt.plot(t3.r[:-1]/t3.r[-2],t3.T[:-1]*t3.T0,'b',lw=2,label=r"$R \gg L_d$")
    plt.plot(t2.r[:-1]/t2.r[-2],t2.T[:-1]*t2.T0,'g',lw=2,label=r"$R \simeq L_d$")
    plt.plot(t1.r[:-1]/t1.r[-2],t1.T[:-1]*t1.T0,'r',lw=2,label=r"$R \ll L_d$")
    
    
    plt.legend(loc="best")
    plt.ylabel(r"T (en K)")
    plt.xlabel(r"r/R")
    #plt.title('Température après : {0:.2}My'.format(t1.t*0.717))
    plt.savefig("graph_sim1_fig2.pdf",bbox_inches='tight')
    
    
def graphes_T_moy():

  
    print("T moy planétésimal rayon fixe 500km")
    dt=0.0001
    n = int(2/(dt*0.717))

    R = 5000
    
    fig = plt.figure(figsize=(8,6))
    
    
    t5 = terre(Ri=R,Ti=300,dt=dt,size=1000,cst=cst)
    T5 = []
    for _ in range(n):
        t5.update_P()
        t5.step()
        t5.fusion()
        T5.append([t5.t*0.717,300*(t5.T*t5.r**2).sum()/(t5.r**2).sum()])
        #print(T5[-1])
        print("t: {:.3} My    ".format(t5.t*0.717),end='\r')
        
    T5 = np.array(T5)
    plt.plot(T5[:,0],T5[:,1],lw=2,label=r'$R = {} km$'.format(R/1000))
        

    #plt.ylim(0,500)
    plt.ylabel(r'T (en K)')
    plt.xlabel(r't (en My)')
    #plt.legend(loc='best')

    plt.savefig("graph_sim2_fig2_4_pla_{}.pdf".format(datetime.now().strftime('%Y%m%d-%H%M%S')),bbox_inches='tight')

#graphe1()
#graphe2()
graphes_T_moy()
plt.show()
