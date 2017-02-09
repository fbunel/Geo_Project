from sim2 import terre
#from sim2_v2 import terre
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axes_grid import Grid
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
    params = {'Ri':5E3,
    'Rf':5E5,
    'beta': 0,
    'ta':1*cst['tau_al']/0.717,#1My
    'Ti':300,
    'dt':1E-4*cst['tau_al'],
    'size':1000,
    }
    n = int(params['ta']/params['dt'])


    params['beta'] = 0
    t1 = terre(**params,cst=cst)
    print("courbe 1")
    for _ in range(n):
        t1.update_P()
        t1.update_m()
        #t.P[:]=0
        t1.step()
        t1.fusion()
        print("t: {:.3} My    ".format(t1.t*0.717),end='\r')
        #print("Energie totale : {}".format((t1.t*(t1.r**2)*t1.R).sum()))


    params['beta'] = 1
    t2 = terre(**params,cst=cst)
    print("courbe 2")
    for _ in range(n):
        t2.update_P()
        t2.update_m()
        #t.P[:]=0
        t2.step()
        t2.fusion()
        print("t: {:.3} My    ".format(t2.t*0.717),end='\r')

    params['beta'] = 2
    t3 = terre(**params,cst=cst)
    print("courbe 3")
    for _ in range(n):
        t3.update_P()
        t3.update_m()
        #t.P[:]=0
        t3.step()
        t3.fusion()
        print("t: {:.3} My    ".format(t3.t*0.717),end='\r')

    """
    params['beta'] = -1
    t4 = terre(**params,cst=cst)
    print("courbe 4")
    for _ in range(n):
        t4.update_P()
        t4.update_m()
        #t.P[:]=0
        t4.step()
        t4.fusion()
        print("t: {:.3} My    ".format(t4.t*0.717),end='\r')

    """
    params['beta'] = -1
    params['Ri'] = params['Rf']
    t5 = terre(**params,cst=cst)
    print("courbe 5")
    for _ in range(n):
        t5.update_P()
        t5.update_m()
        #t.P[:]=0
        t5.step()
        t5.fusion()
        print("t: {:.3} My    ".format(t5.t*0.717),end='\r')


    print()

    plt.figure(figsize=(8,6))
    plt.plot(t1.r[:-1]/t1.r[-2],t1.T[:-1]*t1.T0,'b',label=r"$\beta = 0$")
    plt.plot(t2.r[:-1]/t2.r[-2],t2.T[:-1]*t2.T0,'g',label=r"$\beta = 1$")
    plt.plot(t3.r[:-1]/t3.r[-2],t3.T[:-1]*t3.T0,'r',label=r"$\beta = 2$")
#plt.plot(t4.T[:-1]*t4.T0,t4.r[:-1]/t4.r[-2],'k-.',label=r"$\beta = -1$")
    plt.plot(t5.r[:-1]/t5.r[-2],t5.T[:-1]*t5.T0,'k-.',label=r"$\beta = \emptyset $")
    plt.ylabel(r'T (en K)')
    plt.xlabel(r'r/R')
    #plt.title('Température après : {0:.2}My'.format(t1.t*0.717))
    plt.legend(loc='best')
    from datetime import datetime
    plt.savefig("graph_sim2_{}.pdf".format(datetime.now().strftime('%Y%m%d-%H%M%S')),bbox_inches='tight')
    
def graphe2():
    params = {'Ri':5E3,
    'Rf':5E5,
    'beta': 0,
    'ta':5*cst['tau_al']/0.717,#5My
    'Ti':300,
    'dt':1E-4*cst['tau_al'],
    'size':1000,
    }
    n = int(params['ta']/params['dt'])
#n=10


    params['beta'] = 0
    t1 = terre(**params,cst=cst)
    print("courbe 1")
    for _ in range(n):
        t1.update_P()
        t1.update_m()
        #t.P[:]=0
        t1.step()
        t1.fusion()
        print("t: {:.3} My    ".format(t1.t*0.717),end='\r')
        #print("Energie totale : {}".format((t1.t*(t1.r**2)*t1.R).sum()))


    params['beta'] = 1
    t2 = terre(**params,cst=cst)
    print("courbe 2")
    for _ in range(n):
        t2.update_P()
        t2.update_m()
        #t.P[:]=0
        t2.step()
        t2.fusion()
        print("t: {:.3} My    ".format(t2.t*0.717),end='\r')

    params['beta'] = 2
    t3 = terre(**params,cst=cst)
    print("courbe 3")
    for _ in range(n):
        t3.update_P()
        t3.update_m()
        #t.P[:]=0
        t3.step()
        t3.fusion()
        print("t: {:.3} My    ".format(t3.t*0.717),end='\r')

    """
    params['beta'] = -1
    t4 = terre(**params,cst=cst)
    print("courbe 4")
    for _ in range(n):
        t4.update_P()
        t4.update_m()
        #t.P[:]=0
        t4.step()
        t4.fusion()
        print("t: {:.3} My    ".format(t4.t*0.717),end='\r')

    
    params['beta'] = -1
    params['Ri'] = params['Rf']
    t5 = terre(**params,cst=cst)
    print("courbe 5")
    for _ in range(n):
        t5.update_P()
        t5.update_m()
        #t.P[:]=0
        t5.step()
        t5.fusion()
        print("t: {:.3} My    ".format(t5.t*0.717),end='\r')
    """

    print()

    plt.figure(figsize=(8,6))
    plt.plot(t1.r[:-1]/t1.r[-2],t1.T[:-1]*t1.T0,'b',label=r"$\beta = 0$")
    plt.plot(t2.r[:-1]/t2.r[-2],t2.T[:-1]*t2.T0,'g',label=r"$\beta = 1$")
    plt.plot(t3.r[:-1]/t3.r[-2],t3.T[:-1]*t3.T0,'r',label=r"$\beta = 2$")
    #plt.plot(t4.T[:-1]*t4.T0,t4.r[:-1]/t4.r[-2],'k-.',label=r"$\beta = -1$")
    #plt.plot(t5.T[:-1]*t5.T0,t5.r[:-1]/t5.r[-2],'k-.',label=r"$\beta = \emptyset $")
    plt.ylabel(r'T (en K)')
    plt.xlabel(r'r/R')
    #plt.title('Température après : {0:.3}My'.format(t1.t*0.717))
    plt.legend(loc='best')
    from datetime import datetime
    plt.savefig("graph_sim2_{}.pdf".format(datetime.now().strftime('%Y%m%d-%H%M%S')),bbox_inches='tight')

def graphe_R():
    global R1,R2,t1,t2
    params = {'Ri':5E3,
    'Rf':5E5,
    'beta': 0,
    'ta':1*cst['tau_al']/0.717,#5My
    'Ti':300,
    'dt':1E-3*cst['tau_al']/0.717,
    'size':1000,
    }
    n = int(params['ta']/params['dt'])
#n=10

    print("courbe 1")
    params['beta'] = 0
    t1 = terre(**params,cst=cst)
    t1.d1,t1.d2,t1.d3,t1.u1,t1.u2,t1.u3=[0]*6
    R1 = []
    for _ in range(n):
        t1.update_m()
        t1.t += t1.dt
        R1.append([t1.t*0.717,t1.R/params['Rf']])
        print("t: {:.3} My    ".format(t1.t*0.717),end='\r')
    
    
    print("courbe 2")
    params['beta'] = 1
    t2 = terre(**params,cst=cst)
    t2.d1,t2.d2,t2.d3,t2.u1,t2.u2,t2.u3=[0]*6
    R2 = []
    for _ in range(n):
        t2.update_m()
        t2.t += t2.dt
        R2.append([t2.t*0.717,t2.R/params['Rf']])
        print("t: {:.3} My    ".format(t2.t*0.717),end='\r')

    print("courbe 3")
    params['beta'] = 2
    t3 = terre(**params,cst=cst)
    t3.d1,t3.d2,t3.d3,t3.u1,t3.u2,t3.u3=[0]*6
    R3 = []
    for _ in range(n):
        t3.update_m()
        t3.t += t3.dt
        R3.append([t3.t*0.717,t3.R/params['Rf']])
        print("t: {:.3} My    ".format(t3.t*0.717),end='\r')


    plt.figure(figsize=(8,5))
    plt.plot([0,0.78,1],[0,0.23,1],'k-.',color = '0.5')
    plt.plot([0,0.95,1],[0,0.08,1],'k-.',color = '0.5')
    plt.plot(*list(zip(*R1)),label=r"$\beta = 0$")
    plt.plot(*list(zip(*R2)),label=r"$\beta = 1$")
    plt.plot(*list(zip(*R3)),label=r"$\beta = 2$")


    plt.xlabel(r't (en My)')
    plt.ylabel(r'R(t)/R$_{final}$')
    #plt.title('Évolution du rayon')
    
    plt.legend(loc="best")
    from datetime import datetime
    plt.savefig("graph_sim2_fig2_3_{}.pdf".format(datetime.now().strftime('%Y%m%d-%H%M%S')),bbox_inches='tight')
    
def graphes_T_moy():

  
    print("T moy planétésimal rayon fixe 500km")

    params = {'Ri':5E3,
    'Rf':5E5,
    'beta': 0,
    'ta':1*cst['tau_al']/0.717,#1My
    'Ti':300,
    'dt':1E-4*cst['tau_al'],
    'size':1000,
    }
    n = int(params['ta']/params['dt'])
    fig = plt.figure(figsize=(8,6))
    
    params['beta'] = -1
    params['Ri'] = params['Rf']
    t5 = terre(**params,cst=cst)
    T5 = []
    for _ in range(n):
        t5.update_P()
        t5.update_m()
        t5.step()
        t5.fusion()
        T5.append([t5.t*0.717,300*(t5.T*t5.r**2).sum()/(t5.r**2).sum()])
        #print(T5[-1])
        print("t: {:.3} My    ".format(t5.t*0.717),end='\r')
        
    T5 = np.array(T5)
    #print(T5)
    plt.plot(T5[:,0],T5[:,1],'k--',lw=2,label=r'$\beta = \emptyset$')
        
    print("T moy du planétésimal de réference (R(t) 5km -> 500km)")
    
    
    
    T = np.load("T_pla_1My.npy")
    plt.plot(T[0],T[4],lw=2,label=r'$\beta = 0$')
    plt.plot(T[0],T[5],lw=2,label=r'$\beta = 1$')
    plt.plot(T[0],T[6],lw=2,label=r'$\beta = 2$')
    
    
    plt.ylim(0,2000)
    plt.ylabel(r'T (en K)')
    plt.xlabel(r't (en My)')
    plt.legend(loc='best')

    plt.savefig("graph_sim2_fig2_4_pla_{}.pdf".format(datetime.now().strftime('%Y%m%d-%H%M%S')),bbox_inches='tight')






graphe1()
graphe2()
#graphe_R()
#graphes_T_moy()
plt.show()
