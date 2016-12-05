from sim2 import terre
#from sim2_v2 import terre
import matplotlib.pyplot as plt
import numpy as np


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
	'Cp'    :  1065  , #J K-1 kg-1
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
    print("graphe 1")
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
    print("graphe 2")
    for _ in range(n):
        t2.update_P()
        t2.update_m()
        #t.P[:]=0
        t2.step()
        t2.fusion()
        print("t: {:.3} My    ".format(t2.t*0.717),end='\r')

    params['beta'] = 2
    t3 = terre(**params,cst=cst)
    print("graphe 3")
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
    print("graphe 4")
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
    print("graphe 5")
    for _ in range(n):
        t5.update_P()
        t5.update_m()
        #t.P[:]=0
        t5.step()
        t5.fusion()
        print("t: {:.3} My    ".format(t5.t*0.717),end='\r')


    print()

    plt.figure(figsize=(6,8))
    plt.plot(t1.T[:-1]*t1.T0,t1.r[:-1]/t1.r[-2],'b',label=r"$\beta = 0$")
    plt.plot(t2.T[:-1]*t2.T0,t2.r[:-1]/t2.r[-2],'g',label=r"$\beta = 1$")
    plt.plot(t3.T[:-1]*t3.T0,t3.r[:-1]/t3.r[-2],'r',label=r"$\beta = 2$")
#plt.plot(t4.T[:-1]*t4.T0,t4.r[:-1]/t4.r[-2],'k-.',label=r"$\beta = -1$")
    plt.plot(t5.T[:-1]*t5.T0,t5.r[:-1]/t5.r[-2],'k-.',label=r"$\beta = \varnothing $")
    plt.xlabel('Température en K')
    plt.ylabel('Rayon adimentioné')
    plt.title('Température après : {0:.2}My'.format(t1.t*0.717))
    plt.legend(loc='best')
    from datetime import datetime
    plt.savefig("graph_sim2_{}.png".format(datetime.now().strftime('%Y%m%d-%H%M%S')))
    


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
    params['dt'] = 5E-6*cst['tau_al']
    t1 = terre(**params,cst=cst)
    print("graphe 1")
    for _ in range(n):
        t1.update_P()
        t1.update_m()
        #t.P[:]=0
        t1.step()
        t1.fusion()
        print("t: {:.3} My    ".format(t1.t*0.717),end='\r')
        #print("Energie totale : {}".format((t1.t*(t1.r**2)*t1.R).sum()))


    params['dt'] = 1E-4*cst['tau_al']
    params['beta'] = 1
    t2 = terre(**params,cst=cst)
    print("graphe 2")
    for _ in range(n):
        t2.update_P()
        t2.update_m()
        #t.P[:]=0
        t2.step()
        t2.fusion()
        print("t: {:.3} My    ".format(t2.t*0.717),end='\r')

    params['beta'] = 2
    t3 = terre(**params,cst=cst)
    print("graphe 3")
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
    print("graphe 4")
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
    print("graphe 5")
    for _ in range(n):
        t5.update_P()
        t5.update_m()
        #t.P[:]=0
        t5.step()
        t5.fusion()
        print("t: {:.3} My    ".format(t5.t*0.717),end='\r')


    print()

    plt.figure(figsize=(6,8))
    plt.plot(t1.T[:-1]*t1.T0,t1.r[:-1]/t1.r[-2],'b',label=r"$\beta = 0$")
    plt.plot(t2.T[:-1]*t2.T0,t2.r[:-1]/t2.r[-2],'g',label=r"$\beta = 1$")
    plt.plot(t3.T[:-1]*t3.T0,t3.r[:-1]/t3.r[-2],'r',label=r"$\beta = 2$")
#plt.plot(t4.T[:-1]*t4.T0,t4.r[:-1]/t4.r[-2],'k-.',label=r"$\beta = -1$")
    plt.plot(t5.T[:-1]*t5.T0,t5.r[:-1]/t5.r[-2],'k-.',label=r"$\beta = \varnothing $")
    plt.xlabel('Température en K')
    plt.ylabel('Rayon adimentioné')
    plt.title('Température après : {0:.3}My'.format(t1.t*0.717))
    plt.legend(loc='best')
    from datetime import datetime
    plt.savefig("graph_sim2_{}.png".format(datetime.now().strftime('%Y%m%d-%H%M%S')))



#graphe1()
graphe2()
plt.show()
