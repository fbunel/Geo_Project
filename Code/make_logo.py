from sim2 import terre
#from sim2_v2 import terre
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axes_grid import Grid
import numpy as np
from datetime import datetime

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
    params = {' Ri':5E3,
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

    
def graphe2():
    params = {'Ri':5E3,
    'Rf':5E5,
    'beta': 0,
    'ta':5*cst['tau_al']/0.717,#5My
    'Ti':300,
    'dt':1E-2*cst['tau_al'],
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

    """
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

    quart(t1)
    #quart(t2)
    #quart(t3)

    plt.savefig("logo_{}.eps".format(datetime.now().strftime('%Y%m%d-%H%M%S')))
    plt.show()

def quart(self):
    size = self.T.shape[0]
    y, x = np.mgrid[0:size, 0:size]
    r = np.sqrt(x**2 + y**2).astype(np.uint)
    valid = np.where(r < size)

    d = np.empty((size,size))
    d[valid] = self.T[r[valid]]*self.T0

    plt.figure()
    plt.imshow(np.ma.array(d, mask=(r >= size)),cmap='hot',interpolation='nearest',vmin=0)
    plt.colorbar()

#graphe1()
graphe2()
plt.show()
