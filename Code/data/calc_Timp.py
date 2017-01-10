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

def calc_T_1My():
    global R1,R2,t1,t2
    params = {'Ri':5E3,
    'Rf':5E5,
    'beta': 0,
    'ta':1*cst['tau_al']/0.717,#1My
    'Ti':300,
    'dt':1E-3*cst['tau_al']/0.717,
    'size':1000,
    }
    n = int(params['ta']/params['dt'])
#n=10

    print("courbe 1")
    params['beta'] = 0
    t1 = terre(**params,cst=cst)
    #t1.d1,t1.d2,t1.d3,t1.u1,t1.u2,t1.u3=[0]*6
    T1 = []
    R1 = []
    t = []
    for _ in range(n):
        t1.update_P()
        t1.update_m()
        #t.P[:]=0
        t1.step()
        t1.fusion()
        T1.append((t1.T*(t1.r**2)).sum()/(t1.r**2).sum())
        R1.append(t1.R/params['Rf'])
        t.append(t1.t*0.717)
        print("t: {:.3} My    ".format(t1.t*0.717),end='\r')
    
    print("courbe 2")
    params['beta'] = 1
    t2 = terre(**params,cst=cst)
    T2 = []
    R2 = []
    for _ in range(n):
        t2.update_P()
        t2.update_m()
        t2.step()
        t2.fusion()
        R2.append(t2.R/params['Rf'])
        T2.append((t2.T*(t2.r**2)).sum()/(t2.r**2).sum())
        print("t: {:.3} My    ".format(t2.t*0.717),end='\r')

    print("courbe 3")
    params['beta'] = 2
    t3 = terre(**params,cst=cst)
    T3 = []
    R3 = []
    for _ in range(n):
        t3.update_P()
        t3.update_m()
        #t.P[:]=0
        t3.step()
        t3.fusion()
        T3.append((t3.T*(t3.r**2)).sum()/(t3.r**2).sum())
        R3.append(t3.R/params['Rf'])
        print("t: {:.3} My    ".format(t3.t*0.717),end='\r')
    


    #T_arr = np.array([list(zip(*R1))[0],list(zip(*R1))[1],T1,T2,T3])
    T_arr = np.array([t,R1,R2,R3,T1,T2,T3])
    T_arr[-3:] = T_arr[-3:]*t1.T0
    
    print(T_arr[:,100])
    np.save('T_imp_1My.npy',T_arr)

    plt.plot(T_arr[0],T_arr[4])
    plt.plot(T_arr[0],T_arr[5])
    plt.plot(T_arr[0],T_arr[6])
    plt.show()

def load_T_1My():
    global T_arr
    params = {'Ri':5E3,
    'Rf':5E5,
    'beta': 0,
    'ta':1*cst['tau_al']/0.717,#5My
    'Ti':300,
    'dt':1E-3*cst['tau_al']/0.717,
    'size':1000,
    }
    t1 = terre(**params,cst=cst)

    T_arr = np.load('T_imp_1My.npy')

    plt.plot(T_arr[0],T_arr[4])
    plt.plot(T_arr[0],T_arr[5])
    plt.plot(T_arr[0],T_arr[6])
    plt.show()

def calc_T_5My():
    global R1,R2,t1,t2
    params = {'Ri':1E3,
    'Rf':1E5,
    'beta': 0,
    'ta':5*cst['tau_al']/0.717,#5My
    'Ti':300,
    'dt':1E-3*cst['tau_al']/0.717,
    'size':1000,
    }
    n = int(params['ta']/params['dt'])
#n=10

    print("courbe 1")
    params['beta'] = 0
    t1 = terre(**params,cst=cst)
    #t1.d1,t1.d2,t1.d3,t1.u1,t1.u2,t1.u3=[0]*6
    T1 = []
    R1 = []
    t = []
    for _ in range(n):
        t1.update_P()
        t1.update_m()
        #t.P[:]=0
        t1.step()
        t1.fusion()
        T1.append((t1.T*(t1.r**2)).sum()/(t1.r**2).sum())
        R1.append(t1.R/params['Rf'])
        t.append(t1.t*0.717)
        print("t: {:.3} My    ".format(t1.t*0.717),end='\r')
    
    print("courbe 2")
    params['beta'] = 1
    t2 = terre(**params,cst=cst)
    T2 = []
    R2 = []
    for _ in range(n):
        t2.update_P()
        t2.update_m()
        t2.step()
        t2.fusion()
        R2.append(t2.R/params['Rf'])
        T2.append((t2.T*(t2.r**2)).sum()/(t2.r**2).sum())
        print("t: {:.3} My    ".format(t2.t*0.717),end='\r')

    print("courbe 3")
    params['beta'] = 2
    t3 = terre(**params,cst=cst)
    T3 = []
    R3 = []
    for _ in range(n):
        t3.update_P()
        t3.update_m()
        #t.P[:]=0
        t3.step()
        t3.fusion()
        T3.append((t3.T*(t3.r**2)).sum()/(t3.r**2).sum())
        R3.append(t3.R/params['Rf'])
        print("t: {:.3} My    ".format(t3.t*0.717),end='\r')
    


    #T_arr = np.array([list(zip(*R1))[0],list(zip(*R1))[1],T1,T2,T3])
    T_arr = np.array([t,R1,R2,R3,T1,T2,T3])
    T_arr[-3:] = T_arr[-3:]*t1.T0
    
    print(T_arr)
    np.save('T_imp_5My.npy',T_arr)

    plt.plot(T_arr[0],T_arr[4])
    plt.plot(T_arr[0],T_arr[5])
    plt.plot(T_arr[0],T_arr[6])
    plt.show()

def load_T_5My():
    global T_arr
    params = {'Ri':1E3,
    'Rf':1E5,
    'beta': 0,
    'ta':1*cst['tau_al']/0.717,#5My
    'Ti':300,
    'dt':1E-3*cst['tau_al']/0.717,
    'size':1000,
    }
    t1 = terre(**params,cst=cst)

    T_arr = np.load('T_imp_5My.npy')

    plt.plot(T_arr[0],T_arr[4])
    plt.plot(T_arr[0],T_arr[5])
    plt.plot(T_arr[0],T_arr[6])
    plt.show()

def calc_R():
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



    R_arr = np.array([list(zip(*R1))[0],list(zip(*R1))[1],list(zip(*R2))[1],list(zip(*R3))[1]])
    print(R_arr)
    np.save('R_imp.npy',R_arr)



#calc_R()
calc_T_1My()
load_T_1My()
calc_T_5My()
load_T_5My()
