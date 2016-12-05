import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse  import linalg




class terre:
    """
    On simule une Terre a symmetrie sphérique
    On considère son rayon constant
    On considère que le chauffage provient uniquement de la radioactivité
    de l'Al26 
    """
    def __init__(self, Ri, Ti, dt, size, cst=None):
        """
        R0 : rayon initial de la Terre en m
        T0 : temperature initiale en K
        dt : pas de temps en fraction de tau_al
        size : nombre de pas d'espace
        """

        if cst is None:
            cst = { 'tau_al':  2.26E13,#s
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
                    }
        #cst contient les constantes physiques des materiaux
        self.cst = cst
        #ces constantes permettent de définir les 
        #paramètres d'adimensionnement
        self.t0 = cst['tau_al'] 
        self.r0 = np.sqrt(cst['kT']*cst['tau_al']
                /(cst['Cp']*cst['rho']))
        self.T0 = cst['T_neb'] 
        self.cst['Tafus_m'] = self.cst['Tfus_m']/self.T0
        self.cst['Tafus_s'] = self.cst['Tfus_s']/self.T0

        self.t = 0
        self.dt = dt #!!Le param dt doit etre en fraction de tau_al!!
        self.c0 = self.dt*cst['tau_al']/(cst['rho']*cst['Cp']*cst['T_neb'])

        print("rayon : {} km, dt : {} My, Ti : {} K, size {}".format(Ri*1E-3,dt*0.717,Ti,size))
        
        #l'equation n'est pas definie en r=0 on decale donc le premier 
        # point de dr
        dr = Ri/(self.r0*size)
        r = np.linspace(dr, Ri/self.r0, num=size)
        self.r = r
        self.dr = dr
        self.size = size
        
        
        #condition initiales
        self.T = (Ti/self.T0)*np.ones(size)
        self.P = np.zeros(size)
        self.phi_m = np.zeros(size)
        self.phi_s = np.zeros(size)
        
        
        #quelques matrices de base pour constuire la matrice complete
        ones = np.ones(self.size)
        self.ones = ones
        eye = sparse.eye(size)
        d1 = sparse.dia_matrix(([1*ones,-1*ones],[0,1]),
                shape=(size, size))
        c1 = self.dt*((r+dr/2)/(r*dr))**2 
        # on met c1 dans une diagonale pour la multiplication avec d1
        u1 = sparse.dia_matrix((c1,[0]), shape=(size, size))
        
        d2 = sparse.dia_matrix(([1*ones,-1*ones],[0,-1]),
                shape=(size, size))
        c2 = self.dt*((r-dr/2)/(r*dr))**2 
        u2 = sparse.dia_matrix((c2,[0]), shape=(size, size))

        #condition aux limites  
        d1 = sparse.lil_matrix(d1)
        #d1 = d1.tocsr()
        d1[0,0:3]  = [2,-1,-1]
        d1[-1:,-4:] = 0
        d2 = sparse.lil_matrix(d2)
        #d2 = d2.tocsr()
        d2[0,0:2]  = 0
        d2[-1,-3:] = [0,-2,2]
        #print(d1.toarray())
        #print(d2.toarray())

        self.c1 = c1
        self.c2 = c2
        self.d1 = d1
        self.d2 = d2

        #on construit m, la matrice complete 
        self.m = eye + u1*d1 + u2*d2

        #print(self.m.toarray())
        
        self.m = sparse.csc_matrix(self.m) #format csc plus efficace
        self.lu = linalg.splu(self.m) #on précalcule la décomposition LU
        
    def update_P(self):
        cst = self.cst
        #On part de la chaleur radioactive
        self.P[:] = cst['rho']*cst['H0']*2**(-(self.t+self.dt/2))
        self.P[-1] = 0
        #On complete avec la radiation de corps noir à la surface
        try :
            self.cS
        except AttributeError:
            self.cS = (((self.r[-1]/(self.r[-1]+self.dr))**2)
            *self.dr*self.r0*cst['sigma']*(cst['T_neb']**4)
            /(self.cst['kT']*self.T0))
            #on utilise un point de temperature pour definir le flux radiatif
        #print("T_rad",self.cS*(1-(self.T[-2]**4)))
        self.T[-1] = self.T[-2] + self.cS*(1-(self.T[-2])**4)

        #print("Cs",self.cS)

        # si la température du point fictif est plus petit que 0
        # la formule de la loi de Stephan n'est plus valable et la valeur 
        # calculée explose, il faut donc borner par zero pour que la 
        # diffusion puisse equilbrer 
        if self.T[-1] < 0 :
            self.T[-1] = 0
        #print(self.T[-2:])

    def step(self):
        self.t += self.dt
        self.T = self.lu.solve(self.T + self.c0*self.P)

    def fusion(self):
        """
        Cette fonction calcule la modification du ratio solide/liquide
        en deux etapes : 
        1. Si la condition pour un changement d'état est verifiée
           on utilise toute la chaleur pour changer d'état
        2. Si on a plus de 100% de liquide ou de solide on rebalance
           l'excedant en chaleur
        """
        cst = self.cst 

        fus_m = np.where((self.T>cst['Tafus_m'])*(self.phi_m < 1))
        sol_m = np.where((self.T<cst['Tafus_m'])*(self.phi_m > 0))
        tot_m = [np.concatenate([fus,sol]) for fus,sol in zip(fus_m,sol_m)]

        self.phi_m[tot_m] = self.phi_m[tot_m] + (self.T[tot_m] - cst['Tafus_m']) * self.T0*cst['Cp_m']/(cst['phi']*cst['L_m'])
        self.T[tot_m] = cst['Tafus_m']


        excess_fus = np.where(self.phi_m > 1)
        excess_sol = np.where(self.phi_m < 0)

        self.T[excess_fus] = self.T[excess_fus] + (self.phi_m[excess_fus] - 1) * (cst['phi']*cst['L_m'])/(self.T0*cst['Cp_m'])
        self.phi_m[excess_fus] =  1

        self.T[excess_sol] = self.T[excess_sol] + self.phi_m[excess_sol] * (cst['phi']*cst['L_m'])/(self.T0*cst['Cp_m'])
        self.phi_m[excess_sol] =  0

        fus_s = np.where((self.T>cst['Tafus_s'])*(self.phi_s < 1))
        sol_s = np.where((self.T<cst['Tafus_s'])*(self.phi_s > 0))
        tot_s = [np.concatenate([fus,sol]) for fus,sol in zip(fus_s,sol_s)]

        self.phi_s[tot_s] = self.phi_s[tot_s] + (self.T[tot_s] - cst['Tafus_s']) * self.T0*cst['Cp_s']/((1-cst['phi'])*cst['L_s'])
        self.T[tot_s] = cst['Tafus_s']

        excess_fus = np.where(self.phi_s > 1)
        excess_sol = np.where(self.phi_s < 0)

        self.T[excess_fus] = self.T[excess_fus] + (self.phi_s[excess_fus] - 1) * ((1-cst['phi'])*cst['L_s'])/(self.T0*cst['Cp_s'])
        self.phi_s[excess_fus] =  1

        self.T[excess_sol] = self.T[excess_sol] + self.phi_s[excess_sol] * ((1-cst['phi'])*cst['L_s'])/(self.T0*cst['Cp_s'])
        self.phi_s[excess_sol] =  0

    def store(self):
        self.stored.append({
            't':self.t.copy(),
            'T':self.T.copy(),
            'P':self.P.copy(),})
    
    def show(self):
        size = self.T.shape[0]
        y, x = np.mgrid[0:size, 0:size]
        r = np.sqrt(x**2 + y**2).astype(np.uint)
        valid = np.where(r < size)
        
        d = np.empty((size,size))
        d[valid] = self.T[r[valid]]
        
        plt.imshow(np.ma.array(d, mask=(r >= size)),cmap='plasma',interpolation='nearest',vmin=0,vmax=400)
        plt.colorbar()
        
    def anim(self):
        i=0
        #while True:
        size = self.T.shape[0]
        for _ in range(20):
            if i == 0:
                y, x = np.mgrid[0:size, 0:size]
                r = np.sqrt(x**2 + y**2).astype(np.uint)
                valid = np.where(r < size)

                d = np.empty((size,size))
                d = np.ma.array(d, mask=(r >= size))
                d[valid] = self.T[r[valid]]

                p = plt.imshow(d,cmap='plasma')
                plt.clim()
                plt.colorbar()
                i=1
            else:
                i+=1
                self.step()
                self.BC()
                d[valid] = self.T[r[valid]]
                p.set_data(d)
                plt.title("frame {}, sum : {}".format(i,(self.r*self.T).sum()))

            
            plt.pause(0.5)


def test():
    """
    Test de la classe et comparaison des resultats avec une sol 
    analytique connue
    """
    #on teste vite fait que la classe fonctionne 
    t = terre(Ri=100,Ti=0,dt=0.1,size=4)
    t.update_P()
    t.step()
    t.fusion()


    print("Comparaison solution analytique stationnaire et simulée (la simulation converge en 2 pas ici)")

    r=np.linspace(1E-3,1,1000)
    R=1
    D=1
    P=1

    s = sol_ref2(r,R,D,P)
    plt.plot(r,s,lw=2,label='analytique')

    from collections import defaultdict
    cst = defaultdict(lambda : 1) #dictionnaire avec des valeur par defaut à 1
    size=1000
    t = terre(Ri=1,Ti=0,dt=1,size=1000,cst=cst)
    t.P[:] = 0
    t.P[:] = 1
    t.P[-1] = -t.P.sum()+1
    for _ in range(1,3):
        t.step()
        plt.plot(r,(t.T[:] - t.T[-1]),label='simulé pas {}'.format(_))
	

    plt.legend()
    plt.show()

def test_fusion():
    print("test visuel de la fusion")
    t = terre(Ri=20000,Ti=0,dt=0.1,size=100)
    t.T[:] = np.linspace(0,3000,100)/t.T0
    plt.subplot(121)
    plt.title("temperature")
    plt.plot(t.T*t.T0,label="avant")
    t.fusion()
    plt.plot(t.T*t.T0,label='après')
    plt.legend()

    plt.subplot(122)
    plt.title("pourcentage de fusion")
    plt.plot(t.phi_m,label='phi_m')
    plt.plot(t.phi_s,label='phi_s')
    plt.ylim(-0.5,1.5)
    plt.legend()

    plt.show()

def sol_ref1(r,t,D):
    """
    Solution analytique pour un problème sans bords avec un dirac placé au centre
    """
    return np.exp(-(r**2)/(4*D*t))/np.sqrt(4*np.pi*t*D)**3

def sol_ref2(r,R,D,P):
    """
    Limite à t->infini avec terme source P uniforme et T=0 fixé au bords
    """
    return (R**2 - r**2)*P/(6*D)

if __name__ == '__main__' :
    #test()

    #test_fusion()

    cst = { 'tau_al':  2.2611E13, #s
            'H0'    :  1.5E-7, #W kg-1
            'T_neb' :  300   , #K
            'sigma' :  5.67E-8,#W m-2 K-4
            #'sigma' :  0,#W m-2 K-4
            #'sigma' :  0,#W m-2 K-4
            'phi'   :  0.18  , # sans unité
            #'kT'    :  11.48 , #W K-1 m-1
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
    t = terre(Ri=1000,Ti=300,dt=0.001,size=1000,cst=cst)
    #t.t = 1/0.717
    plt.figure(figsize=(6,8))
    for _ in range(10):
        for __ in range(1000):
            t.update_P()
            #t.P[:]=0
            t.step()
            t.fusion()
            print("t: {} My".format(t.t*0.717))
            """
            if t.t*0.74 > 1.0 :
                break
        if t.t*0.74 > 1.0 :
            break
            """

        plt.plot(t.T[:-1]*t.T0,t.r[:-1]/t.r[-2])
        #plt.plot(t.phi_m,'-')

    plt.xlabel('Température en K')
    plt.ylabel('Rayon adimentioné')
    plt.show()

