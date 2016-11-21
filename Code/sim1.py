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
    def __init__(self,**kwargs):
        #N nombre initial de particules
        
        
        R, r, dr, dt, lmbda, rho, Cp, P, sigma, T0 = [1]*10
        
        try:
            taille = kwargs['taille']
        except:
            taille = 10
            
        try:
            dt = kwargs['dt']
        except:
            dt = 1
        
        try:
            dr = kwargs['dr']
        except:
            dr = 1
            
        try:
            lmbda = kwargs['lmbda']
        except:
            lmbda = 1
            
        try:
            P = kwargs['P']
        except:
            lmbda = np.zeros(taille)        
            
        try:
            rho = kwargs['rho']
        except:
            rho = 1    
        
        #l'equation n'est pas definie en r=0 on decale donc le premier point
        epsilon = 1/taille
        # on fait commencer r après 0 pour l'instant
        r = np.linspace(epsilon, 1, num=taille)
        
        dt = 1
        
        
        dr=1/taille
        lmbda = 1
        sigma = 5.670373E-8 
        T0 = 2.7 #2.7K fond diffus cosmologique
        c = 1
        
        
        
        self.taille   = taille
        self.epsilon  = epsilon  #l'equation n'est pas definie en r=0
        self.R        = R    #rayon de la terre
        self.r        = r    #rayon relatif [0,1], il intervient dans la formule du fait de la sym spherique
        self.dr       = dr   #pas d'espace
        self.dt       = dt    #pas de temps
        self.lmbda    = lmbda #conductivité thermique
        self.rho      = rho  #masse volumique
        self.Cp       = Cp   #capaité calorifique
        self.P        = P    #terme de source (nucleaire)
        self.c        = c
        self.sigma    = sigma
        self.T0       = T0   #temperature du fond diffus cosmologique
        
        #condition initiales
        self.T     = np.zeros(taille)
        #self.T[:] = 300
        #self.rho0  = np.zeros(taille)
        #self.P = np.zeros(taille)
        #self.P[:20] = 1
        
        
        #quelques matrices de base pour constuire la matrice complete
        ones = np.ones(self.taille)
        self.ones = ones
        eye = sparse.eye(taille)
        d1 = sparse.dia_matrix(([1*ones,-1*ones],[0,1]), shape=(taille, taille))
        c1 = ((r+dr/2)**2)/((r*dr)**2)
        u1 = sparse.dia_matrix((c1,[0]), shape=(taille, taille))
        
        d2 = sparse.dia_matrix(([1*ones,-1*ones],[0,-1]), shape=(taille, taille))
        c2 = ((r-dr/2)**2)/((r*dr)**2) 
        # on met c2 dans une diagonale pour la multiplication avec d2
        u2 = sparse.dia_matrix((c2,[0]), shape=(taille, taille))
        #print('u \n',u2.toarray())
        #print('u \n',(u2*d2).toarray())

        #condition aux limites  
        
        d1 = d1.tocsr()
        d1[0,0:3]  = [2,-1,-1]
        d1[-1,-3:] = 0

        
        d2 = d2.tocsr()
        d2[0,0:2]  = 0
        d2[-1,-3:] = [-1,-1,2]
        
        
        #print(d1.toarray())
        #print(d2.toarray())
        #on construit m, la matrice complete 
        c = self.dt*self.lmbda/(self.rho*self.Cp*self.R**2)
        self.m = eye + c*u1*d1 + c*u2*d2
        
        #print(self.m.toarray())
        
        self.m = sparse.csc_matrix(self.m) #format csc plus efficace
        self.lu = linalg.splu(self.m) #on précalcule la décomposition LU
        
        
    def update_P(self,t,T_neb):
        #On part de la chaleur radioactice
        self.P = self.rho*self.H0*np.exp(-np.log(2)*t/self.t_half)
        #On complete avec la radiation de corps noir à la surface
        self.P[-1] -= (self.R**2)*np.pi*4*self.sigma(self.T[-1]**4-T_neb**4)

    def step(self):
        self.T = self.lu.solve(self.T + self.c*self.P)

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
        %matplotlib qt4
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
    t = terre(taille=4)
    t.step()


    print("Comparaison solution analytique stationnaire et simulée (la simulation converge en 2 pas ici)")

    r=np.linspace(1E-3,1,1000)
    R=1
    D=1
    P=1

    s = sol_ref2(r,R,D,P)
    plt.plot(r,s,lw=2,label='analytique')

    P=np.zeros(1000)
    P[0] = 1
    P[-1] = -1
    taille=1000
    t = terre(taille=taille,dt=1,R=1,D=1,P=P,rho=1)
    t.T[:] = 0
    for _ in range(1,3):
	t.step()
	plt.plot(r,1000*(t.T[:] - t.T[-1]),label='simulé pas {}'.format(_))
	

    plt.legend()
    plt.show()

def sol_ref1(r,t,D):
    """
    Solution analytique pour un problème sans bords avec un dirac placé au centre
    """
    return np.exp(-(r**2)/(4*D*t))/np.sqrt(4*np.pi*t*D)**3

def sol_ref2(r,R,D,P):
    """
    Limite à t->infini avec source P au centre et T=0 fixé au bords
    """
    return (R**2 - r**2)*P/(6*D)


