#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sim1_v2
import sim2
import numpy as np
np.set_printoptions(precision=3,linewidth=150)




print(" modèle 1 ".center(40,"#"))
t1 = sim1_v2.terre(Ri=50000,Ti=300,dt=0.1,size=10)#dt en pourcentage de tau_al
#print("cst",t1.cst)
print("dt : {} ua, dr : {} ua".format(t1.dt,t1.dr))
print("c1 : {},\nc2 : {},\nc0 : {} ".format(t1.c1,t1.c2,t1.c0))
print("m : \n",t1.m.toarray())
print()


print(" modèle 2 sans evolution".center(40,"#"))
t2 = sim2.terre(Ri=50000,Rf=100000,Ti=300,dt=0.1*2.26E13,ta=10*2.26E13,size=10,beta=-1)
#print("cst",t2.cst)
print("dt : {} ua, dr : {} ua".format(t2.dt,t2.dr))
#print("c1 {},\n c2 {},\n c0 {},  ".format(t2.c1,t2.c2,t2.c0))
print("c1/R² : {},\nc2/R² : {},\nc0 : {},  ".format(t2.c1/t2.R**2,t2.c2/t2.R**2,t2.c0))
print("m : \n",t2.m.toarray())

print(" modèle 2 avec évolution".center(40,"#"))
t2 = sim2.terre(Ri=50000,Rf=50001,Ti=300,dt=0.1*2.26E13,ta=10*2.26E13,size=10,beta=0)
#print("cst",t2.cst)
print("dt : {} ua, dr : {} ua".format(t2.dt,t2.dr))
#print("c1 {},\n c2 {},\n c0 {},  ".format(t2.c1,t2.c2,t2.c0))
print("c1/R² : {},\nc2/R² : {},\nc0 : {},  ".format(t2.c1/t2.R**2,t2.c2/t2.R**2,t2.c0))
print("m : \n",t2.m.toarray())

print(t2.d3.toarray())

