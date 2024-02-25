# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 12:39:30 2024

@author: VGenot

Solves RH relations with theta_Bn
From "Effect of Oblique Interplanetary Magnetic Field on Shape and Behavior of the Magnetosphere", G. Walters, JGR 1964
 https://doi.org/10.1029/JZ069i009p01769

The aim is to verify that n_Qpara>n_Qperp just behind the shock (for classical Parker IMF) as shown in the Figure 3 (3.76 vs 3.15 for instance)

"""
import math 
import matplotlib.pyplot as plt
import numpy as np

class MonComplexe (): 
    """ modélisation des nombres complexes """ 
    def __init__(self,a=0.0,b=0.0): 
        """le constructeur soit en cartésiennes soit en polaires""" 
        self.x=a 
        self.y=b 
 
    def __str__(self): 
        """ représentation externe pour print et str """ 
        if self.y==0: 
            return str(self.x) 
        if self.x==0 and self.y==1: 
            return "i" 
        if self.x==0 and self.y==-1: 
            return "-i" 
        if self.x==0: 
            return str(self.y)+"i" 
        if self.y==1: 
            return str(self.x)+"+i" 
        if self.y==-1: 
            return str(self.x)+"-i" 
        if self.y>0: 
            return str(self.x)+"+"+str(self.y)+"i" 
        else: 
            return str(self.x)+str(self.y)+"i" 
 
    def __add__(self,other): 
        """somme de deux complexes""" 
        a=self.x+other.x 
        b=self.y+other.y 
        return MonComplexe(a,b) 
 
    def __sub__(self,other): 
        """différence de deux complexes""" 
        a=self.x-other.x 
        b=self.y-other.y 
        return MonComplexe(a,b) 
 
    def __neg__(self): 
        """opposé d'un complexe""" 
        return MonComplexe(-self.x,-self.y) 
 
    def null(self): 
        """test de nullité""" 
        return self.x==0 and self.y==0 
 
    def __mul__(self,other): 
        """produit de deux complexes""" 
        a=self.x*other.x-self.y*other.y 
        b=self.x*other.y+self.y*other.x 
        return MonComplexe(a,b) 
 
    def __truediv__(self,other): 
        """quotient de deux complexes""" 
        return self*(~other) 
 
    def conj(self): 
         """conjugué d'un complexe""" 
         return MonComplexe(self.x,-self.y) 
 
    def module(self): 
        """module d'un complexe""" 
        return math.sqrt(self.x*self.x+self.y*self.y) 
    
    def argument(self): 
        """argument d'un complexe""" 
        if self.x==0 and self.y==0: 
            return 0 
        if self.x==0 and self.y>0: 
            return math.pi/2 
        if self.x==0 and self.y<0: 
            return -math.pi/2 
        if self.x >0: 
            return math.atan(self.y/self.x) 
        return math.pi-math.atan(self.y/(-self.x)) 
 
    def __invert__(self): 
        """inverse d'un complexe""" 
        if self.null(): 
            raise ZeroDivisionError 
        return MonComplexe(self.x/(self.x*self.x+self.y*self.y),-self.y/(self.x*self.x+self.y*self.y)) 
 
    def __pow__(self,n): 
        """puissances d'un complexe""" 
        if n<=0 and self.null(): 
            raise ZeroDivisionError 
        if n >=0: 
            R=MonComplexe(1,0) 
            while n: 
                R=R*self 
                n=n-1 
            return R 
        if n<0: 
            return (~R)**(-n) 
 
    def racines(self,n): 
        """calcule les n racines n-ièmes du nombre""" 
        # on utilise les racines de l'unité 
        return [MonComplexe(self.module()**(1.0/n) *math.cos((k*2*math.pi+self.argument())/n),self.module()**(1.0/n)*math.sin((k*2*math.pi+self.argument())/n)) for k in range(0,n) ] 
 
def Cardan (a,b,c,d): 
    """ a, b, c, d sont les coefficients initiaux de l'équation""" 
    # on commence par mettre sous forme canonique 
    b,c,d=b/a,c/a,d/a 
    p=c-b*b/MonComplexe(3.0,0) 
    q=d-b*c/MonComplexe(3.0,0)-(b**3)/MonComplexe(27.0,0)+(b**3)/MonComplexe(9.0,0) 
    B,C= q,-p*p*p/MonComplexe(27.0,0) 
    D=B*B-MonComplexe(4.0,0)*C 
    R=D.racines(2) 
    U=(-B+R[0])/MonComplexe(2.0,0) 
    roots=U.racines(3) 
    sol1=[u-p/(MonComplexe(3.0,0)*u) for u in roots] 
    sol2=[z-b/MonComplexe(3.0,0) for z in sol1] 
    return sol2 
 
def Resolution (a,b,c,d): 
    """Résout l'équation az^3+bz^2+c^z+d=0""" 
    # les coefficients peuvent être entiers, réels ou complexes 
    # Dans tous les cas on convertit en complexes pour commencer 
    if isinstance(a,float) or isinstance(a,int): 
        a=MonComplexe(float(a),0) 
    if isinstance(b,float) or isinstance (b,int): 
        b=MonComplexe(float(b),0) 
    if isinstance(c,float) or isinstance(c,int): 
        c=MonComplexe(float(c),0) 
    if isinstance(d,float) or isinstance(d,int): 
        d=MonComplexe(float(d),0) 
    Z= Cardan(a,b,c,d) 
    print("Racines") 
    print(Z[0]) 
    print(Z[1]) 
    print(Z[2]) 
    P0=a*Z[0]**3+b*(Z[0]**2)+c*Z[0]+d 
    P1=a*Z[1]**3+b*(Z[1]**2)+c*Z[1]+d 
    P2=a*Z[2]**3+b*(Z[2]**2)+c*Z[2]+d 
    #vérification 
    #print("Vérification") 
    #print(P0) 
    #print(P1) 
    #print(P2) 
 
def Resolution2 (a,b,c,d): 
    """Résout l'équation az^3+bz^2+c^z+d=0""" 
    # les coefficients peuvent être entiers, réels ou complexes 
    # Dans tous les cas on convertit en complexes pour commencer 
    if isinstance(a,float) or isinstance(a,int): 
        a=MonComplexe(float(a),0) 
    if isinstance(b,float) or isinstance (b,int): 
        b=MonComplexe(float(b),0) 
    if isinstance(c,float) or isinstance(c,int): 
        c=MonComplexe(float(c),0) 
    if isinstance(d,float) or isinstance(d,int): 
        d=MonComplexe(float(d),0) 
    Z= Cardan(a,b,c,d)    
    return Z
  
def main(): 
    
    beta=40. # 4pi*rho1*u1**2/H1**2 : ratio of particle-directed flow energy to magnetic energy density in interplanetary medium (trad : Mach Alfvén**2)
    psi=25. # deg : angle between IMF and solar wind
    theta1=90.1 # deg : angle between IMF and normal to the shock front
    theta1Ar=180-np.arange(100)/99*89.9
    A=np.sin(1.5*np.pi-(psi+theta1)*np.pi/180.)**2
    co1=np.cos(theta1*np.pi/180.)
    si1=np.sin(theta1*np.pi/180.)
    Q=2. # H1**2/(8pi*P1) : ratio of magnetic pressure to particle pressure in interplanetary medium (trad: inverse classical beta)
    gamma=5./3.
    S1car=(np.tan(theta1*np.pi/180.))**2
    soluAr=[0]*100

    for i in range(100):

        theta1=theta1Ar[i]
        A=np.sin(1.5*np.pi-(psi+theta1)*np.pi/180.)**2
        co1=np.cos(theta1*np.pi/180.)
        si1=np.sin(theta1*np.pi/180.)        
        S1car=(np.tan(theta1*np.pi/180.))**2

        t3=co1**4+(gamma-1)/gamma*beta*A*Q*(1+S1car)*co1**4
        t2=-2*beta*A*(co1**2+Q*co1**4+Q*co1**2*si1**2)+beta**2*A**2*Q*(si1**2-2*(gamma-1)/gamma*(1+S1car)*co1**2)+(gamma-1)/gamma*beta*A*Q*(1+S1car)*co1**4
        t1=beta**2*A**2*(1+4*Q*co1**2-2*(gamma-1)/gamma*Q*co1**2+Q*si1**2)+(gamma-1)/gamma*beta**3*A**3*Q
        cst=-beta**3*A**3*Q*(gamma+1)/gamma

        Resolution(t3,t2,t1,cst)
    
        solu=Resolution2(t3,t2,t1,cst)
        soluAr[i]=min(MonComplexe.module(solu[0]),MonComplexe.module(solu[1]),MonComplexe.module(solu[2]))
    
        print('alpha1 = ',psi+theta1-180)
    
    print(soluAr)
    plt.plot(180-theta1Ar,soluAr)
    plt.ylabel('Compression ratio')
    plt.xlabel('theta_Bn')
    plt.title('Ma**2=40, beta=0.5, psi=25°')
    
if __name__ == '__main__': 
    main() 