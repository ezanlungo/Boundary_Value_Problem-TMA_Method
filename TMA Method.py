# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 16:19:21 2020

@author: EZEQUIEL
"""

import matplotlib.pyplot as plt
import numpy as np

def pul_cm (p):
    conv = p*2.54
    return conv

a = 5 #radio interno
b = 8 #radio externo

u_a = 0.00387
u_b = 0.00308
N = 10000
h = (b-a)/N
xa = np.linspace(a,b,N-1)
E = 30*(10**6)
v = np.array([0.3,0.5])

def low (h,r):
    return 1-(h/(2*r))

def diag (h,r):
    return -2-h**2/(r**2)

def upp (h,r):
    return 1+(h/(2*r))

A = [low(h,a+i*h) for i in range(N-2)]
B = [diag(h, a+i*h) for i in range(N-1)]
C = [upp(h, a+i*h) for i in range (N-2)]

D = np.zeros(N-1)
D[0] = -A[0]*u_a
D[N-2] = -C[N-3]*u_b

def thomas (a,b,c,d):
    n=len(d)
    w=np.zeros(n-1,float)
    g=np.zeros(n,float)
    p=np.zeros(n,float)
    
    w[0]=c[0]/b[0]
    g[0]=d[0]/b[0]
    
    for i in range(1,n-1):
        w[i]=c[i]/(b[i]-a[i-1]*w[i-1])
        
    for i in range(1,n):
        g[i]=(d[i]-a[i-1]*g[i-1])/(b[i]-a[i-1]*w[i-1])
    
    p[n-1]=g[n-1]
    
    for i in range(n-1,0,-1):
        p[i-1]=g[i-1]-w[i-1]*p[i]
    
    return p

u_thomas=thomas(A,B,C,D)

c1=1.35*10**-4
c2=0.016
u_analitico=np.array([])
for i in range (len(xa)):
    u_an = c1*xa[i]+c2/xa[i]
    u_analitico=np.append(u_analitico,u_an)    

plt.figure(figsize=(8,5))
plt.plot(xa,u_thomas, 'k', label='Thomas')
plt.plot(xa,u_analitico, 'b', label='Analítico')
plt.xlabel(r'$r_i$', fontsize = 15)
plt.ylabel(r'$u(r_k)$', fontsize = 15)
plt.legend(loc='best')
plt.grid(True)
plt.show()

#Punto 2
#SIGMA MAX = E/(1-v^2)*(u/r(a)+v*du/dr(a))

sigma_a=np.array([])
sigma_b=np.array([])
for i in range (len(u_thomas)-2):
    aux=(-u_thomas[i]+u_thomas[i+2])/(2*h)
    uprima=aux
    aux1=(E/(1-v[0]**2))*((u_thomas[i+1])/xa[i+1]+v[0]*uprima)
    aux2=(E/(1-v[1]**2))*((u_thomas[i+1])/xa[i+1]+v[1]*uprima)
    sigma_a=np.append(sigma_a,aux1)
    sigma_b=np.append(sigma_b,aux2)

sigma_analit_a=np.array([])
sigma_analit_b=np.array([])
for i in range (len(u_thomas)-2):
    s1 = (E/(1-v[0]**2))*((c1*xa[i+1]+c2/xa[i+1])/xa[i+1]+v[0]*(c1-c2*(xa[i+1]**-2)))
    s2 = (E/(1-v[1]**2))*((c1*xa[i+1]+c2/xa[i+1])/xa[i+1]+v[1]*(c1-c2*(xa[i+1]**-2)))
    sigma_analit_a=np.append(sigma_analit_a,s1)
    sigma_analit_b=np.append(sigma_analit_b,s2)

xa_sigma = xa.tolist()
xa_sigma.pop(0)
xa_sigma.pop(-1)

xa_s=xa[1:-1]

plt.figure(figsize=(8,5))
plt.plot(xa_s,sigma_a, 'm-', label='Sigma Thomas para v = 0.3')
plt.plot(xa_s,sigma_b, 'g-', label='Sigma Thomas para v = 0.5')
plt.xlabel(r'$r_i$', fontsize = 15)
plt.ylabel(r'$sigma(r_k)$', fontsize = 15)
plt.plot(xa_s,sigma_analit_a, 'm--', label='Sigma analítico para v = 0.3')
plt.plot(xa_s,sigma_analit_b, 'g--', label='Sigma analítico para v = 0.5')
plt.legend(loc='best')
plt.grid(True)
plt.show()

#Punto 3

e_uno=np.array([])
e_dos=np.array([])
for i in range (len(xa_s)):
    Error_uno = np.abs((sigma_analit_a[i]-sigma_a[i])/sigma_analit_a[i])
    e_uno=np.append(e_uno,Error_uno)
    Error_dos = np.abs((sigma_analit_b[i]-sigma_b[i])/sigma_analit_b[i])
    e_dos=np.append(e_dos,Error_dos)
    
plt.figure(figsize=(8,5))
plt.plot(xa_s,e_uno*100, 'm-', label='Error para Poisson 0.3 por Thomas')
plt.plot(xa_s,e_dos*100, 'g-', label='Error para Poisson 0.5 por Thomas')
plt.xlabel(r'$r_i$', fontsize = 15)
plt.ylabel(r'$Error Porcentual$', fontsize = 15)
plt.legend(loc='best')
plt.grid(True)
plt.show()

#Punto 4
#Cambio de variables para pasar a sistema de orden 1
#Z1 = u(r)
#Z2 = Z1' = u'(r)
#Z2' = u''(r)

#EXPRESIÓN : Z2'= -Z2/r + Z1/r2

#Defino las funciones necesarias

coef_a=np.array([u_a,0])
coef_arbitrario=np.array([0,1])

def z1prima(z2):
    return z2

def z2prima(z2,z1,r):
    return (-z2/r)+(z1/r**2)

def rk_solucion(coefa, coefarbitrario, coefa2, coefarbitrario2):
    sol1=np.array([coefa])
    sol2=np.array([coefarbitrario])
    sol21=np.array([coefa2])
    sol22=np.array([coefarbitrario2])
    for i in range(len(xa)):
        k1_1=h*z1prima(coefarbitrario)
        k1_2=h*z2prima(coefarbitrario, coefa, xa[i])
        k2_1=h*z1prima(coefarbitrario+k1_2/2)
        k2_2=h*z2prima(coefarbitrario+k1_2/2, coefa+k1_1/2 ,xa[i]+h/2)
        k3_1=h*z1prima(coefarbitrario+k2_2/2)
        k3_2=h*z2prima(coefarbitrario+k2_2/2, coefa+k2_1/2, xa[i]+h/2)
        k4_1=h*z1prima(coefarbitrario+k3_2)
        k4_2=h*z2prima(coefarbitrario+k3_2, coefa+k3_1, xa[i]+h)
        
        coef1=coefa+(k1_1+2*k2_1+2*k3_1+k4_1)/6
        coef2=coefarbitrario+(k1_2+2*k2_2+2*k3_2+k4_2)/6
        
        coefa=coef1
        coefarbitrario=coef2
        
        k1_1p=h*z1prima(coefarbitrario2)
        k1_2p=h*z2prima(coefarbitrario2, coefa2, xa[i])
        k2_1p=h*z1prima(coefarbitrario2+k1_2p/2)
        k2_2p=h*z2prima(coefarbitrario2+k1_2p/2, coefa2+k1_1p/2 ,xa[i]+h/2)
        k3_1p=h*z1prima(coefarbitrario2+k2_2p/2)
        k3_2p=h*z2prima(coefarbitrario2+k2_2p/2, coefa2+k2_1p/2, xa[i]+h/2)
        k4_1p=h*z1prima(coefarbitrario2+k3_2p)
        k4_2p=h*z2prima(coefarbitrario2+k3_2p, coefa2+k3_1p, xa[i]+h)
        
        coef21=coefa2+(k1_1p+2*k2_1p+2*k3_1p+k4_1p)/6
        coef22=coefarbitrario2+(k1_2p+2*k2_2p+2*k3_2p+k4_2p)/6
        
        coefa2=coef21
        coefarbitrario2=coef22
        
        sol1=np.append(sol1,coef1)
        sol2=np.append(sol2,coef2)
        
        sol21=np.append(sol21,coef21)
        sol22=np.append(sol22,coef22)        
        
    return sol1,sol21

pvi1=rk_solucion(coef_a[0],coef_arbitrario[0],coef_a[1],coef_arbitrario[1])[0]
pvi2=rk_solucion(coef_a[0],coef_arbitrario[0],coef_a[1],coef_arbitrario[1])[1]  
    
w1=u_a
w2=(u_b-pvi1[-1])/(pvi2[-1])

u_disparo=np.array([])
for i in range(len(xa)):
    aux=pvi1[i]+w2*pvi2[i]
    u_disparo=np.append(u_disparo,aux)
    
plt.figure(figsize=(8,5))
plt.grid(True)
plt.plot(xa,u_analitico, 'k',label='Analítico')
plt.plot(xa,u_disparo, 'r',label='Disparo')
plt.xlabel(r'$r_i$', fontsize = 15)
plt.ylabel(r'$u(r_k)$', fontsize = 15)
plt.legend(loc='best')
plt.grid(True)
plt.show()

#Para obtener los sigma se procede a repetir el algoritmo del punto 2 pero para el vector
#De Mus obtenidos mediante el método del disparo

#Punto 2 Disparo

sigma_a_disparo=np.array([])
sigma_b_disparo=np.array([])
for i in range (len(u_thomas)-2):
    aux=(-u_disparo[i]+u_disparo[i+2])/(2*h)
    uprima=aux
    aux1=(E/(1-v[0]**2))*((u_disparo[i+1])/xa[i+1]+v[0]*uprima)
    aux2=(E/(1-v[1]**2))*((u_disparo[i+1])/xa[i+1]+v[1]*uprima)
    sigma_a_disparo=np.append(sigma_a_disparo,aux1)
    sigma_b_disparo=np.append(sigma_b_disparo,aux2)

xa_sigma = xa.tolist()
xa_sigma.pop(0)
xa_sigma.pop(-1)

xa_s=xa[1:-1]

plt.figure(figsize=(8,5))
plt.plot(xa_s,sigma_a_disparo, 'y-', label='Sigma disparo para v = 0.3')
plt.plot(xa_s,sigma_b_disparo, 'c-', label='Sigma disparo para v = 0.5')
plt.xlabel(r'$r_i$', fontsize = 15)
plt.ylabel(r'$sigma(r_k)$', fontsize = 15)
plt.plot(xa_s,sigma_analit_a, 'y--', label='Sigma analítico para v = 0.3')
plt.plot(xa_s,sigma_analit_b, 'c--', label='Sigma analítico para v = 0.5')
plt.legend(loc='best')
plt.grid(True)
plt.show()


#Punto 3 Disparo

e_uno_disparo=np.array([])
e_dos_disparo=np.array([])
for i in range (len(xa_s)):
    Error_uno = np.abs((sigma_analit_a[i]-sigma_a_disparo[i])/sigma_analit_a[i])
    e_uno_disparo=np.append(e_uno_disparo,Error_uno)
    Error_dos = np.abs((sigma_analit_b[i]-sigma_b_disparo[i])/sigma_analit_b[i])
    e_dos_disparo=np.append(e_dos_disparo,Error_dos)
        
plt.figure(figsize=(8,5))
plt.plot(xa_s,e_uno_disparo*100, 'y-', label='Error para Poisson 0.3 por método del Disparo')
plt.plot(xa_s,e_dos_disparo*100, 'c-', label='Error para Poisson 0.5 por método del Disparo')
plt.xlabel(r'$r_i$', fontsize = 15)
plt.ylabel(r'$Error Porcentual$', fontsize = 15)
plt.legend(loc='best')
plt.grid(True)
plt.show()

print('\nError máximo para Esfuerzo Nominal Máx para Poisson 0.3 (Thomas):', np.round(np.max(e_uno)*100,3), '%')
print('Error máximo para Esfuerzo Nominal Máx para Poisson 0.5 (Thomas):', np.round(np.max(e_dos)*100,3), '%')
print('Error máximo para Esfuerzo Nominal Máx para Poisson 0.3 (Disparo):', np.round(np.max(e_uno_disparo)*100,3), '%')
print('Error máximo para Esfuerzo Nominal Máx para Poisson 0.5 (Disparo):', np.round(np.max(e_dos_disparo)*100,3), '%')
