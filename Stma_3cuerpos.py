# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 19:52:42 2018

@author: david
"""

"""
Problema de los 3 cuerpos sin relatividad general

Aviso, todas las unidades están en unidades astronómicas
1AU=1.5e11m; 1año=3.2e7s
En estas unidades, GMs=4*pi**2
"""

#============================================================================== 
#==============================================================================
from numpy import zeros, pi, sqrt
from vpython import canvas, vector, sphere, arrow
from vpython import rate
from timeit import default_timer as timer
from pylab import figure, show, plot

#============================================================================== 
#==============================================================================
start=timer()
Masaestrella=5.7e30                                       #Kg
#Primer planeta
Masaplaneta1=5.2e30                                    #Kg
Semieje_mayor1=-2.                                    #Unidades astronómicas
Excentricidad1=0.056
#Segundo planeta
Masaplaneta2=5.1e30                                     #Kg
Semieje_mayor2=2.                                     #Unidades astronómicas
Excentricidad2=0.056
#Otros datos
Periodoorbital=1                                        #años
Numerodeperiodos=100

#==============================================================================
#==============================================================================
#inicializamos constantes
eps=1e-8
#Masa del sol y masa de los planetas (Kg)
Ms=Masaestrella
Mp=Masaplaneta1
Mp2=Masaplaneta2

#G*Ms y G*Mpl2 (unidades astronómicas)
Gms=4*pi**2
Gmp=Gms*(Mp/Ms)
Gmp2=Gms*(Mp2/Ms)
#Periodo y diferencia finita de tiempos (años)
T=Periodoorbital #en este caso de la tierra
dt=0.00001
#Tiempo inicial y final
t0=0.
tf=T*Numerodeperiodos
#B es el exponente de 1/r más 1 (porque multiplicas y divides por r para meter
#x e y)
B=2
#Número de pasos
N=int(tf//dt)
#Semieje mayor de la elipse
a=Semieje_mayor1
a2=Semieje_mayor2
#Excentricidad
e=Excentricidad1
e2=Excentricidad2
#Velocidad máxima en el afelio
vmin=sqrt(Gms*(1-e)/((1+e)*abs(a)))
vmin2=sqrt(Gms*(1-e2)/((1+e2)*abs(a2)))
#==============================================================================
#==============================================================================
#Condiciones iniciales (suponemos siempre en el afelio)
#Posición inicial en x e y en Ua (la conversión está arriba)
#Posición en x
x0=0
x10=(1+e)*a
x20=(1+e2)*a2
#Posición en y
y0=0
y10=0
y20=0
#Posición en z
z0=0
z10=0
z20=0
#Velocidad en x
vx0=0
vx10=0
vx20=0
#Velocidad en y
vy0=0
vy10=vmin
vy20=-vmin
#Velocidad en z
vz0=0
vz10=0
vz20=0
#Creamos la matriz de la posición de los planetas y el sol
R=zeros((N,3),float)
V=zeros((N,3),float)
R1=zeros((N,3),float)
V1=zeros((N,3),float)
R2=zeros((N,3),float)
V2=zeros((N,3),float)
t=zeros (N,float)

#==============================================================================
#==============================================================================
#Inicializamos posiciones y velocidades
R[0,0]=x0
R[0,1]=y0
R[0,2]=z0
V[0,0]=vx0
V[0,1]=vy0
V[0,2]=vz0
R1[0,0]=x10
R1[0,1]=y10
R1[0,2]=z10
V1[0,0]=vx10
V1[0,1]=vy10
R2[0,0]=x20
R2[0,1]=y20
R2[0,2]=z20
V2[0,0]=vx20
V2[0,1]=vy20
t[0]=t0

#Calculamos las posiciones
for i in range(0,N-1):
    #distancias al sol
    r1=sqrt((R[i,0]-R1[i,0])**2+(R[i,1]-R1[i,1])**2+(R[i,2]-R1[i,2])**2)
    r2=sqrt((R[i,0]-R2[i,0])**2+(R[i,1]-R2[i,1])**2+(R[i,2]-R2[i,2])**2)
    #Distancias al planeta 1
    r12=sqrt((R1[i,0]-R2[i,0])**2+(R1[i,1]-R2[i,1])**2+(R1[i,2]-R2[i,2])**2)
   
    V[i+1,:]=V[i,:]-(Gmp*(R[i,:]-R1[i,:])/r1**3)*dt-(Gmp2*(R[i,:]-R2[i,:])/r2**3)*dt
    V1[i+1,:]=V1[i,:]-(Gms*(R1[i,:]-R[i,:])/r1**3)*dt-(Gmp2*(R1[i,:]-R2[i,:])/r12**3)*dt
    V2[i+1,:]=V2[i,:]-(Gms*(R2[i,:]-R[i,:])/r2**3)*dt-(Gmp*(R2[i,:]-R1[i,:])/r12**3)*dt
    R[i+1,:]=R[i,:]+V[i+1,:]*dt
    R1[i+1,:]=R1[i,:]+V1[i+1,:]*dt
    R2[i+1,:]=R2[i,:]+V2[i+1,:]*dt
    t[i+1]=t[i]+dt
    if r1<2*6.957e8/1.495978707e11 or r2<2*6.957e8/1.495978707e11 or r12<2*6.957e8/1.495978707e11:
        print("Las estrellas han colisionado en t=",dt*i," no se continúan los cálculos")
        break
print('Cálculos terminados')
end=timer()
if (end-start)/60<1:    
    print('Los cálculos han tardado', (end-start),'segundos')
else:
    print('Los cálculos han tardado', (end-start)/60,'minutos')
    
#==============================================================================
#==============================================================================

#Creamos el plot (canvas)
canvas(title = "Sistema 3 cuerpos", x = 500, y = 100, width = 900, 
       height = 900, center = vector(0, 0, 0.4), forward = vector(0, 0, -1),
       background = vector(0, 0, 0))
sol = sphere(pos = vector(R[0,0], R[0,1],R[0,2]), color =vector(1,0.65,0), radius = 0.1,
             make_trail=True)
planeta = sphere(pos = vector(R1[0,0], R1[0,1], R1[0,2]), color = vector(1,1,0), radius = 0.01,
                 make_trail=True)
planeta2 = sphere(pos = vector(R2[0,0], R2[0,1], R2[0,2]), color = vector(1,0,0), radius = 0.01,
                 make_trail=True)

#ploteamos
for i in range(N-1):
    #cuanto mayor sea el número de rate, más rápido va
    rate(100000)
    px = R[i,0]
    py = R[i,1]
    pz = R[i,2]
    px1 = R1[i,0]
    py1 = R1[i,1]
    pz1 = R1[i,2]
    px2 = R2[i,0] 
    py2 = R2[i,1]
    pz2 = R2[i,2]
    sol.pos = vector(px,py,pz)
    planeta.pos = vector(px1, py1, pz1)
    planeta2.pos = vector(px2, py2, pz2)

end=timer()
if (end-start)/60<1:    
    print('Vpython ha terminado su parte en', (end-start),'segundos')
else:
    print('Vpython ha terminado su parte en', (end-start)/60,'minutos')


figure('Sistema solar')
plot(R[:,0],R[:,1],'orange')
plot(R1[:,0],R1[:,1],'y')
plot(R2[:,0],R2[:,1],'r')
show()

end=timer()
if (end-start)/60<1:    
    print('El programa ha tardado', (end-start),'segundos')
else:
    print('El programa ha tardado', (end-start)/60,'minutos')
    
print('Terminado')