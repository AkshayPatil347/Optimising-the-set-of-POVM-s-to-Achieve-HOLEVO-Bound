# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 15:07:03 2022

@author: Vishal
"""
import time
import numpy as np
from numpy.linalg import eig
n=3
rho0=np.reshape([1,0,0,0],(2,2))
j=complex(0,1)
print(j)
rho1=np.reshape([1,np.sqrt(2),np.sqrt(2),2],(2,2))
rho1=rho1/3
rho2=np.reshape([1,np.sqrt(2)*np.exp(-2*np.pi*j/3),np.sqrt(2)*np.exp(2*np.pi*j/3),2],(2,2))
rho2=rho2/3
rho3=np.reshape([1,np.sqrt(2)*np.exp(-4*np.pi*j/3),np.sqrt(2)*np.exp(4*np.pi*j/3),2],(2,2))
rho3=rho3/3
rho=[rho0,rho1,rho2,rho3]
er=1
J={}
stime=time.time()
for g in range(0,360,10):
    for m in range(0,360,10):
        #m is phi
        ph=m
        for r in range(181): #n is theta
            th=r
            e0=np.reshape([1+np.cos(np.deg2rad(r)),np.sin(np.deg2rad(r))*(np.cos(np.deg2rad(ph))-j*np.sin(np.deg2rad(ph))),np.sin(np.deg2rad(r))*(np.cos(np.deg2rad(ph))+j*np.sin(np.deg2rad(ph))),1-np.cos(np.deg2rad(r))],(2,2))/2
            e1=np.reshape([1+np.cos(np.deg2rad(r+g)),np.sin(np.deg2rad(r+g))*(np.cos(np.deg2rad(ph))-j*np.sin(np.deg2rad(ph))),np.sin(np.deg2rad(r+g))*(np.cos(np.deg2rad(ph))+j*np.sin(np.deg2rad(ph))),1-np.cos(np.deg2rad(r+g))],(2,2))/2
            e2=np.eye(2)-e0-e1
            w1 = np.linalg.eigvals(e0)
            w2 = np.linalg.eigvals(e1)
            w3 = np.linalg.eigvals(e2)
            print(w1,w2,w3)
            for z in range(2):
                if w1[z].real>=0 and w1[z].imag==0 and w2[z].real>=0 and w2[z].imag==0 and w3[z].real>=0 and w3[z].imag==0:
                    proceed=True
                else:
                    proceed=False
                    #print(m,r)
            
            if proceed:
                #print("yes")
                em=[e0,e1,e2]
                pxy=np.zeros(shape=(4*n,4*n),dtype=complex)
                px=np.zeros(shape=(4,4),dtype=complex)
                py=np.zeros(shape=(n,n),dtype=complex)
                #print(pxy)
                for j in range(len(em)):
                    pxy[j,j]=np.trace(np.matmul(em[j],rho0))/4
                    pxy[j+n,j+n]=np.trace(np.matmul(em[j],rho1))/4
                    pxy[j+2*n,j+2*n]=np.trace(np.matmul(em[j],rho2))/4
                    pxy[j+3*n,j+3*n]=np.trace(np.matmul(em[j],rho3))/4
               
                
                for i in range(4):
                    tr=0j
                    for k in range(n):
                        tr+=np.trace(np.matmul(em[k],rho[i]))
                    px[i,i]=tr/4
                
                for p in range(n):
                    tr1=0j
                    for q in range(4):
                        tr1+=np.trace(np.matmul(em[p],rho[q]))
                    py[p,p]=tr1/4
                
                lpxy=np.log2(pxy)
                lpxy[lpxy==np.NINF]=0j
                lpx=np.log2(px)
                lpx[lpx==np.NINF]=0j
                lpy=np.log2(py)
                lpy[lpy==np.NINF]=0j
                J1=-1*np.trace(np.matmul(px,lpx))-1*np.trace(np.matmul(py,lpy))+np.trace(np.matmul(pxy,lpxy))
                er=1-J1
                J[(g,m,r)]=J1
                a = [J1.real];
                for k in a:
                    if k>=0 :
                        p = True
                    else:
                        p = False
                if p:
                    print(max(a))
                
               
                if er<=1e-1:
                    print(e0,e1,J1,th,ph)
                
                    a = J.real
                    b = a.max()
                    print(b)
                        
                    break;
            else:
                continue;

    

etime=time.time()
print("Exe time: ",etime-stime) 