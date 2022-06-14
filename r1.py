#!/usr/bin/env python
import numpy as np
import random
import sys
#import time
from energytwopoint_LJpotential_spa_inclu_energyk1 import *
from rotatetwo import *
from rejectcollision_endpoint import *
from rejectcollision_midddlepoint import *
from rejectcollision_twopoints import *

import scipy

def r1(posname,n=5e5,N1=2,l=0.4,psize=50,cubesize_int=2):

    Position3D_raw = np.load(posname)['arr_0']

    for i in range(Position3D_raw.shape[0]):
        if Position3D_raw[i, :].sum() == 0:
            Position3D_raw = Position3D_raw[:i, :] - Position3D_raw[:i, :].min(0)
            break

    Position3D_raw += (psize - Position3D_raw.max(0)) // 2

    Position3D = scipy.array(np.round(Position3D_raw), dtype=int)

    Occupation3D = np.zeros((psize, psize, psize), dtype=int)
    for i,v in enumerate(Position3D):
        Occupation3D[v[0], v[1], v[2]] = 1
    pos_min = 9
    pos_max = psize-10

    Rx = lambda a: np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Ry = lambda a: np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    Rz = lambda a: np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])

    newpos1 = np.zeros(Position3D.shape, dtype=int)


    Occupation3Dold=Occupation3D.copy()
    Position3DB=Position3D.copy()
    t = 0
    t_record=0
    density_cal_save=[]
    j_t=1
    j=1
    RR = [[Rx(b * np.pi / 2), Ry(b * np.pi / 2), Rz(b * np.pi / 2)] for b in range(-1, 2) if b!=0]
    N = Position3D.shape[0]
    while t <= n and j<103:

        #density_var_cal_data=density_var_cal(Occupation3D,6.2,[9,9,9],[40,40,40],4096)
        i = random.randint(0, N - 1)
        angle=[-1,1]
        if i == 0:
            c = random.randint(0, 2)
            b_random = random.randint(0,1)
            b=angle[b_random]
            u = Position3D[0, :] - Position3D[1, :]
            newpos1[i, :] = Position3D[i+1, :] + np.round(np.dot(RR[b_random][c], u))
            if Occupation3D[newpos1[i, 0], newpos1[i, 1], newpos1[i, 2]] >1e-6: continue
            if newpos1[i,0]<pos_min or newpos1[i,0]>pos_max or newpos1[i,1]<pos_min or newpos1[i,1]>pos_max or newpos1[i,2]<pos_min or newpos1[i,2]>pos_max:continue
            if not rejectcollision_endpoint(i,N,N1,b,c,Rx,Ry,Rz,Position3D,Occupation3D,l): continue


            U = energytwopoint_LJpotential_spa_inclu_energyk1(i,Position3D, newpos1,Occupation3D,cubesize_int)
            if not acceptSmp(U): continue

            Occupation3D[Position3D[i, 0], Position3D[i, 1], Position3D[i, 2]] = 0
            Position3D[i, :] = newpos1[i, :].copy()
            Occupation3D[newpos1[i, 0], newpos1[i, 1], newpos1[i, 2]] = 1
            t += 1



        elif i == (N - 1):
            c = random.randint(0, 2)
            b_random = random.randint(0,1)
            b=angle[b_random]
            u = Position3D[i, :] - Position3D[i-1, :]
            newpos1[i, :] = Position3D[i-1, :] + np.round(np.dot(RR[b_random][c], u))
            if Occupation3D[newpos1[i, 0], newpos1[i, 1], newpos1[i, 2]] >1e-6: continue
            if newpos1[i,0]<pos_min or newpos1[i,0]>pos_max or newpos1[i,1]<pos_min or newpos1[i,1]>pos_max or newpos1[i,2]<pos_min or newpos1[i,2]>pos_max:continue
            if not rejectcollision_endpoint(i,N,N1,b,c,Rx,Ry,Rz,Position3D,Occupation3D,l): continue
            U = energytwopoint_LJpotential_spa_inclu_energyk1(i,Position3D, newpos1,Occupation3D,cubesize_int)
            if not acceptSmp(U): continue
            Occupation3D[Position3D[i, 0], Position3D[i, 1], Position3D[i, 2]] = 0
            Position3D[i, :] = newpos1[i, :].copy()
            Occupation3D[newpos1[i, 0], newpos1[i, 1], newpos1[i, 2]] = 1
            t += 1



        elif i >= 2 and i <= N - 3:
            b_random = random.randint(0,1)
            b=angle[b_random]
            (v1,u1,newpos1[i,:],newpos1[i-1,:])=rotatetwo(i,b,Position3D,newpos1)            #to calculate the two points


            if newpos1[i,0]<pos_min or newpos1[i,0]>pos_max or newpos1[i,1]<pos_min or newpos1[i,1]>pos_max or newpos1[i,2]<pos_min or newpos1[i,2]>pos_max:continue
            if newpos1[i-1,0]<pos_min or newpos1[i-1,0]>pos_max or newpos1[i-1,1]<pos_min or newpos1[i-1,1]>pos_max or newpos1[i-1,2]<pos_min or newpos1[i-1,2]>pos_max:continue
            if np.linalg.norm(v1 - u1) < 1e-06:
            #if abs(u1[0]*v1[1]-u1[1]*v1[0])<1e-6 and abs(u1[1]*v1[2]-v1[1]*u1[2])<1e-6:
                a1 = Occupation3D[newpos1[i, 0] , newpos1[i, 1] , newpos1[i, 2] ]
                b1 = Occupation3D[newpos1[i - 1, 0] , newpos1[i - 1, 1] , newpos1[i - 1, 2] ]
                if a1 > 1e-6 or b1> 1e-6:continue

                if not rejectcollision_twopoints(i,N1,b,Position3D,Occupation3D,l,v1): continue
                U_i = energytwopoint_LJpotential_spa_inclu_energyk1(i,Position3D, newpos1,Occupation3D,cubesize_int)
                U_i_1 = energytwopoint_LJpotential_spa_inclu_energyk1(i-1,Position3D, newpos1,Occupation3D,cubesize_int)

                if not acceptSmp(U_i+U_i_1): continue

                Occupation3D[Position3D[i,0],Position3D[i,1],Position3D[i,2]]=0
                Position3D[i,:]=newpos1[i,:].copy()
                Occupation3D[newpos1[i,0],newpos1[i,1],newpos1[i,2]]=1
                Occupation3D[Position3D[i-1,0],Position3D[i-1,1],Position3D[i-1,2]]=0
                Position3D[i-1,:]=newpos1[i-1,:].copy()
                Occupation3D[newpos1[i-1,0],newpos1[i-1,1],newpos1[i-1,2]]=1
                t+=1


            else:
                i0=i+1
                (v1,u1,newpos1[i0,:],newpos1[i0-1,:])=rotatetwo(i0,b,Position3D,newpos1)
                if newpos1[i0,0]<pos_min or newpos1[i0,0]>pos_max or newpos1[i0,1]<pos_min or newpos1[i0,1]>pos_max or newpos1[i0,2]<pos_min or newpos1[i0,2]>pos_max:continue
                if newpos1[i0-1,0]<pos_min or newpos1[i0-1,0]>pos_max or newpos1[i0-1,1]<pos_min or newpos1[i0-1,1]>pos_max or newpos1[i0-1,2]<pos_min or newpos1[i0-1,2]>pos_max:continue
                if np.linalg.norm(v1 - u1) < 1e-06:
                    a1 = Occupation3D[newpos1[i0, 0] , newpos1[i0, 1] , newpos1[i0, 2]]
                    b1 = Occupation3D[newpos1[i0 - 1, 0] , newpos1[i0 - 1, 1] , newpos1[i0 - 1, 2] ]
                    if a1 > 1e-6 or b1> 1e-6:continue
                    if not rejectcollision_twopoints(i0,N1,b,Position3D,Occupation3D,l,v1): continue
                    U_i0 = energytwopoint_LJpotential_spa_inclu_energyk1(i0,Position3D, newpos1,Occupation3D,cubesize_int)
                    U_i0_1 = energytwopoint_LJpotential_spa_inclu_energyk1(i0-1,Position3D, newpos1,Occupation3D,cubesize_int)
                    if not acceptSmp(U_i0+U_i0_1): continue

                    Occupation3D[Position3D[i0,0],Position3D[i0,1],Position3D[i0,2]]=0
                    Position3D[i0,:]=newpos1[i0,:].copy()
                    Occupation3D[newpos1[i0,0],newpos1[i0,1],newpos1[i0,2]]=1
                    Occupation3D[Position3D[i0-1,0],Position3D[i0-1,1],Position3D[i0-1,2]]=0
                    Position3D[i0-1,:]=newpos1[i0-1,:].copy()
                    Occupation3D[newpos1[i0-1,0],newpos1[i0-1,1],newpos1[i0-1,2]]=1
                    t+=1


                else:
                    u3 = Position3D [i,:]- Position3D[i-1,:]
                    u5 = Position3D [i+1,:]- Position3D[i,:]
                    #if u5[0]-u3[0]!=u5[1]-u3[1]:
                    if np.linalg.norm(u3 - u5) > 1e-6:
                        w2 = Position3D [(i+1),:]- Position3D [(i-1),:]
                        w20,w21,w22 = w2/np.linalg.norm(w2)
                        b_random = random.randint(0,1)
                        b=angle[b_random]
                        a = np.pi*b
                        cosa=np.cos(a)
                        sina=np.sin(a)
                        n3 = np.array([[cosa+(1-cosa)*w20**2,(1-cosa)*w20*w21-sina*w22,(1-cosa)*w20*w22+sina*w21],
                                       [(1-cosa)*w20*w21+sina*w22,cosa+(1-cosa)*w21**2,(1-cosa)*w21*w22-sina*w20],
                                       [(1-cosa)*w22*w20-sina*w21,(1-cosa)*w22*w21+sina*w20,cosa+(1-cosa)*w22**2]
                        ])
                        #newpos0[i,:]=np.round(np.dot(n3,u3))
                        newpos1[i,:]=np.round(np.dot(n3,u3))+Position3D[i-1,:]
                        if newpos1[i,0]<pos_min or newpos1[i,0]>pos_max or newpos1[i,1]<pos_min or newpos1[i,1]>pos_max or newpos1[i,2]<pos_min or newpos1[i,2]>pos_max:continue
                        a1 = Occupation3D[newpos1[i, 0] , newpos1[i, 1] , newpos1[i, 2] ]
                        if a1>1e-6:continue
                        if not rejectcollision_middlepoint(i,b,N1,Position3D,Occupation3D,l,w20,w21,w22): continue
                        U = energytwopoint_LJpotential_spa_inclu_energyk1(i,Position3D, newpos1,Occupation3D,cubesize_int)
                        if not acceptSmp(U): continue
                        Occupation3D[Position3D[i,0],Position3D[i,1],Position3D[i,2]]=0
                        Position3D[i,:]=newpos1[i,:].copy()
                        Occupation3D[newpos1[i,0],newpos1[i,1],newpos1[i,2]]=1
                        t+=1


        elif i==1 or i==N-2:
        #abs(i)>1e-6 and abs(i-(N-1))>1e-6:
            i5=2
            if i==1:
                u2=Position3D[i+1,:]-Position3D[i,:]
                v2=Position3D[i+2,:]-Position3D[i-1,:]
                if np.linalg.norm(v2 - u2) > 1e-6:
                    i5=0

            elif i==N-2:
                u1 = Position3D[i,:]- Position3D [i-1,:]
                v1 = Position3D [i+1,:]- Position3D [i-2,:]
                if np.linalg.norm(v1 - u1) > 1e-6:
                    i5=0

            if i5 < 1e-6:
                u3 = Position3D [i,:]- Position3D[i-1,:]
                u5 = Position3D [i+1,:]- Position3D[i,:]
                #if u5[0]-u3[0]!=u5[1]-u3[1]:
                if np.linalg.norm(u3 - u5) > 1e-6:
                    w2 = Position3D [(i+1),:]- Position3D [(i-1),:]
                    w20,w21,w22 = w2/np.linalg.norm(w2)
                    b_random = random.randint(0,1)
                    b=angle[b_random]
                    a = np.pi*b
                    cosa=np.cos(a)
                    sina=np.sin(a)
                    n3 = np.array([[cosa+(1-cosa)*w20**2,(1-cosa)*w20*w21-sina*w22,(1-cosa)*w20*w22+sina*w21],
                                   [(1-cosa)*w20*w21+sina*w22,cosa+(1-cosa)*w21**2,(1-cosa)*w21*w22-sina*w20],
                                   [(1-cosa)*w22*w20-sina*w21,(1-cosa)*w22*w21+sina*w20,cosa+(1-cosa)*w22**2]
                    ])
                    #newpos0[i,:]=np.round(np.dot(n3,u3))
                    newpos1[i,:]=np.round(np.dot(n3,u3))+Position3D[i-1,:]
                    a1 = Occupation3D[newpos1[i, 0] , newpos1[i, 1] , newpos1[i, 2] ]
                    if a1>1e-6:continue
                    if newpos1[i,0]<pos_min or newpos1[i,0]>pos_max or newpos1[i,1]<pos_min or newpos1[i,1]>pos_max or newpos1[i,2]<pos_min or newpos1[i,2]>pos_max:continue
                    if not rejectcollision_middlepoint(i,b,N1,Position3D,Occupation3D,l,w20,w21,w22): continue
                    U = energytwopoint_LJpotential_spa_inclu_energyk1(i,Position3D, newpos1,Occupation3D,cubesize_int)
                    if not acceptSmp(U): continue

                    Occupation3D[Position3D[i,0],Position3D[i,1],Position3D[i,2]]=0
                    Position3D[i,:]=newpos1[i,:].copy()
                    Occupation3D[newpos1[i,0],newpos1[i,1],newpos1[i,2]]=1
                    t+=1
        if t%1e5==0:

            minpos=[9,9,9]
            maxpos=[40,40,40]

            Occupation3Dold=Occupation3D.copy()


    np.savez('%d.npz' % t, [Position3D, Occupation3D])

    for r in range(N-1):
        if np.linalg.norm(Position3D[r,:] - Position3D[r+1,:])>1:
            print 'false'
    for i in range(0,Position3D.shape[0]):
        m0=Position3D[i,0]
        m1=Position3D[i,1]
        m2=Position3D[i,2]
        Occupation3D[m0,m1,m2]=Occupation3D[m0,m1,m2]+1
        if Occupation3D[m0,m1,m2]>2:
            print 'repeat'

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    pp = PdfPages('mpage.pdf')

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(xs=Position3D[:, 0], ys=Position3D[:, 1], zs=Position3D[:, 2], zdir='z', label='ys=0,zdir=z')
    pp.savefig(bbox_inches='tight')
    pp.close()
    plt.show()



if __name__ == '__main__':
    r1(sys.argv[1])





