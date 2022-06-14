#!/usr/bin/env python
import numpy as np

import math, random
#k=1

def acceptSmp(eng):

    return eng >= 0 or math.exp(eng) > random.uniform(0,1)

def energy(dd):
    d2 = dd.dot(dd)
    if d2 < 1e-6: return 0
    nd_1 = 1.0/d2
    #return 0 if (nd_1/0.81) < 2.5**-2 else 4*((nd_1/0.81)**6 - (nd_1/0.81)**3)+1

    #return 0 if nd_1 < 2**-(1.0/12) else 4*(nd_1**6 - nd_1**3)+1
    #return 0 if 1-nd_1>1e-6 else 1
    return 0 if abs(1-d2)>1e-2 else 1

    #return 0 if (nd_1/0.81) < 2**-(1.0/12) else 4*((nd_1/0.81)**6 - (nd_1/0.81)**3)+1


def energytwopoint_LJpotential_spa_inclu_energyk1(i,Position3D, newpos1,Occupation3D,cubesize_int):
    U1 = 0
    U2 = 0
    N=Position3D.shape[0]
    cube_lo=np.array(newpos1[i,:] - cubesize_int, dtype=int)

    cube_hi=np.array(newpos1[i,:] + cubesize_int + 1, dtype=int)

    nearby1=np.zeros((2*cubesize_int+1,3),dtype=int)
    nearby1=[[cube_lo[0],cube_lo[1],cube_lo[2]],
             [cube_lo[0],cube_lo[1],cube_lo[2]+1],
             [cube_lo[0],cube_lo[1],cube_lo[2]+2],
             [cube_lo[0],cube_lo[1],cube_lo[2]+3],
             [cube_lo[0],cube_lo[1],cube_lo[2]+4],
             [cube_lo[0],cube_lo[1]+1,cube_lo[2]],
             [cube_lo[0],cube_lo[1]+1,cube_lo[2]+1],
             [cube_lo[0],cube_lo[1]+1,cube_lo[2]+2],
             [cube_lo[0],cube_lo[1]+1,cube_lo[2]+3],
             [cube_lo[0],cube_lo[1]+1,cube_lo[2]+4],
             [cube_lo[0],cube_lo[1]+2,cube_lo[2]],
             [cube_lo[0],cube_lo[1]+2,cube_lo[2]+1],
             [cube_lo[0],cube_lo[1]+2,cube_lo[2]+2],
             [cube_lo[0],cube_lo[1]+2,cube_lo[2]+3],
             [cube_lo[0],cube_lo[1]+2,cube_lo[2]+4],
             [cube_lo[0],cube_lo[1]+3,cube_lo[2]],
             [cube_lo[0],cube_lo[1]+3,cube_lo[2]+1],
             [cube_lo[0],cube_lo[1]+3,cube_lo[2]+2],
             [cube_lo[0],cube_lo[1]+3,cube_lo[2]+3],
             [cube_lo[0],cube_lo[1]+3,cube_lo[2]+4],
             [cube_lo[0],cube_lo[1]+4,cube_lo[2]],
             [cube_lo[0],cube_lo[1]+4,cube_lo[2]+1],
             [cube_lo[0],cube_lo[1]+4,cube_lo[2]+2],
             [cube_lo[0],cube_lo[1]+4,cube_lo[2]+3],
             [cube_lo[0],cube_lo[1]+4,cube_lo[2]+4],
             [cube_lo[0]+1,cube_lo[1],cube_lo[2]],
             [cube_lo[0]+1,cube_lo[1],cube_lo[2]+1],
             [cube_lo[0]+1,cube_lo[1],cube_lo[2]+2],
             [cube_lo[0]+1,cube_lo[1],cube_lo[2]+3],
             [cube_lo[0]+1,cube_lo[1],cube_lo[2]+4],
             [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]],
             [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+1],
             [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+2],
             [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+3],
             [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+4],
             [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]],
             [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+1],
             [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+2],
             [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+3],
             [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+4],
             [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]],
             [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+1],
             [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+2],
             [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+3],
             [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+4],
             [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]],
             [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+1],
             [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+2],
             [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+3],
             [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+4],
             [cube_lo[0]+2,cube_lo[1],cube_lo[2]],
             [cube_lo[0]+2,cube_lo[1],cube_lo[2]+1],
             [cube_lo[0]+2,cube_lo[1],cube_lo[2]+2],
             [cube_lo[0]+2,cube_lo[1],cube_lo[2]+3],
             [cube_lo[0]+2,cube_lo[1],cube_lo[2]+4],
             [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]],
             [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+1],
             [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+2],
             [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+3],
             [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+4],
             [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]],
             [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]+1],
             [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]+3],
             [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]+4],
             [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]],
             [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+1],
             [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+2],
             [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+3],
             [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+4],
             [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]],
             [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+1],
             [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+2],
             [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+3],
             [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+4],
             [cube_lo[0]+3,cube_lo[1],cube_lo[2]],
             [cube_lo[0]+3,cube_lo[1],cube_lo[2]+1],
             [cube_lo[0]+3,cube_lo[1],cube_lo[2]+2],
             [cube_lo[0]+3,cube_lo[1],cube_lo[2]+3],
             [cube_lo[0]+3,cube_lo[1],cube_lo[2]+4],
             [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]],
             [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+1],
             [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+2],
             [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+3],
             [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+4],
             [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]],
             [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+1],
             [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+2],
             [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+3],
             [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+4],
             [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]],
             [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+1],
             [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+2],
             [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+3],
             [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+4],
             [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]],
             [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+1],
             [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+2],
             [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+3],
             [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+4],
             [cube_lo[0]+4,cube_lo[1],cube_lo[2]],
             [cube_lo[0]+4,cube_lo[1],cube_lo[2]+1],
             [cube_lo[0]+4,cube_lo[1],cube_lo[2]+2],
             [cube_lo[0]+4,cube_lo[1],cube_lo[2]+3],
             [cube_lo[0]+4,cube_lo[1],cube_lo[2]+4],
             [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]],
             [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+1],
             [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+2],
             [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+3],
             [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+4],
             [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]],
             [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+1],
             [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+2],
             [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+3],
             [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+4],
             [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]],
             [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+1],
             [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+2],
             [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+3],
             [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+4],
             [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]],
             [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+1],
             [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+2],
             [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+3],
             [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+4]]
    for v_i,v_near in enumerate(nearby1):

        if Occupation3D[tuple(v_near)]>1e-6:

            U1+=energy(v_near-newpos1[i,:])
            U2+=energy(v_near-Position3D[i,:])

    return U2-U1

def energytwopoint_LJpotential_spa_inclu_energyk1_system(Position3D,Occupation3D,cubesize_int):
    U1 = 0
    U2 = 0
    U2_sys=0
    N=Position3D.shape[0]
    for i in range(N):
        cube_lo=np.array(Position3D[i,:] - cubesize_int, dtype=int)
        cube_hi=np.array(Position3D[i,:] + cubesize_int + 1, dtype=int)
        nearby1=np.zeros((2*cubesize_int+1,3),dtype=float)
        nearby1=[[cube_lo[0],cube_lo[1],cube_lo[2]],
                 [cube_lo[0],cube_lo[1],cube_lo[2]+1],
                 [cube_lo[0],cube_lo[1],cube_lo[2]+2],
                 [cube_lo[0],cube_lo[1],cube_lo[2]+3],
                 [cube_lo[0],cube_lo[1],cube_lo[2]+4],
                 [cube_lo[0],cube_lo[1]+1,cube_lo[2]],
                 [cube_lo[0],cube_lo[1]+1,cube_lo[2]+1],
                 [cube_lo[0],cube_lo[1]+1,cube_lo[2]+2],
                 [cube_lo[0],cube_lo[1]+1,cube_lo[2]+3],
                 [cube_lo[0],cube_lo[1]+1,cube_lo[2]+4],
                 [cube_lo[0],cube_lo[1]+2,cube_lo[2]],
                 [cube_lo[0],cube_lo[1]+2,cube_lo[2]+1],
                 [cube_lo[0],cube_lo[1]+2,cube_lo[2]+2],
                 [cube_lo[0],cube_lo[1]+2,cube_lo[2]+3],
                 [cube_lo[0],cube_lo[1]+2,cube_lo[2]+4],
                 [cube_lo[0],cube_lo[1]+3,cube_lo[2]],
                 [cube_lo[0],cube_lo[1]+3,cube_lo[2]+1],
                 [cube_lo[0],cube_lo[1]+3,cube_lo[2]+2],
                 [cube_lo[0],cube_lo[1]+3,cube_lo[2]+3],
                 [cube_lo[0],cube_lo[1]+3,cube_lo[2]+4],
                 [cube_lo[0],cube_lo[1]+4,cube_lo[2]],
                 [cube_lo[0],cube_lo[1]+4,cube_lo[2]+1],
                 [cube_lo[0],cube_lo[1]+4,cube_lo[2]+2],
                 [cube_lo[0],cube_lo[1]+4,cube_lo[2]+3],
                 [cube_lo[0],cube_lo[1]+4,cube_lo[2]+4],
                 [cube_lo[0]+1,cube_lo[1],cube_lo[2]],
                 [cube_lo[0]+1,cube_lo[1],cube_lo[2]+1],
                 [cube_lo[0]+1,cube_lo[1],cube_lo[2]+2],
                 [cube_lo[0]+1,cube_lo[1],cube_lo[2]+3],
                 [cube_lo[0]+1,cube_lo[1],cube_lo[2]+4],
                 [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]],
                 [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+1],
                 [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+2],
                 [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+3],
                 [cube_lo[0]+1,cube_lo[1]+1,cube_lo[2]+4],
                 [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]],
                 [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+1],
                 [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+2],
                 [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+3],
                 [cube_lo[0]+1,cube_lo[1]+2,cube_lo[2]+4],
                 [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]],
                 [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+1],
                 [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+2],
                 [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+3],
                 [cube_lo[0]+1,cube_lo[1]+3,cube_lo[2]+4],
                 [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]],
                 [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+1],
                 [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+2],
                 [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+3],
                 [cube_lo[0]+1,cube_lo[1]+4,cube_lo[2]+4],
                 [cube_lo[0]+2,cube_lo[1],cube_lo[2]],
                 [cube_lo[0]+2,cube_lo[1],cube_lo[2]+1],
                 [cube_lo[0]+2,cube_lo[1],cube_lo[2]+2],
                 [cube_lo[0]+2,cube_lo[1],cube_lo[2]+3],
                 [cube_lo[0]+2,cube_lo[1],cube_lo[2]+4],
                 [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]],
                 [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+1],
                 [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+2],
                 [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+3],
                 [cube_lo[0]+2,cube_lo[1]+1,cube_lo[2]+4],
                 [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]],
                 [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]+1],
                 [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]+3],
                 [cube_lo[0]+2,cube_lo[1]+2,cube_lo[2]+4],
                 [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]],
                 [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+1],
                 [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+2],
                 [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+3],
                 [cube_lo[0]+2,cube_lo[1]+3,cube_lo[2]+4],
                 [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]],
                 [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+1],
                 [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+2],
                 [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+3],
                 [cube_lo[0]+2,cube_lo[1]+4,cube_lo[2]+4],
                 [cube_lo[0]+3,cube_lo[1],cube_lo[2]],
                 [cube_lo[0]+3,cube_lo[1],cube_lo[2]+1],
                 [cube_lo[0]+3,cube_lo[1],cube_lo[2]+2],
                 [cube_lo[0]+3,cube_lo[1],cube_lo[2]+3],
                 [cube_lo[0]+3,cube_lo[1],cube_lo[2]+4],
                 [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]],
                 [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+1],
                 [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+2],
                 [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+3],
                 [cube_lo[0]+3,cube_lo[1]+1,cube_lo[2]+4],
                 [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]],
                 [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+1],
                 [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+2],
                 [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+3],
                 [cube_lo[0]+3,cube_lo[1]+2,cube_lo[2]+4],
                 [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]],
                 [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+1],
                 [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+2],
                 [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+3],
                 [cube_lo[0]+3,cube_lo[1]+3,cube_lo[2]+4],
                 [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]],
                 [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+1],
                 [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+2],
                 [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+3],
                 [cube_lo[0]+3,cube_lo[1]+4,cube_lo[2]+4],
                 [cube_lo[0]+4,cube_lo[1],cube_lo[2]],
                 [cube_lo[0]+4,cube_lo[1],cube_lo[2]+1],
                 [cube_lo[0]+4,cube_lo[1],cube_lo[2]+2],
                 [cube_lo[0]+4,cube_lo[1],cube_lo[2]+3],
                 [cube_lo[0]+4,cube_lo[1],cube_lo[2]+4],
                 [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]],
                 [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+1],
                 [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+2],
                 [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+3],
                 [cube_lo[0]+4,cube_lo[1]+1,cube_lo[2]+4],
                 [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]],
                 [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+1],
                 [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+2],
                 [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+3],
                 [cube_lo[0]+4,cube_lo[1]+2,cube_lo[2]+4],
                 [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]],
                 [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+1],
                 [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+2],
                 [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+3],
                 [cube_lo[0]+4,cube_lo[1]+3,cube_lo[2]+4],
                 [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]],
                 [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+1],
                 [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+2],
                 [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+3],
                 [cube_lo[0]+4,cube_lo[1]+4,cube_lo[2]+4]]
        for v_i,v_near in enumerate(nearby1):

            if Occupation3D[tuple(v_near)]>1e-6:

                U2+=energy(v_near-Position3D[i,:])
        U2_sys+=U2
    return U2_sys

if __name__ == '__main__':
    energytwopoint_LJpotential_spa_inclu_energyk1()
