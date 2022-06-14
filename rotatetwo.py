#!/usr/bin/env python
import numpy as np
import random
def rotatetwo(i,b,Position3D,newpos1):
    u1 = Position3D [i,:]- Position3D [i-1,:]
    v1 = Position3D [i+1,:]- Position3D[i-2,:]
    w1 = Position3D [i,:]- Position3D[i-2,:]
    w2 = Position3D [i-1,:]- Position3D[i-2,:]

    if np.linalg.norm(v1 - u1) < 1e-06:
        v110,v111,v112=v1/np.linalg.norm(v1)
        a = b*(np.pi/2)
        cosa=np.cos(a)
        sina=np.sin(a)

        n1=np.array([[cosa+(1-cosa)*v110**2,(1-cosa)*v110*v111-sina*v112,(1-cosa)*v110*v112+sina*v111],
        [(1-cosa)*v110*v111+sina*v112,cosa+(1-cosa)*v111**2,(1-cosa)*v111*v112-sina*v110],
        [(1-cosa)*v112*v110-sina*v111,(1-cosa)*v112*v111+sina*v110,cosa+(1-cosa)*v112**2]
        ])


        newpos1[i,:]= np.round(np.dot(n1,w1)) + Position3D[i-2,:]
        newpos1[i-1,:]=np.round(np.dot(n1,w2)) + Position3D[i-2,:]
    return (v1,u1,newpos1[i,:],newpos1[i-1,:])

if __name__ == '__main__':
    rotatetwo()
