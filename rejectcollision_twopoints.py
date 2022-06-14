#!/usr/bin/env python
import numpy as np
def rejectcollision_twopoints(i,N1,b,Position3D,Occupation3D,l,v1):
    N2=1
    v110,v111,v112 = v1/np.linalg.norm(v1)
    a= np.pi/(2*b*N1)
    #if b>0:
     #   a = np.pi/(2*N1)
    #else:
    #    a = np.pi/(-2*N1)
    cosa=np.cos(a)
    sina=np.sin(a)
    n1=np.array([[cosa+(1-cosa)*v110**2,(1-cosa)*v110*v111-sina*v112,(1-cosa)*v110*v112+sina*v111],
    [(1-cosa)*v110*v111+sina*v112,cosa+(1-cosa)*v111**2,(1-cosa)*v111*v112-sina*v110],
    [(1-cosa)*v112*v110-sina*v111,(1-cosa)*v112*v111+sina*v110,cosa+(1-cosa)*v112**2]
    ])

    Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 0
    Occupation3D[Position3D[i-1, 0] , Position3D[i-1, 1] , Position3D[i-1, 2] ] = 0
    newposi = Position3D[i,:].copy()
    newposi1 = Position3D[i-1,:].copy()
    vb = Position3D[i-2,:].copy()
    while N2 <= N1-1:
        newposi= np.dot(n1,newposi-vb)+vb
        newposi1=np.dot(n1,newposi1-vb)+vb

        cube_lo_a=(np.ceil(newposi - l)).astype(int)
        cube_hi_a=(np.floor(newposi + l + 1)).astype(int)
        occu_a = Occupation3D[cube_lo_a[0]:cube_hi_a[0], cube_lo_a[1]:cube_hi_a[1], cube_lo_a[2]:cube_hi_a[2]].sum()
        cube_lo_b=(np.ceil(newposi1 - l)).astype(int)
        cube_hi_b=(np.floor(newposi1 + l + 1)).astype(int)
        occu_b = Occupation3D[cube_lo_b[0]:cube_hi_b[0], cube_lo_b[1]:cube_hi_b[1], cube_lo_b[2]:cube_hi_b[2]].sum()
        if occu_a  <1e-6 and occu_b<1e-6:
            N2 += 1


        else:
            Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 1
            Occupation3D[Position3D[i-1, 0] , Position3D[i-1, 1] , Position3D[i-1, 2] ] = 1
            return False
    Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 1
    Occupation3D[Position3D[i-1, 0] , Position3D[i-1, 1] , Position3D[i-1, 2] ] = 1
    return True

if __name__ == '__main__':
    rejectcollision_twopoints()
