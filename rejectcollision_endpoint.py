#!/usr/bin/env python
import numpy as np
def rejectcollision_endpoint(i,N,N1,b,c,Rx,Ry,Rz,Position3D,Occupation3D,l):
    RR1 = [Rx(np.pi / (2*b*N1)), Ry(np.pi / (2*b*N1)), Rz(np.pi / (2*b*N1))]
    N2 = 1
    Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 0
    newpos = Position3D[i,:].copy()
    if i==0:
        vb = Position3D[i+1,:].copy()
        while N2 <= N1-1:
            newpos = vb + np.dot(RR1[c], newpos - vb)
            cube_lo=(np.ceil(newpos - l)).astype(int)
            cube_hi=(np.floor(newpos + l + 1)).astype(int)
            occu = Occupation3D[cube_lo[0]:cube_hi[0], cube_lo[1]:cube_hi[1], cube_lo[2]:cube_hi[2]].sum()
            if occu <1e-6:
                N2 += 1

            else:
                Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 1
                return False

    if i==N-1:
        vb = Position3D[i-1,:].copy()
        while N2 <= N1-1:
            newpos = vb + np.dot(RR1[c], newpos-vb)
            cube_lo=(np.ceil(newpos - l)).astype(int)
            cube_hi=(np.floor(newpos + l + 1)).astype(int)
            occu = Occupation3D[cube_lo[0]:cube_hi[0], cube_lo[1]:cube_hi[1], cube_lo[2]:cube_hi[2]].sum()
            if occu <1e-6:
                N2 += 1

            else:
                Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 1
                return False
    Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 1
    return True

if __name__ == '__main__':
    rejectcollision_endpoint()
