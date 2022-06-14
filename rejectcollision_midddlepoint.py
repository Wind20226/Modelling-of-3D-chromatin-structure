#!/usr/bin/env python
import numpy as np
def rejectcollision_middlepoint(i,b,N1,Position3D,Occupation3D,l,w20,w21,w22):
    N2=1
    a = np.pi/(2*b*N1)
    cosa=np.cos(a)
    sina=np.sin(a)
    n3 = np.array([[cosa+(1-cosa)*w20**2,(1-cosa)*w20*w21-sina*w22,(1-cosa)*w20*w22+sina*w21],
                   [(1-cosa)*w20*w21+sina*w22,cosa+(1-cosa)*w21**2,(1-cosa)*w21*w22-sina*w20],
                   [(1-cosa)*w22*w20-sina*w21,(1-cosa)*w22*w21+sina*w20,cosa+(1-cosa)*w22**2]
    ])
    Occupation3D[Position3D[i, 0] , Position3D[i, 1] , Position3D[i, 2] ] = 0
    newpos = Position3D[i,:].copy()
    vb = Position3D[i-1,:].copy()
    while N2 <= 2*N1-1:
        newpos = np.dot(n3,newpos-vb) + vb
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
    rejectcollision_middlepoint()
