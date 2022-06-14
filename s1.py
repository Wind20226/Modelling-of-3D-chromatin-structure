#!/usr/bin/env python

import numpy as np
def sawtooth_curve(x0,y0,z0,n):

    Position3D=np.zeros((n*n*(n+1),3))
    i0=0
    j0=0
    k0=0
    v = 1
    u = 1
    i = 1
    Position3D[0,:] = [x0+i0,y0+j0,z0+k0]
    while i<n**2-n/2:
        if 0<=i0+u<n:
            j0+=v
            Position3D[i,:] = [x0+i0,y0+j0,z0+k0]
            i+=1
            i0+=u
            Position3D[i,:] = [x0+i0,y0+j0,z0+k0]
            i+=1
            v=-v
        else:
            u=-u
            v=-v
            j0+=v
            Position3D[i,:] = [x0+i0,y0+j0,z0+k0]
            i+=1
    pos = Position3D[:i,:][::-1,:].copy()
    for z in range(1,n):
        pos[:,2] += 1
        Position3D[z*i:(z+1)*i,:] = pos
        pos = pos[::-1,:]
        if z==n-1:
            pos[:,2]+=1
            Position3D[(z+1)*i:(z+1)*i+n**2/2]=pos[:n**2/2]

    Position3D=Position3D[:n*i+n**2/2,:]

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
    return Position3D

if __name__ == '__main__':
    sawtooth_curve()
