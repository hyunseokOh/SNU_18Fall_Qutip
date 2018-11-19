from qutip import *
import numpy as np
import matplotlib.pyplot as plt

def pulse(amount,ddx,ddz,err0,err1,azimuth=0):

    sx=sigmax()
    sy=sigmay()
    sz=sigmaz()

    dx=(ddx+err1)*np.cos(azimuth)
    dy=(ddx+err1)*np.sin(azimuth)
    dz=ddz+err0

    dxi=ddx*np.cos(azimuth)
    dyi=ddx*np.sin(azimuth)
    dzi=ddz

    omegaideal=np.sqrt(dxi*dxi+dyi*dyi+dzi*dzi)
    omega=np.sqrt(dx*dx+dy*dy+dz*dz)
    prop=qeye(2)

    if omegaideal !=0:
        theta=omega*amount/omegaideal

        prop=qeye(2)*np.cos(theta/2)-1j*(dz*sz+dx*sx+dy*sy)

    return prop


g1=0.0
g2=0.0
amount=np.pi/2
azimuth=0

rmat=np.zeros((9,9))
xlist=np.linspace(-0.25,0.25,9)
ylist=np.linspace(-0.25,0.25,9)
ddx=1/np.sqrt(2)
ddz=1/np.sqrt(2)

op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]]


for i in range(9):
    g1=ddz*xlist[i]

    for j in range(9):
        g2=ddx*ylist[j]

        U1=pulse(amount,ddx,ddz,0,0,azimuth)
        U2=pulse(amount,ddx,ddz,0,0,azimuth)
        chi1=Qobj(qpt(U1,op_basis))
        chi2=Qobj(qpt(U2,op_basis))
        
        #  U1=hadamard_transform()+1
        #  U2=hadamard_transform()+1
        rmat[i,j]=process_fidelity(chi1,chi2,normalize=True)

fig,ax=plt.subplots(figsize=(8,8))
x,y=np.meshgrid(xlist,ylist)
c1=ax.pcolor(x,y,rmat)
fig.colorbar(c1,ax=ax)
plt.show()

