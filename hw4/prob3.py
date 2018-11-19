from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import time
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
        prop=qeye(2)*np.cos(theta/2)-1j*(dz*sz+dx*sx+dy*sy)*np.sin(theta/2)/omega

    return prop


def convertLog(res):
    x,y=res.shape
    scale=np.zeros((x,y))
    for i in range(x):
        for j in range(y):
            if np.abs(1-res[i][j])>0.00001:
                scale[j][i]=-np.log10(1-res[i][j])
            else:
                scale[j][i]=5
    
    return scale


def figure(rmat,title):
    rmat=convertLog(rmat)
    fig,ax=plt.subplots(figsize=(8,8))
    x,y=np.meshgrid(del0,del1)
    c1=ax.pcolor(x,y,rmat,cmap='Spectral_r',vmin=0,vmax=5)
    ax.set_xlabel(r'delta0/w1')
    ax.set_ylabel(r'delta1/w1')
    ax.set_title("Fidelity of "+title)
    fig.colorbar(c1,ax=ax)
    fig.savefig("Prob3_"+title+".png")
    return

start_time=time.time()

g1=0.0
g2=0.0
theta=np.pi/2
w=1

xres=201
yres=201

rmat=np.zeros((xres,yres))
del0=np.linspace(-0.25,0.25,xres)*w
del1=np.linspace(-0.25,0.25,yres)*w


#ideal rotation
xrot=np.cos(theta/2)*qeye(2)-1j*np.sin(theta/2)*sigmax()
ideal=spre(xrot)*spost(xrot.dag())
op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]]
chi2=Qobj(qpt(ideal,op_basis))
rot=np.pi/2

# Rectangular naive pulse
rot=np.pi/2
phi=0

for i in range(xres):
    for j in range(yres):
        U1=pulse(rot,w,0,del0[i],del1[j],phi)
        U1s=spre(U1)*spost(U1.dag())
        chi1=Qobj(qpt(U1s,op_basis))

        rmat[i,j]=process_fidelity(chi1,chi2,normalize=True)

figure(rmat,"Rectangular naive pulse")
print("1/4 finish, %.1fs lapsed" % (time.time()-start_time))

# Five piece SUPCODE
rot2=rot+2*np.pi
t1=1/np.sin(rot2)*(1-2*np.cos(rot2/2)+np.cos(rot2)+np.sqrt(4-8*np.cos(rot2/2)+4*np.cos(rot2)+rot*np.sin(rot2)))
t3=-2*(t1*np.cos(rot2/2)+np.sin(rot2/2))
phi=0

for i in range(xres):
    for j in range(yres):
        U1=pulse(t1,0,0,del0[i],del1[j],phi)
        U2=pulse(rot2/2,w,0,del0[i],del1[j],phi)
        U3=pulse(t3,0,0,del0[i],del1[j],phi)
        U4=pulse(rot2/2,w,0,del0[i],del1[j],phi)
        U5=pulse(t1,0,0,del0[i],del1[j],phi)

        Utot=U5*U4*U3*U2*U1
        Utots=spre(Utot)*spost(Utot.dag())
        chi1=Qobj(qpt(Utots,op_basis))

        rmat[i,j]=process_fidelity(chi1,chi2,normalize=True)

figure(rmat,"Five piece SUPCODE")
print("2/4 finish, %.1fs lapsed" % (time.time()-start_time))

# BB1 pulse
phi=np.arccos(-rot/4/np.pi)

for i in range(xres):
    for j in range(yres):
        U1=pulse(rot/2,w,0,del0[i],del1[j],0)
        U2=pulse(np.pi,w,0,del0[i],del1[j],phi)
        U3=pulse(np.pi*2,w,0,del0[i],del1[j],3*phi)
        U4=pulse(np.pi,w,0,del0[i],del1[j],phi)
        U5=pulse(rot/2,w,0,del0[i],del1[j],0)
        
        Utot=U5*U4*U3*U2*U1
        Utots=spre(Utot)*spost(Utot.dag())
        chi1=Qobj(qpt(Utots,op_basis))

        rmat[i,j]=process_fidelity(chi1,chi2,normalize=True)

figure(rmat,"BB1 pulse")
print("3/4 finish, %.1fs lapsed" % (time.time()-start_time))

# BB1inC sequence
phi=np.arccos(-rot/4/np.pi)
theta1=rot/2-np.arcsin(np.sin(rot/2)/2)
theta2=2*np.pi-2*np.arcsin(np.sin(rot/2)/2)
theta3=2*np.pi-np.arcsin(np.sin(rot/2)/2)

for i in range(xres):
    for j in range(yres):
        U1=pulse(rot/2,w,0,del0[i],del1[j],0)
        U2=pulse(np.pi,w,0,del0[i],del1[j],phi)
        U3=pulse(np.pi*2,w,0,del0[i],del1[j],3*phi)
        U4=pulse(np.pi,w,0,del0[i],del1[j],phi)
        U5=pulse(theta3,w,0,del0[i],del1[j],0)
        U6=pulse(theta2,w,0,del0[i],del1[j],np.pi)
        U7=pulse(theta1,w,0,del0[i],del1[j],0)

        Utot=U7*U6*U5*U4*U3*U2*U1
        Utots=spre(Utot)*spost(Utot.dag())
        chi1=Qobj(qpt(Utots,op_basis))

        rmat[i,j]=process_fidelity(chi1,chi2,normalize=True)

figure(rmat,"BB1inC pulse")
print("4/4 finish, %.1fs lapsed" % (time.time()-start_time))
plt.show()
