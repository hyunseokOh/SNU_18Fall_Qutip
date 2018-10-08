import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from qutip.ui.progressbar import TextProgressBar as ProgressBar

options=Options(nsteps=200000)
use_rwa=True

def calculate():
    rmat = np.zeros((len(frelist),len(tlist)))
    pbar=ProgressBar(len(frelist))

    c_ops=[]
    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())
    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    for i in range(len(frelist)):
        pbar.update(i)
        wmw=frelist[i]
        H0=0
        H1=0
        args={'wmw':wmw}
        if use_rwa:
            H=(w0-wmw)/2.0*sigmaz()+wR/2.0*sigmax()
            output=mesolve(H,rho0,tlist,c_ops,[sigmax(),sigmay(),sigmaz()],options=options)

        else:
            H=[H0,[H1,'np.sin(wmw*t)*np.heaviside(0.5-t,0)']]
            output=mesolve(H,rho0,tlist,c_ops,[sigmax(),sigmay(),sigmaz()],options=options)
                
        rmat[i]=output.expect[2]
    return rmat

w0=10*2*np.pi
wR=2*2*np.pi

tlist=np.linspace(0,10,201)
frelist=np.linspace(6,14,51)*2*np.pi
g2=0.0
g1=0.0

rho0=ket2dm(basis(2,0))

rmat=calculate()


fig,ax=plt.subplots(figsize=(8,8))
t_mat,fre_mat=np.meshgrid(tlist,frelist/(2*np.pi))

ax.pcolor(t_mat,fre_mat,rmat)
ax.set_xlabel(r'time ($\mu$s)')
ax.set_ylabel(r'Frequency $MHz$')

plt.show()
