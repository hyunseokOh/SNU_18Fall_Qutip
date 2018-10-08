from qutip import *

import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=200000)
use_rwa=True

def calculate():
    rmat=np.zeros((len(tlist),3))
    c_ops=[]
    if g1 > 0.0:
        c_ops.append(np.sqrt(g1)*sigmam())
    if g2 > 0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    wmw=freq
    H0=0
    H1=0
    args={'wmw':wmw}

    H=(w0-wmw)/2.0*sigmaz()+wR/2.0*sigmax()
    output=mesolve(H,rho0,tlist,c_ops,[sigmax(),sigmay(),sigmaz()],options=options)

    rmat=[output.expect[0],output.expect[1],output.expect[2]]
    return rmat

w0=10*2*np.pi
wR=20*2*np.pi
freq=8.5*2*np.pi

tlist=np.linspace(0,1,201)

g1=0.0
g2=0.0

rho0=ket2dm(basis(2,0))
rmat=calculate()

b=Bloch()
b.add_points(rmat)
b.show()
