from qutip import *

import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=200000)
use_rwa=True

def calculate(rho0):
    rmat=np.zeros((len(tlist),3))
    c_ops=[]
    if g1 > 0.0:
        c_ops.append(np.sqrt(g1)*sigmam())
    if g2 > 0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    wmw=freq
    args={'wmw':wmw}

    H=(w0-wmw)/2.0*sigmaz()+wR/2.0*sigmax()
    output=mesolve(H,rho0,tlist,c_ops,[sigmax(),sigmay(),sigmaz()],options=options)

    rmat=[output.expect[0],output.expect[1],output.expect[2]]
    return rmat

w0=10*2*np.pi
wR=10*2*np.pi
freq=0

tlist=[0,np.pi/np.sqrt(2)/w0]

g1=0.0
g2=0.0

init0=basis(2,0)
init1=basis(2,1)
rho0=ket2dm(init0)
rho1=ket2dm(init1)

rmat0=calculate(rho0)
rmat1=calculate(rho1)

vec0=[rmat0[0][1],rmat0[1][1],rmat0[2][1]]
vec1=[rmat1[0][1],rmat1[1][1],rmat1[2][1]]

b0=Bloch()
b0.add_states(init0)
b0.add_vectors(vec0)
b0.show()

b1=Bloch()
b1.add_states(init1)
b1.add_vectors(vec1)
b1.show()
