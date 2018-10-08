from qutip import *
import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=20000)

def qubit_integrate(w0,vamp,wmw,g1,g2):
    H0=w0/2.0*sigmaz()
    H1=vamp/2.0*sigmax()
    args={'wmw':wmw}
    H=[H0,[H1,'np.cos(wmw*t)']]
    H_rwa=(w0-wmw)/2.0*sigmaz()+vamp/2.0*sigmax()

    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    e_ops=[sigmax(),sigmay(),sigmaz()]

    output=mesolve(H,psi0,tlist,c_ops,e_ops,args,options=options)
    output2=mesolve(H_rwa,psi0,tlist,c_ops,e_ops,args,options=options)

    return [(output.expect[0],output.expect[1],output.expect[2]),(output2.expect[0],output2.expect[1],output2.expect[2])]

w0=10.0*2*np.pi
vamp=15*2*np.pi
wmw=10.0*2*np.pi
g1=0.0
g2=0.0

psi0=basis(2,0)
tlist=np.linspace(0,1,1000)

[(sx1,sy1,sz1),(sx2,sy2,sz2)]=qubit_integrate(w0,vamp,wmw,g1,g2)

fig,ax=plt.subplots(figsize=(12,6))
ax.plot(tlist,np.real(sz1))
ax.plot(tlist,np.real(sz2))
ax.legend(("non-RWA","RWA"))
ax.set_xlabel('Time')
ax.set_ylabel('Sz expectation value')
plt.show()
