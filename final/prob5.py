from qutip import *
from qutip.ui.progressbar import TextProgressBar as ProgressBar
import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=20000)


#  Implementation of free evolution
def freeEvol(rho0,tau,target):
    wmw1=0
    wmw2=0

    tl=np.linspace(0,tau,2)
    c_ops=[]

    if g1>0.0:
        c_ops.append(tensor(np.sqrt(g1)*sigmam(),qeye(2)))

    if g2>0.0:
        c_ops.append(tensor(np.sqrt(g2)*sigmaz(),qeye(2)))

    sz=sigmaz()
    sx=sigmax()
    sy=sigmay()
    idm=qeye(2)

    H0=(w01-wmw1)/2*tensor(sz,idm)+(w02-wmw2)/2*tensor(idm,sz)+J12*tensor(sz,sz)

    output=mesolve(H0,rho0,tl,c_ops,[],options=options)
    return output.states[-1]


#  Make unitary matrix for random interaction with spin bath
def randomUnitary():
    phi=np.random.random()*2*np.pi
    U=(1j*sigmaz()*phi).expm()
    return tensor(qeye(2),U)


#  Calculate ramsey interferometry without feedback
def ramsey(tau):
    s0=basis(2,0)
    s1=basis(2,1)
    psi0=1/np.sqrt(2)*(s0+s1)
    rho0=ket2dm(tensor(psi0,s0))

    rho1=freeEvol(rho0,tau,1)
    half_flip=tensor(ry(np.pi/2),idm)
    rho2=half_flip.dag()*rho1*half_flip
    return ((rho2*tensor(sigmaz(),qeye(2))).tr()+1)/2


#  Calculate ramsey interferometry with coherent feedback
def ramseyFeedback(tau):
    s0=basis(2,0)
    s1=basis(2,1)
    psi0=1/np.sqrt(2)*(s0+s1)
    rho0=ket2dm(tensor(s0,psi0))

    rho1=freeEvol(rho0,tau,2)
    half_flip=tensor(idm,ry(np.pi/2))
    rho2=half_flip.dag()*rho1*half_flip
    
    hadamard_1=tensor(snot(),idm)
    rho3=hadamard_1.dag()*rho2*hadamard_1
    
    rho4=cnot().dag()*rho3*cnot()
    
    bathInter=randomUnitary()
    rho5=bathInter.dag()*rho4*bathInter

    rho6=cnot().dag()*rho5*cnot()

    rho7=hadamard_1.dag()*rho6*hadamard_1

    rho8=cphase(np.pi).dag()*rho7*cphase(np.pi)

    return ((rho8*tensor(qeye(2),sigmaz())).tr()+1)/2


#  Parameter settings
w01=11*2*np.pi
w02=12*2*np.pi
Bac1=0.2*2*np.pi
Bac2=0.2*2*np.pi
J12=1*2*np.pi
g1=0.0
g2=0.5

#  Basic configuration
sx=sigmax()
sy=sigmay()
sz=sigmaz()
idm=qeye(2)
s0=basis(2,0)
s1=basis(2,1)


tlist=np.linspace(0,2,501)
rmat=np.zeros(len(tlist))
rmat2=np.zeros(len(tlist))

for i in range(len(tlist)):
    tau=tlist[i]
    rmat[i]=ramsey(tau)
    rmat2[i]=ramseyFeedback(tau)

fig,ax=plt.subplots()
ax.plot(tlist,rmat)
ax.plot(tlist,rmat2)
ax.set_title("Comparison of feedback effects on Ramsey interferometry")
ax.set_xlabel(r'Time [$\mu s$]')
ax.set_ylabel("Probability")
ax.legend(loc=0,labels=["No feedback","Coherent feedback"])
fig.savefig("Prob5.png")
plt.show()

