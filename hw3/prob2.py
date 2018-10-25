from qutip import *
from qutip.ui.progressbar import TextProgressBar as ProgressBar
import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=20000)

def calculate():
    rmat=np.zeros((len(frelist),len(tlist)))
    pbar=ProgressBar(len(frelist))
    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    e_ops=[sigmax(),sigmay(),sigmaz()]
    
    for i in range(len(frelist)):
        pbar.update(i)
        wmw=frelist[i]
        args={'wmw':wmw}
        #  H0=(w01/2.0)*tensor(sigmaz(),qeye(2))+(w02/2.0)*tensor(qeye(2),sigmaz())+(J12/2.0)*tensor(sigmaz(),sigmaz())
        #  H1=Bac*tensor(sigmax(),qeye(2))
        #  H=[H0,[H1,'np.cos(wmw*t)']]
        H0=(w01-wmw)/2*tensor(sigmaz(),qeye(2))+(w02-wmw)/2*tensor(qeye(2),sigmaz())+J12/2*tensor(sigmaz(),sigmaz())
        Hint=(Bac1/2)*tensor(sigmax(),qeye(2))+(Bac2/2)*tensor(qeye(2),sigmax())
        H=H0+Hint
        output=mesolve(H,rho0,tlist,c_ops,[tensor(sigmax(),qeye(2)),tensor(sigmay(),qeye(2)),tensor(sigmaz(),qeye(2))],args,options=options)
        rmat[i]=output.expect[2]
    return rmat
       
w01=10*2*np.pi
w02=8*2*np.pi
Bac1=0.2*2*np.pi
Bac2=0
J12=1*2*np.pi
g1=0.0
g2=0.0

psi0=basis(2,0)
rho0=ket2dm(tensor(psi0,psi0))
tlist=np.linspace(0,10,101)
frelist=np.linspace(4,18,41)*2*np.pi
#  frelist=[11*2*np.pi]
rmat=calculate()

wmw=11*2*np.pi
H0=(w01-wmw)/2*tensor(sigmaz(),qeye(2))+(w02-wmw)/2*tensor(qeye(2),sigmaz())+J12/2*tensor(sigmaz(),sigmaz())
Hint=(Bac1/4)*tensor(sigmax(),qeye(2))+(Bac2/4)*tensor(qeye(2),sigmax())

plot_energy_levels([H0,Hint],labels=['non-interacting','interacting'],show_ylabels=True,figsize=(8,4))

fig,ax = plt.subplots(figsize=(8, 8))
t_mat, fre_mat = np.meshgrid(tlist, frelist/(2*np.pi))
ax.pcolor(t_mat, fre_mat, rmat)
ax.set_xlabel(r'time ($\mu$s)')
ax.set_ylabel(r'Frequency $MHz$')
ax.set_title("<Sz1>")
plt.show()

