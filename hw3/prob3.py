from qutip import *
from qutip.ui.progressbar import TextProgressBar as ProgressBar
import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=20000)

#  Find optimal frequency of microwave
def freqFind(rho0):
    rmat=np.zeros((len(frelist),len(tlist)))
    pbar=ProgressBar(len(frelist))
    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    #  H0=(w01/2.0)*tensor(sigmaz(),qeye(2))+(w02/2.0)*tensor(qeye(2),sigmaz())+(J12/2.0)*tensor(sigmaz(),sigmaz())
    #  H1=Bac*tensor(sigmax(),qeye(2))
    #  H=[H0,[H1,'np.cos(wmw*t)']]
    for i in range(len(frelist)):
        pbar.update(i)
        wmw2=frelist[i]
        H0=(w01-wmw1)/2*tensor(sigmaz(),qeye(2))+(w02-wmw2)/2*tensor(qeye(2),sigmaz())+J12*tensor(sigmaz(),sigmaz())
        Hint=(Bac1/4)*tensor(sigmax(),qeye(2))+(Bac2/4)*tensor(qeye(2),sigmax())
        H=H0+Hint
        output=mesolve(H,rho0,tlist,c_ops,[tensor(qeye(2),sigmaz())],options=options)
        rmat[i]=output.expect[0]
    return rmat

#  Parameter settings
w01=9*2*np.pi
w02=10*2*np.pi
wmw1=0*2*np.pi
wmw2=1*2*np.pi
Bac1=0*2*np.pi
Bac2=0.2*2*np.pi
J12=1*2*np.pi
g1=0.0
g2=0.0

fi=6
ff=16
ti=0
tf=10

frelist=np.linspace(fi,ff,101)*2*np.pi
tlist=np.linspace(ti,tf,101)

val=1

rho1=ket2dm(tensor(basis(2,0),basis(2,0)))
rmat1=freqFind(rho1)

rho2=ket2dm(tensor(basis(2,1),basis(2,0)))
rmat2=freqFind(rho2)

opt=np.where(rmat2==np.min(rmat2))
f_opt=frelist[opt[0][0]]
t_opt=tlist[opt[1][0]]

frelist=np.linspace(f_opt-(ff-fi)/25,f_opt+(ff-fi)/25,101)
tlist=np.linspace(t_opt-(tf-ti)/25,t_opt+(tf-ti)/25,101)
tlist=np.insert(tlist,0,[0])

rmat2=freqFind(rho2)
opt=np.where(rmat2==np.min(rmat2))
f_opt=frelist[opt[0][0]]
t_opt=tlist[opt[1][0]]
print("Optimal frequency : %f, Optimal time : %f" % (f_opt/2/np.pi,t_opt))
print()

#  fig,ax=plt.subplots()
#  t_mat,fre_mat=np.meshgrid(tlist[1:],frelist/(2*np.pi))
#  c1=ax.pcolor(t_mat,fre_mat,rmat2[:,1:])
#  fig.colorbar(c1,ax=ax)
#  plt.show()
#
#  fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,8))
#  t_mat,fre_mat=np.meshgrid(tlist,frelist/(2*np.pi))
#  c1=ax1.pcolor(t_mat,fre_mat,rmat1,vmin=np.min(rmat1),vmax=np.max(rmat1))
#  fig.colorbar(c1,ax=ax1)
#  c2=ax2.pcolor(t_mat,fre_mat,rmat2,vmin=np.min(rmat2),vmax=np.max(rmat2))
#  fig.colorbar(c2,ax=ax2)
#
#  ax1.set_title(r'<S\_z2 >, control=0')
#  ax2.set_title(r'<S\_z2 >, control=1')
#  ax1.set_xlabel(r'time ($\mu\s)')
#  ax1.set_ylabel(r'Frequency $MHz$')
#
#  ax2.set_xlabel(r'time ($\mu\s)')
#  ax2.set_ylabel(r'Frequency $MHz$')
#  plt.show()
#

def CNOT_impl(rho0):
    tl=[0,t_opt]
    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    #  H0=(w01/2.0)*tensor(sigmaz(),qeye(2))+(w02/2.0)*tensor(qeye(2),sigmaz())+(J12/2.0)*tensor(sigmaz(),sigmaz())
    #  H1=Bac*tensor(sigmax(),qeye(2))
    #  H=[H0,[H1,'np.cos(wmw*t)']]
    wmw2=f_opt
    H0=(w01-wmw1)/2*tensor(sigmaz(),qeye(2))+(w02-wmw2)/2*tensor(qeye(2),sigmaz())+J12*tensor(sigmaz(),sigmaz())
    Hint=(Bac1/4)*tensor(sigmax(),qeye(2))+(Bac2/4)*tensor(qeye(2),sigmax())
    H=H0+Hint
    output=mesolve(H,rho0,tl,c_ops,[],options=options)
    
    # Rotation without microwave, remove phase shift
    rho1=output.states[-1]
    Hrot=w01/2*tensor(sigmaz(),qeye(2))+w02/2*tensor(qeye(2),sigmaz())+J12*tensor(sigmaz(),sigmaz())
    tlist_rot=np.linspace(0,1.5*np.pi/J12,2)
    output=mesolve(Hrot,rho1,tlist_rot)
    
    return output.states


#  Prob B-1
psi1=basis(2,0)
psi2=basis(2,0)

rho0=ket2dm(tensor(psi1,psi2))
calc=CNOT_impl(rho0)[-1]
ideal=ket2dm(cnot()*tensor(psi1,psi2))

print("Prob B-1")
print(calc)
print("Fidelity for single CNOT : "+str(fidelity(ideal,calc)))
print()


#  Prob B-2
psi1=basis(2,1)
psi2=basis(2,0)

rho0=ket2dm(tensor(psi1,psi2))
calc=CNOT_impl(rho0)[-1]
ideal=ket2dm(cnot()*tensor(psi1,psi2))

print("Prob B-2")
print(calc)
print("Fidelity for single CNOT : "+str(fidelity(ideal,calc)))
print()


#  Prob C
psi2=basis(2,0)
psi1=1/np.sqrt(2)*(basis(2,0)+basis(2,1))

rho0=ket2dm(tensor(psi1,psi2))
calc=CNOT_impl(rho0)[-1]
ideal=ket2dm(cnot()*tensor(psi1,psi2))

print("Prob C")
print(calc)
print("Fidelity for single CNOT : "+str(fidelity(ideal,calc)))
print()


#  Prob D
psi1=1/np.sqrt(2)*(basis(2,0)+basis(2,1))
psi2=basis(2,0)

rho0=ket2dm(tensor(psi1,psi2))
calc=CNOT_impl(rho0)[-1]
calc2=CNOT_impl(calc)[-1]
ideal=ket2dm(cnot()*cnot()*tensor(psi1,psi2))

print("Prob D")
print(calc2)
print("Fidelity for single CNOT : "+str(fidelity(ideal,calc2)))
