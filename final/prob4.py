from qutip import *
from qutip.ui.progressbar import TextProgressBar as ProgressBar
import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=20000)

#  Calculate expectation value of Sz for each spin
def freqCalc(rho0,frelist,tlist,target):
    wmw1=0
    wmw2=0
    wmw3=0

    rmat=np.zeros((len(frelist),len(tlist)))
    pbar=ProgressBar(len(frelist))
    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    sy=sigmay()
    sz=sigmaz()
    idm=qeye(2)
    
    for i in range(len(frelist)):
        #  pbar.update(i)
        
        if target==1:
            e_ops=[tensor(sz,idm,idm)]
            wmw1=frelist[i]
        
        elif target==2:
            e_ops=[tensor(idm,sz,idm)]
            wmw2=frelist[i]

        else:
            e_ops=[tensor(idm,idm,sz)]
            wmw3=frelist[i]

        H0=(w01-wmw1)/2*tensor(sz,idm,idm)+(w02-wmw2)/2*tensor(idm,sz,idm)+(w03-wmw3)/2*tensor(idm,idm,sz)+J12*tensor(sz,sz,idm)+J23*tensor(idm,sz,sz)+J13*tensor(sz,idm,sz)
        Hint=(Bac1/4)*tensor(sy,idm,idm)+(Bac2/4)*tensor(idm,sy,idm)+(Bac3/4)*tensor(idm,idm,sy)
        H=H0+Hint
        
        output=mesolve(H,rho0,tlist,c_ops,e_ops,options=options)
        rmat[i]=output.expect[0]
    return rmat


#  Find optimal frequency and time for each spin
def freqFind(frelist,tlist,control,target):
    s0=basis(2,0)
    s1=basis(2,1)
    
    #  Find the position of <Sz>=1
    rho1=ket2dm(tensor(s0,s0,s0))
    rmat1=freqCalc(rho1,frelist,tlist,target)
    
    avoid=np.where(rmat1==np.min(rmat1))
    f_avoid=frelist[avoid[0][0]]
    t_avoid=tlist[avoid[1][0]]

    #  Find the position of <Sz>=-1
    if control==1:
        rho2=ket2dm(tensor(s1,s0,s0))

    elif control==2:
        rho2=ket2dm(tensor(s0,s1,s0))

    else:
        rho2=ket2dm(tensor(s0,s0,s1))
    
    rmat2=freqCalc(rho2,frelist,tlist,target)
    opt=np.where(rmat2==np.min(rmat2))
   
    f_opt=frelist[opt[0][0]]
    t_opt=tlist[opt[1][0]]
    
    #  Repeat searching with better resolution
    frelist=np.linspace(f_opt-(ff-fi)/25,f_opt+(ff-fi)/25,101)
    tlist=np.linspace(t_opt-(tf-ti)/25,t_opt+(tf-ti)/25,101)
    tlist=np.insert(tlist,0,[0])

    rmat2=freqCalc(rho2,frelist,tlist,target)
    
    opt=np.where(rmat2==np.min(rmat2))
    f_opt=frelist[opt[0][0]]
    t_opt=tlist[opt[1][0]]
   
    print("\nOptimal frequency : %f, Optimal time : %f" % (f_opt/2/np.pi,t_opt))
    print("Avoid frequency : %f, Optimal time : %f" % (f_avoid/2/np.pi,t_avoid))

    return (f_opt,t_opt)


#  Implementation of CNOT gate
def CNOT_impl(rho0,f_opt,t_opt,target):
    wmw1=0
    wmw2=0
    wmw3=0

    #  tl=[0,t_opt]
    tl=np.linspace(0,t_opt,100)
    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())
    
    if target==1:
        wmw1=f_opt
    else:
        wmw3=f_opt
    
    sz=sigmaz()
    sx=sigmax()
    sy=sigmay()
    idm=qeye(2)
    
    H0=(w01-wmw1)/2*tensor(sz,idm,idm)+(w02-wmw2)/2*tensor(idm,sz,idm)+(w03-wmw3)/2*tensor(idm,idm,sz)+J12*tensor(sz,sz,idm)+J23*tensor(idm,sz,sz)+J13*tensor(sz,idm,sz)
    Hint=(Bac1/4)*tensor(sy,idm,idm)+(Bac2/4)*tensor(idm,sy,idm)+(Bac3/4)*tensor(idm,idm,sy)
    H=H0+Hint
    
    output=mesolve(H,rho0,tl,c_ops,[],options=options)
    return output.states


#  Control amplitude of microwaves
def ampControl(target):
    Bac1=0
    Bac2=0
    Bac3=0
    val=0.2*2*np.pi
    if target==1:
        Bac1=val
    elif target==2:
        Bac2=val
    else:
        Bac3=val

    return Bac1,Bac2,Bac3


#  Parameter settings
w01=11*2*np.pi
w02=12*2*np.pi
w03=10*2*np.pi
Bac1=0.2*2*np.pi
Bac2=0.2*2*np.pi
Bac3=0.2*2*np.pi
J12=1*2*np.pi
J23=2*2*np.pi
J13=3*2*np.pi
g1=0.0
g2=0.0
alpha=0

fi=0
ff=20
ti=0
tf=10

frelist=np.linspace(fi,ff,101)*2*np.pi
tlist=np.linspace(ti,tf,101)

#  Find optimal frequency and time for each spin
print("\nSpin 1 frequency finding")
Bac1,Bac2,Bac3=ampControl(1)
f1,t1=freqFind(frelist,tlist,2,1)
print("\nSpin 2 frequency finding")
Bac1,Bac2,Bac3=ampControl(2)
f2,t2=freqFind(frelist,tlist,1,2)
print("\nSpin 3 frequency finding")
Bac1,Bac2,Bac3=ampControl(3)
f3,t3=freqFind(frelist,tlist,1,3)

#  Quantum Circuit
psi1=basis(2,0)
psi2=basis(2,0)
psi3=basis(2,0)

rho0=ket2dm(tensor(psi1,psi2,psi3))

Bac1,Bac2,Bac3=ampControl(2)
rho1=Hadamard2_impl(rho0,f2,t2)[-1]

Bac1,Bac2,Bac3=ampControl(1)
rho2=CNOT_impl(rho1,f1,t1,1)[-1]

Bac1,Bac2,Bac3=ampControl(3)
rhoTime=CNOT_impl(rho2,f3,t3,3)
rho3=rhoTime[-1]

#  Comparison with ideal denstiy matrix
ideal=ghz_state()
print(rho3)
print(ket2dm(ideal))
print("\nFidelity : %f\n" % fidelity(rho3,ket2dm(ideal)))

