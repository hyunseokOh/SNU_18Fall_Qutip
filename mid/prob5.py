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

    sx=sigmax()
    sz=sigmaz()
    idm=qeye(2)
    
    for i in range(len(frelist)):
        pbar.update(i)
        
        if target==1:
            e_ops=[tensor(sz,idm,idm)]
            wmw1=frelist[i]
        
        elif target==2:
            e_ops=[tensor(idm,sz,idm)]
            wmw2=frelist[i]

        else:
            e_ops=[tensor(idm,idm,sz)]
            wmw3=frelist[i]

        H0=(w01-wmw1)/2*tensor(sz,idm,idm)+(w02-wmw2)/2*tensor(idm,sz,idm)+(w03-wmw3)/2*tensor(idm,idm,sz)+J12*tensor(sz,sz,idm)+J23*tensor(idm,sz,sz)
        Hint=(Bac1/4)*tensor(sx,idm,idm)+(Bac2/4)*tensor(idm,sx,idm)+(Bac3/4)*tensor(idm,idm,sx)
        H=H0+Hint
        
        output=mesolve(H,rho0,tlist,c_ops,e_ops,options=options)
        rmat[i]=output.expect[0]
    return rmat


#  Find optimal frequency and time for each spin
def freqFind(frelist,tlist,target):
    rho1=ket2dm(tensor(basis(2,0),basis(2,0),basis(2,0)))
    rmat1=freqCalc(rho1,frelist,tlist,target)
   
    #  Find the position of <Sz>=0
    if target==2:
        rho2=rho1
        rmat2=rmat1
        opt=np.where(np.abs(rmat2)==np.min(np.abs(rmat2)))
    
    #  Find the position of <Sz>=-1
    else:
        rho2=ket2dm(tensor(basis(2,0),basis(2,1),basis(2,0)))
        rmat2=freqCalc(rho2,frelist,tlist,target)
        opt=np.where(rmat2==np.min(rmat2))
    
    f_opt=frelist[opt[0][0]]
    t_opt=tlist[opt[1][0]]
    
    #  Repeat searching with better resolution
    frelist=np.linspace(f_opt-(ff-fi)/25,f_opt+(ff-fi)/25,101)
    tlist=np.linspace(t_opt-(tf-ti)/25,t_opt+(tf-ti)/25,101)
    tlist=np.insert(tlist,0,[0])

    rmat2=freqCalc(rho2,frelist,tlist,target)
    
    if target==2:
        opt=np.where(np.abs(rmat2)==np.min(np.abs(rmat2)))
    else:
        opt=np.where(rmat2==np.min(rmat2))
    
    f_opt=frelist[opt[0][0]]
    t_opt=tlist[opt[1][0]]
    print("\nOptimal frequency : %f, Optimal time : %f\n" % (f_opt/2/np.pi,t_opt))
    
    return (f_opt,t_opt)


#  Implementation of Hadamard gate
def Hadamard2_impl(rho0,f_opt,t_opt):
    wmw1=0
    wmw2=f_opt
    wmw3=0

    sx=sigmax()
    sz=sigmaz()
    idm=qeye(2)
    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())
     
    H0=(w01-wmw1)/2*tensor(sz,idm,idm)+(w02-wmw2)/2*tensor(idm,sz,idm)+(w03-wmw3)/2*tensor(idm,idm,sz)+J12*tensor(sz,sz,idm)+J23*tensor(idm,sz,sz)
    Hint=(Bac1/4)*tensor(sx,idm,idm)+(Bac2/4)*tensor(idm,sx,idm)+(Bac3/4)*tensor(idm,idm,sx)
    H=H0+Hint
    
    tlist_rot=[0,t_opt]
    output=mesolve(H,rho0,tlist_rot,c_ops,[],options=options)
    
     # Rotation without microwave, remove phase shift
    rho1=output.states[-1]
    Hrot=w01/2*tensor(sz,idm,idm)+w02/2*tensor(idm,sz,idm)+w03/2*tensor(idm,idm,sz)+J12*tensor(sz,sz,idm)+J23*tensor(idm,sz,sz)
    
    tlist_rot=np.linspace(0,0.5*np.pi/f_opt,2)

    output=mesolve(Hrot,rho1,tlist_rot,c_ops,[],options=options)
    
    return output.states


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
    
    H0=(w01-wmw1)/2*tensor(sz,idm,idm)+(w02-wmw2)/2*tensor(idm,sz,idm)+(w03-wmw3)/2*tensor(idm,idm,sz)+J12*tensor(sz,sz,idm)+J23*tensor(idm,sz,sz)
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


#  Calculate scala product of two matrices
def scalProduct(m1,m2):
    M=m1*m2
    return M.tr()


#  Calculate mermin product of density matrix
def merminProduct(rho0):
    sx=sigmax()
    sy=sigmay()
    
    xxx=scalProduct(tensor(sx,sx,sx),rho0)
    xyy=scalProduct(tensor(sx,sy,sy),rho0)
    yxy=scalProduct(tensor(sy,sx,sy),rho0)
    xyy=scalProduct(tensor(sx,sy,sy),rho0)

    yyy=scalProduct(tensor(sy,sy,sy),rho0)
    yxx=scalProduct(tensor(sy,sx,sx),rho0)
    xyx=scalProduct(tensor(sx,sy,sx),rho0)
    yxx=scalProduct(tensor(sy,sx,sx),rho0)

    ms1=np.real(xxx-xyy-yxy-xyy)
    ms2=np.real(-yyy+yxx+xyx+yxx)
    mp1=np.real(xxx*xyy*yxy*xyy)
    mp2=np.real(yyy*yxx*xyx*yxx)

    return ms1,ms2,mp1,mp2


#  Parameter settings
w01=11*2*np.pi
w02=12*2*np.pi
w03=10*2*np.pi
Bac1=0.2*2*np.pi
Bac2=0.2*2*np.pi
Bac3=0.2*2*np.pi
J12=1*2*np.pi
J23=2*2*np.pi
g1=0.0
g2=0.0

fi=0
ff=20
ti=0
tf=10

frelist=np.linspace(fi,ff,101)*2*np.pi
tlist=np.linspace(ti,tf,101)

#  Find optimal frequency and time for each spin
print("\nSpin 1 frequency finding")
Bac1,Bac2,Bac3=ampControl(1)
f1,t1=freqFind(frelist,tlist,1)
print("\nSpin 2 frequency finding")
Bac1,Bac2,Bac3=ampControl(2)
f2,t2=freqFind(frelist,tlist,2)
print("\nSpin 3 frequency finding")
Bac1,Bac2,Bac3=ampControl(3)
f3,t3=freqFind(frelist,tlist,3)

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

#  Calculate mermin product of each state
mp1mat=np.zeros(len(rhoTime))
mp2mat=np.zeros(len(rhoTime))
bound=np.ones(len(rhoTime))/(-16)
for i in range(len(rhoTime)):
    ms1,ms2,mp1,mp2=merminProduct(rhoTime[i])
    mp1mat[i]=mp1
    mp2mat[i]=mp2

tl=np.linspace(0,t3,len(rhoTime))

#  Plot evolution of mermin product
fig,ax=plt.subplots()
ax.plot(tl,mp1mat,'r',label='Mp1')
ax.plot(tl,mp2mat,'b',label='Mp2')
ax.plot(tl,bound,'k',label='Threshold')
ax.set_xlabel(r'Time ($\mu s$)')
ax.set_ylabel(r'Mermin product')
ax.set_title('Mermin product evolution')
ax.legend(loc=0)
fig.savefig('Prob5_result.png')
plt.show()

