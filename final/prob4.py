from qutip import *
from qutip.ui.progressbar import TextProgressBar as ProgressBar
import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=20000)

#  Calculate expectation value of Sz for each spin
def freqCalc(rho0, frelist, tlist,target):
    wmw1 = 0
    wmw2 = 0
    wmw3 = 0

    rmat = np.zeros((len(frelist),len(tlist)))
    pbar = ProgressBar(len(frelist))
    c_ops = []

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    sy=sigmay()
    sz=sigmaz()
    idm=qeye(2)
    
    for i in range(len(frelist)):
        #  pbar.update(i)
        
        if target == 1:
            e_ops=[tensor(sz,idm,idm)]
            wmw1=frelist[i]
        
        elif target==2:
            e_ops=[tensor(idm,sz,idm)]
            wmw2=frelist[i]

        else:
            e_ops=[tensor(idm,idm,sz)]
            wmw3=frelist[i]

        H0=(w01-wmw1)/2*tensor(sz,idm,idm)+(w02-wmw2)/2*tensor(idm,sz,idm)+(w03-wmw3)/2*tensor(idm,idm,sz)+J12*tensor(sz,sz,idm)+J23*tensor(idm,sz,sz)
        Hint=(Bac1/4)*tensor(sy,idm,idm)+(Bac2/4)*tensor(idm,sy,idm)+(Bac3/4)*tensor(idm,idm,sy)
        H=H0+Hint
        
        output=mesolve(H,rho0,tlist,c_ops,e_ops,options=options)
        rmat[i]=output.expect[0]
    return rmat


#  Find optimal frequency and time for each spin
def freqFind(frelist,tlist,target):
    rho1=ket2dm(tensor(basis(2,0),basis(2,0),basis(2,0)))
    rmat1=freqCalc(rho1,frelist,tlist,target)
    avoid=np.where(rmat1==np.min(rmat1))

    f_avoid=frelist[avoid[0][0]]
    t_avoid=tlist[avoid[1][0]]
    
    #  Find the position of <Sz>=-1
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
    opt=np.where(rmat2==np.min(rmat2))
    
    f_opt=frelist[opt[0][0]]
    t_opt=tlist[opt[1][0]]
    
    print("Optimal frequency : %f, Optimal time : %f" % (f_opt/2/np.pi,t_opt))
    print("Avoid frequency : %f, Optimal time : %f\n" % (f_avoid/2/np.pi,t_avoid)) 

    return (f_opt,t_opt)


def freqFind2(frelist,tlist,control1,control2,target):
    rho1=ket2dm(tensor(basis(2,0),basis(2,0),basis(2,0)))
    rmat1=freqCalc(rho1,frelist,tlist,target)
    avoid=np.where(rmat1==np.min(rmat1))

    f_avoid=frelist[avoid[0][0]]
    t_avoid=tlist[avoid[1][0]]

    #  Find the position of <Sz>=-1
    
    if control1 == 0:
        pcon1=basis(2,0)
    else:
        pcon1=basis(2,1)

    if control2 == 0:
        pcon2=basis(2,0)
    else:
        pcon2=basis(2,1)


    rho2=ket2dm(tensor(pcon1,basis(2,0),pcon2))
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
  
    print("Optimal frequency : %f, Optimal time : %f" % (f_opt/2/np.pi,t_opt))
    print("Avoid frequency : %f, Optimal time : %f\n" % (f_avoid/2/np.pi,t_avoid)) 
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
    elif target==2:
        wmw2=f_opt
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


#  Implementation of toffoli gate
def toffoli_impl(rho0,f_opt,t_opt,target):
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
    elif target==2:
        wmw2=f_opt
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
def ampControl(tarList):
    Bac=[0,0,0]
    val=0.2*2*np.pi

    for i in tarList:
        Bac[i-1]=val

    return Bac[0],Bac[1],Bac[2]


#  Measuring process
def measure(psi0,op,alpha):

    #  Quantum Circuit
    #  Initial condition
    psi1=basis(2,0)
    psi2=psi0
    psi3=basis(2,0)

    rho0=ket2dm(tensor(psi1,psi2,psi3))

    op_meas=tensor(qeye(2),op,qeye(2))

    #  Encoding
    Bac1,Bac2,Bac3=ampControl([1])
    rho1=CNOT_impl(rho0,f1,t1,1)[-1]

    Bac1,Bac2,Bac3=ampControl([3])
    rho1=CNOT_impl(rho1,f3,t3,3)[-1]

    #  Error
    U=tensor(qeye(2),rx(alpha),qeye(2))
    rho2=U*rho1*U.dag()

    #  Decoding
    Bac1,Bac2,Bac3=ampControl([1])
    rho3=CNOT_impl(rho2,f1,t1,1)[-1]

    Bac1,Bac2,Bac3=ampControl([3])
    rho3=CNOT_impl(rho3,f3,t3,3)[-1]

    exp_err=np.real((rho3*op_meas).tr())/2

    #  Correction
    Bac1,Bac2,Bac3=ampControl([2])
    rho4=toffoli_impl(rho3,f2,t2,2)[-1]
    
    #  Measure
    exp_corr=np.real((rho4*op_meas).tr())/2
    return exp_err,exp_corr


#  Parameter settings
w01=11*2*np.pi
w02=12*2*np.pi
w03=9*2*np.pi
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
print("\nSpin 1 frequency finding\n")
Bac1,Bac2,Bac3=ampControl([1])
f1,t1=freqFind(frelist,tlist,1)
print("\nSpin 2 frequency finding\n")
Bac1,Bac2,Bac3=ampControl([2])
f2,t2=freqFind2(frelist,tlist,1,1,2)
#  freqFind2(frelist,tlist,1,0,2)
#  freqFind2(frelist,tlist,0,1,2)
#  freqFind2(frelist,tlist,0,0,2)
print("\nSpin 3 frequency finding\n")
Bac1,Bac2,Bac3=ampControl([3])
f3,t3=freqFind(frelist,tlist,3)

f1,t1=9*2*np.pi,5.0
f2,t2=6*2*np.pi,5.0
f3,t3=5*2*np.pi,5.0


sx=sigmax()
sy=sigmay()
sz=sigmaz()
s0=basis(2,0)
s1=basis(2,1)

z_ket=s0
mz_ket=s1
x_ket=1/np.sqrt(2)*(s0+s1)
mx_ket=1/np.sqrt(2)*(s0-s1)
y_ket=1/np.sqrt(2)*(s0+1j*s1)
my_ket=1/np.sqrt(2)*(s0-1j*s1)

alpLst=np.linspace(0,np.pi,5)
pLst=(np.sin(alpLst/2))**2
fid_ideal=1-pLst
fid=np.zeros(len(alpLst))
fid2=np.zeros(len(alpLst))

for i in range(len(alpLst)):
    alpha=alpLst[i]
    r_z_z,r_z_z2=measure(z_ket,sz,alpha)
    r_mz_z,r_mz_z2=measure(mz_ket,sz,alpha)
    r_x_x,r_x_x2=measure(x_ket,sx,alpha)
    r_mx_x,r_mx_x2=measure(mx_ket,sx,alpha)
    r_y_y,r_y_y2=measure(y_ket,sy,alpha)
    r_my_y,r_my_y2=measure(my_ket,sy,alpha)

    fid[i]=(1+r_z_z-r_mz_z+r_x_x-r_mx_x+r_y_y-r_my_y)/4
    fid2[i]=(1+r_z_z2-r_mz_z2+r_x_x2-r_mx_x2+r_y_y2-r_my_y2)/4
    print("Calculation ongoing (%d/%d)" % (i+1,len(alpLst)))
    print(fid[i],fid2[i])

fig,ax=plt.subplots()
ax.plot(pLst,fid_ideal)
ax.plot(pLst,fid)
ax.plot(pLst,fid2)
plt.show()
