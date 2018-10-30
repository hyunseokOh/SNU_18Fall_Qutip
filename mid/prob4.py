from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import time

def hamiltonian_t(t,args):
    H0=args[0]
    H1=args[1]

    return H0+t*H1

def compute(tlist):
    evals_mat=np.zeros((len(tlist),2))
    idx=0

    for t in tlist:

        H=-(delta/2.0)*sigmax()-(eps0/2.0)*sigmaz()-t*(A/2.0)*sigmaz()

        evals,ekets=H.eigenstates()
        evals_mat[idx,:]=np.real(evals)
        idx+=1

    return evals_mat

options=Options(nsteps=20000)

def qubit_intergrate(delta,eps0,A,gamma1,gamma2,psi0,tlisti,tp,ep):

    sx=sigmax()
    sz=sigmaz()
    sm=destroy(2)

    H0=delta/2.0*sx+eps0/2.0*sz
    H1=((A+ep*2*np.pi)/2.0)*sz

    c_op_list=[]
    n_th=0.0

    rate=gamma1*(1+n_th)
    if rate>0.0:
        c_op_list.append(np.sqrt(rate)*sm)

    rate=gamma1*n_th
    if rate>0.0:
        c_op_list.append(np.sqrt(rate)*sm.dag())

    rate=gamma2
    if rate>0.0:
        c_op_list.append(np.sqrt(rate)*sm*sm.dag())

    H=[H0,[H1,'np.heaviside(t,0)*np.heaviside('+str(tp)+'-t,0)']]
    output=mesolve(H,psi0,tlist,c_op_list,[(sigmaz()+1)/2],{},options=options)

    return output.expect[0]


# Parameter settings
delta=1*2*np.pi     # qubit sigma_x coefficient
eps0=-10*2*np.pi    # qubit sigma_z coefficient
A=10*2*np.pi         
gamma1=0.01          # relaxation rate
gamma2=0.0          # dephasing rate
psi0=basis(2,0)     # initial state

trange=np.linspace(0,10,10)
erange=np.linspace(-5,5,10)

tlist=np.linspace(0,2,100)

evals_mat=compute(tlist)
fig,ax=plt.subplots(figsize=(12,6))
ax.plot(tlist,(evals_mat)/(2*np.pi),'b')
ax.set_xlabel('Time')
ax.set_ylabel('Eigenenergies')
ax.set_title('Energy spectrum of Landau-Zener Hamiltonian')

res=np.zeros((len(erange),len(trange)))

singletime=2.7
timeRemain=singletime*len(trange)*len(erange)

for i in range(len(trange)):
    ti=trange[i]
    tlist=np.linspace(0,ti,2)

    
    for j in range(len(erange)):
        start_time=time.time()
        ei=erange[j]
        p_ex=np.real(qubit_intergrate(delta,eps0,A,gamma1,gamma2,psi0,tlist,ti,ei))
        res[j,i]=1-p_ex[-1]
        
        tlapse=time.time()-start_time
        timeRemain-=tlapse
        print('time comsumed ='+str(tlapse)+', Time remaining = '+str(timeRemain))

fig,ax=plt.subplots(figsize=(12,8))
t_mat,e_mat=np.meshgrid(trange,erange)
c=ax.pcolor(t_mat,e_mat,res,cmap='RdBu')
ax.set_xlabel(r'time ($\mu$s)')
ax.set_ylabel(r'Energy $eV$')
ax.set_title('Landau-Zener transition')
fig.colorbar(c,ax=ax)
plt.show()

