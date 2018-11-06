from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import time
options=Options(nsteps=20000)

#  Calculate value of Hamiltonian
def hamiltonian_t(t,args):
    H0=args[0]
    H1=args[1]
    return H0+t*H1


#  Calculate eigenenergies of Hamiltonian
def compute(dlist):
    evals_mat=np.zeros((len(dlist),2))
    idx=0

    for e in dlist:
        H=-(delta/2.0)*sigmax()-(eps0/2.0)*sigmaz()-(e+1)*(A/2.0)*sigmaz()
        evals,ekets=H.eigenstates()
        evals_mat[idx,:]=np.real(evals)
        idx+=1

    return evals_mat


#  Calculate probability of excited state when applying pulse
def qubit_intergrate(delta,eps0,A,gamma1,gamma2,psi0,tlist,ep):
    
    sx=sigmax()
    sz=sigmaz()
    sm=destroy(2)
    
    tp=np.max(tlist)
    H0=delta/2.0*sx+eps0/2.0*sz
    H1=(ep+1)*(A/2.0)*sz

    c_op_list=[]
    n_th=0.0

    rate=gamma1*(1+n_th)
    if rate>0.0:
        c_op_list.append(np.sqrt(rate)*sm)

    rate=gamma1*n_th
    if rate>0.0:
        c_op_list.append(np.sqrt(rate)*sm.dag())

    rate=gamma2*ep
    if rate>0.0:
        c_op_list.append(np.sqrt(rate)*sm*sm.dag())

    H=[H0,[H1,'np.heaviside(t,0)*np.heaviside('+str(tp)+'-t,0)']]
    output=mesolve(H,psi0,tlist,c_op_list,[(sigmaz()+1)/2],options=options)
    
    return output.expect[0]


# Parameter settings
delta=1*2*np.pi     # qubit sigma_x coefficient
eps0=-10*2*np.pi    # qubit sigma_z coefficient
A=10*2*np.pi         
gamma1=0.1          # relaxation rate
gamma2=0.2          # dephasing rate
psi0=basis(2,0)     # initial state

trange=np.linspace(0,10,1000)
erange=np.linspace(-0.5,0.5,1000)
dlist=np.linspace(-5,5,1000)

#  Plot energy diagram
evals_mat=compute(dlist)
fig,ax=plt.subplots(figsize=(12,6))
ax.plot(dlist,(evals_mat)/(2*np.pi),'b')
ax.set_xlabel('Detuning')
ax.set_ylabel('Energy')
ax.set_title('Energy spectrum of Landau-Zener Hamiltonian')
fig.savefig('Prob4_energy')

#  Calculate probability for each point varying pulse width and detuning
res=np.zeros((len(erange),len(trange)))
timeLapse=0
counter=len(erange)
maxCount=counter
    
for j in range(len(erange)):
    start_time=time.time()
    ei=erange[j]
    p_ex=np.real(qubit_intergrate(delta,eps0,A,gamma1,gamma2,psi0,trange,ei))
    res[j]=1-p_ex
    
    counter-=1
    tlapse=time.time()-start_time
    timeRemain=tlapse*counter
    timeLapse+=tlapse
    h=timeRemain//3600
    m=(timeRemain-h*3600)//60
    s=(timeRemain-h*3600-m*60)
    h2=timeLapse//3600
    m2=(timeLapse-h2*3600)//60
    s2=(timeLapse-h2*3600-m2*60)

    print('Time consumed = %.2fs, Time lapsed = %.0fh %.0fm %.0fs, Time remaining = %.0fh %.0fm %.0fs, (%d/%d) left' % (tlapse,h2,m2,s2,h,m,s,counter,maxCount))

#  Plot probability 2D map
fig,ax=plt.subplots(figsize=(12,8))
t_mat,e_mat=np.meshgrid(trange,erange)
c=ax.pcolor(t_mat,e_mat,res,cmap='coolwarm')
ax.set_xlabel(r'Pulse width ($ns$)')
ax.set_ylabel(r'Pulse detuning ($\mu eV$)')
ax.set_title('Coherent oscillations')
fig.colorbar(c,ax=ax)
fig.savefig('Prob4_result.png')

plt.show()

