from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import time
options=Options(nsteps=20000)

#  Calculate expectation value of Sz or state when applying pulse
def qubit_intergrate(delta,eps0,A,gamma1,gamma2,psi0,tlist,ep,state):
    
    sx=sigmax()
    sy=sigmay()
    sz=sigmaz()
    sm=destroy(2)
    
    tp=np.max(tlist)
    H0=delta/2.0*sy+eps0/2.0*sz
    H1=(ep+1)*(A/2.0)*sz

    c_op_list=[]
    e_op_list=[]
    
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
    
    if state==0:
        e_op_list=[sz]

    H=[H0,[H1,'np.heaviside(t,0)*np.heaviside('+str(tp)+'-t,0)']]
    output=mesolve(H,psi0,tlist,c_op_list,e_op_list,options=options)
    
    if state==0:
        return output.expect[0]
    else:
        return output.states


#  Implementation of Hadamard gate
def Hadamard_impl(psi0):
    p_ex=np.real(qubit_intergrate(delta,eps0,A,gamma1,gamma2,psi0,trange,ep,0))
    opt=np.where(np.abs(p_ex)==np.min(np.abs(p_ex)))
    t_opt=trange[opt]

    tlist=[0,t_opt]
    calc=qubit_intergrate(delta,eps0,A,gamma1,gamma2,psi0,tlist,ep,1)[-1]
    return calc   


# Parameter settings
delta=1*2*np.pi     # qubit sigma_x coefficient
eps0=-10*2*np.pi    # qubit sigma_z coefficient
A=10*2*np.pi         
gamma1=0.1          # relaxation rate
gamma2=0.2          # dephasing rate
psi0=basis(2,0)     # initial state

trange=np.linspace(0,0.5,1000)
ep=0

#  Calculate fidelity
calc=Hadamard_impl(psi0)
ideal=ket2dm(hadamard_transform()*psi0)

print(calc)
print(ideal)
print("Fidelity : %f" % fidelity(calc,ideal))
