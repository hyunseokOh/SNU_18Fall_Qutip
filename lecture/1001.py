from qutip import *
import numpy as np
import matplotlib.pyplot as plt

def qubit_intedgrate(epsilon, delta, g1, g2, solver):
    H=epsilon/2.0*sigmaz()+delta/2.0*sigmax()

    c_ops=[]
    
    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    output=mesolve(H,psi0,tlist,c_ops,e_ops)
    return output.expect[0],output.expect[1],output.expect[2]
