from qutip import *
import numpy as np
import matplotlib.pyplot as plt

A=1
H=A*sigmaz()
istate=(2*basis(2,0)+basis(2,1))/np.sqrt(5)
idm=ket2dm(istate)

t=np.linspace(0,3,1000)

output=mesolve(H,idm,t,[],[sigmax(),sigmay(),sigmaz()])
point=[output.expect[0],output.expect[1],output.expect[2]]

b=Bloch()
b.add_points(point)
b.show()
