from qutip import *
import numpy as np
import matplotlib.pyplot as plt

def dop(q):
    d=np.zeros([6,6])
    
    for i in range(3):
        for j in range(3):
            d[i+3][j]=clebsch(1,1,1,i-1,q,j-1)

    return Qobj(d).dag()

def pop():
    p=np.zeros([6,6])

    for i in range(3):
        p[i+3][i+3]=1
    
    return Qobj(p)

def dop_dir(direc):
    if direc=='x':
        d=-1.0/np.sqrt(2)*(-dop(1).dag()+dop(-1).dag())
    elif direc=='y':
        d=-1.0j/np.sqrt(2)*(dop(1).dag()+dop(-1).dag())
    else:
        d=dop(0).dag()
    return d.dag()

def ham(w,gamma,direc):
    d=dop_dir(direc)
    h=-w/2*(d+d.dag())
    
    heff=h-1j*gamma/2*pop()
    
    return heff


w=1*2*np.pi
gamma=0.2*2*np.pi
sqgam=np.sqrt(gamma)

Hz=ham(w,gamma,'y')

Utrans=Qobj(np.array([[-1,-np.sqrt(2)*1j,1,0,0,0],[np.sqrt(2),0,np.sqrt(2),0,0,0],[-1,np.sqrt(2)*1j,1,0,0,0],[0,0,0,-1,-np.sqrt(2)*1j,1],[0,0,0,np.sqrt(2),0,np.sqrt(2)],[0,0,0,-1,np.sqrt(2)*1j,1]]))/2
Hy=Utrans*Hz*Utrans.dag()

tlist=np.linspace(0,10,101)


psi0_z=basis(6,0)
psid_z=(basis(6,0)+basis(6,2))/np.sqrt(2)
psi0_y=Utrans*psi0_z
psid_y=Utrans*psid_z

c_ops=[sqgam*dop(-1),sqgam*dop(0),sqgam*dop(1)]
e_ops=[psid_z*psid_z.dag()]
e_ops2=[psid_y*psid_y.dag()]

ntraj=[1,5,25,100,1000]

for n in ntraj:
    output=mcsolve(Hz,psi0_z,tlist,c_ops,e_ops,ntraj=n)
    output2=mcsolve(Hy,psi0_y,tlist,c_ops,e_ops2,ntraj=n)
    fig,ax=plt.subplots()
    ax.plot(tlist,output.expect[0])
    ax.plot(tlist,output2.expect[0])
    ax.set_title("Evolution towards the dark state (N=%d)" % n)
    ax.set_xlabel("Time")
    ax.set_ylabel("Probability")
    ax.legend(loc=0,labels=["z-axis","y-axis"])
    fig.savefig("Result_n=%d.png" % n)


