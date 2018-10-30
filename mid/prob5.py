from qutip import *
from qutip.ui.progressbar import TextProgressBar as ProgressBar
import numpy as np
import matplotlib.pyplot as plt

options=Options(nsteps=20000)

#  Calculate spin states
def calculate(rho0):
    #  rmat=np.zeros((len(tlist)),6)
    c_ops=[]

    if g1>0.0:
        c_ops.append(np.sqrt(g1)*sigmam())

    if g2>0.0:
        c_ops.append(np.sqrt(g2)*sigmaz())

    e_ops=[sigmax(),sigmay(),sigmaz()]
    
    wmw=frelist[0]
    args={'wmw':wmw}
    #  H0=(w01/2.0)*tensor(sigmaz(),qeye(2))+(w02/2.0)*tensor(qeye(2),sigmaz())+(J12/2.0)*tensor(sigmaz(),sigmaz())
    #  H1=Bac*tensor(sigmax(),qeye(2))
    #  H=[H0,[H1,'np.cos(wmw*t)']]
    sx=sigmax()
    sy=sigmay()
    sz=sigmaz()
    idm=qeye(2)
    H0=(w01-wmw)/2*tensor(sz,idm,idm)+(w02-wmw)/2*tensor(idm,sz,idm)+(w03-wmw)/2*tensor(idm,idm,sz)+J12/2*tensor(sz,sz,idm)+J23/2*tensor(idm,idm,sz)
    Hint=(Bac1/4)*tensor(sx,idm,idm)+(Bac2/4)*tensor(idm,sx,idm)+(Bac3/4)*tensor(idm,idm,sx)
    H=H0+Hint
    output=mesolve(H,rho0,tlist,c_ops,[tensor(sx,idm,idm),tensor(sy,idm,idm),tensor(sz,idm,idm),tensor(idm,sx,idm),tensor(idm,sy,idm),tensor(idm,sz,idm),tensor(idm,idm,sx),tensor(idm,idm,sy),tensor(idm,idm,sz)],args,options=options)

    return output.expect

#  Draw Bloch Spheres for each spin
def bloch(rmat):
    b1=Bloch()
    b2=Bloch()
    b3=Bloch()

    for i in range(len(tlist)):
        vec1=[rmat[0][i],rmat[1][i],rmat[2][i]]
        vec2=[rmat[3][i],rmat[4][i],rmat[5][i]]
        vec3=[rmat[6][i],rmat[7][i],rmat[8][i]]
        b1.add_vectors(vec1)
        b2.add_vectors(vec2)
        b3.add_vectors(vec3)    
    b1.show()
    b2.show()
    b3.show()

#  Return final vectors
def mat2vec(rmat):
    #  vec1_i=[rmat[0][0],rmat[1][0],rmat[2][0]]
    vec1_f=[rmat[0][-1],rmat[1][-1],rmat[2][-1]]
    #  vec2_i=[rmat[3][0],rmat[4][0],rmat[5][0]]
    vec2_f=[rmat[3][-1],rmat[4][-1],rmat[5][-1]]
    vec3_f=[rmat[6][-1],rmat[7][-1],rmat[8][-1]]
    
    return vec1_f,vec2_f,vec3_f


#  Transform (x,y,z) coordinates to Qobj
def transform(vec):
    #  Ceiling in case of larger than 1
    for i in range(len(vec)):
        if vec[i]>1:
            vec[i]=1
        if vec[i]<-1:
            vec[i]=-1

    theta=np.arccos(vec[2])
    
    if(np.sin(theta)==0):
        phi=0
    else:
        if vec[0]/np.sin(theta)>1:
            phi=0
        elif vec[0]/np.sin(theta)<-1:
            phi=np.pi
        else :
            phi=np.arccos(vec[0]/np.sin(theta))
    
    if vec[1]<0:
        phi=phi+np.pi
    
    a=np.cos(theta/2)
    b=(np.cos(phi)+1j*np.sin(phi))*np.sin(theta/2)
    
    return a*basis(2,0)+b*basis(2,1)


#  Parameter settings
w01=10*2*np.pi
w02=20*2*np.pi
w03=20*2*np.pi
Bac1=0*2*np.pi
Bac2=1*2*np.pi
Bac3=1*2*np.pi
J12=10*2*np.pi
J23=10*2*np.pi
g1=0.0
g2=0.0
tlist=np.linspace(0,2*np.pi/Bac2,2)
frelist=[10*2*np.pi]


psi1=basis(2,0)
psi2=basis(2,0)
psi3=basis(2,0)

rho0=ket2dm(tensor(psi1,psi2,psi3))

rmat=calculate(rho0)
bloch(rmat)

vec1,vec2,vec3=mat2vec(rmat)
v1=transform(vec1)
v2=transform(vec2)
v3=transform(vec3)

ideal=cnot()*tensor(psi1,psi2,psi3)
calc=ket2dm(tensor(v1,v2,v3))
print(calc)
print("Fidelity for single CNOT : "+str(fidelity(ideal,calc)))



