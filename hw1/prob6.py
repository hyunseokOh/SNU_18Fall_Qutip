from qutip import *
import numpy as np
import matplotlib.pyplot as plt

f=open("prob6_result.txt","w+")

matA=np.array([[1,1],[0,1]])/np.sqrt(3)
matB=np.array([[2,1],[-1,-2]])/np.sqrt(10)
matC=np.array([[2,1,-1],[-2,1,0],[1,0,1]])/np.sqrt(13)

_,a,_=np.linalg.svd(matA,True,True)
_,b,_=np.linalg.svd(matB,True,True)
_,c,_=np.linalg.svd(matC,True,True)

print("Singular values of case A : ",file=f)
print(a,file=f)
print(file=f)
print("Singular values of case B : ",file=f)
print(b,file=f)
print(file=f)
print("Singular values of case C : ",file=f)
print(c,file=f)


