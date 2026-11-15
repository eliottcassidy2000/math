import numpy as np
import sys as _s
tok=open(_s.argv[1] if len(_s.argv)>1 else '../../../04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt').read().split(); pos=0
def nx():
    global pos; v=int(tok[pos]); pos+=1; return v
dim=nx(); nl=nx(); init=np.array([nx() for _ in range(dim)],dtype=object); fin=np.array([nx() for _ in range(dim)],dtype=object)
M={}
for _ in range(nl):
    li=nx(); M[li]=np.array([[nx() for _ in range(dim)] for _ in range(dim)],dtype=object)
def H(word):
    v=init.copy()
    for a in word: v=M[a].dot(v)
    return int(fin.dot(v))
for a in range(8):
    ev=np.max(np.abs(np.linalg.eigvals(np.array(M[a].tolist(),dtype=float))))
    hs=[H((a,)*k) for k in range(1,16)]
    print('letter %d: spectral radius %.5f; H(a^k), k=1..15: %s'%(a,ev,hs))
