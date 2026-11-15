# Exact pair sets of REDUCIBLE tournaments of orders 11..14 from the exact census of orders <= 10:
# an ordinal sum A => B has (H, F2) = (H_A H_B, F_A F_B), F2 = P(2) = 2H + P2.
import numpy as np, glob, sys
V={n:np.loadtxt('pairs_%d.txt'%n,dtype=np.int64).reshape(-1,2) for n in range(2,11)}
V[1]=np.array([[1,0]],dtype=np.int64)
HF={n:np.unique(np.stack([V[n][:,0],2*V[n][:,0]+V[n][:,1]],1),axis=0) for n in V}   # (h, f)
def prod(A,B):
    h=(A[:,0][:,None]*B[:,0][None,:]).ravel(); f=(A[:,1][:,None]*B[:,1][None,:]).ravel()
    return np.unique(np.stack([h,f],1),axis=0)
R={}
for n in range(11,15):
    parts=[]
    for a in range(1,n):
        b=n-a
        A=HF[a] if a<=10 else R.get(a); B=HF[b] if b<=10 else R.get(b)
        if A is None or B is None: continue
        if len(A)*len(B)>4e7: continue
        parts.append(prod(A,B))
    R[n]=np.unique(np.concatenate(parts),axis=0)
    print('order %d: %d exact reducible (h, F2) pairs (blocks of order <= 10 or reducible)'%(n,len(R[n])))
np.save('reducible_hf.npy',np.array([R[n] for n in range(11,15)],dtype=object),allow_pickle=True)
