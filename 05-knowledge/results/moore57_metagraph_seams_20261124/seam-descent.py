#!/usr/bin/env python3
"""Exact seam score identities and a finite descent calibration on degree seven."""
from collections import Counter
from itertools import combinations, permutations
import numpy as np
from scipy.optimize import linear_sum_assignment

def defect(B):
    d,_,n,_=B.shape
    out=np.empty_like(B)
    for i in range(d):
        for j in range(d):
            out[i,j] = sum(B[i,k]@B[k,j] for k in range(d))+B[i,j]-(1 if i!=j else n*np.eye(n,dtype=int))
    return out

def energy(F):
    return int(np.sum(F*F))

def install(B,i,j,Q):
    B[i,j]=Q
    B[j,i]=Q.T

def cost(B,F,i,j):
    d=B.shape[0]
    return F[i,j]+sum(F[i,k]@B[j,k].T+B[i,k]@F[j,k].T
                       for k in range(d) if k not in (i,j))-(2*d-3)*B[i,j]

def best(W,derange):
    n=len(W)
    E=np.eye(n,dtype=int)
    C=W+((2*n*np.max(np.abs(W))+1)*E if derange else 0)
    rng,t=linear_sum_assignment(C)
    assert np.array_equal(rng,np.arange(n))
    assert not derange or np.all(rng!=t)
    return E[t]

def hs():
    V=list(np.ndindex(2,5,5)); N={v:set() for v in V}
    for v,w in combinations(V,2):
        t,i,j=v;s,k,l=w
        if (t==s and i==k and (j-l)%5 in ({1,4} if t==0 else {2,3})) or (t==0 and s==1 and l==(i*k+j)%5):
            N[v].add(w);N[w].add(v)
    r=V[0];parent=sorted(N[r]); L=[sorted(N[x]-{r}) for x in parent]
    for i in range(1,7):
        L[i]=[next(y for y in L[i] if y in N[x]) for x in L[0]]
    B=np.zeros((7,7,6,6),dtype=int)
    for i,j in combinations(range(7),2):
        P=np.array([[int(y in N[x]) for y in L[j]] for x in L[i]])
        assert np.all(P.sum(axis=0)==1) and np.all(P.sum(axis=1)==1)
        install(B,i,j,P)
    return B

def cycle_counts(B):
    d,_,n,_=B.shape; A=B.transpose(0,2,1,3).reshape(d*n,d*n)
    C=A@A
    return int(np.sum(C*A)//6),int(np.sum((C*(C-1)-np.diag(np.diag(C*(C-1))))//1)//8)

def main():
    rng=np.random.default_rng(571)
    known=hs();assert energy(defect(known))==0
    eps=np.eye(6,dtype=int); P=known[1,2].copy(); Q=P[[1,0,2,3,4,5]]
    B=known.copy();install(B,1,2,Q)
    assert energy(defect(B))==88 and cycle_counts(B)==(2,8)
    print('HS single transposition: energy 88, 2 triangles, 8 quadrangles')
    trials=[]
    seams=list(combinations(range(1,7),2))
    P_legal=[eps[np.array(p)] for p in permutations(range(6)) if all(i!=j for i,j in enumerate(p))]
    for r in range(32):
        B=known.copy()
        for i,j in seams:
            install(B,i,j,P_legal[rng.integers(len(P_legal))])
        F=defect(B);e0=energy(F); e=e0;moves=0
        if r==0:
            i,j=seams[0];W=cost(B,F,i,j);P=B[i,j].copy()
            clone=B.copy()
            for p in permutations(range(6)):
                Q=eps[np.array(p)]
                install(clone,i,j,Q)
                assert energy(defect(clone))-e==4*np.sum(W*(Q-P))
            Q=best(W,True)
            assert np.sum(W*Q)==min(np.sum(W*T) for T in P_legal)
        while True:
            choices=[]
            for i,j in seams:
                W=cost(B,F,i,j);Q=best(W,True)
                delta=int(4*np.sum(W*(Q-B[i,j])))
                if delta<0:
                    choices.append((delta,i,j,Q))
            if not choices:break
            delta,i,j,Q=min(choices,key=lambda c:c[:3])
            install(B,i,j,Q)
            F=defect(B)
            assert energy(F)==e+delta
            e+=delta;moves+=1
        trials.append((e0,e,moves))
    print('32 fixed-seed terminal scores:',dict(sorted(Counter(e for _,e,_ in trials).items())))
    print('min terminal',min(e for _,e,_ in trials),'max moves',max(m for _,_,m in trials))
if __name__=='__main__':
    main()
