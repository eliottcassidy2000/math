from pathlib import Path
ROOT = Path(__file__).resolve().parent
from itertools import permutations, combinations
from collections import defaultdict
from math import comb
from fractions import Fraction
import json
import numpy as np
import sympy as sp

def canon_tilings(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    pairs=list(combinations(range(n),2))
    allmask=(1<<len(pairs))-1
    index={p:k for k,p in enumerate(pairs)}
    def code(out,p):
        return sum(((out[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(pairs))
    cache={};keys=[];cm=[]
    for x in range(1<<len(free)):
        out=[0]*n
        for i in range(n-1):out[i]|=1<<(i+1)
        for k,(i,j) in enumerate(free):
            a,b=(i,j) if (x>>k)&1 else (j,i)
            out[a]|=1<<b
        c0=code(out,range(n))
        if c0 not in cache:
            orbit=[code(out,p) for p in permutations(range(n))]
            key=min(min(orbit),allmask-max(orbit))
            for c in orbit:cache[c]=key;cache[allmask^c]=key
        keys.append(cache[c0])
    codes=sorted(set(keys));fi=[codes.index(k) for k in keys]
    return free,codes,fi

def kraw(m,j,d):return sum((-1)**a*comb(d,a)*comb(m-d,j-a) for a in range(max(0,j-(m-d)), min(j,d)+1))

def build(n):
    free,codes,fi=canon_tilings(n);m=len(free);N=1<<m;k=len(codes)
    mass=np.bincount(fi,minlength=k).astype('int64')
    W=[np.zeros((k,k),dtype='int64') for d in range(m+1)]
    for x in range(N):
        for y in range(N):W[(x^y).bit_count()][fi[x],fi[y]]+=1
    Q=[sum(kraw(m,j,d)*W[d] for d in range(m+1)) for j in range(m+1)]
    f=(n-1)//2;q=(m-f)//2
    caps=[(comb(m,j)+sum(comb(f,a)*comb(q,(j-a)//2) for a in range(f+1) if 0<=j-a<=2*q and (j-a)%2==0))//2 for j in range(m+1)]
    # Independently compute Walsh Gram rather than only transform.
    grams=[np.zeros((k,k),dtype='int64') for j in range(m+1)]
    for s in range(N):
        z=np.zeros(k,dtype='int64')
        for x in range(N):z[fi[x]]+=(-1)**((x&s).bit_count())
        grams[s.bit_count()]+=np.outer(z,z)
    assert all(np.array_equal(Q[j],grams[j]) for j in range(m+1))
    ranks=[sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(v.tolist())).rank() for v in Q]
    assert all(r<=c for r,c in zip(ranks,caps))
    assert all(Q[j][u,u]%8==comb(m,j)%8 for j in range(m+1) for u in range(k) if mass[u]%2)
    D=sp.diag(*map(int,mass));wm=[sp.Matrix(v.tolist()) for v in W]
    leak=m*D+2*wm[2]-wm[1]*D.inv()*wm[1] if m>=2 else sp.zeros(k)
    local=np.zeros((N,k),dtype='int64')
    for x in range(N):
        for a in range(m):local[x,fi[x^(1<<a)]]+=1
    z=sp.Matrix(local.tolist());pp=sp.zeros(k)
    for u in range(k):
        rows=[x for x in range(N) if fi[x]==u]
        A=z[rows,:];s=sp.ones(1,len(rows))*A
        pp+=A.T*A-s.T*s/int(mass[u])
    assert pp==leak
    report={'n':n,'m':m,'codes':codes,'mass':mass.tolist(),'fibres':[[x for x in range(N) if fi[x]==u] for u in range(k)],'rank_Q':ranks,'rank_caps':caps,'leak_rank':sp.polys.matrices.DomainMatrix.from_Matrix(leak).rank(),'leak_trace':str(sp.trace(leak)),'leak_matrix':str(leak),'W':[v.tolist() for v in W],'Q':[v.tolist() for v in Q]}
    return report,W,Q,mass,fi

if __name__=='__main__':
    reports=[]
    for n in [3,4,5]:
        report,W,Q,mass,fi=build(n)
        print(n,{a:report[a] for a in ['mass','rank_Q','rank_caps','leak_rank','leak_trace']},flush=True)
        reports.append(report)
    json.dump(reports,open(ROOT/'audit.json','w'),indent=2)
