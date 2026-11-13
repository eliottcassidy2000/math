"""Hostile search for the candidate H-H(delete r) >= indegree(r)*starts(r)."""
import random,json,gzip
import numpy as np
from inflation_probe import zprofile
from fibre_probe import HERE,save

def state(T):
    n=len(T);z=tuple(map(int,zprofile(np.array(T,dtype=np.int64),n-1)))
    d=sum(t>>(n-1)&1 for t in T)
    return z,d,z[0]-z[2]-d*z[1]

def run():
    rng=random.Random(61113);best=[];trials=0;bad=None
    for n in range(8,13):
        record=None
        for seed in range(25):
            T=[0]*n
            for a in range(n):
                for b in range(a+1,n):
                    if rng.randrange(2):T[a]|=1<<b
                    else:T[b]|=1<<a
            z,d,margin=state(T)
            for step in range(100):
                a,b=rng.sample(range(n),2);U=T.copy();U[a]^=1<<b;U[b]^=1<<a
                zz,dd,mm=state(U);trials+=1
                if record is None or mm*record['profile'][0]<record['margin']*zz[0]:record=dict(rows=U,profile=zz,indegree=dd,margin=mm)
                if mm<0:bad=record;break
                if mm*z[0]<=margin*zz[0] or rng.random()<0.015:T,z,d,margin=U,zz,dd,mm
            if bad:break
        best.append(dict(n=n,best=record));print('ORDER',n,'BEST',record,flush=True)
        if bad:break
    save('degree_probe',dict(status='REFUTED' if bad else 'FINITE-EXACT search only; not proved',universe='Deterministic25-start100-step one-arc hill search per order8..12; stops on first negative',evaluations=trials,best=best,counterexample=bad))

if __name__=='__main__':run()
