"""Hostile search for the candidate H-H(delete r) >= indegree(r)*starts(r)."""
import random,json,gzip
import numpy as np
from inflation_probe import zprofile
from fibre_probe import HERE,save,need

def state(T):
    n=len(T);z=tuple(map(int,zprofile(np.array(T,dtype=np.int64),n-1)))
    d=sum(t>>(n-1)&1 for t in T)
    return z,d,z[0]-z[2]-d*z[1]

def run():
    rng=random.Random(61113);best=[];trials=0;bad=None
    rows=json.loads(gzip.decompress((HERE/'fibre_probe.json.gz').read_bytes()))['records']
    small=[]
    for row in rows:
        T=row['rows'];h,s,q,F=row['profile'];d=sum(t>>(len(T)-1)&1 for t in T)
        if h-q-d*s<0:small.append(row)
    hostile=(2,20,25,19,1);hz,hd,hm=state(hostile)
    need(hz==(9,1,5,23) and hd==3 and hm==1,'cyclic-three in-neighbourhood is not a chart')
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
    save('degree_probe',dict(status='REFUTED' if bad or small else 'OPEN; exact bounded and sampled controls only',universe='Complete2618 rooted strong classes through7, then deterministic25-start100-step one-arc hill search per order8..12; stops on first negative',small_classes=len(rows),small_violations=small,evaluations=trials,best=best,counterexample=bad,cyclic_three_hostile=dict(rows=hostile,root=4,profile=hz,indegree=hd,margin=hm)))

if __name__=='__main__':run()
