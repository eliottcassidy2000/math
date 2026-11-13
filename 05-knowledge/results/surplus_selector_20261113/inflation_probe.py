"""Common transitive-module inflations inside exact current/next-target fibres."""
from itertools import permutations
import json,gzip
import numpy as np
from numba import njit
from math import gcd,factorial
from fibre_probe import HERE,save,need

@njit
def end_table(T):
    n=len(T);D=np.zeros((1<<n,n),dtype=np.int64)
    for v in range(n):D[1<<v,v]=1
    for S in range(1,1<<n):
        for v in range(n):
            if not(S>>v&1):continue
            rest=S^(1<<v)
            if rest:
                for u in range(n):
                    if T[u]>>v&1:D[S,v]+=D[rest,u]
    return D

@njit
def zprofile(T,r):
    D=end_table(T);n=len(T);full=(1<<n)-1;H=D.sum(axis=1);H[0]=1
    reverse=np.zeros(n,dtype=np.int64)
    for v in range(n):reverse[v]=full^(1<<v)^T[v]
    s=end_table(reverse)[full,r];F=0
    for S in range(1,1<<n):
        if S>>r&1:F+=D[S,r]*H[full^S]
    return H[full],s,H[full^(1<<r)],F

def relabel(T,order):return tuple(sum(((T[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)

def align(A,B,r):
    n=len(A);best=None
    for p in permutations([v for v in range(n) if v!=r]):
        order=list(p);order.insert(r,r);U=relabel(B,order)
        flips=[(a,b) for a in range(n) for b in range(a+1,n) if ((A[a]^U[a])>>b)&1]
        key=(len(flips),len([e for e in flips if r in e]),tuple(order))
        if best is None or key<best[0]:best=(key,U,flips,order)
    return best[1:]

def inflate(T,w,m):
    vertices=[v for v in range(len(T)) if v!=w];labels=[(v,0) for v in vertices]+[(w,i) for i in range(m)]
    U=[]
    for v,i in labels:
        U.append(sum((i<j if v==u else ((T[v]>>u)&1))<<k for k,(u,j) in enumerate(labels)))
    return tuple(U),vertices

def run():
    source=json.loads(gzip.decompress((HERE/'fibre_probe.json.gz').read_bytes()));records=source['records'];out=[];survivors=[]
    for ai,bi in source['pure_deletion_repairs']:
        a,b=records[ai],records[bi];A=tuple(a['rows']);B,flips,order=align(A,tuple(b['rows']),a['root']);r=a['root'];tests=[]
        for w in range(len(A)):
            if w==r:continue
            states=[];passed=True
            for m in range(1,4):
                U,vs=inflate(A,w,m);V,_=inflate(B,w,m);rr=vs.index(r)
                x=tuple(map(int,zprofile(np.array(U,dtype=np.int64),rr)));y=tuple(map(int,zprofile(np.array(V,dtype=np.int64),rr)))
                states.append(dict(m=m,left=x,right=y,target=x[0]+x[3]-x[1],contents=[gcd(x[0]-x[1],x[2],x[3]),gcd(y[0]-y[1],y[2],y[3])]))
                if x[0]!=y[0] or x[3]-x[1]!=y[3]-y[1]:passed=False;break
            tests.append(dict(vertex=w,passed_three_sizes=passed,states=states))
            if passed:
                survivors.append(dict(ai=ai,bi=bi,vertex=w,left_rows=A,right_rows=B,root=r,states=states))
        out.append(dict(ai=ai,bi=bi,left=A,right=B,root=r,permutation=order,flipped_pairs=flips,tests=tests))
        print('PAIR',ai,bi,'FLIPS',flips,'SURVIVORS',[x['vertex'] for x in tests if x['passed_three_sizes']],flush=True)
    save('inflation_probe',dict(scope='All6 pure-deletion primitive repairs through7; minimum root-preserving arc alignment, then each nonroot vertex cloned at widths1..3 (early failure allowed); three-size survival alone is not a universal identity',alignments=out,survivors=survivors))
    for s in survivors:print('SURVIVOR',s['ai'],s['bi'],s['vertex'],s['states'])

if __name__=='__main__':run()
