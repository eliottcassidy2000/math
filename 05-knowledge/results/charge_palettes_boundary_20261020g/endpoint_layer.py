"""Closed (n-3,2) layer, exact ports, literal three-K ledgers.

Basic tournament operations follow the inherited detours_check.py convention.
No isomorphism, clone contraction or unrecorded relabelling is used.
"""
from itertools import combinations
from collections import Counter
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)
def arc(t,a,b):return (t[a]>>b)&1
def blue(t,p):
    n=len(t)
    return all(arc(t,p[i],p[j])==arc(t,p[n-1-j],p[n-1-i]) for i,j in combinations(range(n),2))
def path(t,p):return len(p)==len(t) and len(set(p))==len(t) and all(arc(t,a,b) for a,b in zip(p,p[1:]))
def move(t,p):
    need(path(t,p),'path directed before K')
    n=len(t);out=[((1<<n)-1)^r^(1<<i) for i,r in enumerate(t)]
    for a,b in zip(p,p[1:]):out[a]^=1<<b;out[b]^=1<<a
    return tuple(out)
def presentations(r):
    free={(a,b) for a in range(r) for b in range(a+2,r)}; orbits=[]
    while free:
        a,b=min(free); o={(a,b),(r-1-b,r-1-a)};free-=o;orbits.append(o)
    for mask in range(1<<len(orbits)):
        t=[0]*r
        for i in range(r-1):t[i]|=1<<(i+1)
        for j,o in enumerate(orbits):
            for a,b in o:
                if mask>>j&1:t[a]|=1<<b
                else:t[b]|=1<<a
        yield tuple(t)
def wrap(v,j):
    r=len(v);z=r+1;t=[sum(1<<i for i in range(1,r+1))]
    t += [(row<<1)|(1<<z) for row in v];t += [1]
    b=r+1-j
    for x,y in ((0,j),(b,z)):t[x]^=1<<y;t[y]^=1<<x
    return tuple(t)
def ledger(t,words):
    now=t
    for p in words:
        if not path(now,p):return 'undirected'
        if blue(now,p):return 'blue'
        now=move(now,p)
    need(now==move(t,tuple(range(len(t)))),'literal marked target')
    return 'paid'
def proposed(v,j):
    r=len(v);a,z=0,r+1;M=tuple(range(1,r+1));b,c=1,r
    if arc(v,0,r-1):
        W=M[1:]
        return 'open',((b,)+W+(z,a),(z,)+W+(b,a),(a,z)+W+(b,))
    W=M[1:-1]
    if not arc(v,0,r-2):
        return 'closed_no_chord',((c,)+W+(b,z,a),(z,c)+W+(b,a),(a,z)+M)
    return 'closed_chord',((c,b,z,a)+W,(z,c,b,a)+W,(a,z)+M)

def main():
    rows=[]; residual=[]; low_rotations=[]
    for r in range(3,9):
        stats=Counter()
        for v in presentations(r):
            for j in range(2,r+1):
                t=wrap(v,j);P=tuple(range(r+2));need(blue(t,P),'marked blue path')
                need(t[0].bit_count()==r-1 and t[-1].bit_count()==2,'endpoint scores')
                cut=P[1:]+P[:1]
                if r>=5:need(not blue(t,cut),'low-degree rotation escape')
                elif blue(t,cut):low_rotations.append({'r':r,'j':j,'rows':t})
                branch,words=proposed(v,j);result=ledger(t,words)
                stats['total']+=1;stats[branch+'_'+result]+=1
                guaranteed=r>=5 and ((branch=='open' and j!=r-1) or (branch=='closed_no_chord' and j!=r))
                if guaranteed:
                    need(result=='paid','new uniform theorem');stats['theorem_covered']+=1
                if result!='paid' and len(residual)<20:residual.append({'r':r,'j':j,'branch':branch,'failure':result,'middle':v})
        rows.append({'r':r,'n':r+2,'counts':dict(stats)})
    print(json.dumps({'status':'FINITE-EXACT presentations; all-order proofs separate',
        'rows':rows,'small_blue_cyclic_cut_hostiles':low_rotations,'retained_residuals':residual},indent=2))
if __name__=='__main__':main()
