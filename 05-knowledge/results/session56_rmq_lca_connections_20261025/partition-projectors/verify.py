#!/usr/bin/env python3
"""Exact finite checks; no assertions, no nonstandard dependencies."""
from fractions import Fraction
from itertools import product
import json, pathlib, hashlib, sys
ROOT=pathlib.Path(__file__).resolve().parent

def need(x,msg):
    if not x: raise RuntimeError(msg)

def partitions(n):
    def rec(a):
        if len(a)==n:
            yield tuple(a); return
        for x in range(max(a,default=-1)+2): yield from rec(a+[x])
    return list(rec([]))

def mask(p):
    return sum(1<<k for k,(i,j) in enumerate((i,j) for i in range(len(p)) for j in range(i+1,len(p))) if p[i]==p[j])

def setup(n):
    ps=sorted(partitions(n),key=lambda p:(mask(p).bit_count(),mask(p)))
    ms=list(map(mask,ps)); ix={m:i for i,m in enumerate(ms)}; N=len(ps)
    le=lambda i,j:ms[i]&ms[j]==ms[i]
    meet=[[ix[a&b] for b in ms] for a in ms]
    mu=[[0]*N for _ in ps]
    for i in range(N):
        mu[i][i]=1
        for j in range(i+1,N):
            if le(i,j): mu[i][j]=-sum(mu[i][k] for k in range(i,j) if le(i,k) and le(k,j))
    def z(f): return [sum(f[j] for j in range(N) if le(i,j)) for i in range(N)]
    def conv(f,g):
        h=[0]*N
        for i,x in enumerate(f):
            if x:
                for j,y in enumerate(g):
                    if y: h[meet[i][j]]+=x*y
        return h
    es=[[mu[i][j] for i in range(N)] for j in range(N)]
    return ps,ms,z,conv,es

def compositions(total,n):
    if n==1: yield (total,); return
    for x in range(total+1):
        for tail in compositions(total-x,n-1): yield (x,)+tail

def main():
    if '--self-test-failure' in sys.argv:
        need(False,'intentional rejection guard'); return
    report={'schema':'partition-projectors-v1','arithmetic':'exact integers and Fraction','partition_checks':[]}
    for n in range(1,6):
        ps,ms,z,conv,es=setup(n); N=len(ps)
        for i,e in enumerate(es):
            need(z(e)==[int(i==j) for j in range(N)],'zeta projector')
            for j,d in enumerate(es): need(conv(e,d)==(e if i==j else [0]*N),'orthogonality')
        need([sum(e[i] for e in es) for i in range(N)]==[0]*(N-1)+[1],'identity sum')
        for i in range(N):
            f=[int(k==i) for k in range(N)]
            for j in range(N):
                g=[int(k==j) for k in range(N)]
                need(z(conv(f,g))==[a*b for a,b in zip(z(f),z(g))],'basis diagonalization')
        report['partition_checks'].append({'n':n,'states':N,'projector_pairs':N*N,'basis_pairs':N*N})
    ps,ms,z,conv,es=setup(3); N=len(ps)
    nonnegative=[]
    for bits in product((0,1),repeat=N):
        f=[sum(bits[j]*es[j][i] for j in range(N)) for i in range(N)]
        need(conv(f,f)==f,'all binary coordinates idempotent')
        if min(f)>=0: nonnegative.append(f)
    need(len(nonnegative)==N+1,'only zero and Dirac positive idempotents')
    report['pi3']={'partitions':ps,'primitive_projectors':es,'idempotents':2**N,'nonnegative_idempotents':nonnegative}
    checked=0; stable=0
    for D in range(1,9):
        for nums in compositions(D,N):
            f=[Fraction(x,D) for x in nums]; h=conv(f,f)
            need(z(h)==[v*v for v in z(f)],'rational diagonalization')
            is_dirac=sum(x!=0 for x in nums)==1
            need((h==f)==is_dirac,'rational probability classification')
            stable+=int(h==f); checked+=1
    report['rational_grid']={'denominators':[1,8],'histograms_counted_with_repetitions':checked,'stable_histograms_counted_with_repetitions':stable}
    a,b=1,2
    f=[Fraction(int(i in (a,b)),2) for i in range(N)]
    h=conv(f,f)
    need(h[0]==Fraction(1,2) and h!=f,'independent mixture hostile')
    correlated=[sum(f[j] for j in range(N) if (ms[j]&ms[j])==ms[i]) for i in range(N)]
    need(correlated==f,'correlated diagonal control')
    report['mixture_hostile']={'p':list(map(str,f)),'independent_meet':list(map(str,h)),'correlated_meet':list(map(str,correlated))}
    # Uniform law on C2, group multiplication (addition mod 2), not meet.
    u=[Fraction(1,2)]*2
    cu=[sum(u[i]*u[j] for i in range(2) for j in range(2) if (i+j)%2==k) for k in range(2)]
    need(cu==u,'group countercontrol')
    report['group_operation_caveat']='uniform C2 is non-Dirac idempotent under group convolution, not under meet'
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
    for item in pins:
        need(hashlib.sha256((ROOT/item['snapshot']).read_bytes()).hexdigest()==item['sha256'],'snapshot source pin')
    report['source_snapshots_verified']=len(pins)
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
