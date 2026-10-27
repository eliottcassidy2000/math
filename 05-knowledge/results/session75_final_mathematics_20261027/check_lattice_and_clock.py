#!/usr/bin/env python3
"""Exact small Paley-field and clocked-occupation checks; no data files needed."""
from itertools import product


def need(ok, message):
    if not ok:
        raise ValueError(message)


# F27 = F3[t]/(t^3-t-1), an irreducible cubic (it has no F3 root).
F = list(product(range(3), repeat=3))

def add(a,b):
    return tuple((x+y)%3 for x,y in zip(a,b))

def neg(a):
    return tuple(-x%3 for x in a)

def mul(a,b):
    c=[0]*5
    for i in range(3):
        for j in range(3): c[i+j]+=a[i]*b[j]
    for k in (4,3):
        c[k-2]+=c[k]; c[k-3]+=c[k]
    return tuple(v%3 for v in c[:3])


def rank(a,p):
    a=[[x%p for x in row] for row in a]; r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]
        c=pow(a[r][j],-1,p); a[r]=[v*c%p for v in a[r]]
        for i in range(len(a)):
            if i!=r:
                c=a[i][j]; a[i]=[(x-c*y)%p for x,y in zip(a[i],a[r])]
        r+=1
    return r


def paley():
    zero=(0,0,0); squares={mul(x,x) for x in F if x!=zero}
    need(len(squares)==13,'field square count')
    def chi(x): return 0 if x==zero else 1 if x in squares else -1
    K=[[chi(add(y,neg(x))) for y in F] for x in F]
    for i in range(27):
        need(sum(K[i])==0,'regularity')
        for j in range(27):
            need(sum(K[i][k]*K[k][j] for k in range(27))==1-27*(i==j),'Paley square')
    # Coordinates e_i-e_26 on the centered lattice.
    M=[[K[i][j]-K[i][26] for j in range(26)] for i in range(26)]
    need(all((M[i][j]-3*(i==j))%2==0 for i in range(26) for j in range(26)),'integrality')
    C=[[(M[i][j]-3*(i==j))//2 for j in range(26)] for i in range(26)]
    need(rank(K,3)==8 and rank(C,3)==7,'full and centered modular ranks')
    print('PASS: F27 Paley identity; full rank8, centered saturation rank7 modulo3.')


def clock():
    # Two sources merge at a shared port, then terminate. All clocks exact.
    successor={0:2,1:2,2:3}; injection=[2,3,0,0]; W=injection[:]
    for source in (0,1):
        q=source
        while q in successor:
            q=successor[q]; W[q]+=injection[source]
    need(W==[2,3,5,5],'occupation')
    for f in product(range(-2,3),repeat=4):
        norm=sum(w*x*x for w,x in zip(W,f))
        pull=sum(W[q]*f[target]**2 for q,target in successor.items())
        need(norm-pull==2*f[0]**2+3*f[1]**2,'energy identity')
    need(sum(W)==15 and W[3]==sum(injection),'mass identities')
    # Interior injection destroys the tempting separate-inequalities shortcut.
    bad=[1,1]; pushed=[0,1]; injected=[1,1]
    need(all(bad[i]>=injected[i] and bad[i]>=pushed[i] for i in range(2)), 'hostile setup')
    need(not all(bad[i]>=injected[i]+pushed[i] for i in range(2)), 'shortcut was not rejected')
    print('PASS:625 clocked energy identities; mass and interior-injection hostile control.')


if __name__=='__main__':
    paley(); clock()
