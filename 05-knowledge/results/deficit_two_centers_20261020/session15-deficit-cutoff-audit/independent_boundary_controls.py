#!/usr/bin/env python3
"""Small exact hostile boundary checks; not a replacement for symbolic proof."""
from itertools import combinations, permutations
from pathlib import Path
import json, math

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def comp(p,q):return tuple(p[q[i]] for i in range(len(p)))
def inv(p):return tuple(p.index(i) for i in range(len(p)))
def order(p):
    z=tuple(range(len(p)));r=z;k=0
    while True:
        k+=1;r=comp(p,r)
        if r==z:return k

def strong(a,V):
    S=sum(1<<i for i in V)
    if not V:return False
    for rev in [False,True]:
        seen=1<<V[0]
        while True:
            new=seen
            for i in V:
                if seen>>i&1:
                    for j in V:
                        if (a[j][i] if rev else a[i][j]):new|=1<<j
            if new==seen:break
            seen=new
        if seen!=S:return False
    return True

def half_support():
    r=tuple(reversed(range(8)));done=0
    for z in permutations(range(4)):
        p=tuple(z)+tuple(range(4,8))
        e=order(p)
        if e%2:continue
        # r_Q*r_P and (q*r_P)^2, respectively.
        a=comp(p,comp(r,comp(inv(p),r)))
        b=comp(comp(p,r),comp(p,r))
        need(order(a)==e and order(b)==e,'bad disjoint support order identity')
        done+=1
    return dict(even_order_permutations_checked=done)

def sharp_connectivity():
    out=[]
    for m in [3,5,7,9]:
        n=2*m;d=(m-1)//2
        a=[[False]*n for _ in range(n)]
        # Labels 0..m-1 for O0; m..2m-1 for O1.
        for i in range(n):
            for j in range(n):
                b,t=divmod(i,m);c,u=divmod(j,m)
                if b==c:a[i][j]=0<(u-t)%m<=d
                elif b==0:a[i][j]=t==u
                else:a[i][j]=t!=u
        need(all(a[i][j]!=a[j][i] for i in range(n) for j in range(i+1,n)),'not tournament')
        R=tuple(m+i if i<m else (i-m+1)%m for i in range(n))
        R2=comp(R,R)
        sigma=tuple((1-(i//m))*m+(-(i%m))%m for i in range(n))
        need(order(R)==n,'not full cycle')
        need(all(a[i][R[i]] for i in range(n)),'not directed Hamilton cycle')
        need(all(a[R2[i]][R2[j]]==a[i][j] for i in range(n) for j in range(n)),'R2 not auto')
        need(all(a[sigma[i]][sigma[j]]==a[j][i] for i in range(n) for j in range(n)),'sigma not anti')
        good=0
        for s in range((m+1)//2):
            for deleted in combinations(range(n),s):
                V=[i for i in range(n) if i not in deleted]
                need(strong(a,V),f'connectivity below bound m={m}, deleted={deleted}')
                good+=1
        D=[i for i in range(n) if a[0][i]]
        need(len(D)==(m+1)//2,'wrong witness size')
        need(not strong(a,[i for i in range(n) if i not in D]),'not a separator')
        out.append(dict(m=m,n=n,connectivity=(m+1)//2,safe_deletions_checked=good,separator=D))
    return out

def c5():
    n=5;a=[[0<(j-i)%5<=2 for j in range(n)] for i in range(n)]
    pp=[p for p in permutations(range(n)) if all(a[p[i]][p[i+1]] for i in range(n-1))]
    def rev(p):
        r=[None]*n
        for i in range(n):r[p[i]]=p[n-1-i]
        return tuple(r)
    need(len(pp)==15,'C5 H changed')
    for p in pp:
        r=rev(p)
        need(all(a[r[i]][r[j]]==a[j][i] for i in range(n) for j in range(n)),'not pure blue')
    for p,q in combinations(pp,2):
        z=[None]*n
        for i in range(n):z[p[i]]=q[i]
        e=order(z)
        if e%2==0:
            S={i for i in range(n) if z[i]!=i};r=rev(p)
            need(bool(S & {r[i] for i in S}),'counterexample to half support lemma')
            return dict(H=15,P=p,Q=q,relative_order=e,support=sorted(S),reflected_support=sorted(r[i] for i in S))
    raise RuntimeError('expected C5 control missing')

out=dict(half_support=half_support(),sharp_q2_connectivity=sharp_connectivity(),overlap_c5=c5())
Path(__file__).with_name('independent_boundary_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
