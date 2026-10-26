"""Construct unimodular lifts, then independently check determinant and residues."""
import itertools
import json
from math import gcd


def need(ok, label):
    if not ok: raise RuntimeError(label)


def ident(n): return [[int(i==j) for j in range(n)] for i in range(n)]


def det(a):
    a = [list(r) for r in a]
    sign, prev, n = 1,1,len(a)
    for k in range(n-1):
        if a[k][k]==0:
            r=next((r for r in range(k+1,n) if a[r][k]),None)
            if r is None: return 0
            a[k],a[r]=a[r],a[k]
            sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=a[i][j]*pivot-a[i][k]*a[k][j]
                need(value%prev==0,'exact Bareiss division')
                a[i][j]=value//prev
            a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]


def primes(m):
    return [p for p in range(2,m+1) if m%p==0 and all(p%d for d in range(2,p))]


def lift_sl(c,m):
    n=len(c)
    a=[[x%m for x in row] for row in c]
    need(det(a)%m==1,'SL input')
    ops=[]
    def rowadd(i,j,t):
        a[i]=[(x+t*y)%m for x,y in zip(a[i],a[j])]
        ops.append((i,j,t))
    for i in range(n-1):
        ps=primes(m)
        desired={p:{j:0 for j in range(i+1,n)} for p in ps}
        for p in ps:
            if a[i][i]%p==0:
                j=next(j for j in range(i+1,n) if a[j][i]%p)
                desired[p][j]=1
        for j in range(i+1,n):
            t=next(t for t in range(m) if all(t%p==desired[p][j] for p in ps))
            if t: rowadd(i,j,t)
        need(gcd(a[i][i],m)==1,'CRT unit pivot')
        u=pow(a[i][i],-1,m)
        v=pow(u,-1,m)
        # Product w(u)w(-1) = diag(u,u^-1), applied rightmost first.
        for r,s,t in [(i,n-1,-1),(n-1,i,1),(i,n-1,-1),
                      (i,n-1,u),(n-1,i,-v),(i,n-1,u)]:
            rowadd(r,s,t)
        for j in range(n):
            if j!=i and a[j][i]: rowadd(j,i,-a[j][i])
    for j in range(n-1):
        if a[j][-1]: rowadd(j,n-1,-a[j][-1])
    need(a==ident(n),'complete modular elimination')
    q=ident(n)
    for i,j,t in reversed(ops):
        q[i]=[x-t*y for x,y in zip(q[i],q[j])]
    need(det(q)==1,'integer determinant one')
    need(all((q[i][j]-c[i][j])%m==0 for i in range(n) for j in range(n)), 'lift residues')
    return q


def lift(c,m):
    d=det(c)%m
    need(d in {1,m-1},'fixed-rank determinant obstruction')
    sign=1 if d==1 else -1
    a=[list(r) for r in c]
    a[0]=[sign*x for x in a[0]]
    q=lift_sl(a,m)
    q[0]=[sign*x for x in q[0]]
    need(abs(det(q))==1,'GL lift')
    need(all((q[i][j]-c[i][j])%m==0 for i in range(len(c)) for j in range(len(c))), 'GL residues')
    return q


def stabilize(c,m):
    n=len(c)
    d=det(c)%m
    need(gcd(d,m)==1,'invertible root characters')
    a=[list(row)+[0] for row in c]+[[0]*n+[pow(d,-1,m)]]
    return lift_sl(a,m)


def main():
    census=[]
    for m in range(2,13):
        total=admissible=blocked=0
        for entries in itertools.product(range(m),repeat=4):
            c=[list(entries[:2]),list(entries[2:])]
            d=(entries[0]*entries[3]-entries[1]*entries[2])%m
            if gcd(d,m)!=1: continue
            total+=1
            if d in {1,m-1}:
                q=lift(c,m)
                need(q[0][0]*q[1][1]-q[0][1]*q[1][0] in {1,-1},'independent 2x2 determinant')
                admissible+=1
            else:
                stabilize(c,m)
                blocked+=1
        census.append({'m':m,'invertible_2_by_2':total,'same_rank_lifts':admissible,
                       'necessary_one_coordinate_repairs':blocked})
    c=ident(6); c[0][0]=2
    q=stabilize(c,7)
    explicit=[[2,7],[7,25]]
    need(det(explicit)==1 and all((explicit[i][j]-[[2,0],[0,4]][i][j])%7==0 for i in range(2) for j in range(2)),'explicit seventh-root repair')
    # Independent finite root-fibre replay in exponent coordinates over cyclic
    # unit groups, retaining the marked extra-root section v=1.
    fibre_cases=0
    for N in range(1,43):
        for b in range(N):
            source={u for u in range(N) if 7*u%N==b}
            mapped={((2*u)%N,(7*u)%N) for u in source}
            target={(w,z) for w in range(N) for z in range(N)
                    if 7*w%N==2*b%N and z==b and (-7*w+2*z)%N==0}
            need(mapped==target,'marked section image and root equations')
            for w,z in target:
                u=(25*w-7*z)%N
                need(u in source and 2*u%N==w,'explicit inverse recovers root')
            fibre_cases+=1
    return {'matrix_census':census,'six_to_seven_mod7_lift':q,
            'two_coordinate_mod7_lift':explicit,'cyclic_unit_group_orders':[1,42],
            'marked_fibre_cases':fibre_cases,
            'universal_unit_sign_moduli_through100':[m for m in range(2,101)
                if all(a in {1,m-1} for a in range(1,m) if gcd(a,m)==1)]}


if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
