#!/usr/bin/env python3
"""Exact self-contained controls for the integral quadratic guard lattice."""
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, prod
import hashlib
import json
from pathlib import Path

CHECKS=0

def require(condition,message):
    global CHECKS
    CHECKS+=1
    if not condition:
        raise RuntimeError(message)


def step(v,a):
    s,u,d=v
    return (u,s+u+d,s) if a=='C' else (s+d,2*u,0)


def state(w):
    v=(1,2,0)
    for a in w:
        v=step(v,a)
    return v


def phi(v):
    s,u,d=v
    return [s*s,u*u,d*d,s*u,s*d,u*d]


def mm(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def determinant(A):
    a=[list(row) for row in A]
    n=len(a)
    if not n:
        return 1
    sign=1;old=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if a[i][k]),None)
        if p is None:
            return 0
        if p!=k:
            a[k],a[p]=a[p],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=pivot*a[i][j]-a[i][k]*a[k][j]
                require(z%old==0,'Bareiss exact division')
                a[i][j]=z//old
            a[i][k]=0
        old=pivot
    return sign*a[-1][-1]


def rank(A,p=None):
    a=[[x%p for x in row] for row in A] if p else [[Fraction(x) for x in row] for row in A]
    if not a:
        return 0
    r=0
    for j in range(len(a[0])):
        k=next((k for k in range(r,len(a)) if a[k][j]),None)
        if k is None:
            continue
        a[r],a[k]=a[k],a[r]
        scale=pow(a[r][j],-1,p) if p else 1/a[r][j]
        a[r]=[(x*scale)%p if p else x*scale for x in a[r]]
        for k in range(r+1,len(a)):
            scale=a[k][j]
            a[k]=[(x-scale*y)%p if p else x-scale*y for x,y in zip(a[k],a[r])]
        r+=1
        if r==len(a):
            break
    return r


W=['','C','D','CC','CCC','DCC']
V=[state(w) for w in W]
E=[phi(v) for v in V]
T_C=[[0,1,1,0,0,1],[1,1,0,1,0,0],[0,1,0,0,0,0],
     [0,2,0,1,1,1],[0,2,0,0,0,1],[0,2,0,1,0,0]]
T_D=[[1,0,0,0,0,0],[0,4,0,0,0,0],[1,0,0,0,0,0],
     [0,0,0,2,0,0],[2,0,0,0,0,0],[0,0,0,2,0,0]]
B_C=[[0,1,0,0,0,0],[0,0,0,1,0,0],[29,30,-3,41,-9,-7],
     [0,0,0,0,1,0],[-12,-6,1,-15,6,3],[3,-1,0,6,2,-1]]
B_D=[[0,0,1,0,0,0],[9,0,0,0,0,0],[25,36,3,54,-12,-9],
     [2,-24,5,-36,8,6],[0,-84,18,-126,28,21],[36,-48,9,-72,16,12]]
B={'C':B_C,'D':B_D}


def literal(w):
    # 0 -> 1 -> 2 -> 0, marked q=0.
    adj=[1<<1,1<<2,1<<0];q=0
    for a in w:
        n=len(adj)
        adj[q]|=1<<n
        adj.append(((1<<n)-1)^(1<<q))
        if a=='C':
            q=n
    return adj,q


def hp(adj,vertices,start=None):
    """Independent endpoint subset DP on a specified induced subtournament."""
    vs=list(vertices);n=len(vs)
    if not n:
        return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i,v in enumerate(vs):
        if start is None or v==start:
            dp[1<<i][i]=1
    for mask in range(1<<n):
        for i,v in enumerate(vs):
            c=dp[mask][i]
            if not c:
                continue
            for j,z in enumerate(vs):
                if not (mask>>j)&1 and (adj[v]>>z)&1:
                    dp[mask|(1<<j)][j]+=c
    return sum(dp[-1])


def strongly_connected(adj):
    n=len(adj)
    for s in range(n):
        seen=1<<s
        while True:
            new=seen
            for i in range(n):
                if (seen>>i)&1:
                    new|=adj[i]
            if new==seen:
                break
            seen=new
        if seen!=(1<<n)-1:
            return False
    return True


def g(n,t):
    vals=[0,1,t]
    while len(vals)<=n:
        vals.append(sum(vals[-3:]))
    return vals[n]


def delta_product(d):
    return prod((2**(j+1)-2**(i+1))**(d-j+1)
                for j in range(1,d+1) for i in range(j))*prod(
                g(n,2**(j+1))**(d-n-j+1)
                for n in range(1,d+1) for j in range(d-n+1))


def factor(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:
            out[str(p)]=out.get(str(p),0)+1;n//=p
        p+=1
    if n>1:
        out[str(n)]=out.get(str(n),0)+1
    return out


def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]


def main():
    require(V==[(1,2,0),(2,3,1),(1,4,0),(3,6,2),(6,11,3),(5,10,4)],'bank states')
    require(determinant(E)==64,'integral index')
    for A,T in [('C',T_C),('D',T_D)]:
        require(mm(E,T)==mm(B[A],E),'exact invariant row lattice '+A)
        for v in product(range(-2,3),repeat=3):
            require(mm([phi(v)],T)[0]==phi(step(v,A)),'quadratic substitution '+A)
    physical=[]
    for w,v in zip(W,V):
        adj,q=literal(w);n=len(adj)
        ins=[r for r in range(n) if (adj[r]>>q)&1]
        require(len(ins)==1,'unique in-neighbor')
        r=ins[0];s=hp(adj,range(n),q);total=hp(adj,range(n))
        d=hp(adj,[i for i in range(n) if i!=q],r)
        require((s,total-s,d)==v,'independent literal Hamilton counts '+w)
        require(strongly_connected(adj),'literal strongness '+w)
        physical.append({'word':w,'vertices':n,'guard':q,'state':list(v),'Hamilton_paths':total})
    divisors=[]
    for k in range(1,7):
        z=0
        for rows in combinations(range(6),k):
            for cols in combinations(range(6),k):
                z=gcd(z,abs(determinant([[E[i][j] for j in cols] for i in rows])))
        divisors.append(z)
    require(divisors==[1,1,2,4,8,64],'Smith determinantal divisors')
    # Independent recurrence propagation controls, beyond the closure proof.
    frontier=[((1,2,0),[1,0,0,0,0,0])];word_count=0
    for depth in range(13):
        nxt=[]
        for v,c in frontier:
            require(mm([c],E)[0]==phi(v),'all-word row propagation')
            word_count+=1
            if depth<12:
                for a in 'CD':
                    nxt.append((step(v,a),mm([c],B[a])[0]))
        frontier=nxt
    # Radius-two obstruction mod 3, and an exact rank check.
    small=[state(''.join(w)) for n in range(3) for w in product('CD',repeat=n)]
    require(rank([phi(v) for v in small],3)==5,'entire radius-two quadratic rank mod 3')
    require(all((s*s-u*u+s*d+u*d)%3==0 for s,u,d in small),'explicit radius-two hostile')
    # All characteristic-two reachable residues, checked as a closed orbit.
    orbit={(1,0,0)};todo=[(1,0,0)]
    while todo:
        v=todo.pop()
        for a in 'CD':
            u=tuple(x%2 for x in step(v,a))
            if u not in orbit:
                orbit.add(u);todo.append(u)
    require(orbit=={(1,0,0),(0,1,1)},'complete characteristic-two orbit')
    require(rank([phi(v) for v in orbit],2)==2,'intrinsic quadratic rank two')
    require(all((s**3*u-s*u**3)%3==0 for s,u,d in product(range(3),repeat=3)),
            'degree-four Frobenius hostile on all prime-field points')
    # Full tangent-relation and wedge-period matrices.
    relations=[[2,-2,-1,1,0,0],[1,-3,-1,0,1,0],[5,-4,-2,0,0,1]]
    require(mm(relations,V)==[[0,0,0]]*3,'exact tangent relations')
    require(determinant(V[:3])==-2,'odd-characteristic spanning')
    period=[]
    basis=[[1,0,0],[0,1,0],[0,0,1]]
    for rel in relations:
        for a in range(3):
            period.append([rel[i]*cross(V[i],basis[j])[a] for i in range(6) for j in range(3)])
    rank_report={}
    for p in [None,3,5,7,11,13,17,19,23,31,101]:
        require(rank(E,p)==6,'quadratic bank field rank')
        require(rank(period,p)==9,'complete wedge-period rank')
        rank_report['Q' if p is None else str(p)]={'evaluation':6,'period':9}
    # Derive common-motion/tangent inclusion directly in the 18-coordinate kernel.
    K=[]
    for a in range(3):
        K.append([int(j==a) for i in range(6) for j in range(3)])
    for i in range(6):
        K.append([V[i][j] if h==i else 0 for h in range(6) for j in range(3)])
    require(mm(period,list(map(list,zip(*K))))==[[0]*9 for _ in range(9)],'common and tangential kernel inclusion')
    for p in [None,3,5,7,11]:
        require(rank(K,p)==9,'kernel dimension by direct basis')
    # All-degree determinant formula controls.
    dets=[]
    for d in range(9):
        mons=[(i,j,d-i-j) for i in range(d+1) for j in range(d+1-i)]
        nodes=[state('D'*j+'C'*k) for k in range(d+1) for j in range(d-k+1)]
        M=[[s**i*u**j*z**k for i,j,k in mons] for s,u,z in nodes]
        exact=abs(determinant(M));predicted=delta_product(d)
        require(len(nodes)==comb(d+2,2),'bank dimension')
        require(exact==predicted,'all-degree determinant product')
        dets.append({'degree':d,'dimension':len(nodes),'determinant':str(exact),'factorization':factor(exact)})
        for p in [3,5,7,11,13,17,19,23,31]:
            criterion=(all((2**(j+1)-2**(i+1))%p for j in range(1,d+1) for i in range(j))
                       and all(g(n,2**(j+1))%p for n in range(1,d+1) for j in range(d-n+1)))
            require(criterion==(exact%p!=0),'exact modular determinant criterion')
    out={'status':'PASS','checks':CHECKS,'bank':physical,'determinant':64,
         'determinantal_divisors':divisors,'smith_factors':[1,1,2,2,2,8],
         'word_controls':word_count,'maximum_control_word_length':12,
         'radius_two_rank_mod_3':5,'orbit_mod_2':[list(v) for v in sorted(orbit)],
         'field_rank_controls':rank_report,'homogeneous_bank_determinants':dets,
         'matrix_certificate_sha256':hashlib.sha256(json.dumps([E,T_C,T_D,B_C,B_D],separators=(',',':')).encode()).hexdigest()}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
