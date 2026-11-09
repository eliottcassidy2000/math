#!/usr/bin/env python3
"""Exact stdlib certificates for the fifteen ternary generators.
Usage: python3 verify_orbit.py [coverage_bound=3000] [--mutate].
--mutate corrupts M1 and must fail, also under python3 -O.
Recovered from retained exact session history; no finite-data cofiniteness claim.
"""
from fractions import Fraction as F
from pathlib import Path
from math import gcd
import argparse

def require(condition,description):
    if not condition:raise RuntimeError('Verification failed: '+description)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('coverage_bound',nargs='?',type=int,default=3000)
parser.add_argument('--mutate',action='store_true')
args=parser.parse_args()
require(args.coverage_bound>=1,'coverage_bound must be positive')
def mm(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mv(A,v):return tuple(sum(x*y for x,y in zip(r,v)) for r in A)
def inverse(A):
    R=[[F(x) for x in A[i]]+[F(i==j) for j in range(3)] for i in range(3)]
    for j in range(3):
        p=next((i for i in range(j,3) if R[i][j]),None)
        if p is None:return None
        R[j],R[p]=R[p],R[j]
        z=R[j][j];R[j]=[x/z for x in R[j]]
        for i in range(3):
            if i!=j:
                z=R[i][j];R[i]=[x-z*y for x,y in zip(R[i],R[j])]
    return tuple(tuple(r[3:]) for r in R)
def rowmul(r,A):return tuple(sum(r[i]*A[i][j] for i in range(3)) for j in range(3))
G=[]
for line in Path(__file__).with_name('generators_base3_15.txt').read_text().splitlines():
    if line.startswith('#') or not line.strip():continue
    z=list(map(int,line.split('|')[0].split()))
    G.append(tuple(tuple(z[3*i:3*i+3]) for i in range(3)))
require(len(G)==15,'len(G)==15')
if args.mutate:G[0]=((G[0][0][0]+1,)+G[0][0][1:],)+G[0][1:]
I=((1,0,0),(0,1,0),(0,0,1));P=((1,1,0),(0,0,1),(0,1,1));Q=inverse(P)
v=(3,1,2);seed=(1,0,1)
print('EXACT ALL-ORDER CERTIFICATES (matrix numbering 1..15)')
for i,A in enumerate(G,1):
    C=mm(mm(Q,A),P)
    require(all(x.denominator==1 and x>=0 for r in C for x in r),'integral nonnegative conjugate')
    bound=rowmul((8,-1,0),mm(A,P))
    require(min(bound)>=0,'upper E bound')
    require(min(rowmul((-1,0,3),A))>=0,'h<=3O')
    low=rowmul((-1,6,0),A)
    if i==1:require(low==(-1,-3,6),'low==(-1,-3,6)')
    elif i==5:require(low==(-1,3,3),'low==(-1,3,3)')
    elif i==6:require(low==(-1,6,6),'low==(-1,6,6)')
    else:require(min(low)>=0,'lower E bound')
    require(A[1][0]+A[1][2]>=1,'positive E after first step')
    require(A[0][0]>=1 and A[0][1]>=3,'growth coefficients')
    require(A[0][0]%3==1 and A[0][1]%3==A[0][2]%3==0,'mod3')
    B=inverse(A)
    if B is not None:
        u=mv(B,v)
        if i==1:require(u==tuple(F(x,2) for x in v),'tower inverse')
        else:require(u[0]<=0 or min(u)<0,'ray preimage obstruction')
        cert='inverse v = '+str(tuple(str(x) for x in u))
    else:
        r=(0,1,-1) if i in (2,3,4) else (1,-2,-1)
        require(rowmul(r,A)==(0,0,0) and sum(x*y for x,y in zip(r,v))!=0,'singular separator')
        cert='left null separator '+str(r)
    print(i,cert,'cone matrix',tuple(tuple(int(x) for x in r) for r in C),'8h-E coefficients',bound)
T=G[0];T3=mm(mm(T,T),T);ell=(1,2,1)
require(T3==tuple(tuple(I[i][j]+v[i]*ell[j] for j in range(3)) for i in range(3)),'tower cube')
require(mv(T,v)==tuple(2*x for x in v),'right eigenvector')
require(rowmul(ell,T)==tuple(2*x for x in ell),'left eigenvector')
require(mv(G[3],seed)==(1,8,8),'sharp upper bound')
print('T^3 = I + v ell^T:',T3,'; sharp E/h=8 witness M4(seed)=(1,8,8)')
def cross(u,v):return (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])
def facets(rays):
    out=[]
    for u,vv in zip(rays,rays[1:]+rays[:1]):
        r=cross(u,vv);d=gcd(gcd(abs(r[0]),abs(r[1])),abs(r[2]))
        out.append(tuple(x//d for x in r))
    return out
def dot(u,vv):return sum(x*y for x,y in zip(u,vv))
def word(w):
    s=seed
    for i in w:s=mv(G[i-1],s)
    return s
for title,words in [('full',[(),(4,6),(4,2),(4,)]),('nonempty',[(4,6,1),(4,6,5),(4,6),(4,2),(4,)])]:
    rays=[word(w) for w in words];fs=facets(rays)
    require(all(dot(f,r)>=0 for f in fs for r in rays),'ray containment')
    for j,r in enumerate(rays):
        require(dot(fs[j],r)==dot(fs[j-1],r)==0,'adjacent facets')
        require(all(dot(f,r)>0 for k,f in enumerate(fs) if k not in (j,(j-1)%len(fs))),'strict polygon convexity')
    require(all(dot(f,mv(A,r))>=0 for A in G for r in rays for f in fs),'cone invariance')
    initials=[seed] if title=='full' else [mv(A,seed) for A in G]
    require(all(dot(f,r)>=0 for f in fs for r in initials),'initial containment')
    print('MINIMAL CONE',title,'attained rays',rays,'facets',fs)
    print('Invariance certificates:',len(G)*len(rays)*len(fs),'exact inequalities PASS')
    if title=='nonempty':
        require(min(F(r[1],r[0]) for r in rays)==F(33,193),'sharp lower E')
        require(min(F(mv(A,r)[0],r[0]) for A in G for r in rays)==F(292,193),'sharp growth')
        print('SHARP nonempty E/h >=33/193 and next h / current h >=292/193')
H=args.coverage_bound;seen={seed};queue=[seed]
for s in queue:
    for A in G:
        t=mv(A,s)
        if t[0]>H:continue
        h,e,o=t
        require(h%3==1 and 0<=e<=o<=h+e and e<=8*h and o<=9*h,'finite cone bounds')
        require(e>=1 and h<=6*e and h<=3*o,'finite lower bounds')
        if s!=seed:require(h>=s[0]+3 and 2*h>=3*s[0],'finite growth')
        require(not(h==3*e and o==2*e),'finite ray exclusion')
        if t not in seen:seen.add(t);queue.append(t)
values={s[0] for s in seen};miss=[h for h in range(1,H+1,3) if h not in values]
print('FINITE-EXACT control: bound',H,'states',len(seen),'misses',miss)
print('PASS. Cofiniteness of HYP-9640(A) remains OPEN.')
