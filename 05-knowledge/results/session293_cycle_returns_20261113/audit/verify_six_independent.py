#!/usr/bin/env python3
"""Independent exact controls; no imports from the primary six-mode lane."""
from fractions import Fraction as Q
from itertools import permutations
from math import lcm

checks = 0
def require(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)

def product(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def unit(n): return [[Q(i==j) for j in range(n)] for i in range(n)]
def mpow(a,n):
    y=unit(len(a))
    for _ in range(n): y=product(y,a)
    return y
def det3(a):
    return sum(a[0][j]*(a[1][(j+1)%3]*a[2][(j+2)%3]-a[1][(j+2)%3]*a[2][(j+1)%3]) for j in range(3))

A=[[Q(0),Q(-3)],[Q(1),Q(3)]]
C=[[Q(0),Q(0),Q(3)],[Q(0),Q(-1),Q(-3)],[Q(1,3),Q(2),Q(3)]]
# Derive the symmetric-square coordinate action directly on basis matrices.
basis=[[[1,0],[0,0]],[[0,1],[1,0]],[[0,0],[0,1]]]
cols=[]
for u in basis:
    out=product(product(A,u),list(map(list,zip(*A))))
    cols.append([out[0][0]/3,out[0][1]/3,out[1][1]/3])
require(list(map(list,zip(*cols)))==C,'derived Sym2 action')
require(mpow(A,6)==[[-27,0],[0,-27]],'projective gain')
require(mpow(C,6)==unit(3),'literal quadratic return')
for k in range(1,6): require(mpow(C,k)!=unit(3),'primitive matrix period')
require(sum(C[i][i] for i in range(3))==2,'characteristic trace')
require(det3(C)==1,'characteristic determinant')
require(sum(mpow(C,2)[i][i] for i in range(3))==0,'characteristic trace square')
# Trace2=0 and trace1=2 give elementary symmetric coefficient2=2,
# hence characteristic polynomial x^3-2x^2+2x-1.
vs=[[Q(3),Q(0),Q(0)]]
for _ in range(5): vs.append([r[0] for r in product(C,[[x] for x in vs[-1]])])
require(vs==[[3,0,0],[0,0,1],[3,-3,3],[9,-6,4],[12,-6,3],[9,-3,1]],'derived six rays')
for v in vs:
    require(v[0]*v[2]==v[1]*v[1] and v[0]>=0 and v[2]>=0,'rank-one PSD')
    require(v[0]+3*v[1]+3*v[2]==3,'positive invariant mass')
require(det3(list(map(list,zip(*vs[:3]))))==9,'span dimension3')
require(product([[1,3,3]],C)==[[1,3,3]],'functional invariance')
require(Q(3)-Q(9,4)>0,'mass form positive definite')
# Rational projective action includes infinity, not division by zero.
def g(x):
    if x is None: return Q(0)
    if x==-3: return None
    return -Q(3)/(x+3)
for x in [None]+[Q(a,b) for a in range(-8,9) for b in range(1,8)]:
    orbit=[x]
    for _ in range(6): orbit.append(g(orbit[-1]))
    require(orbit[6]==x and len(set(orbit[:6]))==6,'rational exact period6')
counts=[]
for d in range(1,7):
    exact=primitive=0
    for perm in permutations(range(d)):
        remaining=set(range(d)); sizes=[]
        while remaining:
            v=next(iter(remaining)); size=0
            while v in remaining:
                remaining.remove(v); size+=1; v=perm[v]
            sizes.append(size)
        exact+=lcm(*sizes)==6
        primitive+=any(s%6==0 for s in sizes)
    require(primitive==(120 if d==6 else 0),'primitive6 permutation lower bound')
    counts.append([d,exact,primitive])
require(counts[4][1]>0,'five states can have abstract order6')
print({'checks':checks,'permutation_counts_n_order6_primitive6':counts,'status':'PASS'})

# General-q cyclotomic image/kernel, independently generated polynomial data.
def trim(a):
    while len(a)>1 and a[-1]==0: a.pop()
    return a

def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return trim(out)

def divrem(a,b):
    a=trim(list(a)); b=trim(list(b)); quotient=[0]*max(1,len(a)-len(b)+1)
    require(b[-1]==1,'monic divisor')
    while len(a)>=len(b) and any(a):
        k=len(a)-len(b); c=a[-1]; quotient[k]=c
        for j,y in enumerate(b): a[j+k]-=c*y
        trim(a)
    return trim(quotient),trim(a)

phis={1:[-1,1]}
for q in range(2,31):
    poly=[-1]+[0]*(q-1)+[1]
    for d in range(1,q):
        if q%d==0:
            poly,rem=divrem(poly,phis[d]); require(rem==[0],'cyclotomic recursive division')
    phis[q]=poly
    fdeg=len(poly)-1; p1=sum(poly)
    require(p1>0,'Phiq1 positive q>=2')
    # Six honest target residues with the correct integral congruence.
    for k in range(-3,4):
        r=[(i+2*k)%5-2 for i in range(fdeg)]
        a=sum(r)+k*p1
        f=r+[0]
        for i,c in enumerate(poly): f[i]+=k*c
        require(len(trim(f))<=q,'constructive target degree fits q states')
        require(sum(f)==a and divrem(f,poly)[1]==trim(r[:]),'constructive signed integral image')
    kernel=mul([-1,1],poly)
    for shift in range(q-fdeg-1):
        f=[0]*shift+kernel
        require(sum(f)==0 and divrem(f,poly)[1]==[0],'kernel generator exact')
    # Source columns obey the lattice congruence, with evaluation on canonical remainder.
    for j in range(q):
        rem=divrem([0]*j+[1],poly)[1]
        require((1-sum(rem))%p1==0,'image congruence for each column')
require({q:len(phis[q]) for q in (6,12,24,30)}=={6:3,12:5,24:9,30:9},'q-to-observed dimension')
print({'total_checks_including_cyclotomic':checks,'q_range':[2,30],'status':'PASS'})
