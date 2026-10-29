#!/usr/bin/env python3
"""Independent exact audit. Stdlib only; no producer modules are imported."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd, isqrt
import json
import sys

checks = []
def require(condition, label):
    if not condition:
        raise RuntimeError('AUDIT FAILURE: ' + label)
    checks.append(label)

def transpose(a):
    return list(map(list, zip(*a)))
def multiply(a, b):
    return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def rank(a, modulus=None):
    a = [[F(x) if modulus is None else int(x)%modulus for x in r] for r in a]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r,len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r],a[pivot] = a[pivot],a[r]
        inv = 1/a[r][c] if modulus is None else pow(a[r][c],-1,modulus)
        a[r] = [x*inv if modulus is None else x*inv%modulus for x in a[r]]
        for i in range(len(a)):
            if i != r:
                q=a[i][c]
                a[i]=[x-q*y if modulus is None else (x-q*y)%modulus for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a): break
    return r

def wedge(a):
    pairs=list(combinations(range(len(a)),2))
    return [[a[i][k]*a[j][l]-a[i][l]*a[j][k] for k,l in pairs] for i,j in pairs]
def ident(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]
def det(a):
    n=len(a)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))*prod(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))
def prod(xs):
    out=1
    for x in xs: out*=x
    return out
# Coefficient lists in ascending order; these prove polynomial identities exactly.
def add(a,b):
    return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]
def scale(a,s): return [x*s for x in a]
def pmul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def charpoly(a):
    n=len(a); out=[0]*(n+1)
    for p in permutations(range(n)):
        term=[(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))]
        for i in range(n): term=pmul(term,[-a[i][p[i]],int(i==p[i])])
        out=add(out,term)
    return out

# Independently encode the printed machine, not the producer implementation.
table=[[0,2,5,5,5],[5,5,0,1,2],[5,5,0,1,4],[3,4,5,5,5],[5,5,3,4,4],[5,5,5,5,5]]
letters=[[[int(table[j][l]==i) for j in range(6)] for i in range(6)] for l in range(5)]
reached={1}
while True:
    new=reached|{table[s][l] for s in reached for l in range(5)}
    if new==reached: break
    reached=new
require(len(reached)==6,'six states reachable')
def minimal_blocks(outputs):
    labels=outputs
    while True:
        keys=[(labels[s],tuple(labels[t] for t in table[s])) for s in range(6)]
        ids={k:i for i,k in enumerate(dict.fromkeys(keys))}
        new=[ids[k] for k in keys]
        if all((labels[i]==labels[j])==(new[i]==new[j]) for i in range(6) for j in range(6)):
            return new
        labels=new
require(len(set(minimal_blocks([0,0,0,1,1,2])))==6,'six-state output minimality')
accept=minimal_blocks([0,0,0,1,1,1])
require(len(set(accept))==4,'four-state acceptance minimality')
require([rank(m) for m in letters]==[3]*5,'six-state letters rank three')
require(rank(multiply(letters[4],letters[4]))==2,'ee rank two')
require({r*(r-1)//2 for r in range(5)}=={0,1,3,6},'dimension-four exterior rank spectrum excludes two')
q=[0,1,2,3,3,3]
for l in range(5):
    targets=[q[table[s][l]] for s in [0,1,2,5]]
    qm=[[int(targets[j]==i) for j in range(4)] for i in range(4)]
    require(rank(qm)==2 and rank(wedge(qm))==1,'four-state quotient ranks letter '+str(l))

# Standard oriented incidence and rational/F2 distinctions.
edges=list(combinations(range(4),2))
b=[[(int(v==j)-int(v==i)) for i,j in edges] for v in range(4)]
require(rank(b)==3 and rank(b,2)==3,'incidence rank three over Q and F2')
lap=multiply(b,transpose(b))
require(lap==[[4*int(i==j)-1 for j in range(4)] for i in range(4)],'K4 Laplacian 4I-J')
h=[[int(i==j+1)-int(i==0) for j in range(3)] for i in range(4)]
cycle=[[h[i][k]*h[j][l]-h[i][l]*h[j][k] for k,l in combinations(range(3),2)] for i,j in edges]
cut=transpose(b)
require(rank(cycle)==3 and multiply(b,cycle)==[[0]*3 for _ in range(4)],'Lambda2H spans cycle kernel')
require(rank([cut[i]+cycle[i] for i in range(6)])==6,'rational cut-cycle direct sum')
require(rank(lap,2)==1 and rank(b,2)-rank(lap,2)==2,'F2 bicycle dimension two')
# *e01=e23, *e02=-e13, *e03=e12; involution.
star=[[0]*6 for _ in range(6)]
for i,j,s in [(0,5,1),(1,4,-1),(2,3,1)]: star[i][j]=star[j][i]=s
require(multiply(star,star)==ident(6),'Hodge involution')
for sign in [1,-1]:
    require(rank([[star[i][j]-sign*int(i==j) for j in range(6)] for i in range(6)])==3,'Hodge eigenspace dimension '+str(sign))
require(multiply(b,multiply(star,cut))==[[0]*4 for _ in range(4)],'star sends cuts to cycles')
require(rank([cut[i]+multiply(star,cycle)[i] for i in range(6)])==3,'star sends cycles to cuts')

# Symbolic Pluecker relation, four indeterminates, full coefficient cancellation.
def dif(i,j):
    ei=tuple(int(t==i) for t in range(4)); ej=tuple(int(t==j) for t in range(4))
    return {ei:-1,ej:1}
p={}
for i,j,k,l,s in [(0,1,2,3,1),(0,2,1,3,-1),(0,3,1,2,1)]:
    for u,a in dif(i,j).items():
        for v,c in dif(k,l).items():
            e=tuple(x+y for x,y in zip(u,v)); p[e]=p.get(e,0)+s*a*c
require(all(v==0 for v in p.values()),'symbolic Pluecker cancellation')
require(1*1-0*0+0*0!=0,'e01+e23 is nondecomposable control')
for x in combinations(range(-4,5),4):
    a=(x[1]-x[0])*(x[3]-x[2]); bb=(x[3]-x[0])*(x[2]-x[1]); c=(x[2]-x[0])*(x[3]-x[1])
    if not (a>0 and bb>0 and c>0 and a+bb==c): raise RuntimeError('ordered positivity')
require(True,'126 ordered potential positivity controls')
# C plus exterior square and the volume-coordinate duality J*Lambda2C=det(C)*C^-T*J.
c=[[1,0,1],[1,1,1],[1,1,0]]; w=wedge(c)
require(det(c)==-1,'determinant C minus one')
require(charpoly(c)==[1,-1,-2,1],'characteristic polynomial C')
require(charpoly(w)==[-1,-2,1,1],'characteristic polynomial exterior square C')
j=[[0,0,1],[0,-1,0],[1,0,0]]
require(multiply(transpose(c),multiply(j,w))==[[-x for x in row] for row in j], 'contragredient volume identity')
l=[[1,0,0,0]]+[[0]+row for row in c]
wl=wedge(l)
require(wl==[[c[i][j] if i<3 and j<3 else w[i-3][j-3] if i>=3 and j>=3 else 0 for j in range(6)] for i in range(6)],'synthetic exterior block identity')

# Integral shifted-dual and cyclotomic-seven bridge (exact, not numeric).
ss=[[0,-1,2],[-1,-2,4],[-1,-3,5]]
d=[[c[i][j]-int(i==j) for j in range(3)] for i in range(3)]
require(det(ss)==1 and multiply(w,ss)==multiply(ss,d),'unimodular shifted exterior dual intertwiner')
require(charpoly(d)==charpoly(w),'shifted counter characteristic polynomial')
def msub(a,b): return [[x-y for x,y in zip(ar,br)] for ar,br in zip(a,b)]
def mscale(a,t): return [[t*x for x in row] for row in a]
i3=ident(3)
f=msub(multiply(d,d),mscale(i3,2))
require(multiply(f,c)==mscale(i3,-1),'root map matrix identity (D^2-2I)C=-I')
f2=msub(multiply(f,f),mscale(i3,2))
f3=msub(multiply(f2,f2),mscale(i3,2))
require(f3==d and f!=d and f2!=d and f2!=f,'root map exact order three')
# z^3*q(z+z^-1) = (z^2+1)^3+z(z^2+1)^2-2z^2(z^2+1)-z^3.
u=[1,0,1]
cyclo=add(add(pmul(pmul(u,u),u),pmul([0,1],pmul(u,u))),add(scale(pmul([0,0,1],u),-2),[0,0,0,-1]))
require(cyclo==[1]*7,'Laurent substitution yields seventh cyclotomic polynomial')
# q(x)=0 implies x(x+1) nonzero; x^2-2 differs from x since gcd(q,x^2-x-2)=1:
# fixed points of x^2-2=x are 2 and -1, neither a root (q(2)=7,q(-1)=1).
require(sum(t*2**i for i,t in enumerate(charpoly(d)))==7 and sum(t*(-1)**i for i,t in enumerate(charpoly(d)))==1 and charpoly(d)[0]==-1,'root map has no poles or fixed roots')
v=[[-1],[-1],[-1]]
bd=[v[i]+multiply(d,v)[i]+multiply(multiply(d,d),v)[i] for i in range(3)]
bw=[v[i]+multiply(w,v)[i]+multiply(multiply(w,w),v)[i] for i in range(3)]
require(det(bd)==-1 and det(bw)==-1 and multiply(ss,bd)==bw,'claimed integral cyclic bases and change of basis')

# Exact all-k coefficient proofs for D1 and D2 from polynomial sums.
k=[0,1]; kp1=[1,1]; twokp1=[1,2]; fourkp1=[1,4]
s1=scale(pmul(k,kp1),F(1,2)); s1_2=pmul(k,twokp1)
s2=scale(pmul(pmul(k,kp1),twokp1),F(1,6)); s2_2=scale(pmul(pmul(k,twokp1),fourkp1),F(1,3))
require(add(scale(s1,2),scale(s1_2,-1))==[0,0,-1],'all-k D1 constant coefficient -k^2')
require(add(scale(s2,2),scale(s2_2,-1))==[0,0,-1,-2],'all-k D2 constant coefficient -k^2(2k+1)')
require(add(k,scale(pmul(k,twokp1),-1))==[0,0,-2],'all-k D2 factor linear coefficient')
require(scale(pmul(pmul(k,k),twokp1),-1)==[0,0,-1,-2],'all-k D2 factor constant coefficient')
def defect(a,k,p,d=1):
    if not isinstance(k,int) or k<1 or d<=0 or p not in [0,1,2]: raise ValueError('domain')
    return sum((a+d*i)**p for i in range(k+1))-sum((a+d*i)**p for i in range(k+1,2*k+1))
for kk,aa,dd in product(range(1,13),range(-30,31),[1,2,5]):
    if defect(aa,kk,1,dd)!=aa-dd*kk**2 or defect(aa,kk,2,dd)!=(aa+dd*kk)*(aa-dd*kk*(2*kk+1)): raise RuntimeError('defect check')
require(True,'2196 direct progression-defect controls')
for kk in range(1,51):
    require(defect(kk*kk,kk,1)==0 and defect(kk*(2*kk+1),kk,2)==0 and defect(-kk,kk,2)==0,'roots k='+str(kk))
    for aa in [F(-5,2),F(0),F(7,3),F(100)]:
        if defect(aa+1,kk,2)-defect(aa,kk,2)!=2*defect(aa,kk,1)+1 or defect(aa+F(1,2),kk,2)-defect(aa-F(1,2),kk,2)!=2*defect(aa,kk,1): raise RuntimeError('finite differences')
require(True,'exact rational forward and centered finite differences')
for args in [(1,0,1),(1,-1,1),(1,1,3),(1,1,2,0)]:
    try: defect(*args)
    except ValueError: pass
    else: raise RuntimeError('illegal argument accepted')
require(True,'invalid-input controls')
# New strengthening: every primitive positive additive triple has this encoding.
for a in range(1,51):
    for bb in range(a,51):
        if gcd(a,bb)!=1: continue
        if a==bb:
            x=(0,2,3,6); factor=6
        else:
            d=bb-a; x=(0,d,2*d,2*bb); factor=2*d
        aa=(x[1]-x[0])*(x[3]-x[2])
        bbb=(x[3]-x[0])*(x[2]-x[1])
        cc=(x[2]-x[0])*(x[3]-x[1])
        if not (list(x)==sorted(set(x)) and (aa,bbb,cc)==(factor*a,factor*bb,factor*(a+bb)) and gcd(gcd(aa,bbb),cc)==factor):
            raise RuntimeError('primitive additive triple lift')
require(True,'primitive additive triple Pluecker realization controls')

def square_graph(n):
    return {v:{u for u in range(1,n+1) if u!=v and isqrt(u+v)**2==u+v} for v in range(1,n+1)}
def components(g):
    unseen=set(g); out=[]
    while unseen:
        seed=min(unseen); found={seed}; pending=[seed]
        while pending:
            v=pending.pop()
            for u in g[v]-found: found.add(u); pending.append(u)
        out.append(found); unseen-=found
    return out
expected=[{1,3,6,8,10},{2,7,9},{4,5,11,12}]
require(components(square_graph(12))==expected,'Q12 three component vertex sets')
require([len(components(square_graph(n))) for n in [12,13,14]]==[3,2,1],'Q12 Q13 Q14 component counts')
require(square_graph(13)[13]=={3,12} and square_graph(14)[14]=={2,11},'merger neighbors at 13 and 14')
require([n for n in range(1,201) if not square_graph(n)[n]]==[1,2,4],'square-sum founders through 200')
for n in range(4,13):
    g=square_graph(n)
    if len(components(g))!=3 or any(len(v)>2 for v in g.values()) or sum(map(len,g.values()))//2!=n-3: raise RuntimeError('three path phase')
require(True,'three disjoint path components for 4 through 12')
require(square_graph(32)[32]=={4,17},'32 is born with two neighbors, not a leaf')
require([m for m in range(1,501) if isqrt(2*m)**2==2*m]==[2*k*k for k in range(1,16)],'excluded square diagonals equal even-balance starting terms')
print(json.dumps({'status':'PASS','optimized':not __debug__,'check_count':len(checks),'checks':checks,'scope':'Independent finite controls plus exact polynomial/matrix identities; source-recovery graph computations verified; provenance is documented by source-recovery workstream; no abc consequence.'},indent=2))
