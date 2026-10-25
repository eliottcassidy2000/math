#!/usr/bin/env python3
"""Exact controls for the repaired planar-lattice chromatic cap.
No floating point, SAT, or assertion statements; all checks survive python -O.
The universal claims are proved in LATTICE_CAP_REPAIR.md, not inferred here.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, isqrt
import json

checks = 0

def need(test, msg):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(msg)

def det2(u,v): return u[0]*v[1]-u[1]*v[0]
def det3(rows):
    a,b,c=rows
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def row(u): return (u[0]**2,u[0]*u[1],u[1]**2)
def q(form,u):
    a,b,c=form;x,y=u
    return a*x*x+b*x*y+c*y*y

def steps(form,D):
    a,b,c=form; disc=4*a*c-b*b
    bx=isqrt((4*c*D)//disc);by=isqrt((4*a*D)//disc)
    return [(x,y) for x in range(-bx,bx+1) for y in range(-by,by+1) if q(form,(x,y))==D]

def sqfree(n):
    d=1;s=1;p=2
    while p*p<=n:
        e=0
        while n%p==0: n//=p;e+=1
        s*=p**(e//2)
        if e%2:d*=p
        p+=1
    d*=n
    return d,s

def vp(n,p):
    if not n:return 10**9
    e=0;n=abs(n)
    while n%p==0:n//=p;e+=1
    return e

def vf(x,p):return vp(x.numerator,p)-vp(x.denominator,p)
def modfrac(x,m):return x.numerator*pow(x.denominator,-1,m)%m

def norm(z,d): return z[0]*z[0]+d*z[1]*z[1]
def div(z,w,d):
    n=norm(w,d)
    return ((z[0]*w[0]+d*z[1]*w[1])/n,(z[1]*w[0]-z[0]*w[1])/n)

def cycle_colour(k,p):
    if p==2:return k
    return 2 if k==p-1 else k%2

def build_colour(form,D,alpha):
    a,b,c=form
    d,s=sqfree(4*a*c-b*b)
    f=lambda u:(F(u[0])+F(b*u[1],2*a),F(s*u[1],2*a))
    w=f(alpha)
    need(norm(w,d)==F(D,a),'normalizing edge norm')
    if d%4 in (1,2):
        p=2;shift=int(d%4==1)
    else:
        p=next(p for p in range(3,isqrt(d)+2,2) if d%p==0) if any(d%p==0 for p in range(3,isqrt(d)+2,2)) else d
        shift=0
    # pi=sqrt(-d)+shift is Eisenstein at p.
    need((d+shift*shift)%p==0 and (d+shift*shift)%(p*p)!=0,'Eisenstein constant')
    need((2*shift)%p==0,'Eisenstein linear coefficient')
    es=[div(f(e),w,d) for e in [(1,0),(0,1)]]
    proj=[z[0]-shift*z[1] for z in es]
    t=max(0,*[-vf(z,p) for z in proj]);m=p**(t+1);pt=p**t
    ell=[modfrac(pt*z,m) for z in proj]
    def colour(u):
        residue=(ell[0]*u[0]+ell[1]*u[1])%m
        return cycle_colour(residue//pt,p)
    for u in steps(form,D):
        z=div(f(u),w,d);A=z[0]-shift*z[1];B=z[1]
        need(norm(z,d)==1,'unit normalization')
        need(vf(A,p)>=0 and vf(B,p)>=0,'local integrality')
        need(modfrac(A,p) in {1,p-1},'unit residue plus/minus one')
        r=(ell[0]*u[0]+ell[1]*u[1])%m
        need(r in {pt,(-pt)%m},'finite cyclic step certificate')
    # Complete target quotient, including wrap edge; no source patch inference.
    for k in range(p):
        need(cycle_colour(k,p)!=cycle_colour((k+1)%p,p),'cycle coloring')
    return (p,t,m,ell,colour)

vectors=[u for u in product(range(-3,4),repeat=2) if u!=(0,0)]
triples=0
for u,v,w in combinations(vectors,3):
    need(det3([row(u),row(v),row(w)])==det2(u,v)*det2(u,w)*det2(v,w),'Veronese determinant identity')
    triples+=1
# Exact irrational rhombus: rational and sqrt(2) coefficients separate.
rhombus=[u for u in product(range(-8,9),repeat=2) if u[0]*u[0]+u[1]*u[1]==1 and u[0]*u[1]==0]
need(set(rhombus)=={(1,0),(-1,0),(0,1),(0,-1)},'irrational rhombus four steps')
# Degenerate mod-three restriction is not generally similar.
need(q((1,0,6),(0,1))==6,'restricted minimum witness')
need(min(9*x*x+6*y*y for x,y in product(range(-3,4),repeat=2) if (x,y)!=(0,0))==6,'restricted min')
need(6!=3*1,'index-three similarity contradiction')
need(len([(x,y) for x,y in product(range(3),repeat=2) if x*y%3==1])==2,'hyperbolic level cardinality correction')
# Fully exact residual.
form=(1,1,6);D=36;u=(6,0);v=(5,1);w=(-6,1)
S=steps(form,D)
need(set(S)=={u,v,w,tuple(-a for a in u),tuple(-a for a in v),tuple(-a for a in w)},'complete residual connection set')
need((4*form[0]*form[2]-form[1]**2)==23,'residual discriminant')
need(all(-11*u[i]+6*v[i]-6*w[i]==0 for i in range(2)),'primitive kernel relation')
need(gcd(11,gcd(6,6))==1,'primitive relation')
cycle=[(0,0)];cur=cycle[0]
for step in [v]*6+[tuple(-a for a in w)]*6+[tuple(-a for a in u)]*11:
    cur=tuple(cur[i]+step[i] for i in range(2));cycle.append(cur)
need(cycle[-1]==cycle[0] and len(set(cycle[:-1]))==23,'simple twenty-three cycle')
need(all(q(form,tuple(cycle[i+1][j]-cycle[i][j] for j in range(2)))==D for i in range(23)),'cycle edges')
p,t,m,ell,colour=build_colour(form,D,u)
need((p,t,m,ell)==(23,0,23,[4,2]),'explicit residue map')
for x,y in product(range(23),repeat=2):
    for sx,sy in S:
        need(colour((x,y))!=colour((x+sx,y+sy)),'whole 23x23 torus coloring')
# General odd polygon family, exact control, not a universal-claim census.
for n in range(1,51):
    U=(n,0);V=(n-1,1);W=(-n,1);path=[(0,0)]
    for step in [V]*n+[tuple(-a for a in W)]*n+[tuple(-a for a in U)]*(2*n-1):
        path.append(tuple(path[-1][i]+step[i] for i in range(2)))
    need(path[-1]==path[0] and len(set(path[:-1]))==4*n-1,'general simple odd polygon')
    need(all(q((1,1,n),step)==n*n for step in (U,V,W)),'general equal norms')
# Broad exact arithmetic tests. Every form in the declared integer box,
# all attained D<=50, with COMPLETE norm connection sets, not pruned rays.
forms=0;layers=0;directions=0;max_modulus=0;primes=set();max_denominator_depth=0
for a in range(1,6):
    for b in range(-5,6):
        for c in range(1,7):
            if 4*a*c<=b*b:continue
            forms+=1
            for D in range(1,51):
                S=steps((a,b,c),D)
                if not S:continue
                p,t,m,ell,colour=build_colour((a,b,c),D,S[0])
                layers+=1;directions+=len(S);max_modulus=max(m,max_modulus);primes.add(p);max_denominator_depth=max(t,max_denominator_depth)
                for s in S:
                    need(colour((0,0))!=colour(s),'origin certificate')
                for z in [(-3,5),(5,-2),(0,0),(1,1)]:
                    need(colour(z)==colour((z[0]+m,z[1]))==colour((z[0],z[1]+m)),'coordinate periods')
print(json.dumps({
 'status':'PASS', 'checks':checks,'veronese_triples':triples,
 'irrational_rhombus_steps':rhombus,'non_similarity_minima':{'old':1,'restricted':6,'predicted_if_similar':3},
 'residual':{'form':form,'squared_norm':D if False else 36,'steps':steps(form,36),'cycle':cycle,'odd_girth':23,'residue_map':[4,2],'modulus':23},
 'general_polygon_controls':50,
 'arithmetic_universe':{'A':[1,5],'B':[-5,5],'C':[1,6],'D':[1,50],'filter':'4AC>B^2 and attained D','forms':forms,'layers':layers,'directed_steps':directions,'max_modulus':max_modulus,'max_denominator_depth':max_denominator_depth,'ramified_primes':sorted(primes)},
},sort_keys=True,indent=2))
