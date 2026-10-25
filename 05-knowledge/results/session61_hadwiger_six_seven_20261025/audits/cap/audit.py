#!/usr/bin/env python3
"""Independent exact audit. Does not import or execute producer code.
Universal conclusions are proved in AUDIT.md; bounded checks are regressions.
"""
from fractions import Fraction as R
from math import gcd, isqrt, lcm
from itertools import combinations, product
from collections import Counter, deque
import json

checks=0
def require(c,why):
    global checks
    checks+=1
    if not c: raise RuntimeError(why)
def ordp(x,p):
    x=R(x)
    if not x:return None
    e=0
    for n,sgn in [(abs(x.numerator),1),(x.denominator,-1)]:
        while n%p==0:n//=p;e+=sgn
    return e
def integral(x,p):return x.denominator%p!=0
def residue(x,m):
    x=R(x);return x.numerator*pow(x.denominator,-1,m)%m
def squarefree_parts(z):
    d=1;s=1;k=2
    while k*k<=z:
        q=0
        while z%k==0:z//=k;q+=1
        d*=k**(q%2);s*=k**(q//2);k+=1
    return d*z,s
def ramified(d):
    out=[]
    if d%4 in (1,2):out.append((2,int(d%4==1)))
    z=d
    if z%2==0:z//=2
    p=3
    while p*p<=z:
        if z%p==0:
            out.append((p,0))
            while z%p==0:z//=p
        p+=2
    if z>2 and z%2:out.append((z,0))
    return out
def chosen(d):return (2,int(d%4==1)) if d%4 in (1,2) else ramified(d)[0]
def cycle_colour(k,p):return k if p==2 else (2 if k==p-1 else k%2)
def q(F,z):
    a,b,c=F;x,y=z;return a*x*x+b*x*y+c*y*y

def complete_steps(F,h):
    # Solve the x-quadratic exactly for each possible y. This is independent
    # of the producer's two-dimensional bounding-box scan.
    a,b,c=F;disc=4*a*c-b*b;by=isqrt((4*a*h)//disc);out=set()
    for y in range(-by,by+1):
        delta=4*a*h-disc*y*y;r=isqrt(delta)
        if r*r!=delta:continue
        for sign in [-1,1]:
            numer=-b*y+sign*r
            if numer%(2*a)==0:out.add((numer//(2*a),y))
    require(all(q(F,z)==h for z in out),'exact quadratic roots')
    return sorted(out)

def certificate(F,h,u,p_shift=None):
    a,b,c=F;d,s=squarefree_parts(4*a*c-b*b)
    p,j=chosen(d) if p_shift is None else p_shift
    # Direct real and imaginary coefficients after division by the edge u.
    ax=R(2*a*u[0]+b*u[1],2*a);ay=R(s*u[1],2*a);N=R(h,a)
    require(ax*ax+d*ay*ay==N,'normalizing norm')
    basis=[]
    for x,y in [(1,0),(0,1)]:
        X=R(2*a*x+b*y,2*a);Y=R(s*y,2*a)
        real=(X*ax+d*Y*ay)/N;imag=(Y*ax-X*ay)/N
        basis.append((real-j*imag,imag))
    depths=[0]+[-ordp(A,p) for A,B in basis if A]
    t=max(depths);pt=p**t;m=p*pt
    ell=[residue(pt*A,m) for A,B in basis]
    for v in complete_steps(F,h):
        A=sum(basis[i][0]*v[i] for i in range(2));B=sum(basis[i][1]*v[i] for i in range(2))
        require(A*A+2*j*A*B+(j*j+d)*B*B==1,'Eisenstein norm equation')
        require(integral(A,p) and integral(B,p),'both edge coefficients integral')
        require(residue(A,p) in {1,p-1},'edge coefficient is plus/minus one')
        require(sum(ell[i]*v[i] for i in range(2))%m in {pt,(-pt)%m},'scaled edge projection')
    return p,t,m,ell

# Independent rational unit-circle parametrization. Check every ramified
# prime, not only the producer-selected prime, on the exact finite universe.
unit_tests=0;ds=[]
for d in range(1,129):
    if squarefree_parts(d)[1]!=1:continue
    ds.append(d)
    for p,j in ramified(d):
        require((d+j*j)%p==0 and (d+j*j)%(p*p)!=0 and 2*j%p==0,'Eisenstein test')
        for a in range(-12,13):
            for b in range(1,13):
                if gcd(a,b)!=1:continue
                denom=b*b+d*a*a;X=R(b*b-d*a*a,denom);Y=R(2*a*b,denom)
                A=X-j*Y;B=Y
                require(X*X+d*Y*Y==1,'rational parametrization')
                require(integral(A,p) and integral(B,p),'local unit parametrization')
                require(residue(A,p) in {1,p-1},'local residue parametrization')
                unit_tests+=1

# Independently count the producer universe and vary the normalizing edge
# over every norm vector, rather than only one selected vector per layer.
forms=layers=directed=normalizers=0;max_t=max_m=0;target_quotients=set()
for a in range(1,6):
    for b in range(-5,6):
        for c in range(1,7):
            if 4*a*c-b*b<=0:continue
            forms+=1
            for h in range(1,51):
                S=complete_steps((a,b,c),h)
                if not S:continue
                layers+=1;directed+=len(S)
                for u in S:
                    p,t,m,ell=certificate((a,b,c),h,u)
                    normalizers+=1;max_t=max(t,max_t);max_m=max(m,max_m)
                    target_quotients.add((p,t))
require((forms,layers,directed)==(280,4798,18848),'independent complete census')

# Deliberately deep p-denominators: cover all scalar target vertices, including
# the wrap edge and negative representatives, not merely a source patch.
deep=[]
for F,p,k in [((1,0,1),2,12),((1,0,2),2,10),((1,1,1),3,6),((1,1,2),7,3),((1,1,6),23,2)]:
    h=p**(2*k);u=(p**k,0);pp,t,m,ell=certificate(F,h,u)
    require(pp==p and t>=k,'deep denominator retained')
    deep.append({'form':F,'h':h,'normalizer':u,'p':p,'t':t,'m':m,'ell':ell,'steps':len(complete_steps(F,h))})
    target_quotients.add((p,t))
scalar_edges=0
for p,t in sorted(target_quotients):
    pt=p**t;m=p*pt
    for r in range(-m,m):
        for e in [-pt,pt]:
            require(cycle_colour((r%m)//pt,p)!=cycle_colour(((r+e)%m)//pt,p),'complete cyclic target colouring')
            scalar_edges+=1

# Residual: complete six generators, exact integer-kernel argument regression,
# and an independent breadth-first odd-cycle detection by parity lift.
F=(1,1,6);S=complete_steps(F,36)
require(set(S)=={(6,0),(-6,0),(5,1),(-5,-1),(-6,1),(6,-1)},'residual connection set')
front={(0,0)};first_odd_return=None
for length in range(1,24):
    front={(x+u,y+v) for x,y in front for u,v in S}
    if length%2 and (0,0) in front:
        first_odd_return=length;break
require(first_odd_return==23,'independent odd return length')
require(certificate(F,36,(6,0))==(23,0,23,[4,2]),'explicit residual certificate')

# General simple polygon, and rational-basis denominator clearing checked
# with three nontrivial rational lattices in each odd-discriminant field.
polygons=embeddings=0
matrices=[((R(2,3),R(1,5)),(R(1,7),R(3,2))),((R(4),R(1)),(R(2),R(3))),((R(1,2),R(-2,3)),(R(3,5),R(7,11)))]
for n in range(1,101):
    U=(n,0);V=(n-1,1);W=(-n,1)
    path=[(0,0)]
    for dx,dy in [V]*n+[(-W[0],-W[1])]*n+[(-U[0],-U[1])]*(2*n-1):
        x,y=path[-1];path.append((x+dx,y+dy))
    require(path[-1]==(0,0) and len(set(path[:-1]))==4*n-1,'simple odd polygon')
    require(all(q((1,1,n),(path[i+1][0]-path[i][0],path[i+1][1]-path[i][1]))==n*n for i in range(4*n-1)),'equal polygon edges')
    # Dividing by n gives unit norm in the same field.
    require(all(q((R(1),R(1),R(n)),(R(x,n),R(y,n)))==1 for x,y in [U,V,W]),'unit normalization of polygon')
    polygons+=1
    for rows in matrices:
        (a,b),(c,d)=rows;det=a*d-b*c
        inv=((d/det,-b/det),(-c/det,a/det));M=lcm(*(x.denominator for row in inv for x in row))
        require(all((M*x).denominator==1 for row in inv for x in row),'scaled order contained in arbitrary rational lattice')
        for x,y in path:
            coeff=[M*(row[0]*x+row[1]*y) for row in inv]
            require(all(z.denominator==1 for z in coeff),'whole polygon in lattice')
        embeddings+=1

# Dense-flower cross-check in an independent basis theta=1+omega,
# theta^2-theta+1=0, N(a+b theta)=a^2+ab+b^2.
def tn(z):return z[0]*z[0]+z[0]*z[1]+z[1]*z[1]
def tm(z,w):
    a,b=z;c,d=w;return (a*c-b*d,a*d+b*c+b*d)
def plus(z,w):return z[0]+w[0],z[1]+w[1]
A=(2,1);B=(3,-1);units=[z for z in product(range(-1,2),repeat=2) if tn(z)==1]
points=sorted({plus(tm(A,u),tm(B,v)) for u in [(0,0),(1,0),(0,1)] for v in [(0,0)]+units})
edges={(i,j) for i,j in combinations(range(len(points)),2) if tn((points[i][0]-points[j][0],points[i][1]-points[j][1]))==7}
triangles=[(i,j,k) for i,j,k in combinations(range(len(points)),3) if {(i,j),(i,k),(j,k)}<=edges]
colors=[(a+2*b)%3 for a,b in points]
hist=dict(sorted(Counter(sum(i in e for e in edges) for i in range(len(points))).items()))
require(len(points)==21 and len(edges)==57 and hist=={5:18,8:3},'dense flower counts')
require(triangles and all(colors[i]!=colors[j] for i,j in edges),'dense flower exact chromatic three')
require(tm(A,B)==(7,0) and plus(A,B)==(5,0),'conjugate ideal norm and sum')
ratio=tuple(R(x,7) for x in tm(A,A))
require(ratio==(R(3,7),R(5,7)) and tn(ratio)==1,'exact relative rotation')
require(2*ratio[0]+ratio[1]==R(11,7),'nonintegral rational trace')
require(tuple(3*(A[i]+B[i])-2*tm(A,B)[i] for i in range(2))==(1,0),'ideal sum Bezout certificate')

result={'status':'PASS','checks':checks,'independent_census':{'forms':forms,'layers':layers,'directed_steps':directed,'normalizing_edges_tested':normalizers,'max_t':max_t,'max_m':max_m},'ramified_unit_tests':unit_tests,'squarefree_d_count':len(ds),'deep_denominator_controls':deep,'complete_scalar_target_edges':scalar_edges,'general_polygons':polygons,'scaled_order_embeddings':embeddings,'residual_odd_girth':first_odd_return,'residual_steps':S,'dense_flower':{'vertices':len(points),'edges':len(edges),'triangles':len(triangles),'degree_histogram':hist,'chi':3,'points_theta':points}}
print(json.dumps(result,sort_keys=True,indent=2))
