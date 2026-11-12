#!/usr/bin/env python3
"""Exact symbolic controls; proofs, not these finite checks, supply all-order scope."""
import json
from pathlib import Path
import sympy as s
z,w=s.symbols('z w')
a,u,v,c,d=s.symbols('a u v c d', real=True)
b=u+s.I*v
R=a*z*w+b*z+s.conjugate(b)*w+c
checks=0

def check(t, label):
    global checks
    if s.simplify(t)!=0:
        raise RuntimeError((label,t))
    checks+=1

def conjugate_z(expr):
    return s.conjugate(expr).subs(s.conjugate(z),w)

hs=[s.Integer(1),z**2,(1+z+z**2)**2,(1+s.I*z)**2,
    (2+(1+s.I)*z+3*z**2-s.I*z**3)**2,s.exp(2*z)]
for k,h in enumerate(hs):
    F=s.integrate(h*(a*z+s.conjugate(b)),z)
    p=d-F-conjugate_z(F)
    q=h*R
    check(s.diff(q,w)+s.diff(p,z),f'equilibrium:{k}')
    check(s.diff(p,z,w),f'compatibility:{k}')
    check(s.diff(q,w,2),f'finite type:{k}')
    check(R-conjugate_z(R.subs(w,s.conjugate(z))),f'reality:{k}')
# Direct Cartesian equation independently from complex conventions.
x,y=s.symbols('x y',real=True)
h=(1+(x+s.I*y)**2)**2
F=s.integrate((1+z**2)**2*(a*z+s.conjugate(b)),z).subs(z,x+s.I*y)
p=s.expand(d-2*s.re(F))
q=s.expand(h*(a*(x*x+y*y)+2*u*x-2*v*y+c))
sxx=s.expand((p+s.re(q))/2)
syy=s.expand((p-s.re(q))/2)
sxy=s.expand(-s.im(q)/2)
check(s.diff(sxx,x)+s.diff(sxy,y),'cartesian equilibrium x')
check(s.diff(sxy,x)+s.diff(syy,y),'cartesian equilibrium y')
check(s.diff(p,x,2)+s.diff(p,y,2),'cartesian compatibility')
# Hermitian determinant and rank-four measurements.
H=s.Matrix([[a,s.conjugate(b)],[b,c]])
xi=s.Matrix([z,1])
check((s.Matrix([[w,1]])*H*xi)[0]-R,'Hermitian carrier')
t,ss=s.symbols('t ss',real=True)
check(H.det().subs({a:t+ss,c:t-ss})-(t*t-ss*ss-u*u-v*v),'Lorentz determinant')
good=[(0,0),(1,0),(0,1),(2,2)]
bad=[(0,0),(1,0),(0,1),(1,1)]
def measurement(points): return s.Matrix([[xx*xx+yy*yy,xx,yy,1] for xx,yy in points])
if measurement(good).rank()!=4 or measurement(bad).rank()!=3:
    raise RuntimeError('measurement ranks')
checks+=2
# Flat five-state connection; Cauchy-Riemann derivatives retained explicitly.
H0,K0,Hx,Kx=s.symbols('H0 K0 Hx Kx')
Ax=s.zeros(5);Ay=s.zeros(5)
Ax[0,1]=1;Ax[1,3]=1;Ax[4,1]=-H0;Ax[4,2]=K0
Ay[0,2]=1;Ay[2,3]=1;Ay[4,1]=K0;Ay[4,2]=H0
curvature=Ax.diff(H0)*(-Kx)+Ax.diff(K0)*Hx-Ay.diff(H0)*Hx-Ay.diff(K0)*Kx+Ax*Ay-Ay*Ax
for value in curvature:check(value,'flat transport')
# Rank-one amplitude -> Pythagorean null vector, distinct from full stress.
m,n=s.symbols('m n',real=True)
check((m*m-n*n)**2+(2*m*n)**2-(m*m+n*n)**2,'rank-one Pythagorean carrier')
# Two-net superposition identity for arbitrary symbolic holomorphic inputs.
AA,BB=s.symbols('AA BB')
cc=s.symbols('cc',real=True)
check(AA*(z+w+cc)+(BB-(z+cc)*AA)-(w*AA+BB),'two-net decomposition')
# Discrete closure matrices: independent rank check on declared rectangles.
grids=[]
for N in range(3,9):
    for M in range(3,9):
        rows=[]
        def row(entries):
            r=[0]*(N*M)
            for xx,yy,coef in entries:r[xx*M+yy]+=coef
            rows.append(r)
        for xx in range(N-1):
            for yy in range(M-1):
                row([(xx+1,yy+1,1),(xx+1,yy,-1),(xx,yy+1,-1),(xx,yy,1)])
        for xx in range(N-2):
            for yy in range(M-2):
                row([(xx+2,yy,1),(xx+1,yy,-2),(xx,yy+2,-1),(xx,yy+1,2)])
        rank=s.Matrix(rows).rank()
        if rank!=N*M-4:raise RuntimeError(('discrete nullity',N,M,rank))
        grids.append([N,M,4]);checks+=1
# Translation and exact coefficient reconstruction as polynomial identities.
A,B,C,D,r,t=s.symbols('A B C D r t')
choose=lambda q:q*(q-1)/2
f=A*(choose(x)+choose(y))+B*x+C*y+D
new=A*(choose(x)+choose(y))+(B+A*r)*x+(C+A*t)*y+D+B*r+C*t+A*(choose(r)+choose(t))
check(f.subs({x:x+r,y:y+t},simultaneous=True)-new,'translation grammar')
check(f.subs({x:2,y:0})-2*f.subs({x:1,y:0})+f.subs({x:0,y:0})-A,'reconstruction curvature')
# Polyanalytic Hermitian carriers m=1..7 have m² real coordinates.
poly=[]
for m in range(1,8):
    basis=[z**j*w**j for j in range(m)]
    for i in range(m):
        for j in range(i):basis += [z**i*w**j+z**j*w**i,s.I*(z**i*w**j-z**j*w**i)]
    if len(basis)!=m*m:raise RuntimeError('Hermitian dimension')
    for term in basis:check(s.diff(term,w,m),'polyanalytic')
    poly.append([m,len(basis)])
print(json.dumps({'status':'PASS','checks':checks,'stress_phases':len(hs),
 'measurement_ranks':[4,3],'grid_nullities':grids,'polyanalytic_dimensions':poly},indent=2))
