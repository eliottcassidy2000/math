#!/usr/bin/env python3
"""Independent polynomial-constraint and flat-connection controls."""
import json
import sympy as S
x,y=S.symbols('x y',real=True); I=S.I
checks={}
# Solve real coefficient equations, not polyanalytic basis-generation.
counts=[]
for m in range(1,5):
    mon=[x**a*y**b for a in range(2*m+1) for b in range(2*m+1-a)]
    cols=[]
    for t in mon:
        for _ in range(m):t=S.expand((S.diff(t,x)+I*S.diff(t,y))/2)
        cols.append(t)
    terms=set()
    for t in cols: terms.update(S.Poly(t,x,y).monoms())
    rows=[]
    for ab in sorted(terms):
        co=[S.Poly(t,x,y).coeff_monomial(ab) for t in cols]
        rows.extend([[S.re(v) for v in co],[S.im(v) for v in co]])
    nullity=len(mon)-S.Matrix(rows).rank()
    if nullity!=m*m:raise RuntimeError(('polynomial kernel',m,nullity))
    counts.append([m,len(mon),nullity])
checks['real_polynomial_kernels']=counts
# Exact flatness for arbitrary symbolic CR-compatible first jets.
H,K,Hx,Hy=S.symbols('H K Hx Hy',real=True)
A=S.zeros(5); B=S.zeros(5)
A[0,1]=1; A[1,3]=1; A[4,1]=-H; A[4,2]=K
B[0,2]=1; B[2,3]=1; B[4,1]=K; B[4,2]=H
Ay=A.diff(H)*Hy+A.diff(K)*Hx
Bx=B.diff(H)*Hx-B.diff(K)*Hy
if Ay+A*B-Bx-B*A!=S.zeros(5):raise RuntimeError('flatness')
checks['symbolic_CR_flatness']=True
# Periods via Laurent coefficient extraction after clearing order-four poles.
z=S.symbols('z'); outputs=[]
for cs,want in [([1],0),([I],1),([1,I],2),([1,I,1+I],3)]:
    g=1+sum(c/(z-j)**2 for j,c in enumerate(cs));h=g*g
    rows=[]
    for j,c in enumerate(cs):
        # Coefficient (z-j)^-1 = third Taylor derivative of cleared quartic /6.
        hc=S.cancel((z-j)**4*h)
        r=S.simplify(S.diff(hc,z,3).subs(z,j)/6)
        t=S.simplify(S.diff(z*hc,z,3).subs(z,j)/6)
        row=[S.im(t),S.im(r),-S.re(r)]
        rows.append(row)
    mat=S.Matrix(rows);rank=mat.rank()
    if rank!=want:raise RuntimeError('period rank')
    outputs.append({'coefficients':list(map(str,cs)),'rank':rank,
                    'determinant':str(mat.det()) if len(cs)==3 else None})
checks['independent_period_ranks']=outputs
print(json.dumps({'status':'PASS','checks':checks},indent=2,sort_keys=True))
