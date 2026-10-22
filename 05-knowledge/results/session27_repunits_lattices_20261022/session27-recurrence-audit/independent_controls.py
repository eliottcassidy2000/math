#!/usr/bin/env python3
"""Independent audit: no imports from reviewed verification script."""
import sympy as s
from itertools import product
x,t,e=s.symbols('x t epsilon')
checks=0
def ck(label,ok):
 global checks
 assert ok,label
 checks+=1
 print('PASS',label,flush=True)
def candidate(d):
 # Chebyshev formula, independently of the author's Lucas-recursion code.
 return s.Poly(s.expand((x-(-1)**d)*s.prod(x*x-(-1)**(d-r)*2*s.chebyshevt(r,t/2)*x+1 for r in range(1,d+1))),x)
def tensor(d,v):
 C=s.Matrix([[0,1,0],[0,0,1],[-1,v-1,v-1]])
 M=s.eye(1)
 for _ in range(d):M=s.kronecker_product(M,C)
 return M
def ev(p,M,reduce=s.expand):
 V=s.zeros(M.rows)
 for a in p.all_coeffs():V=(V*M+a*s.eye(M.rows)).applyfunc(reduce)
 return V
for d in range(4):
 M=tensor(d,t)
 ck(f'full tensor symbolic Z[t] d={d}',ev(candidate(d),M)==s.zeros(M.rows))
 for v in [-2,2]:
  expected=(x-(-1)**d)**(2*d+1) if v==-2 else (x-1)**(d+1)*(x+1)**d
  ck(f'degenerate multiplicities t={v} d={d}',s.expand(candidate(d).as_expr().subs(t,v)-expected)==0)
for d in [1,2,3]:
 M=tensor(d,83)
 cp=M.charpoly(x).as_poly()
 ck(f'generic full tensor squarefree characteristic polynomial d={d}',cp.sqf_part()==s.Poly(candidate(d).as_expr().subs(t,83),x))
for v,m,d in product([e,2+e,-2+e],[4,6],[2,3]):
 def red(a):
  a=s.Poly(s.expand(a),e)
  return (a.nth(0)%m)+(a.nth(1)%m)*e
 M=tensor(d,v).applyfunc(red)
 p=s.Poly(candidate(d).as_expr().subs(t,v),x)
 ck(f'nilpotent ring Z/{m}[epsilon]/epsilon^2 t={v} d={d}',ev(p,M,red)==s.zeros(M.rows))
f=candidate(3).as_expr().subs(t,83)
expected=[1,-564733,-3889880904,322812677460,322812677460,-3889880904,-564733,1]
ck('explicit seventh degree coefficients',s.Poly(f,x).all_coeffs()==expected)
z=s.symbols('z')
seq=[]
for numerator in [1+53*z+9*z*z,2-26*z-12*z*z,2+8*z-10*z*z]:
 expansion=s.Poly(s.series(numerator/(1-82*z-82*z*z+z**3),z,0,7).removeO(),z)
 seq.append([expansion.nth(n) for n in range(7)])
for n in range(7):
 vals=[a[n] for a in seq]
 ck(f'initial residual n={n}',vals[0]**3+vals[1]**3-vals[2]**3==(-1)**n)
 print('TRIPLE',n,vals)
back=tuple(82*a[0]+82*a[1]-a[2] for a in seq)
ck('negative index triple',back==(-9,12,10))
ck('negative index residual',back[0]**3+back[1]**3-back[2]**3==-1)
ck('1729',9**3+10**3==12**3+1==1729)
for v in [-2,2]:
 M=tensor(3,v)
 squarefree=s.Poly(candidate(3).as_expr().subs(t,v),x).sqf_part()
 ck(f'hostile discard multiplicities t={v}',ev(squarefree,M)!=s.zeros(27))
# Direct nonhomogeneous counterexample at d=1: constant 1 + eigenmode (-1)^n.
ck('hostile nonhomogeneous observable',sum(c*(1+(-1)**(3-i)) for i,c in enumerate(s.Poly(candidate(1).as_expr().subs(t,83),x).all_coeffs()))!=0)
# Both published K3 sections: an independent polynomial sanity check, not rank computation.
k=63*(3*t*t-3*t+1)*(t*t+t+1)*(t*t-3*t+3)
for a,b in [(6*t*t-4*t+4,-3*t*t-5*t+5),(4*t*t-4*t+6,5*t*t-5*t-3)]:
 ck('K3 cited section lies on family',s.expand(a**3+b**3-k)==0)
ck('K3 k has six simple roots',s.degree(k,t)==6 and s.degree(s.gcd(k,s.diff(k,t)),t)==0)
print('ALL',checks,'INDEPENDENT CHECKS PASS')
