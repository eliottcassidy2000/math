#!/usr/bin/env python3
"""Exact controls for the universal reciprocal-cubic power annihilator.
Requires SymPy. No inherited script import and no floating-point arithmetic.
"""
import sympy as s
from itertools import product
T, X = s.symbols('t x')
z = s.symbols('z0:3')
checks = 0

def check(label, condition):
    global checks
    assert condition, label
    checks += 1
    print('PASS', label)

def F(t,d):
    L=[s.Integer(2),t]
    for r in range(2,d+1): L.append(s.expand(t*L[-1]-L[-2]))
    f=X-(-1)**d
    for r in range(1,d+1): f*=X**2-(-1)**(d-r)*L[r]*X+1
    return s.Poly(s.expand(f),X)

def companion(t):
    return s.Matrix([[0,1,0],[0,0,1],[-1,t-1,t-1]])

def symmetric(t,d):
    """Rows express degree-d monomials in the coordinates after C_t acts."""
    exps=[(a,b,d-a-b) for a in range(d+1) for b in range(d-a+1)]
    image=(z[1],z[2],-z[0]+(t-1)*z[1]+(t-1)*z[2])
    rows=[]
    for e in exps:
        p=s.Poly(s.prod(v**k for v,k in zip(image,e)),*z)
        rows.append([p.coeff_monomial(s.prod(v**k for v,k in zip(z,q))) for q in exps])
    return s.Matrix(rows)

def evaluate(f,M,mod=None):
    V=s.zeros(M.rows); I=s.eye(M.rows)
    for a in f.all_coeffs():
        V=V*M+a*I
        V=V.applyfunc(lambda v: s.expand(v) if mod is None else v%mod)
    return V

for d in range(4):
    M=symmetric(T,d)
    check(f'symbolic Z[t] symmetric d={d}, dimension={M.rows}', evaluate(F(T,d),M)==s.zeros(M.rows))
for t,d in product([-3,-2,-1,0,1,2,3,83],range(7)):
    M=symmetric(t,d)
    check(f'integer symmetric t={t}, d={d}',evaluate(F(t,d),M)==s.zeros(M.rows))
for t,d,m in product([-2,0,2,3,83],range(1,5),[2,4,6,8,9,25]):
    M=symmetric(t,d)
    check(f'modular symmetric t={t}, d={d}, modulus={m}',evaluate(F(t,d),M,m)==s.zeros(M.rows))
for t,d in product([-2,2,3,83],range(1,4)):
    M=s.ones(1)
    for j in range(d): M=s.kronecker_product(M,companion(t))
    check(f'full tensor t={t}, d={d}, dimension={M.rows}',evaluate(F(t,d),M)==s.zeros(M.rows))
# Independent expansion by multiplying numerator by inverse denominator.
def sequence(num,last=14):
    u=[]
    for n in range(last+1):
        u.append((num[n] if n<len(num) else 0)+sum(c*u[n-j] for j,c in [(1,82),(2,82),(3,-1)] if n>=j))
    return u
A,B,C=[sequence(num) for num in ([1,53,9],[2,-26,-12],[2,8,-10])]
f=F(83,3)
print('F_3(t=83), factored:', s.factor(f.as_expr()))
print('F_3(t=83), coefficients descending:', f.all_coeffs())
check('factor coefficients independently computed',s.factor(f.as_expr())==(X+1)*(X**2-83*X+1)*(X**2+6887*X+1)*(X**2-571538*X+1))
print('Seven initial triples and actual cubic values:')
for n in range(7):
    value=A[n]**3+B[n]**3-C[n]**3
    print(n,A[n],B[n],C[n],value)
    check(f'initial residual n={n}',value-(-1)**n==0)
for name,u in [('A',A),('B',B),('C',C)]:
    cubes=[v**3 for v in u]
    for n in range(len(cubes)-7):
        check(f'nonzero cube sequence {name}, recurrence window {n}',sum(c*cubes[n+7-i] for i,c in enumerate(f.all_coeffs()))==0)
back=tuple(82*u[1]+82*u[0]-u[2] for u in [A,B,C])
check('backward triple n=-1',back==(-9,12,10))
check('backward cubic residual',back[0]**3+back[1]**3-back[2]**3==-1)
check('1729 rearrangement',9**3+10**3==12**3+1==1729)
perturbed=sequence([2,53,9])
check('hostile changed numerator detected',any(perturbed[n]**3+B[n]**3-C[n]**3-(-1)**n !=0 for n in range(7)))
# Dropping a genuine spectral factor must fail at generic t: not only a zero-residual test.
M=symmetric(83,3)
wrong=s.Poly(s.div(f.as_expr(),X+1,X)[0],X)
check('hostile dropped central factor detected',evaluate(wrong,M)!=s.zeros(M.rows))
print('ALL',checks,'EXACT CHECKS PASS')
