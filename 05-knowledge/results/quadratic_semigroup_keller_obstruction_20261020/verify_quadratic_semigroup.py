#!/usr/bin/env python3
"""Exact corroborating controls; the all-m theorem is proved in QUADRATIC_SEMIGROUP_KELLER_OBSTRUCTION.md."""
import json
import math
from pathlib import Path
import sympy as S

x, y, t, s = S.symbols('x y t s')
checks = []
def check(name, condition):
    if not bool(condition):
        raise RuntimeError('FAILED: ' + name)
    checks.append(name)
def zero(v):
    return S.expand(v) == 0

def jac(p, q, z=y):
    return S.expand(S.diff(p,x)*S.diff(q,z)-S.diff(p,z)*S.diff(q,x))

# Arbitrary lower coefficients represented by independent functions of x.
for n in range(1, 9):
    for d in range(1, 9):
        a = S.symbols('a', nonzero=True)
        b = S.Function('b')(x)
        p = a*y**n + sum(S.Function('p'+str(i))(x)*y**i for i in range(n))
        q = b*y**d + sum(S.Function('q'+str(i))(x)*y**i for i in range(d))
        check(f'top_coefficient_{n}_{d}', zero(jac(p,q).coeff(y,n+d-1)+n*a*S.diff(b,x)))

for k in range(1, 11):
    cs = S.symbols('c0:'+str(k))
    p = (y+s)**(2*k)+sum(cs[i]*(y+s)**(2*i) for i in range(k))
    check(f'even_center_coefficient_{k}', zero(S.expand(p).coeff(y,2*k-1)-2*k*s))

for k in range(1, 5):
    coeffs = [S.Function('a'+str(i))(x) for i in range(k)]
    p = t**(2*k)+sum(coeffs[i]*t**(2*i) for i in range(k))
    qcoeffs = [S.Function('b'+str(i))(x) for i in range(2*k+2)]
    q = sum(qcoeffs[i]*t**i for i in range(2*k+2))
    check(f'line_unit_identity_{k}', zero(jac(p,q,t).subs(t,0)-S.diff(coeffs[0],x)*qcoeffs[1]))

# The finite necessary-degree filter is not an embedding realization census.
terminal_candidates = {}
for m in range(3, 42, 2):
    gamma = [d for d in range(1,m+1) if d%2==0 or d==m]
    pairs = [(a,b) for a in gamma for b in gamma if a<b and a%b!=0 and b%a!=0 and m in (a,b)]
    wanted = [(n,m) for n in range(2,m,2)]
    check(f'necessary_terminal_degrees_{m}', pairs==wanted)
    terminal_candidates[str(m)] = pairs

# Even coordinate degrees do not force an imprimitive function field.
p, q = t**4, t**6+t
r = S.expand(q*q-p*p*p)
check('even_even_cancellation', r==2*t**7+t**2)
basis_degrees = [S.degree(z,t) for z in (S.Integer(1),q,r,S.expand(q*r))]
check('even_even_distinct_residues', sorted(int(d%4) for d in basis_degrees)==[0,1,2,3])
check('even_even_primitive_field', math.gcd(4,int(S.degree(r,t)))==1)
check('even_even_delta_five', sum(basis_degrees)-sum(range(4))==20)
gamma_hostile = {4*a+6*b+7*c for a in range(20) for b in range(20) for c in range(20)}
check('even_even_semigroup_gaps', [d for d in range(30) if d not in gamma_hostile]==[1,2,3,5,9])

# Already-known saturated hostiles: their displayed generating pairs are not Keller.
check('saturated_quadratic_hostile', jac(y**3,x*y**2+y)==-3*y**4)
check('rational_center_pole', zero((y**2+y/x).subs(y,t-1/(2*x))-(t**2-1/(4*x*x))))
check('saturated_cubic_hostile', jac(y**4,x*y**3+y**2+y)==-4*y**6)
u=S.Function('u')(x)
c=S.symbols('c', nonzero=True)
check('moving_quadratic_hostile', zero(jac(y**2,y*(y**4-u*y**2+c))-2*S.diff(u,x)*y**4))
P=y**4+x*(4*y*y-2*y+1)
Q=y*y+y+2*x
check('finite_jet_hostile', zero(jac(P,Q)-(1+x*(4-16*y))))

# Positive boundary: an automorphism with an even coordinate after a nonlinear shear.
P=x+(y+x*x)**2
Q=y+x*x
check('nonlinear_shear_positive_jacobian', jac(P,Q)==1)
check('nonlinear_shear_actual_line', zero(P.subs(y,t-x*x)-(x+t*t)))
check('nonlinear_shear_line_injectivity', P.subs(y,-x*x)==x)

result = {
    'status':'PASS',
    'exact_check_count':len(checks),
    'universes':{
        'top_coefficient_degrees':'1<=n,d<=8 with independent symbolic lower-coefficient functions',
        'center_coefficient_even_degrees':'2,4,...,20 with independent symbolic lower even coefficients',
        'line_identity_even_degrees':'2,4,6,8 with independent symbolic coefficient functions',
        'necessary_terminal_degree_filter':'odd m=3,5,...,41; necessary conditions only',
    },
    'checks':checks,
    'even_even_hostile':{'p':'t^4','q':'t^6+t','residual':'2t^7+t^2','basis_degrees':[int(d) for d in basis_degrees],'semigroup':[4,6,7],'gaps':[1,2,3,5,9]},
    'necessary_terminal_candidates':terminal_candidates,
    'scope':'Finite symbolic controls corroborate identities; QUADRATIC_SEMIGROUP_KELLER_OBSTRUCTION.md proves the all-m theorem using cited characteristic bounds and the one-line theorem.'
}
print(json.dumps(result,indent=2,sort_keys=True))
