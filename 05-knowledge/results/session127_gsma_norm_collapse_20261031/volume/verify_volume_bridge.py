#!/usr/bin/env python3
"""Finite exact checks for the Session 127 volume / norm-form bridge."""
import json
from pathlib import Path
import sympy as s

x, p, q = s.symbols('x p q')
OUT = Path(__file__).resolve().parent / 'verification.json'

def check(condition, detail='verification failed'):
    if not condition:
        raise AssertionError(detail)

def same(a, b):
    check(s.expand(a-b) == 0, (a, b))

def modular_irreducibility_certificate(f, prime):
    poly = s.Poly(f, x, modulus=prime)
    degree = poly.degree()
    base = s.Poly(x, x, modulus=prime)
    powers = {0: base}
    for j in range(1, degree+1):
        powers[j] = (powers[j-1] ** prime).rem(poly)
    final = (powers[degree]-base).rem(poly)
    check(final.is_zero)
    gcds = {}
    for ell in s.factorint(degree):
        g = s.gcd(poly, powers[degree//ell]-base).monic()
        check(g.degree() == 0)
        gcds[str(ell)] = str(g.as_expr())
    return {
        'prime': prime,
        'rabin_frobenius_remainders': {str(j):str(v.as_expr()) for j,v in powers.items()},
        'final_minus_x': str(final.as_expr()),
        'gcds_for_prime_divisors_of_degree': gcds,
        'irreducible': True,
    }

families = []
for k, prime, r1, disc in [(3,3,1,-44),(4,2,2,-563),(5,5,1,9584)]:
    f = x**k - sum(x**j for j in range(k))
    F = p**k - sum(p**j*q**(k-j) for j in range(k))
    C = s.zeros(k)
    for j in range(k-1):
        C[j+1,j] = 1
    for i in range(k):
        C[i,k-1] = 1
    same(C.charpoly(x).as_expr(), f)
    same((p*s.eye(k)-q*C).det(), F)
    same(s.resultant(f, p-q*x, x), F)
    same(C.trace(), 1)
    same((C*C).trace(), 3)
    same(C.det(), (-1)**(k+1))
    d1 = s.discriminant(f,x)
    trace_gram = s.Matrix(k,k,lambda i,j:(C**(i+j)).trace())
    d2 = trace_gram.det()
    check(d1 == d2 == disc)
    root_count = s.Poly(f,x).count_roots(-s.oo,s.oo)
    check(root_count == r1)
    check(f.subs(x,1) < 0 < f.subs(x,2))
    irreducibility = modular_irreducibility_certificate(f,prime)
    unit_controls=[]
    e0=s.eye(k)[:,0]
    for m in range(9):
        M = C**(2*m)
        check(M.det() == 1)
        unit_controls.append({'m':m,'theta_power':2*m,'coordinates':list(map(int,M*e0)), 'norm':1})
    # q=0 is the exact equality control for the metric/norm inequality.
    same(F.subs(q,0),p**k)
    V=s.simplify(s.pi**s.Rational(k,2)/s.gamma(s.Rational(k,2)+1))
    constant=s.simplify(2**k/(k**s.Rational(k,2)*V))
    roots=s.nroots(f,n=40,maxsteps=100)
    Snum=sum(abs(complex(root))**2 for root in roots)
    check(Snum>3)
    families.append({
        'k':k,'polynomial':str(f),'binary_norm_form':str(F),
        'multiplication_matrix':C.tolist(),
        'signature':[int(root_count),(k-int(root_count))//2],
        'polynomial_discriminant':int(d1),
        'independent_trace_gram_determinant':int(d2),
        'trace_theta':1,'trace_theta_squared':3,
        'modular_irreducibility_certificate':irreducibility,
        'unit_power_controls':unit_controls,
        'unit_ball_volume':str(V),'ball_norm_bound_constant':str(constant),
        'numerical_S_k_not_an_exact_certificate':Snum,
    })

f3=x**3-x**2-x-1
F3=p**3-p**2*q-p*q**2-q**3
check(F3.subs({p:103,q:56}) == -1)
check(f3.subs(x,s.Rational(9,5)) < 0 < f3.subs(x,2))
# tau in (9/5,2) and S=tau^2+2/tau imply S>81/25+1.
S_lower=s.Rational(106,25)
Q_lower=20291+3136*S_lower
check(Q_lower == s.Rational(839691,25) and Q_lower>180**2)
# Gram determinant in symbols, checked separately from the general prose proof.
K,S,T=s.symbols('K S T', positive=True)
G=s.Matrix([[K,-T],[-T,S]])
same(G.det(),K*S-T**2)
same((s.Matrix([p,q]).T*G*s.Matrix([p,q]))[0],K*p*p-2*T*p*q+S*q*q)
for n in range(16):
    Vn=s.pi**s.Rational(n,2)/s.gamma(s.Rational(n,2)+1)
    Vnext=s.pi**s.Rational(n+2,2)/s.gamma(s.Rational(n+2,2)+1)
    check(s.simplify(Vnext-2*s.pi*Vn/(n+2))==0)
# Unit theta^2 is outside the binary span in all three power bases.
check(all(row['unit_power_controls'][1]['coordinates'][2] == 1 for row in families))

result={
    'status':'PASS',
    'scope':'Exact k=3,4,5 reconstructed k-nacci polynomial family; symbolic p,q; unit powers 0,2,...,16; volume recurrence n=0,...,15.',
    'sympy_version':s.__version__,
    'families':families,
    'tribonacci_hostile':{
        'p':103,'q':56,'exact_norm':-1,
        'exact_Q_expression':'20291+3136*(tau**2+2/tau)',
        'tau_strict_lower_bound':'9/5','tau_strict_upper_bound':'2',
        'Q_strict_lower_bound':str(Q_lower),
        'euclidean_length_strict_lower_bound':180,
    },
    'independent_paths':[
        'Norm: resultant vs multiplication-matrix determinant vs explicit homogeneous polynomial',
        'Discriminant: polynomial discriminant vs trace-power Gram determinant',
        'Irreducibility: exact Rabin modular certificates',
        'Signature: exact real-root count via Sturm method',
    ],
    'not_claimed':['All Thue-Mahler solutions','Effective binary height bounds','Maximal-order discriminants','A recovered historic GSMA axiom system','A completed old normcollapse provenance link'],
}
OUT.write_text(json.dumps(result,indent=2,default=str)+'\n')
print('PASS: exact norm, discriminant, irreducibility, signature, unit and volume controls; '+str(OUT))
