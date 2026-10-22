#!/usr/bin/env python3
"""Session31: exact algebraic controls plus explicitly numerical parameter checks.

No lattice scaling-limit or MRB irrationality assertion is tested here.
Run: python3 check_categorical_parameters.py > categorical_checks.json
Dependencies: sympy, mpmath (already present in the research environment).
"""
import json
import mpmath as mp
import sympy as sp

mp.mp.dps = 90
x = sp.symbols('x')
beta = 2**(-sp.Rational(3,8))*sp.sqrt(2-sp.sqrt(2))
minimal = sp.Poly(sp.minpoly(beta, x), x)
reciprocal = sp.Poly(x**8-8*x**4-2, x)
assert minimal.as_expr() == 2*x**8+8*x**4-1
assert reciprocal.is_irreducible
# Independent exact Eisenstein check at 2.
assert reciprocal.LC() == 1
assert all(int(a) % 2 == 0 for a in reciprocal.all_coeffs()[1:])
assert int(reciprocal.TC()) % 4 != 0

N = sp.zeros(7)
for j in range(6):
    N[j,j+1] = N[j+1,j] = 1
charpoly = sp.factor(N.charpoly(x).as_expr())
assert charpoly == x*(x*x-2)*(x**4-4*x*x+2)
mu8 = sp.sqrt(2+sp.sqrt(2))
assert sp.simplify(-mu8**4+4*mu8**2-2) == 0
# TL_2(0) regular representation: a nonzero nilpotent radical generator.
e = sp.Matrix([[0,0],[1,0]])
assert e != sp.zeros(2) and e*e == sp.zeros(2)
# Finite order alone does not constrain supertrace without parity preservation.
a = sp.symbols('a')
G = sp.Matrix([[a,1],[1-a*a,-a]])
parity = sp.diag(1,-1)
assert sp.simplify(G*G) == sp.eye(2)
assert sp.trace(G) == 0
assert sp.trace(parity*G) == 2*a
assert sp.simplify(parity*G-G*parity) != sp.zeros(2)
# Branch names used in the note apply for m >= 4, not m = 3.
g_first_at3 = 2-sp.Rational(4,3)
g_second_at3 = sp.Rational(4,3)
assert g_first_at3 < 1 < g_second_at3
assert 2-sp.Rational(4,4) == sp.Rational(4,4) == 1

rows=[]
for m in range(4,41):
    assert 0 < sp.Rational(4,m) <= 1 <= 2-sp.Rational(4,m) < 2
    th=mp.pi/m
    mu=2*mp.cos(th)
    n=-2*mp.cos(4*th)
    A=2/mp.pi*mp.power(m,-mp.mpf(1)/m)*mp.sin(th)
    S00=mp.sqrt(mp.mpf(2)/m)*mp.sin(th)
    xc=1/mp.sqrt(2+mp.sqrt(2-n))
    errors={
       'chebyshev_identity': abs(n-(-mu**4+4*mu**2-2)),
       'nienhuis_candidate_identity': abs(xc-1/mu),
       'normalized_s00_identity': abs(A-mp.sqrt(2*m)*mp.power(m,-mp.mpf(1)/m)*S00/mp.pi),
    }
    assert max(errors.values())<mp.mpf('1e-85')
    rows.append({'m':m,'n':mp.nstr(n,40),'A':mp.nstr(A,40),
                 'mu':mp.nstr(mu,40),'S00':mp.nstr(S00,40),
                 'max_identity_error':mp.nstr(max(errors.values()),8)})
qtheta=mp.exp(-mp.pi*mp.sqrt(2))
report={
 'status':'Exact algebra and finite numerical identity checks only',
 'beta_definition':'beta = pi*A(8)',
 'beta_primitive_minimal_polynomial':str(minimal.as_expr()),
 'beta_inverse_polynomial':str(reciprocal.as_expr()),
 'independent_eisenstein_prime':2,
 'beta_is_algebraic_integer':False,
 'A7_characteristic_polynomial':str(charpoly),
 'A7_endpoint_path_counts_n0_through_8':[int((sp.ones(1,7)*N**j*sp.eye(7)[:,0])[0]) for j in range(9)],
 'TL2_0_nonzero_nilpotent':True,
 'supertrace_hypothesis_regression':{
   'matrix':'G(a) = [[a,1],[1-a^2,-a]] on 1|1',
   'square_is_identity':True,
   'ordinary_trace':str(sp.trace(G)),
   'supertrace':str(sp.trace(parity*G)),
   'parity_preserving':False,
   'lesson':'Finite-order supertrace field control requires parity-preserving twining.'
 },
 'branch_label_regression':{
   'valid_named_branch_domain':'m >= 4',
   'branches_meet_at_m4_g':1,
   'm3_2_minus_4_over_m':str(g_first_at3),
   'm3_4_over_m':str(g_second_at3),
   'lesson':'At m=3 the two formulas exchange the dense/dilute ranges.'
 },
 'qtheta_at_i_sqrt2':mp.nstr(qtheta,60),
 'qcharacter_at_i_sqrt2':mp.nstr(qtheta**2,60),
 'm8_A':rows[4]['A'],
 'parameter_checks':rows,
 'limitations':[
   'Nienhuis x_c is evaluated as an algebraic candidate; no criticality claim for general n follows.',
   'The m-grid checks are numerical; general identities are proved in the note.',
   'No assertion about the arithmetic nature of C_MRB or A(8)-C_MRB.',
   'No character or tensor functor realizing the MRB defect was found or constructed.'
 ]
}
print(json.dumps(report,indent=2))
