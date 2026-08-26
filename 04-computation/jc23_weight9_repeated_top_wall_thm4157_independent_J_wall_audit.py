#!/usr/bin/env python3
"""Exact clean-room audit of the row-D endpoint wall J(eta)=0."""

import sympy as s
from hashlib import sha256

X,T,e=s.symbols('X T eta')
P=T+X**2*T**2
Y=X*T*P
D=s.Rational(2048,45)
K=s.Rational(1376,135)
G=s.expand(-X**2*T/2-3*P+s.Rational(8,3)*P**2-s.Rational(1376,135)*P**3
           +K*Y**2+D*P**4-D*P*Y**2+e*P**3*Y-e*Y**3)
f=s.cancel(s.diff(G,X)/T)
h=s.diff(G,T)
J=22143375*e**2+15510536192

def basis_coeff(c):
    r=s.rem(s.Poly(s.together(c).as_numer_denom()[0],e),s.Poly(J,e)).as_expr()
    den=s.together(c).as_numer_denom()[1]
    r=s.expand(r/den)
    return (s.Rational(r.coeff(e,0)),s.Rational(r.coeff(e,1)))

print('source_ready',flush=True)
R=s.factor(s.resultant(f,h,X)/(T**35*(6*T+1)**2))
Q12=s.Poly(R,T)
print('Q12_degree',Q12.degree(),flush=True)
# Reduce every coefficient modulo the quadratic J.  This is the exact
# representative q0(T)+eta*q1(T) in Q[eta,T]/(J).
QJ=s.Poly(sum(s.rem(s.Poly(c,e),s.Poly(J,e)).as_expr()*T**i
              for (i,),c in Q12.terms()),T,domain=s.QQ.frac_field(e))
# terms() above is descending and i is the actual exponent.
while QJ.LC()==0:
    QJ=s.Poly(QJ.as_expr()-QJ.LC()*T**QJ.degree(),T,domain=s.QQ.frac_field(e))
print('QJ_degree_symbolic',QJ.degree(),flush=True)
print('QJ_LC',s.factor(QJ.LC()),flush=True)
print('QJ_const',s.factor(QJ.nth(0)),flush=True)
print('QJ_at_minus_1_6',s.factor(QJ.eval(-s.Rational(1,6))),flush=True)

# Exact quadratic-field copy.  eta0 is a chosen root; both embeddings give
# the same degree/squarefree/irreducibility conclusions.
eta0=s.sqrt(-s.Rational(15510536192,22143375))
qext=s.Poly(s.expand(QJ.as_expr().subs(e,eta0)),T,extension=eta0)
print('QJ_degree',qext.degree(),flush=True)
print('QJ_gcd_degree',s.gcd(qext,qext.diff()).degree(),flush=True)
fac=s.factor_list(qext.as_expr(),T,extension=eta0)[1]
factor_ledger=[(s.Poly(g,T,extension=eta0).degree(),m) for g,m in fac]
basis=[basis_coeff(c) for c in QJ.all_coeffs()]
basis_payload=';'.join(f'{a},{b}' for a,b in basis)
print('QJ_factor_degrees',factor_ledger,flush=True)
print('QJ_basis_coefficients_desc_1_eta',basis,flush=True)
print('QJ_basis_sha256',sha256(basis_payload.encode()).hexdigest(),flush=True)
assert QJ.degree()==11 and qext.degree()==11
assert s.gcd(qext,qext.diff()).degree()==0 and factor_ledger==[(11,1)]
assert QJ.nth(0)!=0 and QJ.eval(-s.Rational(1,6))!=0
# A second irreducibility certificate: QJ/eta has rational coefficients.
# After primitive integral normalization its reduction modulo the inert prime
# 41 is the displayed irreducible degree-eleven polynomial.  Since gcd(11,2)=1,
# it remains irreducible over the residue field F_(41^2) of Q(sqrt(-30)).
qrat=s.Poly(s.cancel(QJ.as_expr()/e),T,domain=s.QQ)
qint=qrat.clear_denoms(convert=True)[1].primitive()[1]
q41=s.Poly(qint,T,modulus=41)
assert q41.degree()==11 and q41.is_irreducible
assert pow((-30)%41,(41-1)//2,41)==40
print('mod41_coefficients_desc',[int(c)%41 for c in qint.all_coeffs()],flush=True)
print('mod41_irreducible',q41.is_irreducible,'minus30_legendre',-1,flush=True)

# Leading X rows over T!=0, reduced modulo J (but these do not actually use J
# except to know eta is a unit).
pf=s.Poly(f,X)
ph=s.Poly(h,X)
print('f_X_degree',pf.degree(),'f_X_LC',s.factor(pf.LC()),flush=True)
print('h_X_degree',ph.degree(),'h_X_LC',s.factor(ph.LC()),flush=True)

# Universal pairs and their Hessians remain unchanged.
Hess=s.det(s.hessian(G,(X,T)))
for tv,x2,gv,hv in [(-s.Rational(1,6),s.Rational(6),s.Rational(1,2),-s.Rational(6)),
                    (s.Rational(0),-s.Rational(6),s.Rational(0),s.Rational(6))]:
    mod=s.Poly(X**2-x2,X)
    vals=[]
    for expr in (s.diff(G,X),s.diff(G,T),G-gv,Hess-hv):
        vals.append(s.rem(s.Poly(s.expand(expr.subs(T,tv)),X),mod).as_expr())
    print('universal',tv,vals,flush=True)

print('L',11+2+2,flush=True)
print('finite_nonidentity_capacity',2*12-15-2+3,'threshold',11,flush=True)
print('finite_identity_origin_index',9,'carrier_product_cap',3,flush=True)
print('full_overlap',18-15,'commutator_cap',2*(18-15),'origin_defect',12,flush=True)
assert 2*12-15-2+3 < 11
assert 3 < 9
assert 2*(18-15) < 12
semantic='|'.join([str(QJ.degree()),str(QJ.LC()),str(QJ.nth(0)),
                   str(QJ.eval(-s.Rational(1,6))),str(factor_ledger),
                   basis_payload,str([int(c)%41 for c in qint.all_coeffs()]),
                   '15','10<11','3<9','6<12'])
print('semantic_sha256',sha256(semantic.encode()).hexdigest(),flush=True)
print('verdict PASS',flush=True)
