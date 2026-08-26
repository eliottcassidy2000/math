#!/usr/bin/env python3
"""Independent exact point-level audit on the repeated-wall H30 stratum.

Unlike the producer certificate, the decisive check evaluates f, h, their
Jacobian, and Hess(G) directly at the repeated-resultant point in the exact
field Q[u,eta]/(h15(u),eta^2-u).  A Rabin certificate modulo 347 proves the
degree-thirty tower polynomial irreducible.
"""

import ast
from hashlib import sha256
from pathlib import Path
import sympy as S
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
FACTOR_OUT = ROOT / '05-knowledge/results/jc23_weight9_repeated_top_wall_factors_thm4157.out'
FACTOR_SHA = '91b91809c99ebeecf26af27fdb7168311fa5124a90f71a8f942c1495133cf483'
X,T,e,u = S.symbols('X T eta u')


def parity_to_u(expr, parity):
    p=S.Poly(expr,e,domain=QQ.frac_field(T))
    assert all(k[0] % 2 == parity for k,_ in p.terms())
    return S.expand(sum(c*u**((k[0]-parity)//2) for k,c in p.terms()))


def powmod_poly(base, exponent, modulus):
    ans=S.Poly(1,base.gens[0],modulus=base.get_modulus())
    while exponent:
        if exponent & 1:
            ans=(ans*base).rem(modulus)
        base=(base*base).rem(modulus)
        exponent >>= 1
    return ans


raw=FACTOR_OUT.read_bytes()
assert sha256(raw).hexdigest()==FACTOR_SHA
lines=raw.decode().splitlines()
i=lines.index('factor=H30;degree_eta=30;degree_u=15')
coeffs=ast.literal_eval(lines[i+1].split('=',1)[1])
h15_integer=S.Poly.from_list(coeffs,u,domain=S.ZZ).primitive()[1]
h15=h15_integer.monic()
H30_integer=S.Poly(h15_integer.as_expr().subs(u,e**2),e,domain=S.ZZ)

# Fully explicit finite-field irreducibility certificate for H30=h15(eta^2).
prime=347
Hmod=S.Poly(H30_integer,e,modulus=prime)
xmod=S.Poly(e,e,modulus=prime)
rabin=[]
for divisor in (2,3,5):
    test=(powmod_poly(xmod,prime**(30//divisor),Hmod)-xmod).rem(Hmod)
    rabin.append(S.gcd(Hmod,test).degree())
frobenius=(powmod_poly(xmod,prime**30,Hmod)-xmod).rem(Hmod)
assert Hmod.degree()==30 and Hmod.LC()!=0
assert rabin==[0,0,0] and frobenius.is_zero

# Direct row-D source; no producer import.
P=T+X**2*T**2
Y=X*T*P
D=S.Rational(2048,45)
K=S.Rational(1376,135)
G=S.expand(-X**2*T/2-3*P+S.Rational(8,3)*P**2-S.Rational(1376,135)*P**3
           +K*Y**2+D*P**4-D*P*Y**2+e*P**3*Y-e*Y**3)
f=S.cancel(S.diff(G,X)/T)
h=S.diff(G,T)
pf=S.Poly(f,X); ph=S.Poly(h,X)
assert pf.degree()==6 and ph.degree()==7
assert S.factor(pf.LC()-7*T**7*e)==0
assert S.factor(ph.LC()-8*T**7*e)==0
detD=S.expand(S.diff(f,X)*S.diff(h,T)-S.diff(f,T)*S.diff(h,X))
hess=S.expand(S.diff(G,X,2)*S.diff(G,T,2)-S.diff(G,X,T)**2)
assert S.cancel(T*detD-hess-f*S.diff(G,X,T))==0

# Primary residual and primitive linear X-subresultant.
R=S.resultant(f,h,X)
Q12=S.Poly(S.cancel(R/(T**35*(6*T+1)**2)),T)
Qu=parity_to_u(Q12.as_expr(),1)
subs=S.subresultants(f,h,X)
lin=S.Poly(subs[-2],X)
content=S.gcd(lin.nth(1),lin.nth(0))
assert S.factor(content-T**35*e*(6*T+1)/6)==0
linp=S.Poly(S.cancel(lin.as_expr()/content),X)
Au=parity_to_u(linp.nth(1),0)
Bu=parity_to_u(linp.nth(0),1)

# Exact degree-15 field K0=Q[u]/h15.
K0=QQ.alg_field_from_poly(h15,alias='u0')
u0=K0.unit

def kelem(expr):
    pp=S.Poly(expr,u,domain=QQ)
    ans=K0.zero
    for c in pp.all_coeffs():
        ans=ans*u0+K0.convert(c)
    return ans

def tpoly(expr):
    pp=S.Poly(expr,T)
    return S.Poly.from_list([kelem(c) for c in pp.all_coeffs()],T,domain=K0)

Qf=tpoly(Qu)
cf=S.gcd(Qf,Qf.diff()).monic()
assert cf.degree()==1
assert S.div(Qf,cf**2)[1].is_zero
assert S.div(Qf,cf**3)[1] != 0
c_anp=cf.rep.to_list()[1]
t0=-c_anp
assert t0!=K0.zero and 6*t0+K0.one!=K0.zero
J0=22143375*u0+K0.convert(15510536192)
assert u0!=K0.zero and J0!=K0.zero

Af=tpoly(Au); Bf=tpoly(Bu)
def eval_anp_poly(poly,value):
    answer=K0.zero
    for coefficient in poly.rep.to_list():
        answer=answer*value+coefficient
    return answer

Aval=eval_anp_poly(Af,t0); Bval=eval_anp_poly(Bf,t0)
assert Aval!=K0.zero
r=-Bval/Aval                 # X=eta*r in the quadratic tower.
assert eval_anp_poly(Qf,t0)==K0.zero
assert eval_anp_poly(cf,t0)==K0.zero

# Evaluate a polynomial at T=t0, X=eta*r, eta^2=u0.  Return its coefficients
# in the basis (1,eta) over K0.
def eval_tower(expr):
    pp=S.Poly(expr,X,T,e,domain=QQ)
    even=K0.zero; odd=K0.zero
    for (ix,it,ie),c in pp.terms():
        power=ix+ie
        term=K0.convert(c)*(r**ix)*(t0**it)*(u0**(power//2))
        if power%2: odd += term
        else: even += term
    return even,odd

point_values={
    'linear':eval_tower(linp.as_expr()),
    'f':eval_tower(f),
    'h':eval_tower(h),
    'detD':eval_tower(detD),
    'hessian':eval_tower(hess),
}
assert all(a==K0.zero and b==K0.zero for a,b in point_values.values())

def anp_digest(value):
    return sha256(','.join(str(q) for q in value.to_list()).encode()).hexdigest()

semantic='\n'.join((
    'p=347;rabin='+str(rabin)+';frobenius=0',
    'gcd_degree=1;multiplicity=2',
    't0='+anp_digest(t0),
    'x_over_eta='+anp_digest(r),
    'Aunit='+anp_digest(Aval),
    'f=h=detD=hessian=0',
))+'\n'
print('INDEPENDENT H30 POINT-LEVEL AUDIT')
print('H30_mod_347_degree=30;Rabin_gcd_degrees='+str(tuple(rabin))+';Frobenius_remainder=0')
print('H30_irreducible_over_Q=yes;h15_irreducible=yes;eta^2-u_irreducible=yes')
print('gcd(Q,Qprime)_degree=1;repeated_multiplicity=2')
print('t0_unit=yes;6t0+1_unit=yes;u_unit=yes;J_unit=yes;A(t0)_unit=yes')
print('X_leading_rows=7*T^7*eta,8*T^7*eta;finite_infinity_loss=0')
print('direct_point_values=subresultant:0,f:0,h:0,detD:0,Hessian:0')
print('t0_sha256='+anp_digest(t0))
print('X_over_eta_sha256='+anp_digest(r))
print('A_at_t0_sha256='+anp_digest(Aval))
print('semantic_sha256='+sha256(semantic.encode()).hexdigest())
print('verdict=PASS')
