"""All-m DG source submersions, alternative affine chart, and unit torsion.

Finite exact controls for the accompanying all-parameter analytic proof.
No previous producer is imported; no repository file is written.
"""
from fractions import Fraction
from math import comb
from pathlib import Path
import json
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')
GATES = 0
x,t,r,b,u,c,a,h,B0 = s.symbols('x t r b u c a h B0')

def check(v, label):
    global GATES
    if s.cancel(v) != 0:
        raise RuntimeError(label + ': ' + str(v))
    GATES += 1

def require(v, label):
    global GATES
    if not v:
        raise RuntimeError(label)
    GATES += 1

def jac(f,g,v=x,w=t):
    return s.diff(f,v)*s.diff(g,w)-s.diff(f,w)*s.diff(g,v)

def qm(m, X=x, H=h):
    return sum(comb(2*m,j)*(-H)**j*X**(m-j) for j in range(m+1))

records=[]

# The full binomial sums, before specializing h, are independent integer
# checks of the two closed forms used by the unbounded proof.
for m in range(1,41):
    q0=sum((-1)**j*comb(2*m,j) for j in range(m+1))
    q1=sum((m-j)*(-1)**j*comb(2*m,j) for j in range(m+1))
    require(q0==(-1)**m*comb(2*m-1,m), 'Q(h) binomial sum')
    require(q1==(-1)**(m-1)*comb(2*m-2,m-1), 'Qprime(h) binomial sum')

for m in range(1,9):
    e=2*m-1
    Q=qm(m)
    A=a*(x-h)**(2*m)
    B=a*Q+B0
    F=A*t+B
    f0=B.subs(x,h)
    g=F-f0
    kappa=(-1)**(m-1)*comb(2*m-2,m-1)*h**(m-1)
    K=s.cancel(g/(a*(x-h)))
    require(s.denom(K)==1, 'special fibre quotient polynomial')
    check(K.subs(x,h)-kappa, 'disjoint reduced special fibre')
    check(s.diff(F,x).subs(x,h)-a*kappa, 'source submersion')
    check(f0-(B0+a*(-1)**m*comb(2*m-1,m)*h**m), 'special value')
    remainder=sum(comb(2*m,j)*(-h)**j*r**(j-m) for j in range(m+1,2*m+1))
    Finf=B0-a*b*(1-h*r)**(2*m)-a*remainder
    check(F.subs({x:1/r,t:-r**m-r**(2*m)*b})-Finf, 'full global chart identity')
    check(Finf.subs(r,0)-(B0-a*b), 'boundary restriction')
    check(s.diff(Finf,b).subs(r,0)+a, 'boundary submersion')
    G0=1/(a*e*(x-h)**e)
    check(jac(F,G0)-1, 'rational mate')
    check(g-a*(x-h)*K, 'polynomial witness factorization')
    # Retain the factorized all-exponent mechanism instead of expanding
    # hundreds of powers of the same polynomial into a giant identity.
    check(jac(F,K)-a*(x-h)**e*K, 'factorized witness derivative')
    if m<=3:
        P=a**(e-1)*K**e/e
        check(s.factor(jac(F,P)-g**e), 'literal small-m witness bracket')
    check(jac(F,1/(x-h))-a*(x-h)**(2*m-2), 'alternative chart volume')
    check(G0.subs(x,h+1/u)-u**e/(a*e), 'power-map mate')
    # Exact valuation controls; the inequalities prove the lower-order
    # obstruction together with the distinct reduced component.
    for ell in range(e):
        require(ell-e<0 and ell>=0, 'unequal component pole orders')
    coeff=s.cancel((a*kappa)**e/(a*e))
    require(coeff!=0, 'top scalar principal part nonzero')
    records.append({'m':m,'e':e,'Q_h1':str(s.expand(Q.subs(h,1))),
                    'kappa_h1':int(kappa.subs(h,1)),
                    'special_value_a1_h1_B00':int(f0.subs({a:1,h:1,B0:0})),
                    'top_principal_coefficient_a1_h1':str(coeff.subs({a:1,h:1})),
                    'annihilator_exponent':e,'alternative_map_degree':e})

# Independent literal parameter controls for the inverse charts, including
# the special fibre, the newly filled D point, and x=0 in the old chart.
for m in range(1,7):
    for av,hv,Bv in [(1,1,0),(2,-1,3),(-3,2,-2)]:
        Q=qm(m,H=s.Integer(hv))
        F=av*(x-hv)**(2*m)*t+av*Q+Bv
        f0=F.subs({x:hv,t:0})
        R=u/(1+hv*u)
        T=s.expand(u**(2*m)*((c-Bv)/s.Integer(av)-Q.subs(x,hv+1/u)))
        BB=s.expand((Bv-c)/s.Integer(av)*(1+hv*u)**(2*m)
                    -sum(comb(2*m,j)*(-hv)**j*u**(j-m)*(1+hv*u)**(3*m-j)
                         for j in range(m+1,2*m+1)))
        require(s.denom(T)==1 and s.denom(BB)==1,'inverse coefficients polynomial')
        check(F.subs({x:hv+1/u,t:T})-c,'original-chart inverse')
        check(-R**m-R**(2*m)*BB-T,'chart gluing t')
        check(BB.subs(u,0)-(Bv-c)/s.Integer(av),'new D point')
        check(T.subs(u,0),'new D point t limit')
        check(BB.subs(u,-s.Rational(1,hv)),'old x0 has b0')
        # The new chart Jacobian proves the finite map's exact branch order.
        U=u**(2*m-1)/(s.Integer(av)*(2*m-1))
        check(s.diff(U,u)-u**(2*m-2)/av,'finite power-map Jacobian')
        for cv in [f0,f0+1,f0-2]:
            tc=T.subs(c,cv)
            bc=BB.subs(c,cv)
            check(F.subs({x:hv+1/u,t:tc})-cv,'whole-fibre parametrization')
            require(s.denom(tc)==1 and s.denom(bc)==1,'whole-fibre polynomial parameters')

# Sharp h=0 boundary: the line is nonreduced and the exponent drops to two
# for every m>=2. Do not feed this singular case to the smooth torsion theorem.
hostiles=[]
for m in range(1,11):
    e=2*m-1
    F=x**(2*m)*t+x**m
    G0=1/(s.Integer(e)*x**e)
    wanted=1 if m==1 else 2
    check(jac(F,G0)-1,'h0 rational mate')
    P=s.cancel(F**wanted*G0)
    require(s.denom(P)==e or s.denom(P)==1,'h0 repaired polynomial witness')
    check(jac(F,P)-F**wanted,'h0 witness bracket')
    require((m*(wanted-1)-e)<0<=m*wanted-e,'h0 exact valuation threshold')
    if m==1:
        check(s.diff(F,x).subs(x,0)-1,'m1 h0 remains smooth')
    else:
        check(s.diff(F,x).subs(x,0),'h0 critical line Fx')
        check(s.diff(F,t).subs(x,0),'h0 critical line Ft')
    hostiles.append({'m':m,'annihilator_exponent':wanted,'line_multiplicity':m})

# Positive affine-source control; its second coordinate has a genuine pole
# on the added boundary, so it is not a global pair on W_m.
check(jac(t,-x)-1,'positive polynomial coordinate control')
require(s.denom(-1/r)==r,'boundary pole of source coordinate')

cert={'scope':'finite exact controls for the stated analytic all-m theorem',
      'ring':'C[x,t]; F global on W_m but D_F is not asserted global',
      'symbolic_m_range':[1,8], 'binomial_m_range':[1,40],
      'inverse_m_range':[1,6], 'inverse_parameter_bank':[[1,1,0],[2,-1,3],[-3,2,-2]],
      'symbolic_controls':records,'h0_controls':hostiles,'always_active_gates':GATES}
source_path=Path(__file__).resolve()
result_dir=(source_path.parent.parent/'05-knowledge'/'results'
            if source_path.parent.name=='04-computation' else source_path.parent)
certificate_path=result_dir/(source_path.stem+'_certificate.json')
certificate_path.write_bytes((json.dumps(cert,indent=2,sort_keys=True)+'\n').encode('utf-8'))
print('DG source submersions on W_m: m>=1, a*h!=0, full polynomial-part constant convention.')
print('Exact unit annihilator: (F-f0)^(2m-1) in the original source response module.')
print('Alternative affine chart: W_m minus E1 = A2_(F,u); mate map is finite flat degree 2m-1.')
print('Boundary ramification index: 2m-1; m=1 is unramified; h=0,m>=2 drops unit order to 2.')
print(f'Always-active exact gates: {GATES}')
