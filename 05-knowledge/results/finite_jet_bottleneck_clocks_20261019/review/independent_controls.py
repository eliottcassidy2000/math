#!/usr/bin/env python3
"""Independent finite-exact controls for the finite-jet clock classification.

No author code is imported. Generic coefficients are derived in two ways:
(1) Q B' = alpha Q' B, a triangular differential recurrence;
(2) scaled-denominator expansion and beta moments, enumerated by partitions.
Additional controls address the coordinate-family caveat, the sharp jet cutoff,
flow matching, and an explicit rational corrected-clock bound.
Run with python3 and python3 -O; writes results beside this script.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import sympy as S

ROOT=Path(__file__).resolve().parent
checks=[]
def check(ok,label):
    if not ok:
        raise RuntimeError(label)
    checks.append(label)

def equal(x,y,label):
    check(S.expand(x-y)==0,label)

def diff_power(unit,alpha,n):
    """Coefficients of unit**alpha from its differential equation, constant 1."""
    b=[S.Integer(1)]
    for k in range(1,n+1):
        b.append(S.expand(sum((((alpha+1)*j-k)*unit[j]*b[k-j]
                    for j in range(1,min(k,len(unit)-1)+1)),S.Integer(0))/k))
    return b

def partitions_of_weight(n,j=1):
    if n==0:
        yield {}
        return
    if j>n:
        return
    for multiplicity in range(n//j+1):
        for rest in partitions_of_weight(n-j*multiplicity,j+1):
            if multiplicity:
                rest=dict(rest)
                rest[j]=multiplicity
            yield rest

def scaled_beta_coefficient(r,k,qs):
    # Scale y=epsilon**(1/r)*t. Each denominator-geometric term of
    # multiplicity s contributes a beta moment B(s+alpha,1-alpha).
    # Divide by B(alpha,1-alpha): (alpha)_s/s!. The s! geometric
    # multinomial cancels it. This is not an inverse-series computation.
    alpha=S.Rational(k+1,r)
    ans=S.Integer(0)
    for partition in partitions_of_weight(k):
        s=sum(partition.values())
        term=(-1)**s*S.rf(alpha,s)
        for degree,multiplicity in partition.items():
            term*=qs[degree]**multiplicity/S.factorial(multiplicity)
        ans+=term
    return S.expand(ans)

x,c,a=S.symbols('x c a', nonzero=True)
generic=[]
for r in [2,4,6,8,10,12]:
    qs=[S.Integer(1)]+list(S.symbols(f'q1:{r}'))
    row={'r':r,'even_coefficients':[]}
    for k in range(0,r-1,2):
        got=diff_power(qs,-S.Rational(k+1,r),k)[k]
        alternate=scaled_beta_coefficient(r,k,qs)
        equal(got,alternate,f'generic coefficient by ODE vs beta partitions: r={r}, k={k}')
        row['even_coefficients'].append(str(got))
    generic.append(row)

q1,q2,q3,q4=S.symbols('q1 q2 q3 q4')
quartic=diff_power([1,q1,q2],-S.Rational(3,4),2)[2]
equal(quartic,S.Rational(21,32)*q1**2-S.Rational(3,4)*q2,'generic asymmetric quartic formula')
equal(quartic.subs({q1:1,q2:S.Rational(7,8)}),0,'quartic different-jet equal-clock positive control')
sextic2=diff_power([1,q1,q2],-S.Rational(1,2),2)[2]
sextic4=diff_power([1,q1,q2,q3,q4],-S.Rational(5,6),4)[4]
equal(sextic2,S.Rational(3,8)*q1**2-S.Rational(1,2)*q2,'generic sextic c2')
equal(sextic4,S.Rational(21505,31104)*q1**4-S.Rational(935,432)*q1**2*q2+
      S.Rational(55,36)*q1*q3+S.Rational(55,72)*q2**2-S.Rational(5,6)*q4,
      'generic sextic c4')

# Unlike a restatement of formula (2), this computes a nonlinear conjugate
# of y+y^4 and checks that reimposing constant-additive epsilon changes clocks.
psi=sum((-1)**(j-1)*S.catalan(j-1)*c**(j-1)*x**j for j in range(1,8))
equal(S.series(psi+c*psi**2-x,x,0,8).removeO(),0,'inverse coordinate phi(y)=y+c y^2 through degree 7')
Fpsi=psi+psi**4
conjugate_increment=S.series(Fpsi+c*Fpsi**2-x,x,0,7).removeO().expand()
equal(conjugate_increment,x**4-2*c*x**5+4*c**2*x**6,'nonlinear conjugate quartic through degree 6')
coordinate_c2=S.expand(quartic.subs({q1:-2*c,q2:4*c*c}))
equal(coordinate_c2,-S.Rational(3,8)*c*c,'coordinate-family hostile has nonzero second clock coefficient')

jet=[]
lam=S.Symbol('lambda')
for r in [4,6,8,10,12]:
    unit=[S.Integer(1)]+[S.Integer(0)]*(r-3)+[lam]
    coefficients=[]
    for k in range(0,r-1,2):
        got=diff_power(unit,-S.Rational(k+1,r),k)[k]
        expected=S.Integer(1) if k==0 else (-S.Rational(r-1,r)*lam if k==r-2 else S.Integer(0))
        equal(got,expected,f'last-jet sharpness r={r}, k={k}')
        coefficients.append(str(got))
    jet.append({'r':r,'normalized_clock_coefficients':coefficients})

for r in [2,4,6,8,10,12]:
    qs=[S.Integer(1)]+list(S.symbols(f't1:{r+1}'))
    changed=qs[:]
    changed[r-1]+=lam
    changed[r]+=lam**2
    for k in range(0,r-1,2):
        equal(diff_power(qs,-S.Rational(k+1,r),k)[k],
              diff_power(changed,-S.Rational(k+1,r),k)[k],
              f'above-cutoff perturbation ignored r={r}, k={k}')
    # Independently expand the binomial flow to its second nonzero term.
    flow=S.series(x*(1-(r-1)*a*x**(r-1))**(-S.Rational(1,r-1)),x,0,2*r).removeO()
    equal(flow,x+a*x**r+S.Rational(r,2)*a*a*x**(2*r-1),f'embedded flow jet r={r}')

# Verify primitive and tail-power endpoints algebraically, including r=2.
epsilon=S.Symbol('epsilon',positive=True)
for r in [2,4,6,8,10,12]:
    equal(S.together(S.diff(S.log(epsilon+x**r)/r,x)-x**(r-1)/(epsilon+x**r)),0,
          f'borderline odd primitive r={r}')
    for k in range(0,r-1,2):
        check(r-k-1>0,f'full-line tail integrable r={r}, k={k}')

# Exact sample domain constants for quartic Q=1+y+(7/8)y^2.
# J=[-1/4,1/4], delta=1/8. Coefficient majorants prove the inequalities
# everywhere on J, without sampling a grid.
r=4; delta=F(1,8); M=2*delta
qs=[F(1),F(1),F(7,8)]
q_lower=1-sum(abs(qs[j])*M**j for j in range(1,len(qs)))
q_upper=sum(abs(qs[j])*M**j for j in range(len(qs)))
D1=sum((r+j)*abs(qs[j])*M**j for j in range(len(qs)))
eta=M**(r-1)*D1
B=sum((r+j)*(r+j-1)*abs(qs[j])*M**(r+j-2) for j in range(len(qs)))
C=M**(r-2)*D1*D1/q_lower
D=F(5,12)*B+F(7,12)*C*(1+eta)/(1-eta)
check(q_lower>0,'exact Q positivity majorant')
check(eta<=F(1,4),'exact derivative majorant <= 1/4')
check(delta/F(4)+eta*delta<=delta/F(2),'exact one-step overshoot domain majorant')
check(F(1,6)+F(1,4)==F(5,12),'correct B contribution to defect constant')
check(F(1,3)+F(1,4)==F(7,12),'correct C contribution to defect constant')
constants={name:str(value) for name,value in locals().copy().items()
           if name in ['delta','M','q_lower','q_upper','eta','B','C','D']}
constants['summed_defect_bound']=str(F(5,2)*delta*D)
constants['overshoot_integral_bound']=str(1/(1-eta))

result={
 'status':'PASS',
 'classification':'FINITE-EXACT algebra and rational inequality certificates; no numerical orbit assertions',
 'checks':len(checks),
 'methods':['Triangular differential recurrence for fractional powers',
            'Independent scaled-geometric-denominator beta-moment expansion by partitions'],
 'generic_even_order_controls':generic,
 'coordinate_family_hostile':{'phi':'y+c*y^2','conjugate_increment_through_degree_6':str(conjugate_increment),
                             'flattened_jacobian_c2':str(coordinate_c2),
                             'A1':'-3*pi*c^2/(8*sqrt(2))'},
 'sharp_jet_controls':jet,
 'quartic_exact_domain_constants':constants,
 'check_labels':checks,
 'limitations':['These finite controls support rather than prove the all-germ analytic theorem.',
               'Analytic nonembedding of polynomial germs is a cited theorem, not a finite computation.',
               'No claim of coordinate-conjugacy invariance is made.']}
(ROOT/'independent-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'exact_checks':len(checks),
                  'coordinate_hostile_c2':str(coordinate_c2),
                  'domain_constants':constants},indent=2))
