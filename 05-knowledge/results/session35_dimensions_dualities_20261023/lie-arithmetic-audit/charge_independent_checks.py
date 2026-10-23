#!/usr/bin/env python3
from fractions import Fraction as Q
from math import comb, gcd
from pathlib import Path
import json
# Independent binomial expansion of (9 + sqrt(-319))^m / 20^m.
traces=[]; charges=[]
for m in range(81):
    real=sum(comb(m,2*j)*9**(m-2*j)*(-319)**j for j in range(m//2+1))
    imag=sum(comb(m,2*j+1)*9**(m-2*j-1)*(-319)**j for j in range((m-1)//2+1))
    a,b=Q(real,20**m),Q(imag,20**m)
    assert a*a+319*b*b==1
    s=2*a; c=13-12*a
    traces.append(s); charges.append(c)
    if m:
        assert 1<c<25 and c.denominator==10**m//2
        assert (c-1)*(c-25)==-319*(12*b)**2
        assert b!=0
        A=s*10**m
        assert A.denominator==1 and gcd(A.numerator,10)==1
assert charges[1]==Q(38,5)
assert len(set(charges))==81
for m in range(1,80):
    assert charges[m+1]==Q(9,10)*charges[m]-charges[m-1]+Q(143,10)
for m in range(1,41):
    assert 25-(13-charges[m])**2/6==charges[2*m]
# Laurent-polynomial constant coefficients certify moment and lag formulas.
def multiply(p,q):
    out={}
    for i,a in p.items():
        for j,b in q.items(): out[i+j]=out.get(i+j,0)+a*b
    return out
p={0:1}
for j in range(21):
    expected=0 if j%2 else 36**(j//2)*comb(j,j//2)
    assert p.get(0,0)==expected
    p=multiply(p,{-1:-6,1:-6})
# Symbolic lag multiplier v: constant in z of 36(z+z^-1)(zv+z^-1v^-1).
terms={}
for z1 in (-1,1):
    for z2 in (-1,1):
        key=(z1+z2,z2)
        terms[key]=terms.get(key,0)+36
assert {v:a for (z,v),a in terms.items() if z==0}=={-1:36,1:36}
# Each lag contribution is 36v and 36/v, hence 36S_r; r=0 variance72.
assert 36*traces[0]==72 and (36*traces[1])/72==Q(9,20)
receipt={'status':'PASS','independent_method':'binomial quadratic-field expansion, not producer recurrence','max_exponent':80,'moments_through_degree':20,'first_charges':[str(c) for c in charges[:6]],'exact_checks':['norm','denominator including m=1','square class via explicit imaginary coordinate','strict interval','distinctness','charge recurrence','doubling identity','Laurent constant moments','lag-one correlation'],'analytic_proof_review':['infinite order','Gelfond-Schneider angle ratio','consecutive-index geometric sums','cosine pushforward arcsine density','fixed-lag covariance'],'not_certified_by_finite_tests':['equidistribution','transcendence','all-index denominator']}
Path(__file__).with_name('charge_independent_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
