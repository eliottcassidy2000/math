#!/usr/bin/env python3
"""Finite-exact regression checks; all-index claims are proved in the note."""
from fractions import Fraction as F
import json
from pathlib import Path

# Exact arithmetic a+b sqrt(-319).
def mul(z,w):
    a,b=z; c,d=w
    return a*c-319*b*d,a*d+b*c

t=(F(9,20),F(1,20))
assert mul(t,(t[0],-t[1]))==(1,0)
N=150
S=[F(2),F(9,10)]
A=[2,9]
for m in range(1,N):
    S.append(F(9,10)*S[-1]-S[-2])
    A.append(9*A[-1]-100*A[-2])
cs=[13-6*s for s in S]
z=(F(1),F(0))
Uprev,Ucur=F(0),F(1) # U_{-1}, U_0 at x=9/20
for m in range(1,N+1):
    z=mul(z,t)
    assert 2*z[0]==S[m]
    assert A[m]==10**m*S[m]
    assert 1<cs[m]<25
    assert cs[m].denominator==10**m//2
    assert (cs[m]-1)*(cs[m]-25)==F(-9*319,25)*Ucur**2
    assert Ucur!=0
    Uprev,Ucur=Ucur,F(9,10)*Ucur-Uprev
    if m>=2:
        assert cs[m]==F(9,10)*cs[m-1]-cs[m-2]+F(143,10)
assert len(set(cs[1:]))==N
P2=lambda c:25-(13-c)**2/F(6)
for m in range(1,N//2+1):
    assert P2(cs[m])==cs[2*m]
assert P2(F(13))==25 # finite-order failure of nonzero square-class preservation
assert cs[:3]==[F(1),F(38,5),F(1007,50)]
receipt={
  'status':'PASS_FINITE_EXACT', 'max_exponent':N,
  'first_charges':[str(c) for c in cs[:8]],
  'checks':['quadratic field norm and trace','integer recurrence','strict open interval',
            'reduced denominator 10^m/2','pairwise distinct finite values',
            'charge affine recurrence','nonzero discriminant U identity',
            'P2 maps exponent m to 2m','hostile finite-order c=13'],
  'proof_only_not_computationally_certified':['infinite multiplicative order',
       'angle/pi transcendence by Gelfond-Schneider','all-integer arcsine limiting law',
       'all-index denominator and square-class statements'],
  'not_claimed':['specific doubling orbit density or equidistribution',
                 'intermediate-E7.5 VOA realization','positive real mirror interpretation'],
}
path=Path(__file__).with_name('checks.json')
path.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
