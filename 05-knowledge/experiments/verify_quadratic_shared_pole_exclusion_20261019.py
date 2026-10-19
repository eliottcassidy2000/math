#!/usr/bin/env python3
"""Standard-library-only exact verification of the elementary shared-pole proof.

All arithmetic is sparse multivariate polynomial arithmetic over Fraction.
No solver, computer-algebra package, modular arithmetic, or removable assertion
is used. The verifier follows the displayed identities in the companion note.
"""
from fractions import Fraction
from collections import Counter
import json

NVAR = 6
ZERO_MON = (0,)*NVAR
checks = Counter()

def check(condition, name):
    checks[name] += 1
    if not condition:
        raise RuntimeError(name)

class Poly:
    def __init__(self, terms=None):
        if isinstance(terms, Poly):
            self.terms = terms.terms.copy()
        elif isinstance(terms, dict):
            self.terms = {m:Fraction(c) for m,c in terms.items() if c}
        elif terms:
            self.terms = {ZERO_MON:Fraction(terms)}
        else:
            self.terms = {}
    def __add__(self, other):
        out=self.terms.copy()
        for m,c in Poly(other).terms.items():
            out[m]=out.get(m,Fraction(0))+c
            if not out[m]:del out[m]
        return Poly(out)
    __radd__=__add__
    def __neg__(self):
        return Poly({m:-c for m,c in self.terms.items()})
    def __sub__(self, other):return self+-Poly(other)
    def __rsub__(self, other):return Poly(other)+-self
    def __mul__(self, other):
        out={}
        for ma,ca in self.terms.items():
            for mb,cb in Poly(other).terms.items():
                m=tuple(a+b for a,b in zip(ma,mb))
                out[m]=out.get(m,Fraction(0))+ca*cb
        return Poly(out)
    __rmul__=__mul__
    def __pow__(self, n):
        if not isinstance(n,int) or n<0:raise ValueError('Nonnegative integer powers only')
        ans=Poly(1)
        for _ in range(n):ans=ans*self
        return ans
    def __eq__(self, other):return self.terms==Poly(other).terms
    def coefficient(self, variable, degree):
        out={}
        for mon,c in self.terms.items():
            if mon[variable]==degree:
                m=list(mon);m[variable]=0;out[tuple(m)]=c
        return Poly(out)
    def specialize(self, replacements):
        out=Poly()
        for mon,c in self.terms.items():
            term=Poly(c)
            for i,n in enumerate(mon):
                term*=Poly(replacements.get(i,variables[i]))**n
            out+=term
        return out

variables=[]
for i in range(NVAR):
    mon=[0]*NVAR;mon[i]=1;variables.append(Poly({tuple(mon):1}))
t,a,b,S,d,eps=variables
D=t*t+1
T=t*(t*t-1)
q=(t-a)**2+b*b
X=t*t-a*t-b
Y=-(b+1)*t+a
check(X*X+Y*Y==D*q,'Gaussian_norm_identity')
E=S*(q*q-X*X+Y*Y)-2*d*X*Y-4*eps*T
F3=a*S-(b+1)*d+2*eps
F2=(5*a*a+3*b*b+4*b+1)*S-2*a*(b+2)*d
F1=(2*a**3+2*a*b*b+2*a*b+a)*S+(-a*a+b*b+b)*d-2*eps
F0=(a**4+2*a*a*b*b+a*a+b**4-b*b)*S+2*a*b*d
for n,wanted in [(3,-2*F3),(2,F2),(1,-2*F1),(0,F0)]:
    check(E.coefficient(0,n)==wanted,'coefficient_identity')
check(all(mon[0]<=3 for mon in E.terms),'leading_term_cancellation')
M=[
    [5*a*a+3*b*b+4*b+1,-2*a*(b+2)],
    [2*a*(a*a+b*b+b+1),-a*a+b*b-1],
    [a**4+2*a*a*b*b+a*a+b**4-b*b,2*a*b],
]
check(F3+F1==M[1][0]*S+M[1][1]*d,'homogeneous_middle_row')
H=a*a+(b+1)**2
expected=[
    H*((4*b+3)*a*a+3*b*b-2*b-1),
    2*a*H*((b+2)*a*a+b**3+b),
    H*(a**4+(2*b+1)*a*a-b*b*(b-1)**2),
]
for (i,j),wanted in zip([(0,1),(0,2),(1,2)],expected):
    check(M[i][0]*M[j][1]-M[j][0]*M[i][1]==wanted,'two_by_two_determinant')
check(3*b*b-2*b-1==(3*b+1)*(b-1),'a_zero_factorization')
check(F2.specialize({1:0,2:1})==8*S,'a_zero_scale_contradiction')
check(q.specialize({1:0,2:-1})==D,'H_zero_excluded_pair')
R=2*b**4+4*b+1
check((4*b+3)*(b**3+b)-(b+2)*(3*b*b-2*b-1)==2*R,'rational_root_eliminant')
for root in [Fraction(1),Fraction(-1),Fraction(1,2),Fraction(-1,2)]:
    check(R.specialize({2:root})!=0,'rational_root_candidate_exclusion')
print(json.dumps({'status':'PASS: elementary shared-pole exclusion verified with exact standard-library Fraction arithmetic',
                  'checks':dict(checks)},sort_keys=True))
