#!/usr/bin/env python3
"""Exact regression controls for the all-parameter proof in README.md."""
from fractions import Fraction as F
from itertools import product
import json
import sympy as s

checks=[]
def require(condition, label):
    if not condition:
        raise RuntimeError(label)
    checks.append(label)

def zero(matrix):
    return all(s.cancel(v)==0 for v in matrix)

w,u,x,g=s.symbols('w u x g')
for n in (1,2,-1,-2):
    # Reduce negative exponents modulo the third-root relation first.
    v=1+w**(n%3)+w**((2*n)%3)
    require(s.rem(v,w*w+w+1,w)==0, 'cyclic_shift_'+str(n))
require(s.expand((1-u)*(2*u+1)-(u-(2*u*u-1)))==0,'sharp_sign_factorization')
q, sigma, gain=s.symbols('q sigma gain')
Ftheta=q*sigma*(s.cos(3*x)/18-s.cos(6*x)/144)
require(s.simplify(s.diff(gain*Ftheta,x)-gain*q*sigma*(-s.sin(3*x)/6+s.sin(6*x)/24))==0,'native_gain_factorization')
a,lam,c,b=s.symbols('a lam c b')
require(s.expand((1-lam)*(c+a*(b-c))+lam*c-(c+a*a*(b-c))).subs(lam,1-a).expand()==0,'branch_induction')
for p in range(1,9):
    cs=s.symbols('c0:'+str(p))
    b0=(1-a)*sum(a**(p-1-k)*cs[k] for k in range(p))/(1-a**p)
    bp=b0
    for ck in cs:
        bp=a*bp+(1-a)*ck
    require(s.cancel(bp-b0)==0,'periodic_formula_p'+str(p))

# Half-open ownership and exact finite universe: 27 abstract templates,
# four rational gains, every word of length 1..4 on three threshold cells.
thresholds=(F(-1,6),F(1,6))
csvals=(F(-1,3),F(0),F(1,3))
lams=(F(1,4),F(1,2),F(3,4),F(1))
def cell(v):
    return sum(h<v for h in thresholds)

def orbit(word,cs,la):
    aa=1-la
    start=la*sum(aa**(len(word)-1-k)*cs[word[k]] for k in range(len(word)))/(1-aa**len(word))
    cur=start
    states=[]
    for i in word:
        states.append(cur)
        cur=aa*cur+la*cs[i]
    if cur!=start:
        raise RuntimeError('periodic closure failed')
    return states
words_checked=valid_words=branch_fixed=0
for cs in product(csvals, repeat=3):
    for la in lams:
        for i,ci in enumerate(cs):
            if cell(ci)==i:
                branch_fixed+=1
                # Endpoints/interiors include both left/right hostile controls.
                for start in (F(-1,3),F(-1,6),F(-1,12),F(0),F(1,12),F(1,6),F(1,3)):
                    if cell(start)!=i:
                        continue
                    cur=start
                    for n in range(1,8):
                        cur=(1-la)*cur+la*ci
                        if cell(cur)!=i or cur!=ci+(1-la)**n*(start-ci):
                            raise RuntimeError('trapping branch failure')
        for p in range(1,5):
            for word in product(range(3), repeat=p):
                states=orbit(word,cs,la)
                valid=all(cell(v)==i for v,i in zip(states,word))
                if valid:
                    valid_words+=1
                    cur=states[0]
                    for i in range(2*p):
                        cur=(1-la)*cur+la*cs[cell(cur)]
                        if cur!=states[(i+1)%p]:
                            raise RuntimeError('realized word failure')
                words_checked+=1
require(words_checked==12960,'declared_word_universe')

# Generic piecewise-affine hostiles; not asserted native-realizable.
for la in lams:
    aa=1-la
    amp=la/(3*(2-la))
    require(aa*(-amp)+la/F(3)==amp and aa*amp-la/F(3)==-amp,'abstract_two_cycle_'+str(la))
require(F(1,2)*F(0)+F(1,2)*F(1,4)==F(1,8),'excluded_endpoint_not_fixed')

# Full symbolic finite-rank identity. No small-gain or truncation.
L=s.diag(-1,-4,-9,-36,-49,-64)
R=L.inv()
U=s.Matrix([[1,0],[0,1],[1,1],[2,-1],[-1,2],[1,3]])
C=s.Matrix([[2,1],[-3,4]])
V=s.Matrix([[1,2,0,-1,0,1],[0,1,1,0,2,-1]])
S=s.eye(2)-g*V*R*U*C
W=R+g*R*U*C*S.inv()*V*R
require(zero((L-g*U*C*V)*W-s.eye(6)),'symbolic_woodbury_full_inverse')
require(s.cancel((L-g*U*C*V).det()/L.det()-S.det())==0,'symbolic_secular_determinant')

# Nonsingular stationary system with divergent naive iteration.
LE=s.diag(-9,-9,-36,-36)
require((LE-18*s.eye(4)).det()!=0 and F(18,-9)==-2,'static_invertible_iteration_unstable')
resonant=LE+9*s.eye(4)
f_bad=s.Matrix([1,0,0,0]); f_good=s.Matrix([0,0,1,0])
require(resonant.row_join(f_bad).rank()>resonant.rank(),'resonant_incompatible_source')
require(resonant.row_join(f_good).rank()==resonant.rank() and len(resonant.nullspace())==2,'resonant_compatible_nonunique')
L2=s.diag(-9,-1); B2=s.Matrix([[0,0],[1,0]])
phi=s.Matrix([-s.Rational(1,9),g/9])
require(zero((L2-g*B2)*phi-s.Matrix([1,0])),'exact_mode_leakage')

print(json.dumps({'status':'PASS','arithmetic':'exact rational and symbolic',
 'symbolic_checks':len(checks),'abstract_templates':27,'gains':4,
 'periodic_words_checked':words_checked,'admissible_words':valid_words,
 'branch_fixed_instances':branch_fixed,'secular_polynomial':str(s.factor(S.det())),
 'native_parity_bound':'1/3','native_fixed_memory_interval':['-1/3','1/3'],
 'native_sharpness_scope':'complete two-harmonic carrier, not prescribed event realization'},sort_keys=True,indent=2))
