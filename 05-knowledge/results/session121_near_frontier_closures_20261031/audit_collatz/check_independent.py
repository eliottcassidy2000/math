#!/usr/bin/env python3
"""Independent exact-arithmetic audit; no producer imports or floor sums."""
from fractions import Fraction as F
import json

def require(test, message):
    if not test:
        raise RuntimeError(message)

eta = F(183105426, 48828125)
D12, D14, q = 8*5**11, 8*5**13, 5**14
T11, T12 = 2*5**10, 2*5**11
lam = lambda e: F(5,8)*(1-F(1,5**e))
d = lam(14)-lam(12)
require(d == F(3,5**13), 'slope increment')
require(F(21,8) < eta-1 < F(22,8), 'minimum eighth-score is 22')
cap = lambda j: D12*(F(j,8)-eta+1)
require([cap(j) for j in (22,23,24)] == [342,48828467,97656592], 'score caps')
forms = [('bare_high', F(1), 0, 0, 101725616),
         ('T11', F(1), F(31,125), T11, 814158),
         ('T12', F(1), F(6,25), T12, 4069366),
         ('bare_low', F(4,5), 0, 0, 20345408)]
for name, target, padding_fee, padding_clock, N in forms:
    bound = (target-eta+3-padding_fee)/d
    require(N < bound <= N+1, name+' strict integer cutoff')
    if name == 'bare_high':
        require(F(25-16,8)-F(N,D14) >= target, 'bare high excludes j>=25')
    else:
        shift = 16 if padding_clock == 0 else 14
        require(F(23-shift,8)-F(N+padding_clock,D14) >= target,
                name+' excludes j>=23')
# All-one H_e at its own level: exact valuation sum for n<=T_e.
# v5((3/2)^n-1)=1+v5(n/2) for even n and zero for odd n.
for e,T,expected in [(11,T11,F(31,125)),(12,T12,F(6,25))]:
    S = sum(T//(2*5**(k-1)) for k in range(1,e+1))
    require(F(5*T,8)-S == F(1,4), 'all-one eighth-score')
    require(lam(14)*T-S == expected, 'all-one price')
# Native endpoint tangent -v=(1,1), so multiplier 1+5^13.
a11 = pow(3*pow(2,-1,q)%q, T11, q)
a12 = pow(3*pow(2,-1,q)%q, T12, q)
targets = {'bare_high':1, 'bare_low':1+5**13,
           'T11':pow(a11,-1,q), 'T12':pow(a12,-1,q)}
# For identity target, inverse (a,a-1) is (a^-1,a^-1-1).
def compose(f,g): # chronological: f first, g second
    m,t=f; n,u=g
    return n*m%q, (n*t+u)%q
for a in (a11,a12):
    b=pow(a,-1,q)
    require(compose((a,a-1),(b,b-1)) == (1,0), 'prefix endpoint')
    require(compose((b,b-1),(a,a-1)) == (1,0), 'suffix endpoint')
p3=[pow(3,L,q) for L in range(343)]
hits={name:[] for name in targets}
pairs=0
for A in range(1,343):
    rhs={name:t*pow(2,A,q)%q for name,t in targets.items()}
    for L in range(1,A):
        pairs+=1
        for name,value in rhs.items():
            if p3[L] == value:
                hits[name].append([A,L])
require(all(not v for v in hits.values()), 'small-clock multiplier obstruction')
# Stronger universal route: no small non-one11 core can have identity multiplier.
q11=5**11
native11_hits=[]
for A in range(1,343):
    rhs=pow(2,A,q11)
    for L in range(1,A):
        if pow(3,L,q11)==rhs:
            native11_hits.append([A,L])
require(not native11_hits, 'universal primitive11 small-clock obstruction')
r=F(24,25)*F(23,8)+(eta-1)/25-2
require(r==F(1062011676,1220703125), 'uniform native13 margin')
require(r>F(4,5), 'low-tangent closure')
require(r+F(31,125)>1 and r+F(6,25)>1, 'uniform padding exclusion')
require(2*(1+r)-1>1, 'two-core exclusion')
require(F(24,25)*F(25,8)+(eta-1)/25-2>1, 'score>=25 exclusion without grammar')
# (2+r)C-B-1-r = (1+r)(C-1)+(C-B), for C>=B>=1.
# Both coefficients are nonnegative; this is an all-count algebraic proof.
require(1+r>0, 'all-count monotonicity')
print(json.dumps({'status':'PASS','conditional_on':'pinned E10 native11 fee and native-cut identities; no E10 replay',
                  'uniform_native13_margin':str(r),'native11_small_clock_hits':native11_hits,'native11_small_clock_pairs':pairs,'D12':D12,'D14':D14,'q':q,'score_clock_caps':{str(j):int(cap(j)) for j in (22,23,24)},
                  'multipliers':targets,'literal_clock_length_pairs':pairs,'hits':hits,
                  'conclusion':'low tangent closed; four padded high-tangent forms excluded; bare high only j=23,24'},
                 sort_keys=True,indent=2))
