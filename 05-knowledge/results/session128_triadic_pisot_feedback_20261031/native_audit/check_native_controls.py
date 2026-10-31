#!/usr/bin/env python3
"""Exact finite controls for the native feedback audit.
All gates use rational/integer arithmetic and survive python -O.
Threshold examples are ABSTRACT controls, not certified Amanda realizations.
"""
from fractions import Fraction as F
from itertools import product
import json

checks=0
def require(ok,label):
    global checks
    checks+=1
    if not ok: raise RuntimeError(label)
def sgn(x): return (x>0)-(x<0)
def seed(h,b): return 1 if h>=b else -1
def cell(h,b): return sum(x<b for x in h)
def update(b,lam,c): return (1-lam)*b+lam*c

# Exact algebraic three-shift certificate: for frequencies 1 and 2 in x=3 theta,
# sum_{j=0}^2 z^{j*k} vanishes modulo z^2+z+1.
# Polynomial remainders evaluated in basis (1,z), where z^2=-z-1.
def zpow(n):
    a,b=1,0
    for _ in range(n): a,b=-b,a-b
    return a,b
for k in (1,2,-1,-2):
    terms=[zpow((j*k)%3) for j in range(3)]
    require(tuple(sum(t[i] for t in terms) for i in (0,1))==(0,0),'three-shift Fourier cancellation')
for a,b in product(range(-50,51),repeat=2):
    triple=(a,b,-a-b)
    require(abs(sum(map(sgn,triple)))<=1,'three-signed-zero-sum bound')
require(sum(1 if v>=0 else -1 for v in (0,0,0))==3,'zero-field hostile to using sgn_plus')

# Exact rational unit-circle sample points, with positive and negative gain controls.
gains=[F(1,1000),F(1,3),F(1),F(7,2),F(10**6)]
sample_count=0
for coeff in product(range(-1,2),repeat=4):
    for t in [F(n,d) for n in range(-8,9) for d in (1,2,3)]:
        x=(1-t*t)/(1+t*t); y=2*t/(1+t*t)
        vals=(y,x,2*x*y,x*x-y*y)
        g=sum(c*v for c,v in zip(coeff,vals))
        for q in gains: require(sgn(q*g)==sgn(g),'positive common gain invariance')
        require(sgn(-g)==-sgn(g),'parity oddness pointwise')
        if g: require(sgn(0*g)!=sgn(g),'zero gain outside invariance scope')
        sample_count+=1

# Bounds for arbitrary allowed parity sequences are stronger than native-only tests.
for lam in (F(1,9),F(1,3),F(1,2),F(2,3),F(1)):
    for cs in product((F(-1,3),F(0),F(1,3)),repeat=6):
        b=F(0)
        for n,c in enumerate(cs,1):
            b=update(b,lam,c)
            envelope=(1-(1-lam)**n)/3
            require(abs(b)<=envelope,'finite-time exact memory envelope')
            require(seed(F(1,3),b)==1,'positive boundary freeze')
            if lam<1: require(seed(F(-1,3),b)==-1,'negative finite-time boundary freeze lambda<1')
require(seed(F(-1,3),F(-1,3))==1,'negative boundary tie lambda=1')

# Threshold tie handling, including duplicate thresholds: cell index counts EVENTS
# strictly below b, not events below or equal; realizable patterns group duplicates.
for hs in ((F(-1,4),F(0),F(1,4)),(F(0),F(0),F(1,4))):
    for b in (F(-1,3),F(-1,4),F(0),F(1,8),F(1,4),F(1,3)):
        pattern=tuple(seed(h,b) for h in hs)
        i=cell(hs,b)
        require(pattern==(-1,)*i+(1,)*(len(hs)-i),'right-closed cell pattern')
    points=[min(hs)-1]+list(hs)+[max(hs)+1]
    patterns=set(tuple(seed(h,b) for h in hs) for b in points)
    require(len(patterns)==len(set(hs))+1,'distinct-threshold cell count')

# Synthetic fixed-map controls: fixed-point iff c belongs to its own cell.
hs=(F(-1,6),F(1,6)); possible=(F(-1,3),F(-1,6),F(0),F(1,6),F(1,3))
for cs in product(possible,repeat=3):
    for lam in (F(1,3),F(1,2),F(1)):
        for i,c in enumerate(cs):
            actual=update(c,lam,cs[cell(hs,c)])==c
            # Coincident c values can make the same fixed point appear in more than
            # one candidate list, so per-cell membership is the exact criterion.
            if cell(hs,c)==i: require(actual,'membership sufficient for fixed point')
        for b in possible:
            i=cell(hs,b)
            require((update(b,lam,cs[i])==b)==(b==cs[i]),'fixed-point equation')
    # lambda=1: after one map application, every orbit lies in at most M+1 images.
    for b0 in possible:
        b=b0; seen={}; seq=[]
        for step in range(5):
            b=cs[cell(hs,b)]
            require(b in cs,'lambda-one finite image')
            if b in seen:
                require(step-seen[b]<=3,'lambda-one period bound')
                break
            seen[b]=step;seq.append(b)
        else: raise RuntimeError('lambda-one failed to repeat')

# An explicit abstract discontinuity control. It obeys the parity range and the
# all-plus/all-minus oddness c_2=-c_0, but native realizability is NOT asserted.
# b0=0 -> 1/8 -> 3/16; thereafter b_n=1/8+2^(-n-2), tending
# to the excluded left endpoint of the middle cell, which is not a fixed point.
hs=(F(1,8),F(1,4)); cs=(F(1,4),F(1,8),F(-1,4));lam=F(1,2)
b=F(0); orbit=[b]
for n in range(1,21):
    b=update(b,lam,cs[cell(hs,b)]);orbit.append(b)
    if n>=2: require(b==F(1,8)+F(1,2**(n+2)),'excluded-endpoint asymptotic control')
require(update(F(1,8),lam,cs[cell(hs,F(1,8))])!=F(1,8),'limit can be nonfixed at discontinuity')

print(json.dumps({'status':'FINITE-EXACT controls passed','checks':checks,
 'rational_carrier_samples':sample_count,
 'scope':'Proof audit controls; synthetic threshold examples do not certify native realizability.',
 'boundary_limit_control':{'thresholds':list(map(str,hs)),'cell_parities':list(map(str,cs)),
 'lambda':str(lam),'first_states':list(map(str,orbit[:7])),'limit':'1/8','limit_is_fixed':False}},indent=2,sort_keys=True))
