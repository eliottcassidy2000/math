"""Exact observer controls for nextforest_20260926_moments.md.
No imports from another experiment; explicit gates survive python -O.
"""
from fractions import Fraction
from math import comb
from collections import Counter

COUNTS = Counter()
def check(ok, label):
    COUNTS[label] += 1
    if not ok:
        raise RuntimeError(label)

def T(n):
    return (3*n+1)//2 if n & 1 else n//2

def positions(n):
    out=[]
    while n:
        bit=n & -n
        out.append(bit.bit_length()-1)
        n-=bit
    return out

def digits(n):
    return {s:1 for s in positions(n)}

def carries(n, initial=1):
    out={}; c=initial; s=0
    while n or c:
        c=(3*(n&1)+c)//2
        s+=1
        if c: out[s]=c
        n >>= 1
    return out

def jet(poly, r):
    return tuple(sum(a*comb(s,k) for s,a in poly.items() if s>=k) for k in range(r))

def evaluate(poly,z):
    return sum((a*z**s for s,a in poly.items()),Fraction())

def valuation(n,p):
    if n==0: raise ValueError('valuation at zero')
    k=0
    while n%p==0: k+=1;n//=p
    return k

def polys(n):
    return (n,n+1,n*n+1)

def newton_poly(n):
    ps=positions(n); S=len(ps)
    powers=[sum(s**k for s in ps) for k in range(S+1)]
    elementary=[1]
    for k in range(1,S+1):
        num=sum((-1)**(i-1)*elementary[k-i]*powers[i] for i in range(1,k+1))
        check(num%k==0,'Newton integrality')
        elementary.append(num//k)
    coeff=[(-1)**k*elementary[k] for k in range(S+1)]
    return coeff

def horner(coeff,x):
    v=0
    for a in coeff:v=v*x+a
    return v

def rational_decode(x,p,q):
    # Positive rational p/q, p>=2, gcd(p,q)=1; input is an actual digit value.
    out=0;s=0
    while x:
        b=(x.numerator*pow(x.denominator,-1,p))%p
        check(b in (0,1),'rational digit admissibility')
        out |= b<<s
        x=Fraction(q,p)*(x-b)
        s+=1
        if s>10000:raise RuntimeError('decoder termination')
    return out

print('FINITE-EXACT: observer tests; no Collatz termination assertion')
# Exhaustive small universe: every positive n below 1024, including dense and sparse controls.
for n in range(1,1024):
    coeff=newton_poly(n)
    ps=positions(n)
    check(all(horner(coeff,s)==0 for s in ps),'Newton source roots')
    check(len(coeff)-1==len(ps),'Newton root count')
    M=sum(ps)
    check(n.bit_length()-1<=M or n==1,'first moment proper')
    for p,q in [(3,2),(2,3),(5,2),(2,1),(3,1)]:
        z=Fraction(p,q); val=evaluate(digits(n),z)
        check(val.denominator==q**(n.bit_length()-1),'rational exact denominator')
        check(rational_decode(val,p,q)==n,'rational decode')
# Deliberately hostile noninjective algebraic observer: z^2=z+1.
check(digits(4)=={2:1} and digits(3)=={0:1,1:1},'golden ratio collision support')

largest_bits=0; rows=[]
for horizon in [0,3,8,16]:
    base=[27]
    for _ in range(horizon):base.append(T(base[-1]))
    primes=(2,3,5)
    exps={p:1+max(valuation(v,p) for a in base for v in polys(a)) for p in primes}
    Q=2**(horizon+exps[2])*3**exps[3]*5**exps[5]
    while Q<=max(base):Q*=2
    Qs=[Q]
    for a in base[:-1]:Qs.append(Qs[-1]*(3 if a&1 else 1)//2)
    minimum=max((3*a+1).bit_length() for a in base)
    minimum=max(minimum,max((3*q).bit_length() for q in Qs))+2
    for r in range(1,6):
      even=[j for j in range(2**r) if j.bit_count()%2==0]
      odd=[j for j in range(2**r) if j.bit_count()%2==1]
      for gap in [minimum,minimum+7,minimum+19]:
        sums=[sum(1<<(gap*(j+1)) for j in group) for group in [even,odd]]
        pair=[27+Q*s for s in sums]
        largest_bits=max(largest_bits,max(pair).bit_length())
        check(max(pair)>2**(gap-1)*min(pair),'Prouhet unbounded ratio bound')
        for t in range(horizon+1):
            check(pair==[base[t]+Qs[t]*s for s in sums],'actual orbit affine blocks')
            check(jet(digits(pair[0]),r)==jet(digits(pair[1]),r),'paired digit jets')
            check(jet(carries(pair[0]),r)==jet(carries(pair[1]),r),'paired carry jets')
            # Independent block assembly for all carry coefficients.
            for ix,group in enumerate([even,odd]):
                block=dict(carries(base[t]))
                for j in group:
                    shift=gap*(j+1)
                    for s,a in carries(Qs[t],0).items():
                        check(s+shift not in block,'carry block separation')
                        block[s+shift]=a
                check(block==carries(pair[ix]),'literal carry block formula')
                for p in primes:
                    check(tuple(valuation(v,p) for v in polys(pair[ix]))==tuple(valuation(v,p) for v in polys(base[t])), 'orbit valuation fibre')
            if t<horizon:pair=[T(n) for n in pair]
        if r in [1,3,5] and gap==minimum+19:
            rows.append((horizon,r,gap,largest_bits))
print('Prouhet horizon/r/gap controls:',rows)
print('largest tested source bit length:',largest_bits)

# Same S at an actual growing U edge; first-moment gain diverges.
for R in range(1,81):
    n=51+sum(150<<(10*j) for j in range(1,R+1))
    m=77+sum(225<<(10*j) for j in range(1,R+1))
    check(T(n)==m and m&1,'growing actual odd edge')
    check(n.bit_count()==m.bit_count()==4+4*R,'same digit-count fibre')
    check(sum(positions(m))-sum(positions(n))==1+4*R,'unbounded first moment gain')
    check(m>n,'growing magnitude')
    if R==80:
        for d in range(1,7):
            w=lambda s:s**d-100*(s**(d-1))
            delta=sum(w(s) for s in positions(m))-sum(w(s) for s in positions(n))
            check(delta>0,'signed polynomial leading term control')
# Hasse versus ordinary moments are intentionally different coordinates.
check(jet(digits(51),2)==(4,10) and jet(digits(77),2)==(4,11),'minimal low block moment witness')
print('block family R=80: digit count324; first moment gain321')
# Full integer coefficient expansions, independently built from source bits
# and the valuation-normalized actual target, rather than a sampled X grid.
def root_polynomial(roots):
    coeff=[1]
    for root in roots:
        nxt=[0]*(len(coeff)+1)
        for k,c in enumerate(coeff):
            nxt[k]-=root*c
            nxt[k+1]+=c
        coeff=nxt
    return coeff
for n in range(1,1024,2):
    raw=3*n+1
    a=valuation(raw,2); m=raw>>a
    cs=carries(n)
    carry_roots=[s for s,c in cs.items() for _ in range(c)]
    left=3*positions(n)+[0]+carry_roots
    right=[s+a for s in positions(m)]+[s-1 for s in carry_roots for _ in range(2)]
    check(root_polynomial(left)==root_polynomial(right),'nonlinear support full coefficients')
    check(len(left)==len(right),'nonlinear support degrees')
for n in range(2,1024,2):
    check(root_polynomial(positions(n//2))==root_polynomial([s-1 for s in positions(n)]),'even support translation')
print('support transport:512 odd and511 even exact polynomial identities')
print('gates:',dict(sorted(COUNTS.items())))
print('TOTAL',sum(COUNTS.values()),'PASS')
