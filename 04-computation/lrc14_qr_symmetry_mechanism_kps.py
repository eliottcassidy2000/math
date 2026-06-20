#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WHY does sum_c D7(c) vanish on every #QR-pattern class?  Find the exact symmetry.

Candidate mechanisms (Euler even-factor cancellation <-> quadratic char):
  (S1) Galois: D7(a*c) = sigma_a( D7(c) ) for a in F_7^* acting as zeta->zeta^a.
       Summing over the F_7^* orbit of c => sum_a sigma_a(D7(c)) = TRACE -> rational.
  (S2) Multiplicativity / dilation: scaling all coords by a in F_7^*.
  (S3) coordinate-permutation antisymmetry (the (-1)^|T| alternation).
  (S4) D7 vanishes UNLESS the coords are "balanced" -> support of the sum.

We test:
  (A) the Galois relation D7(a c) = a-Frobenius of D7(c) EXACTLY in Z[zeta_7];
  (B) the per-coordinate structure: does D7 factor as prod_j d(c_j) * (alternating)?
  (C) sum over a full F_7^* dilation orbit {a*c : a in F_7^*}; is it 0 or trace-rational?
  (D) decompose D7 in the basis {1=trace part} vs {Gauss-sum part}: write D7(c) =
      A(c) + B(c)*g with g the Gauss sum; A rational, B rational? (i.e. D7 in Q(g)?)
"""
import sys, itertools
from fractions import Fraction
import cmath, math
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# exact Z[zeta_7]
ZERO=(Fraction(0),)*6
def zexp(e):
    e%=7
    if e<=5:
        v=[Fraction(0)]*6; v[e]=Fraction(1); return tuple(v)
    return (Fraction(-1),)*6
def zadd(a,b): return tuple(x+y for x,y in zip(a,b))
def zsub(a,b): return tuple(x-y for x,y in zip(a,b))
def zscale(a,s): return tuple(x*s for x in a)
def zmul(a,b):
    out=[Fraction(0)]*7
    for i in range(6):
        if a[i]==0: continue
        for j in range(6):
            if b[j]==0: continue
            out[(i+j)%7]+=a[i]*b[j]
    c6=out[6]
    return tuple(out[i]-c6 for i in range(6))
Z=cmath.exp(2j*math.pi/7)
def znum(a): return sum(a[k]*(Z**k) for k in range(6))
def galois(a, x):
    """apply zeta -> zeta^a to x in Z[zeta_7]."""
    r=ZERO
    for k in range(6):
        if x[k]!=0:
            r=zadd(r, zscale(zexp((a*k)%7), x[k]))
    return r

Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:Fraction((-1)**len(T)) for T in Tlist}
def sigma_T(T,m):
    r=ZERO
    for t in T: r=zadd(r, zexp((-m*t)%7))
    return r
SIG={(T,m):sigma_T(T,m) for T in Tlist for m in range(1,7)}
PREF={m:zsub(zexp(0), zexp((-m)%7)) for m in range(1,7)}
def D7(c):
    pref=zexp(0)
    for cj in c: pref=zmul(pref,PREF[cj])
    acc=ZERO
    for T in Tlist:
        p=zexp(0)
        for cj in c: p=zmul(p, SIG[(T,cj)])
        acc=zadd(acc, zscale(p, SGN[T]))
    return zmul(pref,acc)

QR={1,2,4}; NQR={3,5,6}
def chi(m):
    m%=7
    return 0 if m==0 else (1 if m in QR else -1)
def gsum():
    g=ZERO
    for r in range(1,7): g=zadd(g, zscale(zexp(r), Fraction(chi(r))))
    return g

if __name__=="__main__":
    g=gsum()
    print(f"g num={znum(g):.5f}")

    # (A) Galois relation: D7(a*c) =? galois_a(D7(c))
    print("\n(A) Galois dilation test  D7(a c) ?= sigma_a(D7(c)):")
    import random
    ok=True
    for _ in range(8):
        c=tuple(random.choice([1,2,3,4,5,6]) for _ in range(6))
        for a in range(1,7):
            lhs=D7(tuple((a*cj)%7 for cj in c))
            rhs=galois(a, D7(c))
            if lhs!=rhs: ok=False; print(f"    FAIL c={c} a={a}"); break
    print(f"   D7(a c) = sigma_a(D7(c)) for all tested: {ok}")
    print("   => D7 is GALOIS-EQUIVARIANT under F_7^* dilation of coords.")

    # (C) sum over dilation orbit {a c} = sum_a sigma_a(D7(c)) = Trace_{Q(z7)/Q}(D7(c))  (RATIONAL!)
    print("\n(C) Dilation-orbit sum = field trace (rational):")
    for c in [(1,1,1,1,1,1),(1,2,4,3,5,6),(1,1,1,1,1,2),(1,2,3,1,2,3)]:
        orb=ZERO
        for a in range(1,7):
            orb=zadd(orb, D7(tuple((a*cj)%7 for cj in c)))
        # trace should be rational => coeffs in basis are all equal? trace of zeta^k = -1 (k!=0), trace(1)=6
        # element r is rational iff r = q*(1) with the basis being... actually rational element = (q,0,0,0,0,0)? No:
        # 1 has coeffs (1,0,0,0,0,0). A rational q is (q,0,0,0,0,0). check:
        israt = all(orb[k]==0 for k in range(1,6))
        print(f"   c={c}: orbit-sum coeffs={[str(x) for x in orb]} rational={israt} num={znum(orb):.4f}")

    # (D) Is D7(c) in Q(g) = Q(sqrt(-7))?  i.e. D7 = A + B g with A,B rational.
    #     Q(g) is the unique quadratic subfield. Element x is in Q(g) iff fixed by the
    #     index-2 subgroup of Gal = <2> (since <2>={1,2,4}=QR fixes Q(sqrt-7)).
    print("\n(D) Is D7(c) in the quadratic field Q(sqrt(-7)) = fixed field of QR-Galois <2>?")
    cntin=0; tot=0
    import random
    sample=[tuple(random.choice([1,2,3,4,5,6]) for _ in range(6)) for _ in range(400)]
    sample += [(1,1,1,1,1,1),(1,2,4,3,5,6),(3,3,3,3,3,3),(1,2,3,4,5,6)]
    for c in sample:
        tot+=1
        d=D7(c)
        if galois(2,d)==d: cntin+=1
    print(f"   #cosets (sample) with D7(c) in Q(sqrt-7) [fixed by sigma_2]: {cntin} / {tot}")
    # also: is D7 fixed by the WHOLE Galois (i.e. rational)?  check sigma_3 (a non-residue generator)
    cntrat=0
    for c in sample:
        d=D7(c)
        if galois(3,d)==d: cntrat+=1
    print(f"   #cosets (sample) fixed by sigma_3 (NQR) too [=> rational]: {cntrat} / {tot}")
