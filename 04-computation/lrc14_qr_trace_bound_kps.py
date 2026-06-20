#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DOES the Galois-trace collapse give a BOUND on the correction sum K(n)?

ESTABLISHED EXACTLY (prior scripts):
  K(n) = arch(n) * D7(n mod 7),   arch(n) = prod_j 1/(2 pi i n_j)  (REAL up to i^6=-1),
  D7(a c) = sigma_a(D7(c))  [Galois-equivariant],
  sum_{a in F7*} D7(a c) = Tr(D7(c)) in Z  [rational, the cleanest collapse].

THE LATTICE QUESTION.  Group relations n by the dilation action on residues.
For a FIXED integer relation n with residue c = n mod 7, the six "dilates" are
NOT integer multiples of n (dilation is mod 7 on residues only).  But the relation
lattice Lambda(E) = {n : sum n_j e_j = 0} is a genuine sublattice of Z^6 that is
INVARIANT under n -> -n (so c -> -c = 6c, a=6=NQR).  Key facts we test:

 (T1) arch(n) is REAL (i^6 = -1 => arch(n) = -1/( (2pi)^6 prod n_j ), a real signed number).
      So K(n) = (real weight w(n)) * D7(c).  The IMAGINARY part of the correction must
      cancel coordinate-wise, and the cyclotomic D7(c) carries ALL the oscillation.

 (T2) PAIRING n <-> -n: arch(-n)=arch(n) (even # of sign flips: 6 coords => (-1)^6=1),
      and D7(-c) = D7(6c) = sigma_6(D7(c)) = conj(D7(c)) (since 6 = -1 mod 7, sigma_{-1}=conjugation).
      THEREFORE  K(n)+K(-n) = w(n)[D7(c)+conj D7(c)] = 2 w(n) Re D7(c)  -- REAL, as required,
      and the imaginary parts cancel EXACTLY pairwise.  This is the first telescoping.

 (T3) The QR/NQR magnitude isometry: |D7(c)| is invariant under c -> 3c (3 = NQR generator)?
      NO (Galois changes value).  But |D7| is invariant under c -> a c for ALL a (|sigma_a x|...
      actually |sigma_a(x)| != |x| in general). Test what |D7| symmetry actually holds:
      the scan showed mean|D7| symmetric under #QR <-> 6-#QR. Pin the exact involution.

 (T4) THE BOUND: Re D7(c) summed over the F7* orbit's REAL parts.  Since
      sum_a D7(ac) = Tr(D7(c)) (rational), and pairing gives Re, we get
      sum_{a} Re D7(ac) = Tr(D7(c))  (because Tr is real = its own Re, and the 6 terms
      come in conj pairs a<->6a).  So the orbit-real-sum is an INTEGER.  Bound |Tr|.
"""
import sys, itertools
from fractions import Fraction
import cmath, math
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

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
def galois(a,x):
    r=ZERO
    for k in range(6):
        if x[k]!=0: r=zadd(r, zscale(zexp((a*k)%7), x[k]))
    return r
def ztrace(x):
    """Tr_{Q(z7)/Q}: trace(1)=6, trace(z^k)=-1 (k=1..6)."""
    # x = x0 + sum_{k>=1} xk z^k ; Tr = 6 x0 + sum_{k=1..5} xk * (-1)  but careful basis only to z^5
    # Tr(z^0)=6, Tr(z^k)=-1 for k=1..6. In basis z^0..z^5: Tr = 6 x0 - (x1+..+x5).
    return Fraction(6)*x[0] - sum(x[1:], Fraction(0))

Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:Fraction((-1)**len(T)) for T in Tlist}
def _sig(T,m):
    r=ZERO
    for t in T: r=zadd(r, zexp((-m*t)%7))
    return r
SIG={(T,m):_sig(T,m) for T in Tlist for m in range(1,7)}
PREF={m:zsub(zexp(0),zexp((-m)%7)) for m in range(1,7)}
def D7(c):
    pref=zexp(0)
    for cj in c: pref=zmul(pref,PREF[cj])
    acc=ZERO
    for T in Tlist:
        p=zexp(0)
        for cj in c: p=zmul(p,SIG[(T,cj)])
        acc=zadd(acc, zscale(p,SGN[T]))
    return zmul(pref,acc)

QR={1,2,4}; NQR={3,5,6}
def chi(m): return 0 if m%7==0 else (1 if m%7 in QR else -1)

if __name__=="__main__":
    print("LRC(14) Galois-trace collapse -> bound on the correction (kps)")

    # T2: D7(-c) = conj(D7(c)); K(n)+K(-n) real.
    print("\n(T2) Pairing n<->-n:  D7(6c) ?= conj(D7(c)) and K-pair is real:")
    import random
    okpair=True
    for _ in range(6):
        c=tuple(random.choice([1,2,3,4,5,6]) for _ in range(6))
        d=D7(c); dneg=D7(tuple((6*cj)%7 for cj in c))
        # conj in Z[z7] = sigma_6
        if dneg!=galois(6,d): okpair=False
    print(f"   D7(6c)=sigma_6(D7(c))=conj for all tested: {okpair}")
    print("   => imaginary part of correction cancels EXACTLY pairwise n<->-n. Correction is REAL.")

    # T4: orbit-real-sum = Tr(D7(c)) is an INTEGER. tabulate Tr over residue classes,
    # grouped by QR/NQR pattern, and bound it.
    print("\n(T4) Field trace Tr(D7(c)) (integer) by #QR-coords pattern:")
    from collections import defaultdict
    byq=defaultdict(list)
    allT=[]
    # full scan of traces is feasible (rational arithmetic) but 46656 * D7 is slow;
    # sample heavily + cover structured points. Use orbit reps to cut by 6.
    seen=set(); reps=[]
    for c in itertools.product(range(1,7),repeat=6):
        # orbit under dilation a in F7*: canonical rep = min over a of (a c)
        key=min(tuple((a*cj)%7 for cj in c) for a in range(1,7))
        if key in seen: continue
        seen.add(key); reps.append(c)
    print(f"   #dilation-orbit reps = {len(reps)} (of 46656; ~/6)")
    for c in reps:
        tr=ztrace(D7(c))
        nq=sum(1 for x in c if x in QR)
        byq[nq].append(int(tr))
        allT.append(int(tr))
    print(f"   total #orbit-reps traced: {len(allT)}")
    print(f"   max |Tr| over all reps = {max(abs(t) for t in allT)}")
    print(f"   sum of Tr over all orbit reps = {sum(allT)}")
    for nq in range(7):
        L=byq[nq]
        if not L: continue
        print(f"     #QR={nq}: count={len(L)}  min={min(L)} max={max(L)} sum={sum(L)} mean={sum(L)/len(L):.2f}")

    # The cleanest statement: per orbit, sum_a Re D7(ac) = Tr(D7(c)). Since |Re D7| <= |D7|,
    # and the orbit has 6 members in 3 conj-pairs, Tr controls the *signed* real content.
    print("\n  INTERPRETATION: within each F7* dilation orbit, the SIGNED real content of D7")
    print("  collapses to the integer Tr(D7(c)). The cyclotomic oscillation is gone after")
    print("  averaging over the multiplicative group -- the Euler/QR even-factor cancellation.")
