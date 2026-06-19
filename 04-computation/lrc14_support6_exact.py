#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXACT check of the THM-538 support-6 floor, using cyclotomic (Q(zeta_7)) algebra
so there is NO float error.  Represent every Fourier coefficient as an exact
element of the group ring Z[1/7][zeta_7] = Fraction-combination of zeta_7^a.

chat_T(n) for n != 0, n not mult 7:
   chat_T(n) = -(1/(2 pi i n)) sum_{j in T} (e^{-2pi i n(j+1)/7} - e^{-2pi i n j/7})
Let w = e^{-2pi i n /7} = zeta_7^{-n mod 7} (a primitive-ish 7th root).
Then e^{-2pi i n j/7} = w^j.  And -1/(2 pi i n) sum_j (w^{j+1}-w^j)
   = -1/(2 pi i n) (w-1) sum_{j in T} w^j.
The factor 1/(2 pi i n) is NOT in Q(zeta_7); but K(n) is built from PRODUCTS over
coords and a SUM over T, and we only need to test whether the *coefficient*
structure forces K(n)=0 for low support.  The cleanest exact route:

K(n) = integral_0^1...0^1 ??? -- instead use the DIRECT definition:
K(n) = sum_{T} (-1)^|T| prod_j chat_T(n_j).
Each chat_T(n_j) = (1 - |T|/7) if n_j=0, else  alpha(n_j) * S_T(n_j)
where alpha(n_j) = (w_j - 1)/(-2 pi i n_j) is a coord-only scalar and
S_T(m) = sum_{j in T} w_{m}^j  with w_m = zeta_7^{(-m) mod 7}.

The pi/n scalars factor OUT of the T-sum per coordinate, so
K(n) = (prod_j alpha(n_j)) * sum_T (-1)^|T| prod_j S_T(n_j)   [for all coords nonzero]
The T-SUM  Q(n) := sum_{T subseteq {1..6}} (-1)^|T| prod_j S_T(n_j)
is an EXACT element of Z[zeta_7].  We compute Q(n) exactly.

CLAIM (support-6 floor): Q(n) = 0 whenever fewer than 6 coords are active.
If Q(n)=0 then K(n)=0 regardless of the (nonzero) pi/n prefactors. We test Q(n)
exactly in Z[zeta_7] for the same low-support cases. (For coords that are 0 or
mult-of-7, chat_T = (1-|T|/7) constant, handled separately.)
"""
import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# Element of Z[zeta_7] as a length-7 Fraction vector over basis zeta^0..zeta^6,
# reduced mod (1+z+...+z^6 = 0) by subtracting the all-ones part is optional;
# we keep degree<7 and compare via the relation that zeta^7=1 and
# 1+z+..+z^6=0.  To get a canonical zero-test, reduce using zeta^6 = -(1+..+z^5).
def reduce6(v):
    # v: list length>=7 of Fractions for powers; fold zeta^7=1 then zeta^6 relation
    # first fold exponents mod 7
    w = [Fraction(0)]*7
    for e,c in enumerate(v):
        w[e % 7] += c
    # apply zeta^6 = -(1+z+z^2+z^3+z^4+z^5)
    c6 = w[6]
    w[6] = Fraction(0)
    for e in range(6):
        w[e] -= c6
    return w[:6]

def zmul(a, b):
    # multiply two length-7 (or 6) poly reps in Z[zeta_7]
    a = a + [Fraction(0)]*(7-len(a)) if len(a)<7 else a
    b = b + [Fraction(0)]*(7-len(b)) if len(b)<7 else b
    out = [Fraction(0)]*13
    for i in range(7):
        if a[i]==0: continue
        for j in range(7):
            if b[j]==0: continue
            out[i+j]+=a[i]*b[j]
    return reduce6(out)

def zadd(a,b):
    a = list(a)+[Fraction(0)]*(6-len(a))
    b = list(b)+[Fraction(0)]*(6-len(b))
    return [a[i]+b[i] for i in range(6)]

def zzero():
    return [Fraction(0)]*6

def zfromexp(e):
    # zeta^e as reduced element
    v=[Fraction(0)]*7
    v[e%7]=Fraction(1)
    return reduce6(v)

def S_T(T, m):
    # sum_{j in T} zeta_7^{((-m) mod 7) * j}
    base = (-m) % 7
    acc = zzero()
    for j in T:
        acc = zadd(acc, zfromexp((base*j) % 7))
    return acc

SUBS = [frozenset(c) for r in range(7) for c in itertools.combinations(range(1,7), r)]

def Q_exact(active_freqs):
    """active_freqs: list of nonzero, not-mult-7 frequencies (the active coords).
       Returns Q(n) = sum_T (-1)^|T| prod_{coord} S_T(freq) in Z[zeta_7]."""
    total = zzero()
    for T in SUBS:
        prod = [Fraction(1)] + [Fraction(0)]*5  # = 1
        for m in active_freqs:
            prod = zmul(prod, S_T(T, m))
        sign = (-1)**len(T)
        prod = [sign*c for c in prod]
        total = zadd(total, prod)
    return reduce6(total + [Fraction(0)])

def is_zero(v):
    return all(c==0 for c in v[:6])

if __name__ == "__main__":
    print("EXACT Q(n) in Z[zeta_7] for low-support active-coord patterns:")
    print("(support-6 floor claims Q=0 for <6 active coords)\n")
    cases = {
        "1 active":   [1],
        "2 active":   [1,-1],
        "3 active":   [1,1,-1],
        "4 active":   [1,1,1,-1],
        "5 active":   [1,1,1,1,-1],
        "6 active (=)": [1,1,1,1,1,-1],
    }
    for name, freqs in cases.items():
        q = Q_exact(freqs)
        z = is_zero(q)
        print(f"  {name:14s} active={len(freqs)}  Q==0 exactly: {z}")
        if not z:
            print(f"      Q = {[str(c) for c in q]}")
