#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_signed_cut_explicit_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

We established (inequality script): a VALID+TIGHT+CONSEC-MAX signed cut C=sum_A lambda_A a[A]
exists using PROPER atoms only (no constant a[empty]).  NOW we EXTRACT THE EXPLICIT cut and
answer the three task questions:
  (1) SPARSEST signed combination (min support);
  (2) does the SIGN pattern follow the residue/Mobius parity (-1)^|A| ?
  (3) is the cut UNIFORM in k (same Boolean-type combination, scaled)?

THE STRUCTURAL (manifestly-nonneg) FORM -- the cleanest non-circular certificate.
Because  measS7(E) = sum_A (-1)^|A| a[A](E)   (inclusion-exclusion, EXACT),
and every cumulative atom  a[A](E) >= 0  (it is a measure), ANY upper cut
   C(E) = sum_A lambda_A a[A](E)
satisfies  C(E) - measS7(E) = sum_A (lambda_A - (-1)^|A|) a[A](E).
=> If  beta_A := lambda_A - (-1)^|A| >= 0 for all A,  then C(E) >= measS7(E) STRUCTURALLY
   (no per-shape checking; this is the even-Bonferroni / Mobius-truncation principle).
So the genuine certificate is: choose a set S of atoms; set lambda_A = (-1)^|A| for A in S and
require the OMITTED atoms to be exactly the negative-parity ones (so truncation stays an UPPER
bound) -- i.e. C = even-truncated inclusion-exclusion.  Then consec-max needs (M): the truncated
sum is consec-maximal.  We test BOTH:
  (E1) the EVEN-BONFERRONI truncations B_R = sum_{|A|<=R, R even} (-1)^|A| a[A] (R=0,2,4,6):
       valid upper bounds by Bonferroni; is B_R consec-MAXIMAL over the stratum? at which R?
  (E2) the LP-optimal sparse cut from the inequality script, made exact, with its support & signs.

If an EVEN-BONFERRONI truncation B_R is consec-maximal at finite R uniformly in k, THAT is the
explicit signed cut: support = all A with |A|<=R, signs = (-1)^|A| (pure Mobius parity), uniform
in k (same combination for every k).  That CLOSES consec-max on the stratum modulo proving B_R
consec-max (a finite, structural statement).  Honest report of the smallest working R.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)
INNER = list(range(1, 7)); SUBMASKS = list(range(64))
def msize(m): return bin(m).count("1")
def mset(m): return tuple(s for s in INNER if (m >> (s-1)) & 1)

def exact_mask_atoms(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(7*e+1): bps.add(F(m, 7*e))
    bps = sorted(bps); q = defaultdict(F)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2; hit = {int(7*e*mid) % 7 for e in E}; mask = 0
        for s in range(1, 7):
            if s not in hit: mask |= 1 << (s-1)
        q[mask] += b - a
    return dict(q)
def cont(q): return {A: sum(v for M, v in q.items() if (M & A) == A) for A in SUBMASKS}
def prim(E): return reduce(gcd, [e for e in E if e], 0) == 1
def fullres(E): return len({e % 7 for e in E}) == 7
def measS7(q): return q.get(0, F(0))
def fmt(x): return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"

WIN = 14
def build(k):
    recs = []; consec = None
    for rest in itertools.combinations(range(1, WIN+1), k-1):
        E = (0,) + rest
        if not prim(E) or not fullres(E): continue
        q = exact_mask_atoms(E); a = cont(q); s7 = measS7(q)
        recs.append((E, a, s7))
        if E == tuple(range(k)): consec = (E, a, s7)
    return recs, consec

# even-Bonferroni truncation B_R(E) = sum_{|A|<=R} (-1)^|A| a[A](E), R even => upper bound on measS7
def B_R(a, R):
    return sum(((-1)**msize(A)) * a[A] for A in SUBMASKS if msize(A) <= R)

print("="*100)
print("THREAD B EXPLICIT CUT: even-Bonferroni / Mobius-truncation certificate B_R = sum_{|A|<=R}(-1)^|A| a[A]")
print("  B_R >= measS7 STRUCTURALLY for even R (Bonferroni). Need: B_R consec-MAXIMAL over stratum.")
print("  Smallest even R with consec=argmax(B_R) and B_R(consec) closing under cap = the signed cut.")
print("="*100)

caps = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

for k in [8, 9, 10]:
    recs, consec = build(k)
    Ec, ac, s7c = consec
    print(f"\nk={k}: #shapes={len(recs)} measS7(consec)={fmt(s7c)}={float(s7c):.5f} cap={fmt(caps[k])}={float(caps[k]):.5f}")
    for R in [0, 2, 4, 6]:
        bc = B_R(ac, R)
        # is consec the argmax of B_R over stratum? count beaters; find argmax value
        beaters = 0; maxv = bc; argmaxE = Ec
        for (E, a, s7) in recs:
            v = B_R(a, R)
            if v > bc: beaters += 1
            if v > maxv: maxv = v; argmaxE = E
        gap = bc - s7c   # B_R(consec) - measS7(consec) >=0; =0 iff R tight at consec
        is_argmax = (beaters == 0)
        print(f"  R={R}: B_R(consec)={float(bc):.5f}  consec=argmax? {is_argmax} (beaters={beaters})  "
              f"slack@consec={float(gap):.5f}  argmaxE={argmaxE if not is_argmax else 'consec'}")

print("\n" + "="*100)
print("READING: the smallest even R with consec=argmax => explicit signed cut B_R (Mobius parity signs,")
print("  support = all A with |A|<=R, UNIFORM in k by construction).  R where it FIRST works is the level.")
print("  If NO even R<=6 makes consec the argmax of the structural B_R => the manifestly-nonneg cut FAILS;")
print("  consec-max needs the LP-tuned (non-(-1)^|A|) magnitudes => report that as 'signs follow parity")
print("  but MAGNITUDES must be tuned (not pure truncation)'.")
print("="*100)
print("\nDONE.")
