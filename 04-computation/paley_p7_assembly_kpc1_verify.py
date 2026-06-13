#!/usr/bin/env python3
"""
paley_p7_assembly_kpc1_verify.py
ADVERSARIAL VERIFIER, kind-pasteur-2026-06-10-S1, thread E-R31-prediction.

Verifies claim E3 end-to-end at p=7 by FULL brute force over all 7! = 5040
orderings of F_7 (no DP, no engine -- a method independent of everything the
worker ran):

  (a) H(T_7) = number of orderings whose consecutive differences are all QR.
  (b) For EVERY subset S of the 6 edge slots:
        T(S) = sum over orderings sigma of prod_{k in S} chi(a_{k+1}-a_k)
      and check the identity  T(S) = A_joint(runs(S)) * (7-n)!  where
      A_joint is computed by separate injective-tuple brute force and
      n = sum (L_i + 1) over the maximal runs of S.
  (c) R(7) = sum_S T(S) / 7!  ==  H * 2^6 / 7!  ==  12/5  exactly, and the
      claimed decomposition  1 + 1 + 2/5 + 1/5 - 1/5  by run-class:
        empty -> 1, {2} -> 1, {4} -> 2/5, {6} -> 1/5, {2,2} -> -1/5,
        all classes containing an odd run -> 0 (E2 lemma check).
"""
from fractions import Fraction
from itertools import permutations
from math import factorial
from collections import defaultdict

p = 7
chi = [0] * p
qr = set((x * x) % p for x in range(1, p))
for d in range(1, p):
    chi[d] = 1 if d in qr else -1

# ---- (a) + (b): one pass over all orderings -------------------------------
H = 0
T = [0] * 64                      # T[S] over subsets of 6 slots
for sigma in permutations(range(p)):
    cv = [chi[(sigma[k + 1] - sigma[k]) % p] for k in range(6)]
    if all(c == 1 for c in cv):
        H += 1
    # accumulate T(S) for all 64 subsets via subset products
    prod = [1] * 64
    for S in range(64):
        w = 1
        for k in range(6):
            if (S >> k) & 1:
                w *= cv[k]
        T[S] += w

print(f"H(T_7) by full enumeration over 5040 orderings = {H}  "
      f"(claimed 189: {'CONFIRMED' if H == 189 else 'REFUTED'})")

# run decomposition of a slot subset
def runs_of(S):
    rr, l = [], 0
    for k in range(6):
        if (S >> k) & 1:
            l += 1
        else:
            if l:
                rr.append(l)
            l = 0
    if l:
        rr.append(l)
    return tuple(sorted(rr))

# joint integrals by separate injective-tuple brute force
def joint(runlens):
    if not runlens:
        return 1
    npts = sum(l + 1 for l in runlens)
    if npts > p:
        return 0
    prev, pos = [], 0
    for l in runlens:
        prev.append(-1)
        for i in range(l):
            prev.append(pos + i)
        pos += l + 1
    tot = 0
    for tup in permutations(range(p), npts):
        w = 1
        for i in range(npts):
            if prev[i] >= 0:
                w *= chi[(tup[i] - tup[prev[i]]) % p]
        tot += w
    return tot

JT = {}
ok_identity = True
class_T = defaultdict(int)
for S in range(64):
    lam = runs_of(S)
    if lam not in JT:
        JT[lam] = joint(lam)
    n = sum(l + 1 for l in lam)
    expect = JT[lam] * factorial(p - n) if n <= p else 0
    if T[S] != expect:
        ok_identity = False
        print(f"  IDENTITY FAILS at S={S:06b} runs={lam}: T={T[S]} expected={expect}")
    class_T[lam] += T[S]
print(f"identity T(S) = A_joint(runs) * (7-n)!  for ALL 64 subsets: "
      f"{'CONFIRMED' if ok_identity else 'REFUTED'}")

# ---- (c) exact assembly ----------------------------------------------------
print()
print("per-run-class contributions to R(7) (exact Fractions):")
total = Fraction(0)
for lam in sorted(class_T, key=lambda t: (len(t), t)):
    fr = Fraction(class_T[lam], factorial(p))
    total += fr
    odd = any(l % 2 == 1 for l in lam)
    print(f"  runs={str(lam):>10}:  sum T = {class_T[lam]:>8d}   contributes {fr}"
          f"{'   (odd run -> 0 by E2 lemma: ' + ('OK' if fr == 0 else 'VIOLATED') + ')' if odd else ''}")
R_direct = Fraction(H * 2 ** 6, factorial(p))
print(f"\nsum over all classes = {total}")
print(f"H * 2^6 / 7!         = {R_direct}")
print(f"equal and == 12/5:    "
      f"{'CONFIRMED' if total == R_direct == Fraction(12, 5) else 'REFUTED'}")

dec = {(): Fraction(1), (2,): Fraction(1), (4,): Fraction(2, 5),
       (6,): Fraction(1, 5), (2, 2): Fraction(-1, 5)}
ok_dec = all(Fraction(class_T[lam], factorial(p)) == v for lam, v in dec.items())
others_zero = all(Fraction(class_T[lam], factorial(p)) == 0
                  for lam in class_T if lam not in dec)
print(f"claimed decomposition 1 + 1 + 2/5 + 1/5 - 1/5 (and all else 0): "
      f"{'CONFIRMED' if ok_dec and others_zero else 'REFUTED'}")
