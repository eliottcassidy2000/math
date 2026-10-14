#!/usr/bin/env python3
"""Audit A (independent), THM-4640 (8) / HYP-9270 (b): Haar coverage of odd 2-adic sources by downward uniform
switches at -1, recomputed exactly by a different method, split into its parts, and compared with the root switch.

Source law: U-word letters i.i.d. P(a) = 2^-a; the word is 1^r u, u reduced; P(r) = 2^-(r+1); within u the first
letter a_1 >= 2 has P = 2^(1-a_1) and later letters 2^-a.
Switch families (THM-4555 (iv)) from a canonical primitive pair (s, t), |s| < |t|, D = |t| - |s| (even >= 2):
  canonical: cylinder 1^r head(s) [next letter >= cmin], r >= D,       cmin = max(1, 1 + e_t - e_s)
  twin     : cylinder 1^r head((1).s) [next letter >= cmin], r >= D - 1  (head starts with 2: reset-2 sources)
  root     : the reset switch (a) ~ (2, a-2): cylinder 1^r [a_1 >= 3], r >= 1 (NOT included in the audited c(E))
Exact union measure: sum over heads h of the family set that are 'first covering' (no proper prefix h' in the set
with h[len h'] >= cmin(h')) of P(prefix h) P(next >= cmin(h)) = 2^(2 - A_h - cmin(h)); integer arithmetic.
Usage: python3 a7_coverage.py [EMAX=24]"""
import sys, time, math
from fractions import Fraction as Fr
import numpy as np

EMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
t0 = time.time()


def canon_reps(K, emax):
    out = []

    def back(K, ehi, tail):
        if K == -1:
            out.append(tuple(reversed(tail)))
            return
        for e in range(2, ehi + 1):
            r = K - (1 << e)
            if r % 3:
                continue
            P = r // 3
            if P == -1 or 1 <= P <= 3 ** (e - 1) - 2 ** e:
                back(P, e - 1, tail + [e])

    back(K, emax, [])
    return out


def head_word(rep):
    if not rep:
        return ()
    return (rep[0] + 1,) + tuple(rep[j] - rep[j - 1] for j in range(1, len(rep)))


# primitive values per level (set recursion) with their birth level
S = np.array([-1], dtype=np.int64); P = np.array([], dtype=np.int64)
births = {}
for e in range(1, EMAX + 1):
    T = 3 * S + (1 << e)
    idx = np.searchsorted(S, T); idx[idx == len(S)] = len(S) - 1
    hit = T[S[idx] == T]
    for k in np.setdiff1d(hit, P, assume_unique=True).tolist():
        births[k] = e
    P = np.union1d(P, hit)
    S = np.union1d(S, T)
del S
# all canonical primitive pairs with exponents <= EMAX, tagged by their max exponent
pairs = []
for k, b in births.items():
    if k == -1:
        continue
    reps = canon_reps(k, EMAX)
    for i in range(len(reps)):
        for j in range(i + 1, len(reps)):
            s, t = reps[i], reps[j]
            if s[-1] == t[-1]:
                continue
            if len(s) > len(t):
                s, t = t, s
            pairs.append((max(s[-1], t[-1]), s, t))
print(f"{len(pairs)} canonical primitive pairs with exponents <= {EMAX}  [{time.time()-t0:.1f}s]", flush=True)


def union_measure(fams):
    """fams: dict head -> cmin; exact measure (as Fraction) under the reduced-word law"""
    tot = 0
    SH = 200
    for h, c in fams.items():
        okfirst = True
        for L in range(1, len(h)):
            c2 = fams.get(h[:L])
            if c2 is not None and h[L] >= c2:
                okfirst = False
                break
        if okfirst:
            tot += 1 << (SH + 2 - sum(h) - c)
    return Fr(tot, 1 << SH)


def coverage(E, kinds):
    fam_by_D = {}
    for (mx, s, t) in pairs:
        if mx > E:
            continue
        D = len(t) - len(s)
        cmin = max(1, 1 + t[-1] - s[-1])
        if "canonical" in kinds:
            fam_by_D.setdefault(D, []).append((head_word(s), cmin))
        if "twin" in kinds:
            fam_by_D.setdefault(D - 1, []).append((head_word((1,) + s), cmin))
    if not fam_by_D:
        return Fr(0), []
    R0 = max(fam_by_D)
    total = Fr(0); cur = {}
    for r in range(0, R0 + 1):
        for (h, c) in fam_by_D.get(r, []):
            cur[h] = min(cur.get(h, 10 ** 9), c)
        u = union_measure(cur) if cur else Fr(0)
        weight = Fr(1, 2 ** (r + 1)) if r < R0 else Fr(1, 2 ** R0)    # r >= R0 all have the same union
        total += weight * u
    return total, sorted(fam_by_D)


print(" E   c(E) [canonical+twin]   canonical part   twin part   c/(3/8)   twin/(1/4)   root-included 1/4+twin "
      " E*(c(E)-c(E-2))", flush=True)
prev = None
for E in range(8, EMAX + 1, 2):
    cc, Ds = coverage(E, ("canonical", "twin"))
    ca, _ = coverage(E, ("canonical",))
    tw, _ = coverage(E, ("twin",))
    assert cc == ca + tw and ca <= Fr(1, 8) and tw <= Fr(1, 4)
    inc = f"{(E - 1) * float(cc - prev):.4f}" if prev is not None else "-"
    print(f"{E:2d}   {float(cc):.6f}              {float(ca):.6f}         {float(tw):.6f}    {float(cc / Fr(3, 8)):.4f}    "
          f"{float(tw * 4):.4f}       {float(Fr(1, 4) + tw):.6f}        {inc}   (switch distances {Ds})  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    prev = cc
