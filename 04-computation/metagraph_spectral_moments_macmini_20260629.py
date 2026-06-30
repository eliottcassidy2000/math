#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Metagraph spectral MOMENTS -- a small new invariant extending klein's signed cycle index
P_n (THM-587), and the Burnside<->Siegel moment analogy for the LRC. (mac-mini-2026-06-29-S19)

klein THM-587: the metagraph eigenvalue multiplicities are the coeffs of
   P_n(x) = (1/n!) sum_{sigma in S_n} prod_{cycles c of sigma's signed action on the C(n,2) pairs}
            (1 + s_c x^{ell_c}),
with P_n(1)=A000568(n) (count), P_n(-1)=SC(n) (antipodal Euler #).

NEW (here): the LEVEL-MOMENTS of the spectrum. P_n(x)=sum_k mult(k) x^k, eigenvalue lambda=d-2k.
  - mean level kbar = P_n'(1)/P_n(1);  Var(k) = P_n''(1)/P_n(1) + kbar - kbar^2.
  - mean metagraph eigenvalue = d - 2 kbar;  spectral variance = 4 Var(k).
These are computable from n! permutations (past the 2^{C(n,2)} wall) -- a cheap new invariant.

Inspired by arXiv:2507.05905 (Han-Lee, Siegel-transform 1st/2nd MOMENTS with CONGRUENCE
conditions): P_n is the FINITE Burnside-average analog of the Siegel transform (SL_n(Z)-average);
P_n(1) is the 'first moment' (Siegel: mean count), and the spectral 2nd moment is the Burnside
analog of Siegel's SECOND moment -- the quantity the LRC floor (THM-579, sheet-count variance
CV(N_R)) needs. The congruence conditions of the paper match the COVERING (divisibility) structure.
"""
from __future__ import annotations
import functools, itertools, math
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def signed_cycle_poly(n):
    """Return P_n as a dict {level: mult} (Fraction mult, summed/n!)."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {p: t for t, p in enumerate(pairs)}
    m = len(pairs)
    # P_n(x) = (1/n!) sum_sigma prod_cycles (1 + s_c x^{ell_c}); accumulate coeffs
    from collections import defaultdict
    total = defaultdict(int)   # level -> count (integer, before /n!)
    for sigma in itertools.permutations(range(n)):
        # signed action on pairs: pair (i<j) -> {sigma[i],sigma[j]} with sign +1 if sorted-order kept
        seen = [False] * m
        factors = []   # list of (sign s_c, length ell_c)
        for start in range(m):
            if seen[start]:
                continue
            # follow the cycle
            cur = pairs[start]; length = 0; sign = 1
            while True:
                t = idx[tuple(sorted(cur))]
                if seen[t]:
                    break
                seen[t] = True; length += 1
                a, b = cur
                na, nb = sigma[a], sigma[b]
                if na > nb:
                    sign = -sign
                    na, nb = nb, na
                cur = (na, nb)
            factors.append((sign, length))
        # expand prod (1 + s_c x^{ell_c}) into level->coeff
        poly = {0: 1}
        for (s, l) in factors:
            new = defaultdict(int)
            for lev, co in poly.items():
                new[lev] += co            # the '1'
                new[lev + l] += co * s     # the s_c x^{ell_c}
            poly = new
        for lev, co in poly.items():
            total[lev] += co
    fact = math.factorial(n)
    return {lev: F(co, fact) for lev, co in total.items() if co != 0}, m


def main():
    print("=" * 80)
    print("Metagraph spectral moments (extending klein THM-587 P_n) + Siegel/LRC analogy")
    print("=" * 80)
    A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
    SC = {3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176}
    print(f"\n{'n':>2} {'d=C(n,2)':>8} {'P_n(1)=A000568':>14} {'P_n(-1)=SC':>11} "
          f"{'mean k':>8} {'Var(k)':>9} {'mean eig':>9} {'spec.var':>9}")
    for n in range(3, 8):
        poly, d = signed_cycle_poly(n)
        P1 = sum(poly.values())                       # = A000568 (Fraction -> int)
        Pm1 = sum(co * (-1)**lev for lev, co in poly.items())
        # moments of level k (weighted by mult)
        m0 = P1
        m1 = sum(F(lev) * co for lev, co in poly.items())
        m2 = sum(F(lev * lev) * co for lev, co in poly.items())
        kbar = m1 / m0
        vark = m2 / m0 - kbar * kbar
        mean_eig = d - 2 * kbar
        spec_var = 4 * vark
        ok1 = "OK" if int(P1) == A000568[n] else "MISMATCH"
        okm = "OK" if int(Pm1) == SC[n] else "MISMATCH"
        print(f"{n:>2} {d:>8} {str(int(P1)):>11}({ok1}) {str(int(Pm1)):>6}({okm}) "
              f"{float(kbar):>8.3f} {float(vark):>9.3f} {float(mean_eig):>9.3f} {float(spec_var):>9.3f}")
    print("\n  mean k -> d/2 (spectrum centered); the level VARIANCE is the new invariant = how")
    print("  spread the metagraph spectrum is = a 2nd moment, computable from n! (past the wall).")

    # H-distribution over iso classes (the LRC-relevant 'size' moment), small n exact
    print("\n  H-distribution over iso classes (H = #Ham paths; the metagraph vertex 'size'):")
    print(f"  {'n':>2} {'#classes':>8} {'mean H':>9} {'Var H':>11} {'CV(H)=sd/mean':>13}")
    for n in range(3, 7):
        classes = enumerate_iso_H(n)
        Hs = list(classes.values())
        meanH = sum(Hs) / len(Hs)
        varH = sum((h - meanH)**2 for h in Hs) / len(Hs)
        cv = math.sqrt(varH) / meanH if meanH else 0
        print(f"  {n:>2} {len(Hs):>8} {meanH:>9.2f} {varH:>11.2f} {cv:>13.4f}")
    print("\n  CV(H) over iso classes is the Burnside analog of Siegel's 2nd-moment CV; the LRC floor")
    print("  THM-579 needs CV(N_R) (sheet count). arXiv:2507.05905 computes such 2nd moments with")
    print("  CONGRUENCE conditions = the LRC COVERING (divisibility) structure. The bridge:")
    print("  metagraph P_n = Burnside avg over S_n  <->  Siegel transform = avg over SL_n(Z);")
    print("  both give 1st (count/mean) and 2nd (variance) moment formulas.")
    print("=" * 80)


def enumerate_iso_H(n):
    """Map each tournament iso class to its H (#Hamiltonian paths)."""
    import itertools as it
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    perms = list(it.permutations(range(n)))
    seen = {}
    for bits in range(1 << len(pairs)):
        arc = [[False]*n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if (bits >> b) & 1: arc[i][j] = True
            else: arc[j][i] = True
        # canonical form under relabeling
        canon = None
        for sg in perms:
            key = tuple(1 if arc[sg[i]][sg[j]] else 0 for i in range(n) for j in range(n) if i != j)
            if canon is None or key < canon: canon = key
        if canon in seen: continue
        # H = #Hamiltonian directed paths
        H = sum(1 for p in perms if all(arc[p[k]][p[k+1]] for k in range(n-1)))
        seen[canon] = H
    return seen


if __name__ == "__main__":
    main()
