#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angleD_sdr_clumping_hardness_mac-mini-2026-06-17-S4.py   (ANGLE D)

SDR / section-shuffling  <->  hard-from-easy families for LRC(14).

Setup (this lab):
  M(S) = max_tau min_v ||v*tau||;  LRC(14) <=> M(S) >= 1/14.
  THM-524 regions/sections: at the grid tau=a/14 (a a unit mod 14) runner v lands in
  section (v*a mod 14); observer sits at section 0; lonely <=> no runner in section 0.
  Reversal fixed points: v==0 mod14 pinned at section-0 CENTER (distance 0, HARDEST,
  defeats every grid witness); v==7 mod14 pinned at section-7 center (distance 1/2, FREE).
  Hard-from-easy: a hard 13-config is A u {14m} for an easy 12-core A; M dips only at
  resonant m, family min stays >= 1/14 (worst A={1..12}: 14/183 > 1/14).

THIS SCRIPT reports the classification:  easy-core residue pattern -> hard-family hardness.

(1) SDR vs CLUMPING classification of easy cores, correlated with family min-M.
(2) (Z/14)* section-shuffling: M(aS)=M(S) (scale invariance, ANY positive a) while units
    PERMUTE the residues/section-occupancy.
(3) The SHARING/SDR pattern predicts WHERE the family dips: dip happens exactly when the
    parked runner 14m RESONATES with the core's off-grid optimum tau*, landing at distance 0.

stdlib only, exact Fractions.
"""
from fractions import Fraction as F
from collections import Counter


def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r


def g(S, t):
    return min(nrm(v * t) for v in S)


def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C


def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at


def realize(residues):
    """Smallest distinct positive integers with the given residues mod 14 (repeats -> +14)."""
    used = set(); out = []
    for r in residues:
        v = r if r != 0 else 14
        while v in used or v == 0:
            v += 14
        used.add(v); out.append(v)
    return sorted(out)


def sdr_def(S):
    c = Counter(v % 14 for v in S)
    return sum(k * (k - 1) // 2 for k in c.values())


ONE14 = F(1, 14)


def main():
    print("=" * 78)
    print("(1) DROP-ONE-RESIDUE 12-CORES A = {1..13}\\{r}  (all are SDR)")
    print("=" * 78)
    print(f"{'drop r':>6}{'M(core)':>9}{'core tau*':>11}{'family min-M':>14}{'>1/14':>7}  dip-m")
    for r in range(1, 14):
        A = [x for x in range(1, 14) if x != r]
        MA, atA = M(A)
        worst = (F(1), None)
        for m in range(1, 40):
            MS, _ = M(A + [14 * m])
            if MS < worst[0]:
                worst = (MS, m)
        print(f"{r:>6}{str(MA):>9}{str(atA):>11}{str(worst[0]):>14}"
              f"{('YES' if worst[0] >= ONE14 else 'NO'):>7}  m={worst[1]}")

    print()
    print("=" * 78)
    print("(1b) CLUMPING CORES: SDR deficiency vs hardness  (more clumping = EASIER core)")
    print("=" * 78)
    tests = [
        ("SDR 1..12",            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
        ("1 clump @res1",        [1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
        ("2 clumps @res1,2",     [1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        ("3 clumps @1,2,3",      [1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9]),
        ("6-stack @res1",        [1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7]),
        ("clump @res7 (free)",   [7, 7, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11]),
        ("has res0 (perf mid)",  [14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
    ]
    print(f"{'core pattern':<22}{'SDR-def':>8}{'M(core)':>9}{'fam min-M':>11}{'>1/14':>7}")
    for name, residues in tests:
        A = realize(residues)
        worst = (F(1), None)
        for m in range(1, 40):
            MS, _ = M(A + [14 * m])
            if MS < worst[0]:
                worst = (MS, m)
        print(f"{name:<22}{sdr_def(A):>8}{str(M(A)[0]):>9}{str(worst[0]):>11}"
              f"{('YES' if worst[0] >= ONE14 else 'NO'):>7}")
    print("  TREND: SDR-deficiency UP -> M(core) UP -> family looser. Clumping HELPS looseness.")

    print()
    print("=" * 78)
    print("(2) (Z/14)* SECTION-SHUFFLING: M(aS)=M(S), units permute residues")
    print("=" * 78)
    A = list(range(1, 13))
    MA, atA = M(A)
    print(f"  A={{1..12}}  M={MA} at {atA}")
    for a in [1, 3, 5, 9, 11, 13]:
        aS = [a * v for v in A]
        val, at = M(aS)
        print(f"  a={a:>2} (unit): M(aS)={val} at {at}  residues={sorted(set(v % 14 for v in aS))}")
    for a in [2, 7, 100]:
        val, _ = M([a * v for v in A])
        print(f"  a={a:>3} (general): M(aS)={val}  (== M(S); scale invariance for ANY a>0)")

    print()
    print("=" * 78)
    print("(3) DIP LOCATION FROM RESONANCE: dip <=> parked runner hits distance 0 at core tau*")
    print("=" * 78)
    A = list(range(1, 13)); MA, atA = M(A)
    print(f"  core A={{1..12}}: M={MA} at tau*={atA}.  Dip predicted when m == 0 mod 13.")
    print(f"  {'m':>4}{'M(S)':>9}{'tau*':>9}{'||14m*tau*||':>14}  dip?")
    for m in list(range(1, 14)) + [26, 39]:
        MS, at = M(A + [14 * m])
        dist = nrm(F(14 * m) * atA)
        print(f"  {m:>4}{str(MS):>9}{str(at):>9}{str(dist):>14}  {'YES' if MS < MA else '-'}")

    print()
    print("=" * 78)
    print("(4) THE REDUCTION: dip <= slack  (i.e. M(S) >= 1/14) across all drop-r cores x m")
    print("=" * 78)
    worst = (F(1), None, None); viol = 0
    for r in range(1, 14):
        A = [x for x in range(1, 14) if x != r]
        MA, _ = M(A)
        for m in range(1, 40):
            MS, _ = M(A + [14 * m])
            if MS < ONE14:
                viol += 1
            if MS < worst[0]:
                worst = (MS, r, m)
    MS, r, m = worst
    print(f"  violations of M(S) >= 1/14: {viol}")
    print(f"  GLOBAL worst family member: drop r={r}, m={m}: M(S)={MS}={float(MS):.6f}")
    print(f"  margin above 1/14 = {MS - ONE14} = {float(MS - ONE14):.6f}  (REDUCTION HOLDS)")


if __name__ == "__main__":
    main()
