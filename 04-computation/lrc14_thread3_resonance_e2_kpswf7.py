#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- ANGLE 2: quasimodular E_2 / theta-correlation structure of the RESONANCE ERROR.

The WIDE bound is closed up to the SIGNED resonance error
   err(E) = p0(E) - p0_decorr(E),
where p0 = true far-coverage, p0_decorr = the moment-dual/Riesz product baseline (independent
clocks). |err| <= 0.012 << margin 0.132 over commensurable far ratios (kps-S24). To CLOSE it we
need a PROVEN bound on err.

ANGLE 2 claim: err is an OFF-DIAGONAL correlation sum. The decorrelated baseline = the DIAGONAL
(self) terms of a Walsh/Fourier square; err = the cross terms = sum over nonzero resonances
n in Lambda(E) (the resonance lattice {n : sum n_i e_i = 0}) of a Fourier coefficient. For a
commensurable far pair (Cq, Cp) the resonances are at lattice points -> a sum like
   sum_{(a,b)!=0, aq=bp... } 1/(stuff)  ~ Eisenstein E_2-type off-diagonal sum.
If err is dominated by the LEADING resonance (lowest harmonic), a Kronecker-limit / E_2 bound
gives a closed-form cap on |err|.

5-min PROBE:
  (1) compute err(E) exactly for commensurable far pairs (Cq, Cp), C large, over a grid of
      ratios p/q and bounded bases; confirm |err| <= 0.012 and find the worst ratio.
  (2) test the HARMONIC decay: does |err(p/q)| ~ const / (lowest resonance harmonic)? i.e. plot
      |err| vs 1/lcm(p,q) or vs the leading Fourier mode. If err ~ C/N with small C, an E_2
      tail bound (sum 1/n^2 type) closes it.
  (3) test SIGN structure: is err mostly NEGATIVE (resonance REDUCES coverage below decorrelated,
      so decorrelated is already an UPPER bound and err only helps)? If err<=0 generically, then
      p0 <= p0_decorr <= Q(k-1) with NO error term needed for the upper bound!
"""
import sys, random
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm
sys.stdout.reconfigure(line_buffering=True)

def covered6(E, x):
    # cover all 6 nonzero residues mod 7
    s = set()
    for e in E:
        v = F(e) * x; v = v - F(v.numerator // v.denominator)
        s.add((v.numerator * 7) // v.denominator)
    return s >= set(range(1, 7)) or len(s & set(range(1, 7))) == 6

def p0(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in E:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        if len(s & set(range(1, 7))) == 6: tot += hi - lo
    return tot

# decorrelated baseline: treat the far pair as 2 INDEPENDENT uniform clocks on the base coverage.
# p0_decorr(B u {f1,f2}) = sum over base-x of P(both far clocks land in the "needed" residues).
# Use the moment-dual closed form proxy: integrate over base the prob each missing residue is hit
# independently by each far clock (uniform). Equivalent to averaging the base-conditional coverage
# with INDEPENDENT far residues. Compute by Monte-Carlo over far residues (exact in the limit) OR
# the simple surrogate: p0_decorr = E_x[ prod over missing residues (1 - (1 - 1/7)^{#far}) ]...
# Simplest faithful surrogate: replace the two far speeds by INDEPENDENT random residues and average.
def p0_decorr(B, nfar=2, samples=4000, seed=0):
    # average over base x of the probability that nfar independent uniform Z/7 clocks complete the cover
    Bs = sorted(set(int(e) for e in B))
    bps = {F(0), F(1)}
    for e in Bs:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    rng = random.Random(seed)
    tot = 0.0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in Bs:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        need = set(range(1, 7)) - s
        if not need:
            tot += float(hi - lo); continue
        # prob nfar independent uniform residues (in 0..6) cover 'need'
        hit = 0
        for _ in range(samples):
            got = set(rng.randrange(7) for _ in range(nfar))
            if need <= got: hit += 1
        tot += float(hi - lo) * hit / samples
    return tot

def main():
    print("=== ANGLE 2: resonance error structure (sign, harmonic decay) ===\n")
    C = 350
    ratios = []
    for q in range(1, 7):
        for p in range(q + 1, 3 * q + 1):
            if gcd(p, q) == 1 and 1 < p / q <= 2.2:
                ratios.append((p, q))
    ratios = sorted(set(ratios), key=lambda t: t[0] / t[1])
    print(f"ratios p/q: {[f'{p}/{q}' for p,q in ratios]}\n")

    for k in (8, 9):
        base = list(range(k - 2))  # consec base
        bdec = p0_decorr(base, nfar=2, samples=6000, seed=k)
        print(f"k={k} consec base={base}  p0_decorr={bdec:.4f}")
        rows = []
        for (p, q) in ratios:
            E = sorted(set(base + [C * q, C * p]))
            if len(E) != k or reduce(gcd, E) != 1: continue
            true = float(p0(E))
            err = true - bdec
            N = lcm(p, q)
            rows.append((p, q, true, err, N))
        rows.sort(key=lambda r: r[3])
        worst_neg = min(r[3] for r in rows); worst_pos = max(r[3] for r in rows)
        print(f"   err range: [{worst_neg:+.4f}, {worst_pos:+.4f}]   (|err| max = {max(abs(worst_neg),abs(worst_pos)):.4f})")
        npos = sum(1 for r in rows if r[3] > 1e-9); nneg = sum(1 for r in rows if r[3] < -1e-9)
        print(f"   sign: {nneg} negative, {npos} positive (out of {len(rows)})  "
              f"{'=> err mostly <=0 => decorr is UPPER bound!' if npos==0 else '=> mixed sign'}")
        print("   harmonic-decay check (|err| vs 1/lcm(p,q)):")
        for p, q, true, err, N in sorted(rows, key=lambda r: -abs(r[3]))[:6]:
            print(f"     {p}/{q}: true={true:.4f} err={err:+.4f} lcm={N} |err|*lcm={abs(err)*N:.3f}")
        print()
    print("VERDICT: if err<=0 always => decorrelated IS an upper bound (error term unneeded!).")
    print("         if |err|*lcm bounded => E_2/harmonic tail bound closes it.")

if __name__ == "__main__":
    main()
