#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S27 -- HYP-4562: THE CROSS-n Dx<D (Selberg spacing) closure.

kps-S36's general ladder M(B+x)=mu*x/(x+rho) samples x at resonances spaced D
apart; a rung lands in the gap (1/(k+1), 2/(2k+1)) iff some jD sits in the
x-window (LO*rho/(mu-LO), HI*rho/(mu-HI)) of width Dx.  My assigned lane (kps):
make "Dx < D, no rung in gap" COMPREHENSIVE + UNIFORM at n=12, and understand
the CROSS-n pattern (why small k catches, k=12 skips).

For each k=5..13: sweep many (k-1)-speed bases + outlier x; report catching
bases (rung in gap) and the crossover.  Also n=12 comprehensive (many bases).
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def float_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = 0.0
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

def bases_for_k(k, rng):
    """many (k-1)-speed bases: defected {1..k}, AP-cores, near-AP, dilated."""
    out = set()
    full = list(range(1, k + 1))  # {1..k} is k speeds; base is k-1
    # single-defect {1..k}\{i}
    for i in range(1, k + 1):
        out.add(tuple(v for v in full if v != i))
    # {1..k}\{i} with one element shifted (defected, rho>=2)
    for i in range(1, k + 1):
        B = [v for v in full if v != i]
        for j in range(len(B)):
            for d in (1, 2, -1):
                C = list(B); C[j] += d
                if len(set(C)) == k - 1 and min(C) >= 1:
                    out.add(primitive(tuple(sorted(C))))
    # random near-AP (k-1 speeds, height <= 2k)
    for _ in range(400):
        B = primitive(tuple(sorted(rng.sample(range(1, 2 * k + 2), k - 1))))
        if len(B) == k - 1:
            out.add(B)
    return [B for B in out if len(B) == k - 1]

def main():
    print("=" * 80)
    print("CROSS-n: does ANY (k-1)-speed base + outlier catch the gap (1/(k+1),2/(2k+1))?")
    print("=" * 80)
    rng = random.Random(27)
    print(f"  {'k(=n-1)':>7} {'gap':>18} {'#bases':>7} {'catching bases':>15} {'min in-gap M':>14}")
    for k in range(5, 14):
        lo, hi = F(1, k + 1), F(2, 2 * k + 1)
        flo, fhi = float(lo), float(hi)
        bases = bases_for_k(k, rng)
        catchers = []
        min_gap = None
        xmax = 5 * k + 10
        for B in bases:
            mx = max(B)
            for x in range(mx + 1, xmax):
                W = primitive(tuple(sorted(set(B) | {x})))
                if len(W) != k:
                    continue
                fm = float_M(W)
                if flo - 1e-6 < fm < fhi + 1e-6:
                    M = exact_M(W)
                    if lo < M < hi:
                        catchers.append((M, B, x, W))
                        if min_gap is None or M < min_gap:
                            min_gap = M
        cbases = len(set(c[1] for c in catchers))
        print(f"  {k:>7} {f'(1/{k+1},2/{2*k+1})':>18} {len(bases):>7} {cbases:>15} "
              f"{(str(min_gap) if min_gap else 'EMPTY'):>14}")
        # show a catcher example for small k (validation)
        if catchers and k <= 8:
            M, B, x, W = sorted(catchers)[0]
            print(f"         e.g. base {list(B)} + {x} = {list(W)} -> M={M}")
    print()
    print("  => the CROSSOVER k* = largest k with a catching base.  For k > k*,")
    print("     NO single-outlier family lands in the gap (uniform Dx<D).  If k=12")
    print("     is EMPTY (> k*), the single-outlier bounded case CLOSES at LRC(14).")
    print("     Small k catch (denser ladder j/(kj+r), rung inside); large k skip")
    print("     (Farey rungs 1/(k+1),2/(2k+1) bracket the gap).")

if __name__ == "__main__":
    main()
