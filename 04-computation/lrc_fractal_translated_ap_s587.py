#!/usr/bin/env python3
"""
claudebox-2026-06-03-S587 : the recursive fractal of translated APs.

S585 found that translating an AP above its own diameter destroys all 3-term relations (sum-free)
while keeping additive energy. RECURSE this: use a sum-free translated-AP as the DIGIT set of a
base-b construction with b > 2·max(digit), so there are no carries:

    S_d = { Σ_{i<d} a_i · b^i : a_i ∈ digits }    (a self-similar / Cantor-like set, |digits|^d speeds)

No carries ⇒ the additive structure factorizes across scales:
  * character polynomial is a RIESZ PRODUCT:  p_{S_d}(x) = Π_{i<d} D(x^{b^i}),  D(y)=Σ_{a∈digits} y^a
  * additive energy is MULTIPLICATIVE:  E(S_d) = E(digits)^d
  * 3-term relations are digit-wise ⇒ sum-free digits ⇒ sum-free at EVERY scale (fractal circuit-free)

Two fractal families:
  - UPPER-half digits {h+1,…,2h}  → recursively sum-free → Lemma A (safe, margin>0) at all scales
  - LOWER digits {1,…,k}          → recursively 3-term-rich → Lemma B (hard) at all scales

This file builds the fractal sets, verifies the factorizations, and tracks the gap across scales.
Rep-theory lens: Riesz product / lacunary Fourier.  Haskell lens: S_d = base-b eval of `replicateM d digits`.
"""

import itertools, math
from functools import reduce

def frac_dist(x):
    y = x - math.floor(x)
    return min(y, 1 - y)

def gap(S, N=24000):
    best = 0.0
    for a in range(1, N):
        t = a / N
        m = min(frac_dist(v * t) for v in S)
        if m > best:
            best = m
    return best

def good_measure(S, delta, N=24000):
    cnt = sum(1 for a in range(N) if all(frac_dist(v * (a / N)) >= delta for v in S))
    return cnt / N

def three_term_count(S):
    Sset = set(S); c = 0
    for a, b in itertools.combinations(S, 2):
        if a + b in Sset: c += 1
    for a in S:
        if 2 * a in Sset: c += 1
    return c

def additive_energy(S):
    from collections import Counter
    sums = Counter()
    for a in S:
        for b in S:
            sums[a + b] += 1
    return sum(v * v for v in sums.values())

def fractal_set(digits, base, d):
    """S_d = { Σ a_i base^i : a_i ∈ digits } — the level-d fractal, |digits|^d speeds."""
    out = []
    for combo in itertools.product(digits, repeat=d):
        out.append(sum(a * base ** i for i, a in enumerate(combo)))
    return sorted(out)

def riesz_check(digits, base, d):
    """Verify p_{S_d}(x) = Π_{i<d} D(x^{b^i}) as multisets of exponents (no-carry factorization)."""
    # LHS exponents = the fractal set itself; RHS = sumset of {a·b^i} blocks — identical by construction.
    S = fractal_set(digits, base, d)
    # RHS via iterated sumset of scaled digit blocks
    blocks = [[a * base ** i for a in digits] for i in range(d)]
    rhs = [0]
    for blk in blocks:
        rhs = [r + x for r in rhs for x in blk]
    return sorted(S) == sorted(rhs), len(S)

def main():
    print("=" * 78)
    print("S587  the recursive fractal of translated APs  (rep-theory: Riesz products)")
    print("=" * 78)

    # --- [1] the factorization laws: energy multiplicative, 3-term digit-wise ----
    print("\n[1] FACTORIZATION across scales (no-carry base):  E(S_d)=E(digits)^d, sum-free preserved")
    for name, digits in [("upper-half {4,5,6} (sum-free)", [4, 5, 6]),
                         ("lower {1,2,3} (has 1+2=3)", [1, 2, 3])]:
        base = 2 * max(digits) + 1
        Ed = additive_energy(digits)
        print(f"  digits={digits}  base={base}  E(digits)={Ed}  3term(digits)={three_term_count(digits)}")
        for d in range(1, 5):
            S = fractal_set(digits, base, d)
            E = additive_energy(S); tt = three_term_count(S)
            ok, sz = riesz_check(digits, base, d)
            print(f"    d={d}: |S|={sz:3d}  E(S_d)={E:7d}  E(digits)^d={Ed**d:7d}  "
                  f"match={E==Ed**d}  3term={tt:4d}  Riesz-factorization={ok}")

    # --- [2] the two fractal regimes: gap across scales --------------------------
    print("\n[2] THE TWO REGIMES — gap G vs δ=1/(|S|+1) across fractal scales")
    for name, digits in [("UPPER-half {4,5,6}: fractal SUM-FREE → Lemma A (safe)", [4, 5, 6]),
                         ("LOWER {1,2,3}: fractal 3-term-rich → Lemma B (hard)", [1, 2, 3])]:
        base = 2 * max(digits) + 1
        print(f"  {name}")
        print(f"    d | |S| | 3term | G       | δ        | margin G-δ | good-meas")
        for d in range(1, 4):
            S = fractal_set(digits, base, d)
            k = len(S); delta = 1 / (k + 1)
            tt = three_term_count(S)
            G = gap(S); gm = good_measure(S, delta)
            print(f"    {d} | {k:3d} |  {tt:4d} | {G:.5f} | {delta:.5f} | {G-delta:+.5f}   | {gm:.4f}")

    # --- [3] does the gap itself have a recursive/self-similar law? --------------
    print("\n[3] SELF-SIMILARITY of the gap: is G(S_d) governed by the single-digit gap?")
    print("    (G is dilation-invariant; the no-carry blocks are dilated copies at scale b^i)")
    digits = [4, 5, 6]; base = 2 * max(digits) + 1
    for d in range(1, 4):
        S = fractal_set(digits, base, d)
        G = gap(S)
        # single-block (digit) gap, and the coarsest-block gap
        Gdig = gap(digits)
        print(f"    d={d}: G(S_d)={G:.5f}   G(digits)={Gdig:.5f}   "
              f"ratio G_d/G_dig={G/Gdig:.4f}")
    print("    => if G(S_d) stabilizes, the fractal has a scale-invariant gap (a 'fractal dimension')")

    # --- [4] HASKELL lens: the fractal as base-b eval of replicateM -------------
    print("\n[4] HASKELL lens: S_d = map (base-b eval) (replicateM d digits) — a d-fold fold;")
    print("    REP lens: p_{S_d} = Π_{i<d} D(x^{b^i}) (Riesz product / lacunary spectrum).")
    digits = [4, 5, 6]; base = 13
    for d in [1, 2, 3]:
        S = fractal_set(digits, base, d)
        diffs = len(set(a - b for a in S for b in S))
        print(f"    d={d}: |S|={len(S)}, |S−S|={diffs}, |S−S|/|S|²={diffs/len(S)**2:.3f} "
              f"(lacunary ⇒ near-Sidon, difference set spreads)")

if __name__ == "__main__":
    main()
