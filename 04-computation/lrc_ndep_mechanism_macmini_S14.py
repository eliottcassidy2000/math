#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S14 -- HYP-4422 (mechanism): WHY is the second gap nonempty
at n=7,8 but (believed) empty at n=13?

VERIFIED: n=7 (1,5,6,11,16,17) M=5/33 and n=8 M=3/23 are gap members, and the
n=7 one is a LIFTED TRANSVERSAL (mod 13 residues hit all 6 shells; 16=3+13,
17=4+13 are lifts).  So lifting a transversal LOWERS M into the gap at n=7 --
the M-minimizer property (my S11) FAILS at n=7 but HOLDS at n=13.

This isolates the n-dependence: for each n, search lifted transversals (the
n=7 mechanism) and find the minimal gap member (or its absence).  Characterize
the ESCAPE: at n=7, a 2-shell lift resonates at a small denominator (33=3*11)
BELOW 2/(2n-1); at n=13 does the analogous lift resonate in the gap, or does
25=5^2's structure block it?
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

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

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

def shells(C):
    """unit antipodal shells {a, C-a} mod C."""
    return [(a, C - a) for a in range(1, C // 2 + 1) if gcd(a, C) == 1 and a < C - a]

def search_lifted_transversals(n, n_iter=120000, kmax=4, seed=0):
    """search transversals mod C=2n-1 with some shells lifted by +C*k for a
    gap member.  Returns minimal gap member + nearest-miss."""
    C = 2 * n - 1
    k = n - 1
    lo, hi = F(1, n), F(2, 2 * n - 1)
    flo, fhi = float(lo), float(hi)
    sh = shells(C)
    if len(sh) < k:
        return None, None, C, len(sh)  # not enough unit shells (composite C)
    rng = random.Random(seed)
    best_gap = None
    nearest = (1.0, None)   # (float dist above hi that we got closest to gap from)
    tested = 0
    for _ in range(n_iter):
        # pick k shells (one element each), lift some by +C
        chosen_shells = rng.sample(sh, k)
        W = []
        for (a, b) in chosen_shells:
            base = rng.choice([a, b])
            base += C * rng.randint(0, kmax)
            W.append(base)
        W = primitive(tuple(sorted(set(W))))
        if len(W) != k:
            continue
        tested += 1
        fm = float_M(W)
        # track nearest approach to the gap from above
        if fm >= fhi and fm - fhi < nearest[0]:
            nearest = (fm - fhi, W)
        if flo - 1e-6 < fm < fhi + 1e-6:
            M = exact_M(W)
            if lo < M < hi:
                if best_gap is None or M < best_gap[0]:
                    best_gap = (M, W)
    return best_gap, nearest, C, len(sh)

def main():
    print("=" * 88)
    print("THE n-DEPENDENCE MECHANISM: lifted-transversal gap members per n")
    print("=" * 88)
    print(f"  {'n':>3} {'2n-1':>5} {'fact':>8} {'#unit shells':>12} {'k=n-1':>6} "
          f"{'gap member?':>28} {'nearest-miss above 2/(2n-1)':>28}")
    for n in range(6, 14):
        C = 2 * n - 1
        from sympy import factorint
        f = factorint(C)
        fs = "*".join(f"{p}^{e}" if e > 1 else str(p) for p, e in f.items())
        ni = 150000 if n <= 9 else (80000 if n <= 11 else 40000)
        bg, nm, C, nsh = search_lifted_transversals(n, n_iter=ni, kmax=4, seed=n)
        gapstr = f"M={bg[0]}={float(bg[0]):.4f}" if bg else "NONE FOUND"
        nmstr = f"+{nm[0]:.5f}" if nm and nm[1] else "--"
        note = "" if nsh >= n - 1 else f" (only {nsh}<{n-1} unit shells!)"
        print(f"  {n:>3} {C:>5} {fs:>8} {nsh:>12} {n-1:>6} {gapstr:>28} {nmstr:>28}{note}")
    print()
    print("  READING: at n where a lifted transversal lands in the gap, (G) FAILS")
    print("  (n=7,8 confirmed).  Where 'NONE FOUND' + nearest-miss is bounded away,")
    print("  the lift MECHANISM is blocked -- the arithmetic of 2n-1 prevents the")
    print("  small-denominator resonance.  n=13 (2n-1=25=5^2): the (G) target.")
    print()
    print("  KEY COMPARISON -- the n=7 member's escape vs n=13 analog:")
    # n=7 member and a direct n=13 analog (lift 2 shells of the mod-25 AP transversal)
    m7 = (1, 5, 6, 11, 16, 17)
    print(f"    n=7  {list(m7)}: M={exact_M(m7)} (gap (1/7,2/13)); "
          f"lifted shells 3->16, 4->17; resonance denom {exact_M(m7).denominator}")
    # n=13 analog: AP {1..12} with 2 elements lifted by +25 (mimicking +C lift)
    for lift_pair in [(3, 4), (5, 7), (2, 9)]:
        W = list(range(1, 13))
        for x in lift_pair:
            W[W.index(x)] = x + 25
        W = primitive(tuple(sorted(W)))
        M = exact_M(W)
        lo, hi = F(1, 13), F(2, 25)
        tag = "IN-GAP!" if lo < M < hi else ("clears" if M >= hi else "tight/below")
        print(f"    n=13 AP lift {lift_pair}->+25: {list(W)} M={M}={float(M):.4f} [{tag}]")

if __name__ == "__main__":
    main()
