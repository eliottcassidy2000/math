#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Perfect numbers ARE arc-hypercube dimensions at Forcade orders; arc-flip (Q_d-edge)
invariance of the per-type Ky Fan parity; the Mersenne-DRT doubling gauge.
(mac-mini-2026-06-29-S13)

THREE THREADS the owner named -- symmetry gauge, Mersenne DRT, perfect numbers -- merged
into the OPEN-Q-059 / Ky-Fan thread (per-type oriented-path parity, HYP-3545; klein THM-587
signed cycle index; klein THM-584 complement=antipodal of the arc-hypercube Q_d, d=C(n,2)).

(A) PERFECT-NUMBER = ARC-HYPERCUBE DIMENSION at FORCADE orders.
    Forcade 1973: at n=2^k EVERY oriented type has ODD count (verified HYP-3545).
    The arc-hypercube Q_d has dimension d = C(n,2) = #arcs. For n=2^p:
        d = C(2^p,2) = 2^{p-1}(2^p-1) = a PERFECT NUMBER  <=>  2^p-1 is a Mersenne PRIME.
    So the Forcade orders n=4,8,32,128 (p=2,3,5,7 Mersenne-prime) are EXACTLY where the
    Ky-Fan count is maximally degenerate (all types odd) AND d is perfect.

(B) ARC-FLIP (Q_d-EDGE) INVARIANCE of the per-type parity: N_tau(T) mod 2 is invariant
    under flipping ONE arc => constant on the connected hypercube Q_d => Forcade. This is
    the discrete realization: Forcade = Q_d-edge-invariance = the discrete Ky Fan.

(C) MERSENNE-DRT doubling gauge: B_0(T_{2m-1}) = T_{m-1} (out-nbhd of vertex 0 = previous
    level). Mersenne orders 3,7,15,31. n=7 (Paley/DRT) -> B_0 = T_3. Per-type structure of
    the DRT (the Ham-Sandwich score-balanced point).
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)


def is_prime(m):
    if m < 2: return False
    i = 2
    while i * i <= m:
        if m % i == 0: return False
        i += 1
    return True


def sigma(m):
    s = 0
    for d in range(1, m + 1):
        if m % d == 0: s += d
    return s


# ---------- (A) perfect numbers = arc-hypercube dimension at Forcade orders ----------
def thread_A():
    print("=" * 78)
    print("(A) PERFECT NUMBER = arc-hypercube dim d=C(n,2) at FORCADE orders n=2^p")
    print("=" * 78)
    print(f"{'p':>2} {'n=2^p':>6} {'d=C(n,2)':>9} {'=2^(p-1)(2^p-1)':>16} "
          f"{'2^p-1 prime?':>12} {'perfect?':>9}")
    for p in range(2, 8):
        n = 2 ** p
        d = n * (n - 1) // 2
        mers = 2 ** p - 1
        mp = is_prime(mers)
        perfect = (sigma(d) == 2 * d)
        print(f"{p:>2} {n:>6} {d:>9} {2**(p-1)*(2**p-1):>16} "
              f"{('PRIME' if mp else 'comp'):>12} {('PERFECT' if perfect else '-'):>9}")
    print("\n  => d=C(2^p,2) is PERFECT exactly when 2^p-1 is a Mersenne prime (p=2,3,5,7).")
    print("  These n=4,8,32,128 are the Forcade all-types-odd orders.")
    print("  LRC(14): apex prime 7=2^3-1 (Mersenne) -> Forcade order 8 -> d=28=2*14 (perfect).")
    print("  n=8 = Havet-Thomasse/Rosenfeld threshold = arXiv:2512 arc-deletion threshold.\n")


# ---------- (B) arc-flip (Q_d-edge) invariance of per-type parity ----------
def type_of(perm, arc, n):
    return tuple(1 if arc[perm[k]][perm[k + 1]] else 0 for k in range(n - 1))


def counts_by_type(arc, n):
    from collections import Counter
    c = Counter()
    for perm in itertools.permutations(range(n)):
        c[type_of(perm, arc, n)] += 1
    return c


def thread_B(n=5):
    print("=" * 78)
    print(f"(B) ARC-FLIP (Q_d-edge) invariance of N_tau mod 2  (n={n}, all tournaments)")
    print("=" * 78)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    bad = 0
    checked = 0
    # for each tournament and each single arc-flip, check N_tau mod2 unchanged for ALL tau
    for bits in range(1 << m):
        arc = [[False] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if (bits >> b) & 1: arc[i][j] = True
            else: arc[j][i] = True
        base = counts_by_type(arc, n)
        for b, (i, j) in enumerate(pairs):
            arc[i][j] = not arc[i][j]; arc[j][i] = not arc[j][i]   # flip arc {i,j}
            flipped = counts_by_type(arc, n)
            allt = set(base) | set(flipped)
            if any((base.get(t, 0) - flipped.get(t, 0)) % 2 for t in allt):
                bad += 1
            checked += 1
            arc[i][j] = not arc[i][j]; arc[j][i] = not arc[j][i]   # restore
    print(f"  checked {checked} (tournament, arc-flip) pairs; parity-CHANGING flips: {bad}")
    print(f"  => single arc-flip preserves N_tau mod 2 for ALL types: "
          f"{'CONFIRMED' if bad == 0 else 'FAILS'}")
    print("  Q_d is connected by arc-flips => N_tau mod 2 is CONSTANT = Forcade.")
    print("  i.e. Forcade's theorem = Q_d-EDGE-invariance of the graded Ky Fan count.\n")


# ---------- (C) Mersenne-DRT doubling gauge: the Paley/DRT tournament on 7 ----------
def paley(q):
    """Paley tournament on q vertices (q prime ==3 mod 4): i->j iff (j-i) is a QR mod q."""
    qrs = set((x * x) % q for x in range(1, q))
    arc = [[False] * q for _ in range(q)]
    for i in range(q):
        for j in range(q):
            if i != j and ((j - i) % q) in qrs:
                arc[i][j] = True
    return arc


def thread_C():
    print("=" * 78)
    print("(C) MERSENNE-DRT doubling gauge: Paley T_7 (=DRT on 7=2^3-1) and B_0(T_7)")
    print("=" * 78)
    q = 7
    arc = paley(q)
    # out-neighborhood of vertex 0
    out0 = [v for v in range(q) if arc[0][v]]
    print(f"  Paley T_7: out-nbhd of vertex 0 = {out0} (size {len(out0)} = (7-1)/2 = 3).")
    # the induced sub-tournament on out0 should be ~ T_3 (the 3-cycle) -- the doubling B_0
    sub = [[arc[a][b] for b in out0] for a in out0]
    # count its directed-3cycle / check it is the 3-cycle (each vertex out-deg 1)
    outdeg = [sum(1 for b in range(3) if sub[a][b]) for a in range(3)]
    print(f"  induced sub-tournament on out-nbhd: out-degrees {outdeg} "
          f"=> {'3-CYCLE (=T_3, the Mersenne doubling B_0(T_7)=T_3)' if sorted(outdeg)==[1,1,1] else 'other'}")
    # per-type counts of T_7 (the DRT) -- Ham-Sandwich balanced point
    c = counts_by_type(arc, q)
    odd = sum(1 for t in c if c[t] % 2 == 1)
    tot = len(c)
    from statistics import pstdev, mean
    vals = list(c.values())
    print(f"  Paley T_7 per-type counts: {tot} types present, {odd} odd; "
          f"mean={mean(vals):.1f}, stdev={pstdev(vals):.1f} (DRT = maximally balanced gauge).")
    print(f"  score sequence (all equal = regular = Ham-Sandwich score-bisection fixed locus): "
          f"{[sum(1 for b in range(q) if arc[a][b]) for a in range(q)]}\n")


if __name__ == "__main__":
    thread_A()
    thread_B(5)
    thread_C()
    print("=" * 78)
    print("MERGE: the Forcade/Ky-Fan all-odd orders n=2^p are the PERFECT-NUMBER arc-hypercube")
    print("dimensions d=C(n,2)=2^{p-1}(2^p-1) (Mersenne-prime p); arc-flip = Q_d edge realizes")
    print("Forcade as discrete Ky Fan; the DRT is the Ham-Sandwich-balanced gauge; the Mersenne")
    print("doubling B_0 descends it. LRC apex 7 -> order 8 -> perfect 28 = 2*14.")
    print("=" * 78)
