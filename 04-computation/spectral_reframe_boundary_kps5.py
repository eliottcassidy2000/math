#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Where does the spectral reframe STOP? Is H (the OCF Hamiltonian-path count)
determined by the spectrum, like the cycle counts c_k = tr(A^k)/k are?
kind-pasteur-2026-06-13-S5.

This session (THM-498/HYP-2492): cycle-count gaps are SKEW-SPECTRAL (power-sum)
exclusions, because c_k = tr(A^k)/k for k<=5.  HYP-2492 conjectured H in {7,21}
are also power-sum exclusions.  But H = I(Omega,2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
depends on the DISJOINTNESS structure (alpha_2 = # vertex-disjoint odd-cycle pairs),
which is NOT a pure spectral quantity.  So the prediction is that H is NOT
spectrally determined -- and finding COSPECTRAL tournaments with DIFFERENT H pins
exactly where the spectral reframe stops (cycle counts: yes; H: no, needs the
conflict graph).

The spectral signature (exact, integer): sig(T) = (tr A, tr A^2, ..., tr A^n)
<=> the characteristic polynomial of A <=> the spectrum.  (Note tr A = 0, tr A^2 = 0
since no self-loops/2-cycles; tr A^k = k*c_k for k=3,4,5.)

TESTS (exhaustive small n):
 (A) group tournaments by spectral signature; does any spectral class carry >1 H value?
     (cospectral but different H => H NOT spectrally determined; report the witness pair.)
 (B) group by the cycle-count vector (c3,c5,c7); same question (does H need MORE than
     the odd-cycle counts -- i.e. the disjointness alpha_2?).
 (C) the converse sanity: confirm c3,c5 ARE spectrally determined (must be, =tr/k).
"""

import sys, itertools
from collections import defaultdict
import numpy as np
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def all_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        yield A


def spectral_sig(A, n):
    sig = []
    P = np.eye(n, dtype=np.int64)
    for k in range(1, n + 1):
        P = P @ A
        sig.append(int(np.trace(P)))
    return tuple(sig)


def ham_path_count(A, n):
    """H(T) = number of Hamiltonian paths (orderings v1..vn with vi->v(i+1))."""
    # DP over subsets: dp[mask][v] = # paths covering 'mask' ending at v
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v, w]:
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))


def cycle_counts(A, n):
    c3 = int(np.trace(A @ A @ A)) // 3
    A2 = A @ A; A4 = A2 @ A2
    c5 = int(np.trace(A4 @ A)) // 5
    # c7 = tr(A^7) minus degeneracy (k>=6 corrections) -- skip exact c7; use c3,c5
    return (c3, c5)


def main():
    for n in (5, 6):
        by_sig = defaultdict(set)        # sig -> set of H values
        by_cyc = defaultdict(set)        # (c3,c5) -> set of H
        sig_examples = defaultdict(list)
        count = 0
        for A in all_tournaments(n):
            count += 1
            sig = spectral_sig(A, n)
            H = ham_path_count(A, n)
            cyc = cycle_counts(A, n)
            by_sig[sig].add(H)
            by_cyc[cyc].add(H)
            if len(sig_examples[sig]) < 2 and (not sig_examples[sig] or sig_examples[sig][0][1] != H):
                sig_examples[sig].append((A.copy(), H))
        split_sig = {s: hs for s, hs in by_sig.items() if len(hs) > 1}
        split_cyc = {c: hs for c, hs in by_cyc.items() if len(hs) > 1}
        print(f"=== n={n} ({count} tournaments) ===", flush=True)
        print(f"   distinct spectral classes: {len(by_sig)}; classes with >1 H value (cospectral, diff H): {len(split_sig)}", flush=True)
        if split_sig:
            print(f"   => H is NOT spectrally determined at n={n}. Example H-splits within a spectral class:", flush=True)
            shown = 0
            for s, hs in sorted(split_sig.items(), key=lambda x: -len(x[1]))[:4]:
                print(f"      sig(tr A^1..A^{n})={s}: H values {sorted(hs)}", flush=True)
                shown += 1
        else:
            print(f"   => every spectral class has a single H (H IS spectrally determined at n={n}).", flush=True)
        print(f"   distinct (c3,c5) classes: {len(by_cyc)}; with >1 H: {len(split_cyc)} "
              f"(H needs MORE than (c3,c5) -- the disjointness alpha_2 -- iff >0)", flush=True)
        # forbidden-H reframe check: are 7,21 the only missing H values, and are they spectral?
        allH = sorted(set().union(*by_sig.values()))
        print(f"   realized H values: {allH[:25]}{'...' if len(allH)>25 else ''}", flush=True)
        print(flush=True)

    print("VERDICT: cycle counts c3,c5 ARE spectral (=tr/k); if H splits within spectral", flush=True)
    print("classes, then H's forbidden values {7,21} are NOT power-sum exclusions -- they", flush=True)
    print("live in the conflict-graph disjointness layer (alpha_2), a DIFFERENT kind of", flush=True)
    print("obstruction than the c5=10 spectral gap. That boundary is the refined reframe.", flush=True)


if __name__ == "__main__":
    main()
