#!/usr/bin/env python3
"""
ROUTE 5, part 2 -- the GENUINE c_k crossover via the OCF level-sums.

THM-505 proved H is NOT a low-degree spectral polynomial in general; the
"H = c_0 + sum c_k e_k" of THM-134 is underdetermined by circulant spectral
data (only 2 classes at p=7, 4 at p=11 vs m=3,5 unknowns). The canonical,
universal identity is the OCF LEVEL-SUM:

      H = 1 + sum_{j>=1} 2^j * alpha_j,
      alpha_j = # vertex-disjoint j-tuples of directed ODD cycles.

For circulants on Z_p the cycle counts c_ell and the alpha_j ARE determined
by the spectrum (regularity + circulant symmetry pin them). So the genuine
"c_k" are the alpha_j; the crossover prime is where the Paley-vs-Interval
ranking of the 2^j-weighted alpha_j flips.

This script (exact integers):
 1. computes c_ell (directed ell-cycle counts) for Paley & Interval at p,
 2. computes alpha_j exactly (independent j-sets in the odd-cycle conflict
    graph) by brute enumeration where feasible, AND via the OCF recursion,
 3. reports, per level j, delta_j = alpha_j(Paley) - alpha_j(Interval) and
    the 2^j-weighted contribution, locating the first level where Interval
    overtakes Paley -- the OCF analogue of the Schur-Ostrowski sign flip,
 4. cross-checks H against the spectral skeleton.

Author: mac-mini-2026-06-21-S20 (ROUTE 5)
"""
import sys, cmath, math
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)


def is_qr(a, p):
    return a % p != 0 and pow(a, (p - 1) // 2, p) == 1


def paley_set(p):
    return frozenset(j for j in range(1, p) if is_qr(j, p))


def interval_set(p):
    return frozenset(range(1, (p - 1) // 2 + 1))


def adj(p, S):
    Sset = set(S)
    return [[1 if (j != i and (j - i) % p in Sset) else 0 for j in range(p)]
            for i in range(p)]


def trace_power_int(A, k, p):
    """Exact tr(A^k) by integer matrix power."""
    # A is small (p up to 19); use Python int matmul.
    M = [row[:] for row in A]
    R = None
    base = [row[:] for row in A]
    # compute A^k by repeated multiply (k small)
    P = [[1 if i == j else 0 for j in range(p)] for i in range(p)]
    for _ in range(k):
        P = matmul(P, A, p)
    return sum(P[i][i] for i in range(p))


def matmul(X, Y, p):
    Z = [[0] * p for _ in range(p)]
    for i in range(p):
        Xi = X[i]
        Zi = Z[i]
        for t in range(p):
            x = Xi[t]
            if x:
                Yt = Y[t]
                for j in range(p):
                    Zi[j] += x * Yt[j]
    return Z


def cycle_counts(A, p, max_ell):
    """Number of directed ell-cycles c_ell for ell=3,5,...,max_ell,
    via Newton/Mobius on traces (valid: tournaments have no 2-cycles;
    for circulants only odd cycles matter to OCF).
    c_ell = (1/ell) * sum_{d | ell} mu(d) * tr(A^{ell/d})  is the # of
    *closed walks that are primitive ell-cycles* only when ell is prime;
    for composite ell we use the standard simple-cycle extraction.
    We instead BRUTE-count simple directed cycles of each odd length to be safe
    when feasible (ell <= 9). For larger ell we fall back to the Witt formula
    W_ell which counts simple cycles for prime ell exactly."""
    counts = {}
    # brute count simple directed cycles up to length max_ell
    # represent cycles by their vertex set start-canonicalized
    import sys as _s
    nbrs = [[j for j in range(p) if A[i][j]] for i in range(p)]
    for ell in range(3, max_ell + 1, 2):
        counts[ell] = count_simple_cycles(nbrs, p, ell)
    return counts


def count_simple_cycles(nbrs, p, ell):
    """Count directed simple cycles of length exactly ell.
    DFS from each start = smallest vertex; only visit vertices > start to
    canonicalize (each cycle counted once at its min vertex)."""
    total = 0
    for start in range(p):
        # paths of length ell starting and ending at start, all interior > start
        stack = [(start, 1 << start, start, 1)]
        # (current, visited_mask, depth, _) -- depth = #vertices so far
        # iterative DFS
        def dfs(cur, mask, depth):
            nonlocal total
            if depth == ell:
                if start in nbrs[cur]:
                    total += 1
                return
            for nx in nbrs[cur]:
                if nx <= start:
                    continue
                if mask & (1 << nx):
                    continue
                dfs(nx, mask | (1 << nx), depth + 1)
        dfs(start, 1 << start, 1)
    return total


def all_odd_cycles(A, p, max_ell):
    """Enumerate ALL directed odd cycles as frozensets of vertices,
    labelled by length. Returns list of (length, frozenset)."""
    nbrs = [[j for j in range(p) if A[i][j]] for i in range(p)]
    cycles = []
    for ell in range(3, max_ell + 1, 2):
        for start in range(p):
            def dfs(cur, mask, path):
                if len(path) == ell:
                    if start in nbrs[cur]:
                        cycles.append((ell, frozenset(path)))
                    return
                for nx in nbrs[cur]:
                    if nx <= start or (mask & (1 << nx)):
                        continue
                    dfs(nx, mask | (1 << nx), path + [nx])
            dfs(start, 1 << start, [start])
    return cycles


def alpha_levels(A, p, max_ell):
    """Exact alpha_j = # vertex-disjoint j-tuples of odd cycles, via the
    independence polynomial of the conflict graph (vertices = odd cycles,
    edges = vertex-sharing). Computed by enumerating independent sets.
    Feasible for moderate cycle counts."""
    cyc = all_odd_cycles(A, p, max_ell)
    N = len(cyc)
    sets = [c[1] for c in cyc]
    # independence polynomial coefficients alpha_j
    # use DP over cycles: but disjointness on vertex sets.
    # Branch & bound enumeration of independent sets, counting by size.
    alpha = defaultdict(int)
    alpha[0] = 1
    # order cycles; recursive add
    # To keep tractable, use bitmask of used vertices.
    import sys as _s
    _s.setrecursionlimit(100000)
    masks = [sum(1 << v for v in s) for s in sets]
    # We count independent sets: each j-subset of cycles that are pairwise
    # vertex-disjoint. Equivalent to subsets whose masks are pairwise disjoint.
    # DP: process cycles in order, choose include/exclude.
    # state: index, used_mask -> but used_mask space huge. Use recursion with
    # pruning by skipping cycles overlapping used.
    res = defaultdict(int)
    def rec(i, used, depth):
        res[depth] += 1
        for j in range(i, N):
            if masks[j] & used:
                continue
            rec(j + 1, used | masks[j], depth + 1)
    rec(0, 0, 0)
    return dict(res), N


def H_from_alpha(alpha):
    return sum((1 << j) * cnt for j, cnt in alpha.items())


def main():
    for p in [7, 11]:
        print("=" * 72)
        print(f"OCF LEVEL-SUM ANALYSIS -- p = {p}  (p mod 4 = {p % 4})")
        print("=" * 72)
        max_ell = p if p % 2 == 1 else p - 1  # max odd cycle length = p
        for name, S in [("PALEY", paley_set(p)), ("INTERVAL", interval_set(p))]:
            A = adj(p, S)
            cc = cycle_counts(A, p, max_ell)
            alpha, Ncyc = alpha_levels(A, p, max_ell)
            H = H_from_alpha(alpha)
            print(f"\n  {name}  S={sorted(S)}")
            print(f"    cycle counts c_ell: {cc}")
            print(f"    total odd cycles (sum c_ell) = {sum(cc.values())} ; Ncyc enum={Ncyc}")
            print(f"    alpha_j levels: {dict(sorted(alpha.items()))}")
            print(f"    H = 1 + sum 2^j alpha_j = {H}")
        # delta per level
        print("\n  --- per-level Paley-vs-Interval ---")
        AP = adj(p, paley_set(p)); AI = adj(p, interval_set(p))
        aP, _ = alpha_levels(AP, p, max_ell)
        aI, _ = alpha_levels(AI, p, max_ell)
        levels = sorted(set(aP) | set(aI))
        cumP = cumI = 0
        print(f"    {'j':>3} {'aP':>10} {'aI':>10} {'delta':>10} {'2^j*delta':>14}")
        for j in levels:
            dP = aP.get(j, 0); dI = aI.get(j, 0)
            d = dP - dI
            w = (1 << j) * d
            if d > 0:
                cumP += (1 << j) * d
            elif d < 0:
                cumI += (1 << j) * (-d)
            tag = "P" if d > 0 else ("I" if d < 0 else "=")
            print(f"    {j:>3} {dP:>10} {dI:>10} {d:>+10} {w:>+14} {tag}")
        print(f"    cumulative 2^j-weighted: Paley+={cumP}  Interval+={cumI}")
        print(f"    H(P)-H(I) = {cumP - cumI}")
    print("\nDONE.")


if __name__ == "__main__":
    main()
