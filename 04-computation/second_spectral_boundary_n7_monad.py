#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THE SECOND SPECTRAL BOUNDARY (monad-explorer-2026-06-13).

Builds directly on THM-499 (kind-pasteur-S5):
  * At n<=6: H = 1 + 2(c3+c5) + 4D, D=alpha_2 = #vertex-disjoint triangle pairs.
  * c3,c5 are SPECTRAL (=tr A^k / k, THM-118); D=alpha_2 is the FIRST non-spectral
    OCF ingredient; the spectral-reframe boundary is the onset of alpha_2 at n=6.

THM-499 showed H splits within cospectral classes at n=6 because alpha_2 does.
But the ODD-CYCLE COUNT alpha_1 = c3+c5 is STILL spectral at n=6.

QUESTION (the second boundary): at n=7 the odd cycles are 3-,5-,7-cycles, and
  alpha_1 = c3 + c5 + c7,   c7 = #directed Hamiltonian (7-)cycles.
c3,c5 are spectral, but is c7? tr(A^7) is spectral, but c7 = tr(A^7)/7 ONLY if no
shorter closed walk pads to length 7. The first padding is a triangle + a 4-cycle
sharing a vertex (3+4=7). So c7 may be NON-spectral -- in which case alpha_1 itself
becomes non-spectral from n=7: the spectrum can no longer even COUNT the odd cycles.

Plan:
 (A) Verify the natural n=7 extension of the seed formula on a sample:
       H = 1 + 2(c3+c5+c7) + 4*DTP,  DTP = #vertex-disjoint triangle pairs.
     (At n=7 the only disjoint odd-cycle pair is triangle+triangle (6 vtx);
      alpha_3 needs 9 vtx => 0. So this is the exact OCF truncation at n=7.)
 (B) Bucket sampled tournaments by EXACT integer spectral signature
       sig = (tr A^1, ..., tr A^7)  <=>  char poly  <=>  spectrum.
     Within each cospectral bucket, measure the spread of c7, alpha_1, DTP, H.
     A bucket with >1 c7 value is an explicit witness that c7 (hence alpha_1) is
     NOT spectrally determined -- the SECOND boundary.
 (C) Sanity: c3,c5 must be constant within every bucket (they are tr/k).
"""

import sys, itertools, random
from collections import defaultdict
import numpy as np

random.seed(20260613)


def random_tournament(n):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if random.getrandbits(1):
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A


def spectral_sig(A, n):
    sig = []
    P = np.eye(n, dtype=np.int64)
    for k in range(1, n + 1):
        P = P @ A
        sig.append(int(np.trace(P)))
    return tuple(sig)


def out_adj(A, n):
    return [[j for j in range(n) if A[i, j]] for i in range(n)]


def count_ham_cycles(adj, n):
    """Count directed Hamiltonian cycles. Fix start=0 to avoid n-fold overcount;
    each directed cycle is counted once (one rotation pinned at 0)."""
    cnt = 0
    start = 0
    visited = 1 << start
    stack = [(start, visited, 1)]
    # iterative DFS
    def dfs(v, vis, depth):
        nonlocal cnt
        if depth == n:
            if 0 in adj[v]:   # close the cycle v -> start(0)
                cnt += 1
            return
        for w in adj[v]:
            b = 1 << w
            if not (vis & b):
                dfs(w, vis | b, depth + 1)
    dfs(start, visited, 1)
    return cnt


def count_triangles(adj, A, n):
    """#directed 3-cycles."""
    return int(np.trace(A @ A @ A)) // 3


def count_5cycles_via_trace(A, n):
    A2 = A @ A
    A4 = A2 @ A2
    return int(np.trace(A4 @ A)) // 5


def disjoint_triangle_pairs(A, n):
    """#unordered pairs of vertex-disjoint directed triangles."""
    tris = []
    for i in range(n):
        for j in range(n):
            if not A[i, j]:
                continue
            for k in range(n):
                if A[j, k] and A[k, i] and i < j and i < k:
                    # canonical: i is min, ensure directed i->j->k->i; count each triangle once
                    if i < j and i < k:
                        tris.append(frozenset((i, j, k)))
    tris = list(set(tris))
    cnt = 0
    for a in range(len(tris)):
        for b in range(a + 1, len(tris)):
            if tris[a].isdisjoint(tris[b]):
                cnt += 1
    return cnt, len(tris)


def ham_path_count(A, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            row = A[v]
            for w in range(n):
                if mask & (1 << w):
                    continue
                if row[w]:
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))


def main():
    n = 7
    NSAMPLE = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    NFORMULA = 4000     # exhaustive-formula check subsample (also computes H)

    # ---- (A) formula check: H = 1 + 2(c3+c5+c7) + 4*DTP at n=7 ----
    print(f"=== (A) n=7 OCF formula extension check on {NFORMULA} random tournaments ===", flush=True)
    bad = 0
    for _ in range(NFORMULA):
        A = random_tournament(n)
        adj = out_adj(A, n)
        c3 = count_triangles(adj, A, n)
        c5 = count_5cycles_via_trace(A, n)
        c7 = count_ham_cycles(adj, n)
        dtp, ntri = disjoint_triangle_pairs(A, n)
        H = ham_path_count(A, n)
        pred = 1 + 2 * (c3 + c5 + c7) + 4 * dtp
        if pred != H:
            bad += 1
            if bad <= 5:
                print(f"   MISMATCH H={H} pred={pred} c3={c3} c5={c5} c7={c7} dtp={dtp}", flush=True)
    if bad == 0:
        print(f"   VERIFIED: H = 1 + 2(c3+c5+c7) + 4*DTP for all {NFORMULA}/{NFORMULA} (0 exceptions).", flush=True)
        print(f"   => the seed's n<=6 formula H=1+2(c3+c5)+4D extends to n=7 with c7 and DTP=triangle-pairs.", flush=True)
    else:
        print(f"   {bad}/{NFORMULA} MISMATCHES -- formula does NOT extend naively; investigate.", flush=True)
    print(flush=True)

    # ---- (B) cospectral bucketing: does c7 split within a spectral class? ----
    print(f"=== (B) bucket {NSAMPLE} random n=7 tournaments by exact spectral signature ===", flush=True)
    by_sig_c7 = defaultdict(set)
    by_sig_c3c5 = defaultdict(set)
    by_sig_dtp = defaultdict(set)
    by_sig_a1 = defaultdict(set)
    witness_examples = {}   # sig -> dict c7 -> example (sig vector kept, plus invariants)
    for s in range(NSAMPLE):
        A = random_tournament(n)
        sig = spectral_sig(A, n)
        adj = out_adj(A, n)
        c3 = sig[2] // 3
        c5 = sig[4] // 5
        c7 = count_ham_cycles(adj, n)
        dtp, _ = disjoint_triangle_pairs(A, n)
        a1 = c3 + c5 + c7
        by_sig_c7[sig].add(c7)
        by_sig_c3c5[sig].add((c3, c5))
        by_sig_dtp[sig].add(dtp)
        by_sig_a1[sig].add(a1)
        if sig not in witness_examples:
            witness_examples[sig] = {}
        if c7 not in witness_examples[sig]:
            witness_examples[sig][c7] = dict(c3=c3, c5=c5, c7=c7, dtp=dtp, a1=a1,
                                             H=1 + 2 * a1 + 4 * dtp, A=A.copy())

    nbuckets = len(by_sig_c7)
    multi = {s: v for s, v in by_sig_c7.items() if len(v) > 1}
    c3c5_split = {s: v for s, v in by_sig_c3c5.items() if len(v) > 1}
    a1_split = {s: v for s, v in by_sig_a1.items() if len(v) > 1}
    dtp_split = {s: v for s, v in by_sig_dtp.items() if len(v) > 1}
    print(f"   distinct spectral classes hit: {nbuckets}", flush=True)
    print(f"   classes with >1 c7 value (cospectral, DIFFERENT Hamiltonian-cycle count): {len(multi)}", flush=True)
    print(f"   classes with >1 (c3,c5)  [MUST be 0 -- spectral sanity]: {len(c3c5_split)}", flush=True)
    print(f"   classes with >1 alpha_1=c3+c5+c7 (odd-cycle count NON-spectral?): {len(a1_split)}", flush=True)
    print(f"   classes with >1 DTP (disjoint triangle pairs, non-spectral): {len(dtp_split)}", flush=True)
    print(flush=True)

    if multi:
        print("   *** SECOND BOUNDARY CONFIRMED: c7 (Hamiltonian-cycle count) is NOT spectral at n=7. ***", flush=True)
        # show the witness with the largest c7 spread
        best = max(multi.items(), key=lambda kv: max(kv[1]) - min(kv[1]))
        sig, c7set = best
        print(f"   widest c7-split spectral class: sig={sig}", flush=True)
        print(f"     c7 values realized: {sorted(c7set)}", flush=True)
        for c7v in sorted(witness_examples[sig]):
            ex = witness_examples[sig][c7v]
            print(f"     -> c7={ex['c7']}: c3={ex['c3']} c5={ex['c5']} alpha_1={ex['a1']} DTP={ex['dtp']} H={ex['H']}", flush=True)
        # print explicit adjacency of two split witnesses
        c7vals = sorted(witness_examples[sig])
        if len(c7vals) >= 2:
            print("   explicit cospectral witnesses (upper-triangle arc orientations):", flush=True)
            for c7v in c7vals[:2]:
                A = witness_examples[sig][c7v]['A']
                arcs = [(i, j) if A[i, j] else (j, i) for i in range(n) for j in range(i + 1, n)]
                print(f"     c7={c7v}: arcs={arcs}", flush=True)
    else:
        print("   No c7 split found in this sample -- c7 may be spectral at n=7, OR sample too small.", flush=True)
        print("   (Increase NSAMPLE; cospectral-but-non-isomorphic pairs are rare under uniform sampling.)", flush=True)

    print(flush=True)
    print("INTERPRETATION:", flush=True)
    print(" - n=6 boundary (THM-499): H non-spectral via alpha_2 (disjoint odd-cycle PAIR).", flush=True)
    print("   But alpha_1 = c3+c5 is still spectral there.", flush=True)
    print(" - n=7 SECOND boundary: if c7 splits, then alpha_1 = c3+c5+c7 is non-spectral too.", flush=True)
    print("   The spectrum then cannot even COUNT the odd cycles. The two non-spectral", flush=True)
    print("   sources at n=7 are c7 (a single 7-cycle that revisiting walks of length 7", flush=True)
    print("   from triangle+4-cycle cannot be untangled spectrally) and DTP (disjointness).", flush=True)


if __name__ == "__main__":
    main()
