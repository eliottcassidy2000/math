#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Efficiency becomes proof: the trace formula c_k = tr(A^k)/k for k=3,4,5 in
tournaments, the O(n^3) speedup it gives, and the spectral reframing of the
c5=10 gap (THM-498).  kind-pasteur-2026-06-13-S4.

THE SPEEDUP-THEOREM.  A tournament is 2-cycle-free (A[i][j]A[j][i]=0).  A closed
directed k-walk with a repeated vertex decomposes into closed sub-walks whose
lengths sum to k; the smallest closed sub-walk in a 2-cycle-free digraph has
length 3 (no loops, no 2-cycles).  For k in {3,4,5} there is NO way to write k as
a sum of parts all >= 3 with >= 2 parts (3,4,5 < 3+3=6).  Hence EVERY closed
k-walk (k<=5) is a genuine k-cycle, counted k times (rotations), so

        c_k(T) = trace(A^k) / k     for k = 3, 4, 5   (any tournament).

This FAILS first at k=6 (figure-8 of two triangles at a shared vertex: 6=3+3).
So c5 = tr(A^5)/5 EXACTLY -- an O(n^3) computation, vs the O(C(n,5)*5!) explicit
cycle enumeration used in THM-498.  The speedup IS a theorem; its correctness
proof is the cycle-decomposition identity.

USES:
 (A) VERIFY the formula against the explicit counter (the THM-498 method), n<=6.
 (B) the speedup makes the n=7 c5 spectrum (2^21 tournaments) FEASIBLE -- compute
     it exhaustively (was infeasible last session with the O(n^5) counter), extend
     THM-498: does c5=10 persist as a gap? what is the n=7 forbidden set?
 (C) SPECTRAL REFRAMING of the gap: c5=10 <=> tr(A^5)=50.  Since A+A^T=J-I, work
     out which tr(A^5) values are achievable and why 50 is excluded at n=6.
"""

import sys, itertools, time
import numpy as np
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def c5_explicit(adj, n):
    """THM-498's method: per-5-subset directed-Ham-cycle enumeration."""
    c = 0
    for S in itertools.combinations(range(n), 5):
        s0 = S[0]; rest = S[1:]
        for perm in itertools.permutations(rest):
            if adj[s0][perm[0]] and all(adj[perm[x]][perm[x+1]] for x in range(3)) and adj[perm[3]][s0]:
                c += 1
    return c


def trace_A5(A):
    A2 = A @ A
    A4 = A2 @ A2
    A5 = A4 @ A
    return int(np.trace(A5))


def edges_of(bits, pairs, n):
    A = np.zeros((n, n), dtype=np.int64)
    adj = [[False]*n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i][j] = 1; adj[i][j] = True
        else:
            A[j][i] = 1; adj[j][i] = True
    return A, adj


def part_A_verify():
    print("=== (A) verify c5 = tr(A^5)/5 vs explicit counter (n=4,5,6) ===", flush=True)
    for n in (4, 5, 6):
        pairs = list(itertools.combinations(range(n), 2))
        mism = 0; total = 0
        rng = np.random.default_rng(0)
        # exhaustive for n<=5, sample for n=6 (cross-check)
        space = range(1 << len(pairs)) if n <= 5 else (int(rng.integers(0, 1<<len(pairs))) for _ in range(4000))
        for bits in space:
            A, adj = edges_of(bits, pairs, n)
            t = trace_A5(A)
            if t % 5 != 0:
                mism += 1
            c_fast = t // 5
            c_slow = c5_explicit(adj, n)
            if c_fast != c_slow:
                mism += 1
            total += 1
        print(f"   n={n}: {total} tournaments, tr(A^5) divisible by 5 & matches explicit: "
              f"{'ALL OK' if mism==0 else str(mism)+' MISMATCH'}", flush=True)


def c5_spectrum_exhaustive(n, tlimit=480):
    """exhaustive c5 spectrum via tr(A^5)/5; returns set, or None if timed out."""
    pairs = list(itertools.combinations(range(n), 2))
    m = len(pairs)
    vals = set()
    t0 = time.time()
    # precompute pair index -> (i,j)
    for bits in range(1 << m):
        if (bits & 0xFFFF) == 0 and time.time() - t0 > tlimit:
            return None, bits  # timeout, partial
        A = np.zeros((n, n), dtype=np.int64)
        for k in range(m):
            i, j = pairs[k]
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        A2 = A @ A; A4 = A2 @ A2
        vals.add(int(np.trace(A4 @ A)) // 5)
    return vals, None


def gaps(vals):
    vmax = max(vals)
    return [v for v in range(vmax + 1) if v not in vals]


def main():
    t0 = time.time()
    part_A_verify()

    print("\n=== (B) c5 spectrum via the speedup ===", flush=True)
    for n in (5, 6):
        vals, _ = c5_spectrum_exhaustive(n)
        g = gaps(vals)
        print(f"   n={n} (exhaustive): c5 in [0,{max(vals)}], gaps={g if g else 'NONE'}", flush=True)

    print("   n=7 (2^21 tournaments) — feasible ONLY via the O(n^3) trace speedup:", flush=True)
    vals7, partial = c5_spectrum_exhaustive(7, tlimit=460)
    if vals7 is None:
        print(f"      timed out (reached {partial}/{1<<21}); falling back to sampling", flush=True)
        pairs = list(itertools.combinations(range(7), 2)); m = len(pairs)
        rng = np.random.default_rng(7); vals7 = set()
        for _ in range(600000):
            bits = int(rng.integers(0, 1 << m))
            A = np.zeros((7,7), dtype=np.int64)
            for k in range(m):
                i,j = pairs[k]
                if (bits>>k)&1: A[i,j]=1
                else: A[j,i]=1
            A2=A@A; A4=A2@A2; vals7.add(int(np.trace(A4@A))//5)
        g7 = gaps(vals7)
        print(f"      SAMPLED c5 in [0,{max(vals7)}], missing-in-range={g7 if len(g7)<20 else str(g7[:20])+'...'}", flush=True)
    else:
        g7 = gaps(vals7)
        print(f"      EXHAUSTIVE c5 in [0,{max(vals7)}], #values={len(vals7)}, "
              f"GAPS={g7 if g7 else 'NONE (gap-free at n=7!)'}", flush=True)

    print("\n=== (C) spectral reframing of the c5=10 gap at n=6 ===", flush=True)
    print("   c5=10 <=> tr(A^5)=50. Achievable tr(A^5) values at n=6:", flush=True)
    vals6, _ = c5_spectrum_exhaustive(6)
    tr_vals = sorted(5*v for v in vals6)
    print(f"      tr(A^5) achievable = {tr_vals}", flush=True)
    print(f"      50 in achievable tr(A^5)? {'YES' if 50 in tr_vals else 'NO — the gap is tr(A^5) skips 50'}", flush=True)
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
