#!/usr/bin/env python3
"""
score_perspective_s113.py — The Score Sequence Perspective (n ≤ 7)
opus-2026-03-21-S113

Two views of the same data:

VIEW 1 (S112): Group by ISOMORPHISM CLASS.
  n=6: 56 classes, 22 score classes, gap = 34 (61%).
  This measures: how many structurally distinct tournaments exist?

VIEW 2 (this session): Group by SCORE SEQUENCE.
  For each score, count: #labeled tournaments, #iso classes within,
  H distribution, OCR contribution.
  This measures: how much does each score contribute to the shadow?

At n=7: 2^21 = 2,097,152 labeled tournaments, 59 score sequences.
Full iso classification is too expensive, but score-based analysis is fast.

THE SCORE PERSPECTIVE answers:
  - Which scores have the most labeled tournaments? (= most "common")
  - Which scores have the most H values? (= most "ambiguous")
  - Which scores contribute most to Var(H)? (= most "important")
  - Which scores contribute most to Var(H|scores)? (= the PoS classes)
"""

from math import comb, factorial
from fractions import Fraction
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

print("=" * 76)
print("  THE SCORE SEQUENCE PERSPECTIVE (n ≤ 7)")
print("  opus-2026-03-21-S113")
print("=" * 76)
print()

# =============================================================================
# HELPER: Fast H computation
# =============================================================================

def H_dp_flat(adj_flat, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n): dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp[S * n + v]
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj_flat[v * n + u]:
                    dp[(S | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

# =============================================================================
# SCORE-BASED ANALYSIS FOR n=3..7
# =============================================================================

for n in range(3, 8):
    m = comb(n, 2)
    N = 1 << m

    print(f"{'='*76}")
    print(f"  n = {n}: {N:,} labeled tournaments, C(n,2) = {m} arcs")
    print(f"{'='*76}")
    print()

    t0 = time.time()

    # For each score class: collect H values, count tournaments
    score_data = defaultdict(lambda: {'count': 0, 'H_sum': 0, 'H2_sum': 0,
                                       'H_vals': Counter(), 'H_min': 10**9, 'H_max': 0})

    for mask in range(N):
        scores = [0] * n
        adj_flat = [0] * (n * n)
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (mask >> idx) & 1:
                    adj_flat[i * n + j] = 1
                    scores[i] += 1
                else:
                    adj_flat[j * n + i] = 1
                    scores[j] += 1
                idx += 1

        sc = tuple(sorted(scores))
        H = H_dp_flat(adj_flat, n)

        d = score_data[sc]
        d['count'] += 1
        d['H_sum'] += H
        d['H2_sum'] += H * H
        d['H_vals'][H] += 1
        d['H_min'] = min(d['H_min'], H)
        d['H_max'] = max(d['H_max'], H)

    elapsed = time.time() - t0
    n_scores = len(score_data)

    # Global statistics
    total_H = sum(d['H_sum'] for d in score_data.values())
    total_H2 = sum(d['H2_sum'] for d in score_data.values())
    E_H = Fraction(total_H, N)
    E_H2 = Fraction(total_H2, N)
    Var_H = E_H2 - E_H * E_H

    # Within-class variance
    Var_within = Fraction(0)
    for sc, d in score_data.items():
        ns = d['count']
        E_H_s = Fraction(d['H_sum'], ns)
        E_H2_s = Fraction(d['H2_sum'], ns)
        Var_s = E_H2_s - E_H_s * E_H_s
        Var_within += Fraction(ns, N) * Var_s

    OCR = 1 - Var_within / Var_H if Var_H > 0 else Fraction(1)

    print(f"  Computed in {elapsed:.1f}s")
    print(f"  Score sequences: {n_scores}")
    print(f"  E[H] = {E_H} = {float(E_H):.4f}")
    print(f"  Var(H) = {Var_H} = {float(Var_H):.4f}")
    print(f"  Var(H|scores) = {Var_within} = {float(Var_within):.6f}")
    print(f"  OCR = {OCR} = {float(OCR):.10f}")
    print()

    # Classify scores
    S2_vals = {}
    for sc in score_data:
        S2_vals[sc] = sum((s - (n - 1) / 2) ** 2 for s in sc)

    # Self-complementary check
    SC_flags = {}
    for sc in score_data:
        comp = tuple(sorted(n - 1 - s for s in sc))
        SC_flags[sc] = (sc == comp)

    # Sort by S₂ (ascending = most regular first)
    sorted_scores = sorted(score_data.keys(), key=lambda sc: S2_vals[sc])

    # SCORE PERSPECTIVE TABLE
    print(f"  THE SCORE PERSPECTIVE TABLE:")
    print()
    header = (f"  {'Score':>25} | {'Count':>7} | {'%':>5} | {'S₂':>5} | {'SC':>2} | "
              f"{'#H vals':>6} | {'H range':>12} | {'E[H|s]':>8} | {'Var(H|s)':>10} | "
              f"{'Contrib%':>8}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for sc in sorted_scores:
        d = score_data[sc]
        ns = d['count']
        pct = ns / N * 100
        S2 = S2_vals[sc]
        is_sc = SC_flags[sc]
        n_H = len(d['H_vals'])
        H_range = f"[{d['H_min']},{d['H_max']}]"
        E_H_s = d['H_sum'] / ns
        Var_s = d['H2_sum'] / ns - E_H_s ** 2
        contrib = (ns / N * Var_s) / float(Var_within) * 100 if float(Var_within) > 0 else 0

        # Only show full detail for ambiguous classes or if n<=6
        if n <= 6 or n_H > 1 or S2 <= 4:
            print(f"  {str(sc):>25} | {ns:7,} | {pct:4.1f}% | {S2:5.1f} | {'✓' if is_sc else ' ':>2} | "
                  f"{n_H:6d} | {H_range:>12} | {E_H_s:8.2f} | {Var_s:10.4f} | "
                  f"{contrib:7.1f}%{'★' if contrib > 5 else ''}")

    if n >= 7:
        # Show summary of omitted rows
        n_shown = sum(1 for sc in sorted_scores if len(score_data[sc]['H_vals']) > 1 or S2_vals[sc] <= 4)
        n_omitted = n_scores - n_shown
        print(f"  ... ({n_omitted} unambiguous high-S₂ classes omitted)")

    print()

    # ==========================================================================
    # TOP CONTRIBUTORS TO THE RESIDUAL
    # ==========================================================================

    if float(Var_within) > 0:
        print(f"  TOP CONTRIBUTORS TO Var(H|scores) (the OCR residual):")
        print()

        # Compute contribution of each score class
        contributions = []
        for sc in sorted_scores:
            d = score_data[sc]
            ns = d['count']
            E_H_s = d['H_sum'] / ns
            Var_s = d['H2_sum'] / ns - E_H_s ** 2
            if Var_s > 0:
                weight = ns / N
                contrib_abs = weight * Var_s
                contributions.append((sc, ns, Var_s, contrib_abs))

        contributions.sort(key=lambda x: -x[3])
        cumulative = 0

        print(f"  {'Rank':>4} | {'Score':>25} | {'Count':>7} | {'Var(H|s)':>10} | {'Contrib':>10} | {'Cum%':>6}")
        print(f"  " + "-" * 75)

        for rank, (sc, ns, var_s, contrib) in enumerate(contributions[:15], 1):
            cumulative += contrib
            cum_pct = cumulative / float(Var_within) * 100
            print(f"  {rank:4d} | {str(sc):>25} | {ns:7,} | {var_s:10.4f} | {contrib:10.6f} | {cum_pct:5.1f}%")

        print()
        print(f"  Top 1 class: {contributions[0][3]/float(Var_within)*100:.1f}% of residual")
        print(f"  Top 3 classes: {sum(c[3] for c in contributions[:3])/float(Var_within)*100:.1f}% of residual")
        if len(contributions) >= 5:
            print(f"  Top 5 classes: {sum(c[3] for c in contributions[:5])/float(Var_within)*100:.1f}% of residual")

    print()

    # ==========================================================================
    # SIZE DISTRIBUTION
    # ==========================================================================

    print(f"  SCORE CLASS SIZE DISTRIBUTION:")
    sizes = sorted([(d['count'], sc) for sc, d in score_data.items()], reverse=True)
    print(f"    Largest class:  {sizes[0][1]} with {sizes[0][0]:,} tournaments ({sizes[0][0]/N*100:.1f}%)")
    print(f"    Smallest class: {sizes[-1][1]} with {sizes[-1][0]:,} tournaments ({sizes[-1][0]/N*100:.2f}%)")
    print(f"    Median class size: {sizes[len(sizes)//2][0]:,}")
    print()

    # Self-complementary class statistics
    sc_classes = [sc for sc in score_data if SC_flags[sc]]
    non_sc_classes = [sc for sc in score_data if not SC_flags[sc]]
    sc_total = sum(score_data[sc]['count'] for sc in sc_classes)
    non_sc_total = sum(score_data[sc]['count'] for sc in non_sc_classes)
    print(f"  SELF-COMPLEMENTARY ANALYSIS:")
    print(f"    SC score classes: {len(sc_classes)} ({len(sc_classes)/n_scores*100:.0f}% of classes)")
    print(f"    SC tournaments: {sc_total:,} ({sc_total/N*100:.1f}% of tournaments)")
    if sc_classes:
        avg_sc = sc_total / len(sc_classes)
    else:
        avg_sc = 0
    if non_sc_classes:
        avg_non_sc = non_sc_total / len(non_sc_classes)
    else:
        avg_non_sc = 0
    if avg_non_sc > 0:
        print(f"    Avg SC class size: {avg_sc:,.0f} vs non-SC: {avg_non_sc:,.0f} (ratio: {avg_sc/avg_non_sc:.2f}×)")
    print()

# =============================================================================
# SYNTHESIS
# =============================================================================

print("=" * 76)
print("  SYNTHESIS: THE SCORE SEQUENCE PERSPECTIVE")
print("=" * 76)
print()

print("""
  THE TWO PERSPECTIVES:

  ISO CLASS PERSPECTIVE (S112):
    "How many structurally distinct tournaments exist?"
    Gap sequence: 0, 0, 3, 34, ? — grows exponentially
    At n=6: 61% of iso classes are "hidden" from scores

  SCORE SEQUENCE PERSPECTIVE (this session):
    "How much of the observable (H) is explained by scores?"
    OCR sequence: 1, 1, 0.970, 0.959, ? — decreases but slowly
    At n=6: only 4% of H VARIANCE is "hidden" from scores

  WHY THE DISCREPANCY:
    The gap counts UNWEIGHTED iso classes.
    The OCR counts WEIGHTED (by labeled tournament count).

    Most labeled tournaments are in LOW-S₂ classes where scores suffice.
    Most iso classes are in HIGH-S₂ classes where scores fail.

    Think of it like income distribution:
    Most INDIVIDUALS are middle-class (well-determined by score).
    Most CATEGORIES are at extremes (poorly determined by score).

  THE REPRESENTATION THEORY INTERPRETATION:
    The score sequence is the CHARACTER of the "score representation."
    The OCR measures the TRACE of this character over the regular representation.
    The gap measures the DIMENSION of the kernel of the character.

    The character has high trace (OCR ~96%) but large kernel (gap ~61%).
    This is typical of "nearly faithful" representations:
    they capture almost all the weight but miss many small subrepresentations.
""")

print("opus-2026-03-21-S113: THE SCORE SEQUENCE PERSPECTIVE")
