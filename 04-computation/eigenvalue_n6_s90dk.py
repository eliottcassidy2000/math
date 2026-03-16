#!/usr/bin/env python3
"""
eigenvalue_n6_s90dk.py — Transition matrix eigenvalues of the tournament flip chain at n=6
opus-2026-03-16-S90dk

Computes:
1. All 2^15 = 32768 tournaments on 6 vertices
2. Classifies into 56 isomorphism classes under S_6
3. Builds the 56x56 flip count matrix F[ci][cj]
4. Transition matrix T = D^{-1} F where D_ii = size_i * C(6,2)
5. Eigenvalues of T (exact fractions)
6. Tests conjecture: all eigenvalues are multiples of 1/5
7. Checks spectral gap = 2/5

Key conjecture: the denominator of all eigenvalues = largest prime <= n = 5.
"""

import numpy as np
from itertools import permutations
from math import comb, factorial
from fractions import Fraction
from collections import Counter
import time
import sys

N = 6
NUM_ARCS = comb(N, 2)  # 15
NUM_TOURNAMENTS = 1 << NUM_ARCS  # 32768

# Precompute arc index mapping: pair (i,j) with i<j -> index
ARC_INDEX = {}
idx = 0
for i in range(N):
    for j in range(i+1, N):
        ARC_INDEX[(i, j)] = idx
        idx += 1

def bits_to_adj(bits):
    """Convert a bitmask to adjacency matrix."""
    A = [[0]*N for _ in range(N)]
    idx = 0
    for i in range(N):
        for j in range(i+1, N):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def adj_to_bits(A):
    """Convert adjacency matrix back to bitmask."""
    bits = 0
    idx = 0
    for i in range(N):
        for j in range(i+1, N):
            if A[i][j] == 1:
                bits |= (1 << idx)
            idx += 1
    return bits

def canonical_bits(bits, perms_as_lists):
    """
    Compute canonical form of a tournament given as bitmask.
    Uses precomputed permutation action on arc indices.
    Returns the minimum bitmask over all relabelings.
    """
    best = bits
    for perm_action in perms_as_lists:
        new_bits = 0
        for arc_idx, (orig_val_mask, new_arc_idx, flipped) in enumerate(perm_action):
            # orig_val_mask = 1 << arc_idx (the bit to check in original)
            has_arc = (bits >> arc_idx) & 1
            if flipped:
                has_arc = 1 - has_arc
            if has_arc:
                new_bits |= (1 << new_arc_idx)
        if new_bits < best:
            best = new_bits
    return best

def precompute_perm_actions():
    """
    For each permutation p in S_6, precompute the action on arc indices.

    Arc (i,j) with i<j maps to (p[i], p[j]). If p[i] < p[j], it maps to
    arc index of (p[i], p[j]) with same orientation. If p[i] > p[j], it
    maps to arc index of (p[j], p[i]) with FLIPPED orientation.

    Returns list of actions. Each action is a list of tuples:
      (orig_arc_idx, new_arc_idx, flipped)
    where flipped means the orientation reverses.
    """
    actions = []
    for perm in permutations(range(N)):
        action = []
        idx = 0
        for i in range(N):
            for j in range(i+1, N):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    new_idx = ARC_INDEX[(pi, pj)]
                    flipped = False
                else:
                    new_idx = ARC_INDEX[(pj, pi)]
                    flipped = True
                action.append((idx, new_idx, flipped))
                idx += 1
        actions.append(action)
    return actions

def canonical_bits_fast(bits, perm_actions):
    """
    Compute canonical form using precomputed permutation actions.
    Returns min bitmask over all permutations.
    """
    best = bits
    for action in perm_actions:
        new_bits = 0
        for orig_idx, new_idx, flipped in action:
            has_arc = (bits >> orig_idx) & 1
            if flipped:
                has_arc = 1 - has_arc
            if has_arc:
                new_bits |= (1 << new_idx)
        if new_bits < best:
            best = new_bits
    return best

def flip_arc(bits, arc_idx):
    """Flip the arc at position arc_idx in the bitmask."""
    return bits ^ (1 << arc_idx)

def main():
    print("TOURNAMENT FLIP CHAIN — TRANSITION MATRIX EIGENVALUES AT n=6")
    print("=" * 70)
    print(f"n = {N}, arcs = {NUM_ARCS}, tournaments = {NUM_TOURNAMENTS}")
    print()

    # Step 1: Precompute permutation actions
    t0 = time.time()
    print("Step 1: Precomputing S_6 permutation actions on arc indices...")
    perm_actions = precompute_perm_actions()
    print(f"  {len(perm_actions)} permutations, {NUM_ARCS} arcs each")
    print(f"  Time: {time.time()-t0:.2f}s")
    print()

    # Step 2: Classify all tournaments into isomorphism classes
    t1 = time.time()
    print("Step 2: Classifying 32768 tournaments into iso classes...")
    canon_map = {}   # canonical_bits -> class_index
    classes = []     # list of canonical bitmasks
    t_class = [0] * NUM_TOURNAMENTS  # tournament -> class index

    for bits in range(NUM_TOURNAMENTS):
        cb = canonical_bits_fast(bits, perm_actions)
        if cb not in canon_map:
            canon_map[cb] = len(classes)
            classes.append(cb)
        t_class[bits] = canon_map[cb]

    num_classes = len(classes)
    print(f"  Found {num_classes} isomorphism classes")

    # Compute class sizes
    sizes = [0] * num_classes
    for ci in t_class:
        sizes[ci] += 1
    print(f"  Class sizes: min={min(sizes)}, max={max(sizes)}, sum={sum(sizes)}")
    assert sum(sizes) == NUM_TOURNAMENTS
    print(f"  Time: {time.time()-t1:.2f}s")
    print()

    # Step 3: Build flip count matrix F[ci][cj]
    # F[ci][cj] = number of (tournament T in class ci, arc e) pairs such that
    #             flipping e in T gives a tournament in class cj
    t2 = time.time()
    print("Step 3: Building 56x56 flip count matrix F...")
    F = np.zeros((num_classes, num_classes), dtype=np.int64)

    for bits in range(NUM_TOURNAMENTS):
        ci = t_class[bits]
        for arc_idx in range(NUM_ARCS):
            flipped_bits = flip_arc(bits, arc_idx)
            cj = t_class[flipped_bits]
            F[ci][cj] += 1

    # Verify: each row of F sums to sizes[ci] * NUM_ARCS
    for ci in range(num_classes):
        row_sum = F[ci].sum()
        expected = sizes[ci] * NUM_ARCS
        assert row_sum == expected, f"Row {ci}: sum={row_sum}, expected={expected}"

    print(f"  F matrix built and verified.")
    print(f"  Row sums = sizes[ci] * 15 (verified for all 56 classes)")
    print(f"  F diagonal (self-flips): {[F[i][i] for i in range(min(10, num_classes))]}")
    print(f"  Time: {time.time()-t2:.2f}s")
    print()

    # Step 4: Transition matrix T = D^{-1} F
    # T[ci][cj] = F[ci][cj] / (sizes[ci] * C(6,2))
    # T is a stochastic matrix (rows sum to 1)
    t3 = time.time()
    print("Step 4: Building transition matrix T = D^{-1} F...")
    T = np.zeros((num_classes, num_classes), dtype=np.float64)
    for ci in range(num_classes):
        denom = sizes[ci] * NUM_ARCS
        for cj in range(num_classes):
            T[ci][cj] = F[ci][cj] / denom

    # Verify stochastic
    for ci in range(num_classes):
        rs = T[ci].sum()
        assert abs(rs - 1.0) < 1e-12, f"Row {ci} sum = {rs}"
    print(f"  T is stochastic (all row sums = 1, verified)")
    print(f"  Time: {time.time()-t3:.2f}s")
    print()

    # Step 5: Compute eigenvalues
    t4 = time.time()
    print("Step 5: Computing eigenvalues of T...")
    eigenvalues = np.linalg.eigvals(T)

    # Check that eigenvalues are real (T might not be symmetric, but for
    # reversible Markov chains it should have real eigenvalues)
    max_imag = max(abs(e.imag) for e in eigenvalues)
    print(f"  Max imaginary part: {max_imag:.2e}")
    if max_imag > 1e-8:
        print("  WARNING: Non-negligible imaginary parts detected!")
        print("  Complex eigenvalues:")
        for e in sorted(eigenvalues, key=lambda x: -abs(x)):
            if abs(e.imag) > 1e-10:
                print(f"    {e}")

    eigs_real = sorted([e.real for e in eigenvalues], reverse=True)
    print(f"  Eigenvalues (real parts, descending):")
    for i, e in enumerate(eigs_real):
        print(f"    λ_{i:2d} = {e:+.12f}")
    print(f"  Time: {time.time()-t4:.2f}s")
    print()

    # Step 6: Round to exact fractions
    print("Step 6: Rounding to exact fractions...")
    print("  Testing denominators up to 1000...")

    eigs_frac = []
    for e in eigs_real:
        best_frac = None
        best_err = 1.0
        for d in range(1, 1001):
            num = round(e * d)
            err = abs(e - num/d)
            if err < best_err:
                best_err = err
                best_frac = Fraction(num, d)
            if err < 1e-10:
                break
        eigs_frac.append(best_frac)

    print(f"\n  Exact eigenvalues as fractions:")
    for i, f in enumerate(eigs_frac):
        print(f"    λ_{i:2d} = {str(f):>8}  ({float(f):+.12f})")

    # Multiplicities
    eig_mult = Counter(eigs_frac)
    print(f"\n  Distinct eigenvalues with multiplicities:")
    for val, mult in sorted(eig_mult.items(), key=lambda x: -float(x[0])):
        print(f"    {str(val):>8}: multiplicity {mult}")
    print(f"  Total: {sum(eig_mult.values())} eigenvalues, {len(eig_mult)} distinct")
    print()

    # Step 7: Test the conjecture — are all eigenvalues multiples of 1/5?
    print("=" * 70)
    print("Step 7: TESTING CONJECTURE — All eigenvalues multiples of 1/5?")
    print("=" * 70)
    print(f"  Largest prime <= {N} is 5.")
    print(f"  Conjecture: every eigenvalue = k/5 for some integer k.")
    print()

    all_mult_of_fifth = True
    for f in eigs_frac:
        if f == 0:
            continue
        # Check if f = k/5 for some integer k
        # f = p/q in lowest terms. f is a multiple of 1/5 iff q divides 5.
        if f.denominator not in (1, 5):
            all_mult_of_fifth = False
            print(f"  FAILS: {f} has denominator {f.denominator}")

    if all_mult_of_fifth:
        print(f"  *** CONJECTURE CONFIRMED: All eigenvalues are multiples of 1/5! ***")
    else:
        print(f"  *** CONJECTURE FAILS ***")
        print(f"  Denominators found: {sorted(set(f.denominator for f in eigs_frac if f != 0))}")

    # Also check denominators
    denoms = set()
    for f in eigs_frac:
        if f != 0:
            denoms.add(f.denominator)
    print(f"  All denominators (in lowest terms): {sorted(denoms)}")
    print()

    # Step 8: Spectral gap
    print("=" * 70)
    print("Step 8: SPECTRAL GAP ANALYSIS")
    print("=" * 70)

    lambda_1 = eigs_frac[0]
    lambda_2 = eigs_frac[1]  # second largest eigenvalue
    lambda_min = eigs_frac[-1]  # most negative eigenvalue

    spectral_gap = lambda_1 - lambda_2
    abs_spectral_gap = lambda_1 - max(abs(lambda_2), abs(lambda_min))

    print(f"  λ_1 (largest)       = {lambda_1}")
    print(f"  λ_2 (second largest)= {lambda_2}")
    print(f"  λ_min (most neg.)   = {lambda_min}")
    print(f"  Spectral gap (1 - λ_2)        = {spectral_gap}")
    print(f"  Absolute spectral gap          = {abs_spectral_gap}")
    print(f"  (uses 1 - max(|λ_2|, |λ_min|))")
    print()

    # Conjecture: spectral gap = 2/5
    expected_gap = Fraction(2, 5)
    print(f"  Conjectured spectral gap: 2/5 = {expected_gap}")
    print(f"  Actual spectral gap:      {spectral_gap}")
    print(f"  Match (1 - λ_2 = 2/5)?   {'YES' if spectral_gap == expected_gap else 'NO'}")
    print()

    # Additional analysis: eigenvalue structure
    print("=" * 70)
    print("ADDITIONAL ANALYSIS")
    print("=" * 70)

    # Check for ± pairing
    print("\n  Positive-negative pairing of eigenvalues:")
    pos_eigs = sorted([f for f in eig_mult if f > 0], reverse=True)
    neg_eigs = sorted([f for f in eig_mult if f < 0])
    zero_mult = eig_mult.get(Fraction(0), 0)
    print(f"  Positive eigenvalues: {[str(f) for f in pos_eigs]}")
    print(f"  Negative eigenvalues: {[str(f) for f in neg_eigs]}")
    print(f"  Zero multiplicity: {zero_mult}")

    paired = True
    for f in pos_eigs:
        if f == Fraction(1):
            continue
        neg_f = -f
        if neg_f not in eig_mult:
            print(f"  No pair for +{f}")
            paired = False
        elif eig_mult[f] != eig_mult[neg_f]:
            print(f"  Multiplicity mismatch: {f} has mult {eig_mult[f]}, {neg_f} has mult {eig_mult[neg_f]}")
            paired = False
    print(f"  All non-trivial eigenvalues paired? {'YES' if paired else 'NO'}")

    # Trace check
    trace_from_eigs = sum(eigs_frac)
    trace_from_T = sum(Fraction(int(F[ci][ci]), sizes[ci] * NUM_ARCS) for ci in range(num_classes))
    print(f"\n  Trace from eigenvalues: {trace_from_eigs} = {float(trace_from_eigs):.6f}")
    print(f"  Trace from T diagonal:  {trace_from_T} = {float(trace_from_T):.6f}")
    print(f"  Match? {'YES' if trace_from_eigs == trace_from_T else 'CLOSE' if abs(float(trace_from_eigs) - float(trace_from_T)) < 1e-8 else 'NO'}")

    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  n = {N}")
    print(f"  Tournaments: {NUM_TOURNAMENTS}")
    print(f"  Isomorphism classes: {num_classes}")
    print(f"  Distinct eigenvalues: {len(eig_mult)}")
    print(f"  All multiples of 1/5: {'YES' if all_mult_of_fifth else 'NO'}")
    print(f"  Spectral gap = {spectral_gap} {'(= 2/5 as conjectured)' if spectral_gap == expected_gap else ''}")
    print(f"  Total time: {time.time()-t0:.2f}s")

    # Comparison with known results for n=3,4,5
    print(f"""
{'='*70}
COMPARISON WITH KNOWN RESULTS
{'='*70}
  n=3: eigenvalues {{1, -1/3}}           — denominator 3, gap 2/3 = 2/3
  n=4: eigenvalues {{1, 1/3, 0, -1/3}}   — denominator 3, gap 2/3 = 2/3
  n=5: eigenvalues {{1, 3/5, 2/5, 1/5(x4), 0, -1/5(x3), -3/5}}
                                        — denominator 5, gap 2/5 = 2/5
  n=6: {num_classes} eigenvalues         — denominator {5 if all_mult_of_fifth else '?'}, gap {spectral_gap}

  The largest prime <= n pattern:
    n=3: p=3, denominator=3, gap=2/3 ✓
    n=4: p=3, denominator=3, gap=2/3 ✓
    n=5: p=5, denominator=5, gap=2/5 ✓
    n=6: p=5, denominator={'5 ✓' if all_mult_of_fifth else 'FAILS'}, gap={'2/5 ✓' if spectral_gap == expected_gap else str(spectral_gap) + ' FAILS'}
""")

    print(f"opus-2026-03-16-S90dk: EIGENVALUE COMPUTATION COMPLETE")

if __name__ == "__main__":
    main()
