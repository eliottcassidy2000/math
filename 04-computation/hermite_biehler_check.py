#!/usr/bin/env python3
"""
hermite_biehler_check.py — oracle-2026-05-21-S1

CRITICAL: Test the Hermite-Biehler condition for TRRT.

I(Omega, x) = A(x) + x * B(x) where:
  A(x) = I(Omega \ C*, x)  [delete C*, keep all others]
  B(x) = I(Omega - N[C*], x) [delete C* and all conflicting cycles]

Hermite-Biehler theorem: A + xB is real-rooted
  iff BOTH A and B are real-rooted AND B interlaces A.

(More precisely: if deg(A) = deg(B) + 1 = n, then A+xB is Hurwitz stable
iff B interlaces A, meaning the roots alternate.)

TESTS:
1. Does B(x) interlace A(x) when degree matches?
2. What happens when degree of B < degree of A - 1?
3. Systematic check across all n=6 tournaments.
4. The CONNECTION to unit distance: "split primes" (degree-1-drop cases)
   vs "inert primes" (degree-drop > 1 cases).

Also investigate: the STABILITY of the polynomial I(Omega,x) under the
deletion-contraction recursion — does it always produce a "stable" sum?
"""

import sys, os, random, time
import numpy as np
from math import comb, sqrt
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (find_odd_cycles, hamiltonian_path_count,
                             tournament_from_bits, random_tournament)

def ip_coeffs(cycles, n):
    m = len(cycles)
    if m == 0: return [1]
    vsets = [frozenset(c) for c in cycles]
    adj = [0]*m
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a] & vsets[b]: adj[a] |= 1<<b; adj[b] |= 1<<a
    max_d = n//3; co = [0]*(max_d+2); co[0]=1; co[1]=m
    pairs = [(a,b) for a in range(m) for b in range(a+1,m) if not(adj[a]>>b&1)]
    co[2] = len(pairs)
    if max_d >= 3:
        trips = [(a,b,c) for a,b in pairs for c in range(b+1,m)
                 if not(adj[a]>>c&1) and not(adj[b]>>c&1)]
        co[3] = len(trips)
    while len(co) > 1 and co[-1] == 0: co.pop()
    return co, adj

def ip_from_subset(cycles, mask, n):
    sub = [cycles[i] for i in range(len(cycles)) if mask>>i&1]
    if not sub: return [1], []
    return ip_coeffs(sub, n)

def real_roots(co):
    if len(co) <= 1: return []
    rs = np.roots(list(reversed(co)))
    return sorted([-r.real for r in rs if abs(r.imag) < 1e-7 and r.real < -1e-10], reverse=True)

def interlaces(P, Q):
    """P interlaces Q: deg(P) = deg(Q) + 1, roots alternate P_1>=Q_1>=P_2>=..."""
    if len(P) != len(Q) + 1: return None
    for i in range(len(Q)):
        if P[i] < Q[i] - 1e-9: return False
        if Q[i] < P[i+1] - 1e-9: return False
    return True

def poly_prod(A, B):
    """Coefficient product of polynomials A and B."""
    d = len(A) + len(B) - 1
    C = [0]*d
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            C[i+j] += a*b
    return C

def eval_poly(co, x):
    return sum(co[k]*x**k for k in range(len(co)))

def add_polys(A, xB):
    """Compute A(x) + x * B(x). A and B are coefficient lists."""
    n = max(len(A), len(xB)+1)
    C = [0]*n
    for i, a in enumerate(A): C[i] += a
    for i, b in enumerate(xB): C[i+1] += b
    while len(C) > 1 and C[-1] == 0: C.pop()
    return C

# ─────────────────────────────────────────────────────────────
# Main Hermite-Biehler check
# ─────────────────────────────────────────────────────────────

def run_hermite_biehler(n, num_samples):
    print(f"{'='*70}")
    print(f"HERMITE-BIEHLER CHECK: n={n}, {num_samples} tournaments")
    print(f"{'='*70}")
    print(f"""
I(Omega,x) = A(x) + x*B(x) where:
  A = I(Omega\\C*) [delete C* from conflict graph]
  B = I(Omega-N[C*]) [keep only cycles DISJOINT from C*]

Hermite-Biehler: I is real-rooted if A,B both real-rooted AND B interlaces A.
For deg(B) = deg(A) - 1: need alternating roots B_1>=A_1>=B_2>=...>=B_{d-1}>=A_d.
""")

    stats = defaultdict(int)
    hb_total = 0  # total (T, C*) pairs where deg(A)=deg(B)+1 [the Hermite-Biehler case]
    hb_holds = 0  # pairs where B interlaces A

    degree_drop_stats = defaultdict(int)

    t0 = time.time()
    for _ in range(num_samples):
        if n == 6:
            T = tournament_from_bits(n, random.randint(0, 2**15))
        else:
            T = random_tournament(n)

        cycles = find_odd_cycles(T)
        m = len(cycles)
        if m < 2: continue

        co_full, adj = ip_coeffs(cycles, n)
        d_full = len(co_full) - 1
        if d_full < 1: continue

        for star in range(min(m, 6)):  # test first 6 cycles
            # A = I(Omega \ C*) — delete C*, keep all others
            mask_A = (1<<m) - 1 - (1<<star)
            co_A, _ = ip_from_subset(cycles, mask_A, n)
            d_A = len(co_A) - 1

            # B = I(Omega - N[C*]) — keep only cycles DISJOINT from C*
            disjoint_mask = 0
            for i in range(m):
                if i != star and not (adj[star]>>i&1):
                    disjoint_mask |= 1<<i
            co_B, _ = ip_from_subset(cycles, disjoint_mask, n)
            d_B = len(co_B) - 1

            degree_drop = d_full - d_A
            degree_drop_stats[degree_drop] += 1

            # Verify: A + xB should equal I(Omega)
            co_check = add_polys(co_A, co_B)
            if co_check != co_full:
                stats['recursion_violation'] += 1
                if stats['recursion_violation'] <= 3:
                    print(f"  RECURSION VIOLATION: I={co_full}, A+xB={co_check}")

            # Get real roots
            roots_full = real_roots(co_full)
            roots_A = real_roots(co_A)
            roots_B = real_roots(co_B)

            # Check: A interlaces I(Omega)?
            r_AI = interlaces(roots_full, roots_A)
            if r_AI is not None:
                if r_AI:
                    stats['A_interlaces_I'] += 1
                else:
                    stats['A_not_interlaces_I'] += 1

            # MAIN CHECK: B interlaces A? (Hermite-Biehler condition)
            if d_A == d_B + 1:
                hb_total += 1
                r_BA = interlaces(roots_A, roots_B)
                if r_BA is True:
                    hb_holds += 1
                    stats['HB_holds'] += 1
                elif r_BA is False:
                    stats['HB_fails'] += 1
                    if stats['HB_fails'] <= 5:
                        H = hamiltonian_path_count(T)
                        print(f"  HB FAIL: H={H}, A_roots={[f'{r:.3f}' for r in roots_A]}, "
                              f"B_roots={[f'{r:.3f}' for r in roots_B]}")

            # Alternative: B interlaces A with degree B < A-1?
            elif d_B <= d_A - 2:
                # B has degree at least 2 less than A — weaker condition
                # Check: all roots of B lie between extremes of A
                if roots_A and roots_B:
                    all_in_range = all(roots_A[-1] <= r <= roots_A[0] for r in roots_B)
                    if all_in_range:
                        stats['B_in_A_range'] += 1
                    else:
                        stats['B_out_A_range'] += 1

    elapsed = time.time() - t0
    print(f"\nResults ({elapsed:.1f}s):")
    print(f"\nDegree-drop distribution: {dict(sorted(degree_drop_stats.items()))}")
    print(f"\nA interlaces I(Omega): {stats['A_interlaces_I']} holds, {stats['A_not_interlaces_I']} fails")
    print(f"\nHermite-Biehler (B interlaces A, deg(A)=deg(B)+1):")
    print(f"  Total HB cases: {hb_total}")
    print(f"  HB holds: {hb_holds} ({100*hb_holds//max(hb_total,1)}%)")
    print(f"  HB fails: {stats['HB_fails']}")
    if stats['HB_fails'] == 0:
        print(f"  ✓ HERMITE-BIEHLER HOLDS in all {hb_total} testable cases!")
    else:
        print(f"  ✗ Hermite-Biehler FAILS in {stats['HB_fails']} cases!")

    print(f"\nDegree-drop > 1 cases: B in A's range: {stats['B_in_A_range']}, out: {stats['B_out_A_range']}")
    print(f"Recursion violations: {stats['recursion_violation']}")

    return stats

# ─────────────────────────────────────────────────────────────
# Deep check: all n=6 degree-2 polynomials
# ─────────────────────────────────────────────────────────────

def n6_deep_check():
    n = 6
    print(f"\n{'='*70}")
    print(f"N=6 DEEP CHECK: Hermite-Biehler for all degree-2 cases")
    print(f"{'='*70}")

    hb_total = 0; hb_holds = 0; hb_fails = []

    for bits in range(0, 2**15, 7):  # stride 7 → ~4700 samples
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        m = len(cycles)
        if m < 2: continue

        co_full, adj = ip_coeffs(cycles, n)
        d_full = len(co_full) - 1
        if d_full != 2: continue  # only degree-2

        roots_full = real_roots(co_full)
        if len(roots_full) != 2: continue

        for star in range(m):
            mask_A = (1<<m) - 1 - (1<<star)
            co_A, _ = ip_from_subset(cycles, mask_A, n)
            d_A = len(co_A) - 1

            disjoint_mask = 0
            for i in range(m):
                if i != star and not (adj[star]>>i&1):
                    disjoint_mask |= 1<<i
            co_B, _ = ip_from_subset(cycles, disjoint_mask, n)
            d_B = len(co_B) - 1

            if d_A != 1 or d_B != 0:
                continue  # only check the simplest HB case: deg(A)=1, deg(B)=0

            roots_A = real_roots(co_A)
            roots_B = real_roots(co_B)  # empty for degree 0

            hb_total += 1
            # For deg(A)=1 and deg(B)=0=constant: B=1, trivially interlaces A
            # (no roots to check)
            hb_holds += 1  # trivially true

    print(f"  Degree-2 with deg(A)=1, deg(B)=0: {hb_total} cases, all trivially HB ✓")

# ─────────────────────────────────────────────────────────────
# NEW INSIGHT: Verify the recursion I = A + xB directly
# ─────────────────────────────────────────────────────────────

def verify_recursion():
    print(f"\n{'='*70}")
    print(f"RECURSION VERIFICATION: I(Omega) = A(x) + x*B(x)")
    print(f"{'='*70}")
    print(f"""
Checking: I(Omega, x) = I(Omega\\C*, x) + x * I(Omega-N[C*], x)
for all cycles C* in Omega.

This is the standard deletion-contraction for independence polynomials:
  I(G, x) = I(G-v, x) + x * I(G - N[v], x)
where G-v removes vertex v and G-N[v] removes v and all its neighbors.
""")

    n = 6; violations = 0; checks = 0
    for bits in range(0, 2**15, 23):
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        m = len(cycles)
        if m == 0: continue
        co_full, adj = ip_coeffs(cycles, n)

        for star in range(min(m, 4)):
            mask_A = (1<<m) - 1 - (1<<star)
            co_A, _ = ip_from_subset(cycles, mask_A, n)

            disjoint_mask = 0
            for i in range(m):
                if i != star and not (adj[star]>>i&1):
                    disjoint_mask |= 1<<i
            co_B, _ = ip_from_subset(cycles, disjoint_mask, n)

            co_check = add_polys(co_A, co_B)
            # Compare
            d = max(len(co_full), len(co_check))
            match = all(
                (co_full[k] if k < len(co_full) else 0) == (co_check[k] if k < len(co_check) else 0)
                for k in range(d)
            )
            checks += 1
            if not match:
                violations += 1
                if violations <= 3:
                    print(f"  VIOLATION at bits={bits}, star={star}:")
                    print(f"    I(Omega) = {co_full}")
                    print(f"    A+xB    = {co_check}")

    print(f"Checked {checks} (T, C*) pairs: {violations} violations "
          f"({'✓ all correct' if violations==0 else '✗ FOUND!'})")

# ─────────────────────────────────────────────────────────────
# Key experiment: when degree(B) = degree(A), is there interlacing?
# ─────────────────────────────────────────────────────────────

def equal_degree_check():
    print(f"\n{'='*70}")
    print(f"EQUAL-DEGREE HERMITE-BIEHLER: deg(B) = deg(A)")
    print(f"{'='*70}")
    print(f"""
When deg(A) = deg(B) (both have same degree d):
  I = A + xB has degree d+1.
  Hermite-Biehler: I real-rooted iff A and B are real-rooted
  and their roots alternate: A_1 >= B_1 >= A_2 >= ... >= B_d >= A_{d+1}.
  (Here B "interlaces" A in the sense that B's roots fit between A's.)
""")

    n = 7
    equal_deg_total = 0
    equal_deg_holds = 0

    for _ in range(500):
        T = random_tournament(n)
        cycles = find_odd_cycles(T)
        m = len(cycles)
        if m < 3: continue

        co_full, adj = ip_coeffs(cycles, n)
        d_full = len(co_full) - 1
        if d_full < 2: continue

        for star in range(min(m, 4)):
            mask_A = (1<<m) - 1 - (1<<star)
            co_A, _ = ip_from_subset(cycles, mask_A, n)
            d_A = len(co_A) - 1

            disjoint_mask = 0
            for i in range(m):
                if i != star and not (adj[star]>>i&1):
                    disjoint_mask |= 1<<i
            co_B, _ = ip_from_subset(cycles, disjoint_mask, n)
            d_B = len(co_B) - 1

            if d_A != d_B: continue  # only equal degree

            equal_deg_total += 1
            roots_A = real_roots(co_A)
            roots_B = real_roots(co_B)

            if len(roots_A) != d_A or len(roots_B) != d_B: continue

            # Check B interlaces A with equal degrees:
            # B_1 >= A_1 >= B_2 >= A_2 >= ... >= B_d >= A_d  (or reversed)
            # Or: all B roots fit between consecutive A roots
            ok = True
            for i in range(min(len(roots_A), len(roots_B))):
                # Each B root should be "between" two A roots
                pass  # complex check — just check numerically

            # Simple check: for same-degree case, interlacing means
            # roots_A and roots_B alternate when merged and sorted
            merged = sorted(roots_A + roots_B, reverse=True)
            # Check alternation: A,B,A,B,... or B,A,B,A,...
            from_A = [i%2 for i, r in enumerate(merged) if r in roots_A]
            # This doesn't work for repeated roots, use approximate check
            eps = 1e-7
            idx_A = []
            idx_B = []
            for i, r in enumerate(merged):
                if any(abs(r-ra) < eps for ra in roots_A): idx_A.append(i)
                elif any(abs(r-rb) < eps for rb in roots_B): idx_B.append(i)

            # Check strict interleaving
            all_idx = sorted(idx_A + idx_B)
            parity_ok = all(
                (all_idx[k] in idx_A) != (all_idx[k+1] in idx_A)
                for k in range(len(all_idx)-1)
            ) if len(all_idx) >= 2 else True

            if parity_ok:
                equal_deg_holds += 1

    print(f"Equal-degree (A,B same degree) at n=7: {equal_deg_total} cases, "
          f"{equal_deg_holds} interlace ({100*equal_deg_holds//max(equal_deg_total,1)}%)")


if __name__ == '__main__':
    random.seed(42)
    np.random.seed(42)

    verify_recursion()
    run_hermite_biehler(6, 300)
    run_hermite_biehler(7, 200)
    equal_degree_check()

    print(f"\n{'='*70}\nDONE\n{'='*70}")
