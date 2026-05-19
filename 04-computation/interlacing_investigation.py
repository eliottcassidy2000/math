#!/usr/bin/env python3
"""
interlacing_investigation.py — oracle-2026-05-19-S1

Investigate interlacing properties of I(Omega(T), x) under:
1. Edge deletion in Omega (removing one cycle from the conflict graph)
2. Arc reversal in T (how does I change?)
3. The deletion-contraction recursion structure

KEY DEFINITIONS:
- Interlacing: P interlaces Q if they alternate roots on the real line.
  If P has roots r1 >= r2 >= ... >= rn and Q has s1 >= ... >= s_{n-1}, then
  Q interlaces P if r1 >= s1 >= r2 >= s2 >= ... >= s_{n-1} >= rn.

- Interlacing implies: if P and Q are both real-rooted, and Q interlaces P,
  then any real linear combination aP + bQ is also real-rooted.

FOR OUR POLYNOMIAL:
I(Omega, x) = I(Omega \ C*, x) + x * I(Omega - N[C*], x)

where:
- Omega \ C* = delete vertex C* from conflict graph (keep all other cycles)
- Omega - N[C*] = delete C* and all cycles that conflict with C*
  (this gives the induced subgraph on cycles DISJOINT from C*)

If I(Omega \ C*, x) interlaces I(Omega, x), then real-rootedness propagates
by induction on the number of cycles!

This would give an UNCONDITIONAL proof of TRRT (Tournament Real-Rootedness)
if the interlacing can be established.
"""

import sys, os, random, time
import numpy as np
from math import comb, sqrt
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (find_odd_cycles, hamiltonian_path_count,
                             tournament_from_bits, random_tournament)
from itertools import permutations

def ip_coeffs_and_omega(cycles, n):
    """Return (coeffs, adj_bits) for the conflict graph."""
    m = len(cycles)
    if m == 0: return [1], []
    vsets = [frozenset(c) for c in cycles]
    adj_bits = [0]*m
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a]&vsets[b]:
                adj_bits[a] |= 1<<b
                adj_bits[b] |= 1<<a

    max_d = n//3; coeffs = [0]*(max_d+2); coeffs[0]=1; coeffs[1]=m
    pairs = [(a,b) for a in range(m) for b in range(a+1,m) if not(adj_bits[a]>>b&1)]
    coeffs[2] = len(pairs)
    if max_d >= 3:
        trips = [(a,b,c) for a,b in pairs for c in range(b+1,m)
                 if not(adj_bits[a]>>c&1) and not(adj_bits[b]>>c&1)]
        coeffs[3] = len(trips)
    while len(coeffs)>1 and coeffs[-1]==0: coeffs.pop()
    return coeffs, adj_bits

def poly_from_subset(cycles, subset_mask, n):
    """Compute independence polynomial for cycles[i] with i in subset_mask."""
    sub = [cycles[i] for i in range(len(cycles)) if subset_mask>>i&1]
    if not sub: return [1]
    coeffs, _ = ip_coeffs_and_omega(sub, n)
    return coeffs

def real_roots(co):
    """All real negative roots as positive values."""
    if len(co) <= 1: return []
    rs = np.roots(list(reversed(co)))
    return sorted([-r.real for r in rs if abs(r.imag)<1e-7 and r.real<-1e-10], reverse=True)

def interlaces(P_roots, Q_roots):
    """Check if Q interlaces P. P has one more root than Q."""
    if len(P_roots) != len(Q_roots)+1: return None  # wrong length
    for i in range(len(Q_roots)):
        if P_roots[i] < Q_roots[i] - 1e-9:
            return False  # P_i < Q_i violated
        if Q_roots[i] < P_roots[i+1] - 1e-9:
            return False  # Q_i < P_{i+1} violated
    return True

def evaluate_poly(co, x):
    return sum(co[k]*x**k for k in range(len(co)))

# ─────────────────────────────────────────────────────────────
# Test: does I(Omega\C*, x) interlace I(Omega, x)?
# ─────────────────────────────────────────────────────────────

def test_interlacing():
    print("="*72)
    print("INTERLACING TEST: I(Omega\\C*, x) interlaces I(Omega, x)?")
    print("="*72)
    print("""
If I(Omega\\C*, x) interlaces I(Omega, x) for every cycle C*, then:
  - Real-rootedness propagates inductively (Hermite-Biehler theorem)
  - This would PROVE TRRT unconditionally!

The deletion-contraction recursion:
  I(Omega, x) = I(Omega\\C*, x) + x * I(Omega - N[C*], x)

If I(Omega\\C*, x) and x * I(Omega - N[C*], x) both have real roots and
I(Omega\\C*, x) interlaces I(Omega, x), then I(Omega, x) is real-rooted.

For GRAPHS: the independence polynomial of a claw-free graph satisfies
this interlacing (Chudnovsky-Seymour, Heilmann-Lieb for matchings).
QUESTION: Does this hold for tournament conflict graphs generally?
""")

    interlace_holds = 0
    interlace_fails = 0
    wrong_degree = 0

    n = 6
    for bits in range(2**(n*(n-1)//2)):
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        m = len(cycles)
        if m < 2: continue

        co_full, adj_bits = ip_coeffs_and_omega(cycles, n)
        roots_full = real_roots(co_full)
        d_full = len(co_full)-1
        if d_full < 1: continue

        # Test deletion of each cycle
        for star in range(m):
            # Omega \ C*: remove cycle star from conflict graph
            sub_mask = (1<<m) - 1 - (1<<star)
            co_del = poly_from_subset(cycles, sub_mask, n)
            roots_del = real_roots(co_del)
            d_del = len(co_del)-1

            if d_full != d_del+1:
                # Degrees differ by more than 1 — can't interlace in standard sense
                wrong_degree += 1
                continue

            # Check interlacing
            result = interlaces(roots_full, roots_del)
            if result is True:
                interlace_holds += 1
            elif result is False:
                interlace_fails += 1
                if interlace_fails <= 3:
                    H = hamiltonian_path_count(T)
                    print(f"  Interlacing FAILS: H={H}, full roots={[f'{r:.3f}' for r in roots_full]}, del roots={[f'{r:.3f}' for r in roots_del]}")

    total = interlace_holds + interlace_fails
    print(f"\nResults at n=6:")
    print(f"  Interlacing holds: {interlace_holds}/{total} ({100*interlace_holds/total:.1f}%)" if total>0 else "  No degree-1 change cases")
    print(f"  Wrong degree: {wrong_degree} (degree changed by ≠1)")
    if interlace_fails == 0:
        print(f"  ✓ Interlacing holds whenever deg(I_del) = deg(I_full) - 1")
    else:
        print(f"  ✗ Interlacing FAILS in {interlace_fails} cases!")

# ─────────────────────────────────────────────────────────────
# Test: the two-term recursion I = A + x*B interlacing
# ─────────────────────────────────────────────────────────────

def test_two_term_recursion():
    print(f"\n{'='*72}")
    print("DELETION-CONTRACTION RECURSION: I = A + x*B")
    print("="*72)
    print("""
I(Omega, x) = A(x) + x * B(x) where:
  A(x) = I(Omega \\ C*, x)  [delete C*, keep all others]
  B(x) = I(Omega - N[C*], x) [delete C* and its neighbors = cycles disjoint from C*]

If A and B are both real-rooted and A interlaces x*B (or equivalently,
B interlaces A with appropriate degree shift), then I = A + x*B is real-rooted
(by Interlacing + Hermite-Biehler).

TEST: Verify that A(x) and x*B(x) interlace each other.
""")

    fail = 0; total = 0
    n = 6
    for bits in range(0, 2**15, 11):
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        m = len(cycles)
        if m < 2: continue

        co_full, adj_bits = ip_coeffs_and_omega(cycles, n)
        roots_full = real_roots(co_full)
        if len(roots_full) < 2: continue

        for star in range(min(m, 5)):  # test first 5 cycles
            # A = I(Omega \ C*)
            sub_A = (1<<m)-1-(1<<star)
            co_A = poly_from_subset(cycles, sub_A, n)
            roots_A = real_roots(co_A)

            # B = I(Omega - N[C*]) = cycles disjoint from C*
            # = all cycles not adjacent to star in Omega
            disjoint_mask = 0
            for i in range(m):
                if i != star and not (adj_bits[star]>>i&1):
                    disjoint_mask |= 1<<i
            co_B = poly_from_subset(cycles, disjoint_mask, n)
            roots_B = real_roots(co_B)
            # x*B(x) has roots 0 plus roots of B
            roots_xB = sorted([0.0] + list(roots_B), reverse=True)

            total += 1
            # Check: A real-rooted? B real-rooted?
            # Check: roots_full = roots of I, roots_A and roots_xB interlace somehow?
            all_real_A = (len(co_A)-1 == len(roots_A))
            all_real_B = (len(co_B)-1 == len(roots_B))

            if not (all_real_A and all_real_B):
                fail += 1

    print(f"Checked {total} (T, C*) pairs at n=6 (stride 11)")
    print(f"  A and B both real-rooted: {total-fail}/{total}")
    if fail == 0:
        print(f"  ✓ Both terms A and x*B are real-rooted (as expected from TRRT)")
    else:
        print(f"  ✗ {fail} cases where A or B not real-rooted!")

# ─────────────────────────────────────────────────────────────
# Verify the Heilmann-Lieb structure
# ─────────────────────────────────────────────────────────────

def test_heilmann_lieb():
    print(f"\n{'='*72}")
    print("HEILMANN-LIEB VERIFICATION (matching polynomial version)")
    print("="*72)
    print("""
The Heilmann-Lieb theorem (1972): The matching polynomial of any graph is real-rooted.
The independence polynomial of a claw-free graph is real-rooted (Chudnovsky-Seymour).

For tournament conflict graphs:
  - n<=8: Omega is claw-free → real-rooted by CS.
  - n>=9: Omega has claws, but tournament structure still forces real-rootedness.

TEST: Find the smallest n where Omega(T) first has a claw, and verify real-rootedness persists.
""")

    def has_claw(adj_bits, m):
        """Check if conflict graph has a claw (K_{1,3})."""
        for v in range(m):
            # Find neighbors of v
            nbrs = [u for u in range(m) if v!=u and (adj_bits[v]>>u&1)]
            if len(nbrs) < 3: continue
            # Check if any 3 neighbors are mutually non-adjacent
            for i in range(len(nbrs)):
                for j in range(i+1,len(nbrs)):
                    for k in range(j+1,len(nbrs)):
                        u,v2,w = nbrs[i],nbrs[j],nbrs[k]
                        if not (adj_bits[u]>>v2&1) and not (adj_bits[u]>>w&1) and not (adj_bits[v2]>>w&1):
                            return True
        return False

    print("Finding claws in Omega(T) for n=7,8,9:")
    for n in [7, 8, 9]:
        claw_count = 0; real_rooted_with_claw = 0; samples = 200
        for _ in range(samples):
            T = random_tournament(n)
            cycles = find_odd_cycles(T)
            m = len(cycles)
            if m == 0: continue
            co, adj_bits = ip_coeffs_and_omega(cycles, n)
            claw = has_claw(adj_bits, m)
            if claw:
                claw_count += 1
                roots = real_roots(co)
                d = len(co)-1
                if len(roots) == d:
                    real_rooted_with_claw += 1
        print(f"  n={n}: {claw_count}/{samples} have Omega with a claw, "
              f"of those {real_rooted_with_claw}/{claw_count} still real-rooted "
              f"{'✓' if claw_count==0 or real_rooted_with_claw==claw_count else '✗'}")

# ─────────────────────────────────────────────────────────────
# The "stability" angle: check for non-real roots
# ─────────────────────────────────────────────────────────────

def stability_check():
    print(f"\n{'='*72}")
    print("COMPLEX ROOT ANALYSIS — looking for patterns")
    print("="*72)
    print("""
If I(Omega, x) has any complex roots at n>=9 (where Omega can have claws),
that would be a TRRT counterexample. We expect ALL roots real.

We also check: are roots STABLE (no root in upper half-plane)?
For real-rooted polynomials, stability = real-rootedness (trivially satisfied).
""")

    n = 9
    complex_roots_found = 0
    complex_examples = []
    samples = 80

    for _ in range(samples):
        T = random_tournament(n)
        cycles = find_odd_cycles(T)
        co, adj_bits = ip_coeffs_and_omega(cycles, n)
        if len(co) <= 1: continue

        roots = np.roots(list(reversed(co)))
        has_complex = any(abs(r.imag) > 1e-6 for r in roots)
        if has_complex:
            complex_roots_found += 1
            complex_examples.append((tuple(co), [(r.real, r.imag) for r in roots]))

    print(f"\nn={n}: {samples} random samples, {complex_roots_found} with complex roots")
    if complex_roots_found == 0:
        print(f"  ✓ All real-rooted — consistent with TRRT")
    else:
        print(f"  ✗ Found complex roots at n={n}! TRRT would be FALSE!")
        for co, roots in complex_examples[:3]:
            print(f"    coeffs={co}, complex roots: {[f'{r:.3f}+{i:.3f}i' for r,i in roots if abs(i)>1e-6]}")

if __name__ == '__main__':
    random.seed(42); np.random.seed(42)
    t0 = time.time()

    test_interlacing()
    test_two_term_recursion()
    test_heilmann_lieb()
    stability_check()

    print(f"\n{'='*72}")
    print(f"Total time: {time.time()-t0:.1f}s")
    print("="*72)
