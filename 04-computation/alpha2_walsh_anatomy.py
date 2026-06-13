#!/usr/bin/env python3
"""
alpha2_walsh_anatomy.py -- Anatomy of the OCF cancellation at p=7

At p=7:
- H = 1 + 2*alpha_1 + 4*alpha_2 (alpha_j=0 for j>=3, since max 2 disjoint cycles)
- c_3 = 14 = CONSTANT (no Walsh content)
- But the WHICH 14 triples form 3-cycles CHANGES with orientation
- alpha_2 = # disjoint 3-3 pairs, varies by orientation
- h_hat_H = 2*h_hat_alpha_1 + 4*h_hat_alpha_2 = 2*(-15.75) + 4*(8.75) = 3.5

KEY QUESTION: What determines sign(h_hat_alpha_2) = chi(ab)?
Is it the change in WHICH triples are 3-cycles, or the change in disjointness?

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def classify_resonance(a, b, p):
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def get_3cycle_sets(A, p):
    """Get all 3-vertex sets that form directed 3-cycles."""
    c3_sets = []
    for triple in combinations(range(p), 3):
        a, b, c = triple
        # Check if this triple has a directed 3-cycle
        if (A[a][b] and A[b][c] and A[c][a]) or \
           (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset(triple))
    return c3_sets


def get_all_odd_cycle_sets(A, p, max_k=None):
    """Get all odd cycle vertex sets."""
    if max_k is None:
        max_k = p
    all_cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_ham_cycles(A, verts)
            for _ in range(n_cyc):
                all_cycles.append(frozenset(subset))
    return all_cycles


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def count_disjoint_pairs(cycles):
    """Count pairs of cycles with disjoint vertex sets."""
    n = len(cycles)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            if not (cycles[i] & cycles[j]):
                count += 1
    return count


def main():
    p = 7
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    half = 1 << m

    print("=" * 70)
    print(f"ALPHA_2 WALSH ANATOMY at p={p}")
    print("=" * 70)

    # For each orientation, compute:
    # - which 14 triples are 3-cycles
    # - all odd cycle sets
    # - alpha_1, alpha_2
    # - H

    data = {}  # bits -> (c3_sets, alpha_1, alpha_2, H)

    for bits in range(half):
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                   for i in range(m))
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]

        A = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in S:
                A[v][(v + s) % p] = 1

        c3_sets = get_3cycle_sets(A, p)
        all_cycles = get_all_odd_cycle_sets(A, p)

        alpha_1 = len(all_cycles)
        alpha_2 = count_disjoint_pairs(all_cycles)
        H = 1 + 2 * alpha_1 + 4 * alpha_2

        # Count cycles by length
        by_len = defaultdict(int)
        for fs in all_cycles:
            by_len[len(fs)] += 1

        data[bits] = {
            'sigma': sigma,
            'S': S,
            'c3_sets': set(c3_sets),
            'alpha_1': alpha_1,
            'alpha_2': alpha_2,
            'H': H,
            'by_len': dict(by_len),
            'all_cycles': all_cycles
        }

        print(f"\n  bits={bits:03b}, sigma={sigma}, S={S}")
        print(f"    c3_sets: {len(c3_sets)}")
        print(f"    by_len: {dict(sorted(by_len.items()))}")
        print(f"    alpha_1={alpha_1}, alpha_2={alpha_2}, H={H}")

    # Which triples are 3-cycles for each orientation?
    print(f"\n\n{'='*60}")
    print("3-CYCLE VERTEX SETS BY ORIENTATION")
    print("=" * 60)

    all_triples = list(combinations(range(p), 3))
    triple_to_idx = {frozenset(t): i for i, t in enumerate(all_triples)}

    # For each triple, list which orientations have it as a 3-cycle
    triple_orientations = defaultdict(set)
    for bits in range(half):
        for fs in data[bits]['c3_sets']:
            triple_orientations[fs].add(bits)

    print(f"\n  Total triples: {len(all_triples)}")
    print(f"  Always 3-cycle: {sum(1 for fs in all_triples if len(triple_orientations[frozenset(fs)]) == half)}")
    print(f"  Never 3-cycle: {sum(1 for fs in all_triples if len(triple_orientations.get(frozenset(fs), set())) == 0)}")
    print(f"  Sometimes 3-cycle: {sum(1 for fs in all_triples if 0 < len(triple_orientations.get(frozenset(fs), set())) < half)}")

    # Detail: for each triple, what fraction of orientations have it
    for triple in sorted(all_triples):
        fs = frozenset(triple)
        count = len(triple_orientations.get(fs, set()))
        if count not in [0, half]:
            print(f"    {triple}: {count}/{half} orientations")

    # Walsh decomposition of alpha_1, alpha_2, and H
    print(f"\n\n{'='*60}")
    print("WALSH DECOMPOSITION OF alpha_1, alpha_2, H")
    print("=" * 60)

    for name, vals_fn in [("alpha_1", lambda d: d['alpha_1']),
                          ("alpha_2", lambda d: d['alpha_2']),
                          ("c_3", lambda d: len(d['c3_sets'])),
                          ("c_5", lambda d: d['by_len'].get(5, 0)),
                          ("c_7", lambda d: d['by_len'].get(7, 0)),
                          ("H", lambda d: d['H'])]:
        print(f"\n  {name}:")
        print(f"    Values: {[vals_fn(data[bits]) for bits in range(half)]}")

        # Degree-0 (mean)
        mean = sum(vals_fn(data[bits]) for bits in range(half)) / half
        print(f"    Mean: {mean:.4f}")

        # Degree-2 Walsh
        for a in range(m):
            for b in range(a + 1, m):
                total = 0
                for bits in range(half):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += vals_fn(data[bits]) * sa * sb
                h = total / half
                if abs(h) > 0.001:
                    gap_a, gap_b = a + 1, b + 1
                    chi_ab = legendre(gap_a * gap_b, p)
                    print(f"    h_hat[{{{a},{b}}}] = {h:>10.4f}, "
                          f"gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}, "
                          f"sign {'OK' if (h>0)==(chi_ab>0) else 'FAIL'}")

    # OCF verification
    print(f"\n\n{'='*60}")
    print("OCF DECOMPOSITION: h_hat_H = 2*h_hat_alpha1 + 4*h_hat_alpha2")
    print("=" * 60)

    for a in range(m):
        for b in range(a + 1, m):
            h_a1 = sum((1 if bits & (1<<a) else -1) * (1 if bits & (1<<b) else -1) * data[bits]['alpha_1'] for bits in range(half)) / half
            h_a2 = sum((1 if bits & (1<<a) else -1) * (1 if bits & (1<<b) else -1) * data[bits]['alpha_2'] for bits in range(half)) / half
            h_H = sum((1 if bits & (1<<a) else -1) * (1 if bits & (1<<b) else -1) * data[bits]['H'] for bits in range(half)) / half

            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)
            print(f"\n  Pair ({a},{b}), gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}:")
            print(f"    h_hat_alpha1 = {h_a1:>10.4f} (sign {'+' if h_a1>0 else '-'})")
            print(f"    h_hat_alpha2 = {h_a2:>10.4f} (sign {'+' if h_a2>0 else '-'})")
            print(f"    2*h_hat_a1   = {2*h_a1:>10.4f}")
            print(f"    4*h_hat_a2   = {4*h_a2:>10.4f}")
            print(f"    h_hat_H      = {h_H:>10.4f} = {2*h_a1:>10.4f} + {4*h_a2:>10.4f}")
            print(f"    Cancellation: {abs(2*h_a1)/(abs(2*h_a1)+abs(4*h_a2))*100:.1f}% / "
                  f"{abs(4*h_a2)/(abs(2*h_a1)+abs(4*h_a2))*100:.1f}%")

    # PART 2: What makes alpha_2 follow the product law?
    print(f"\n\n{'='*60}")
    print("PART 2: alpha_2 VARIATION ANATOMY")
    print("=" * 60)

    # For each orientation, list the disjoint cycle pairs
    for bits in range(half):
        sigma = data[bits]['sigma']
        S = data[bits]['S']
        all_cyc = data[bits]['all_cycles']
        alpha_2 = data[bits]['alpha_2']

        # Find the actual disjoint pairs
        disjoint = []
        for i in range(len(all_cyc)):
            for j in range(i+1, len(all_cyc)):
                if not (all_cyc[i] & all_cyc[j]):
                    disjoint.append((all_cyc[i], all_cyc[j]))

        print(f"\n  bits={bits:03b}, S={S}, alpha_2={alpha_2}:")
        for c1, c2 in disjoint:
            free = set(range(p)) - (c1 | c2)
            print(f"    {sorted(c1)} + {sorted(c2)}, free={sorted(free)}, "
                  f"lens=({len(c1)},{len(c2)})")

    # PART 3: Symmetric difference analysis
    print(f"\n\n{'='*60}")
    print("PART 3: 3-CYCLE SET CHANGES BETWEEN ORIENTATIONS")
    print("=" * 60)

    # For pairs of orientations that differ in exactly one chord, what changes?
    for a in range(m):
        for bits1 in range(half):
            bits2 = bits1 ^ (1 << a)
            if bits1 > bits2:
                continue
            gained = data[bits2]['c3_sets'] - data[bits1]['c3_sets']
            lost = data[bits1]['c3_sets'] - data[bits2]['c3_sets']
            if gained or lost:
                print(f"\n    Flip chord {a} (gap {a+1}<->{p-a-1}): "
                      f"bits {bits1:03b}->{bits2:03b}")
                print(f"      Lost: {[sorted(fs) for fs in lost]}")
                print(f"      Gained: {[sorted(fs) for fs in gained]}")
                print(f"      c3: {len(data[bits1]['c3_sets'])}->{len(data[bits2]['c3_sets'])}")
                break  # just one example per chord

    # For pairs that differ in exactly TWO chords
    print(f"\n  Two-chord flips and alpha_2 change:")
    for a in range(m):
        for b in range(a + 1, m):
            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)
            # Compare bits=0 with bits=(1<<a)|(1<<b)
            bits1 = 0
            bits2 = (1 << a) | (1 << b)

            da1 = data[bits2]['alpha_1'] - data[bits1]['alpha_1']
            da2 = data[bits2]['alpha_2'] - data[bits1]['alpha_2']
            dH = data[bits2]['H'] - data[bits1]['H']

            # c3 set changes
            gained = data[bits2]['c3_sets'] - data[bits1]['c3_sets']
            lost = data[bits1]['c3_sets'] - data[bits2]['c3_sets']

            print(f"\n    Flip ({a},{b}), gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}:")
            print(f"      delta_alpha1={da1:+d}, delta_alpha2={da2:+d}, delta_H={dH:+d}")
            print(f"      3-cycles gained: {len(gained)}, lost: {len(lost)}")
            print(f"        gained: {[sorted(fs) for fs in gained]}")
            print(f"        lost: {[sorted(fs) for fs in lost]}")


if __name__ == '__main__':
    main()
