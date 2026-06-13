#!/usr/bin/env python3
"""
Investigate the pairing structure of f(S) for a general proof.

H(T) - H(T') = sum_{S subset V\{i,j}} f(S)

where f(S) = h_end(S+{i}, i) * h_start({j}+R, j)
            - h_end(S+{j}, j) * h_start({i}+R, i)
and R = V\{i,j}\S.

Key question: does f(S) + f(S^c) have a universal form?
At n=4: f(S) + f(S^c) = -sum(s_x) for all S, giving -2*sum(s_x) total.

Instance: kind-pasteur-2026-03-05-S7
"""

from itertools import permutations, combinations
from collections import defaultdict


def h_end_v(T, verts, v):
    """Count Ham paths on verts ending at v."""
    if len(verts) == 1:
        return 1 if v == verts[0] else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        if all(T(p[k], p[k+1]) for k in range(len(p)-1)):
            count += 1
    return count


def h_start_v(T, verts, v):
    """Count Ham paths on verts starting at v."""
    if len(verts) == 1:
        return 1 if v == verts[0] else 0
    count = 0
    for p in permutations(verts):
        if p[0] != v:
            continue
        if all(T(p[k], p[k+1]) for k in range(len(p)-1)):
            count += 1
    return count


def analyze_pairing(n):
    """Analyze f(S) + f(S^c) structure."""
    print(f"=== Pairing Structure at n={n} ===\n")

    I, J = 0, 1
    others = list(range(2, n))

    # Set up variables
    interface = [(x, I) for x in others] + [(J, x) for x in others]
    internal = [(a, b) for a in others for b in others if a < b]
    all_vars = interface + internal
    n_vars = len(all_vars)

    print(f"Variables: {n_vars}, Cases: {1 << n_vars}")

    # For each variable assignment, compute f(S) for all S
    # and check pairing structure

    num_others = len(others)
    all_subsets = list(range(1 << num_others))  # bitmask over others

    def subset_from_mask(smask):
        return [others[bit] for bit in range(num_others) if smask & (1 << bit)]

    # Statistics: for each pairing (S, S^c), what is f(S)+f(S^c)?
    full_mask = (1 << num_others) - 1

    # Track: for each S-size, what is avg f(S)?
    # And for each pairing, what is f(S)+f(S^c)?
    pair_stats = defaultdict(list)  # indexed by |S| (for |S| <= |others|/2)

    # Also check: is there a universal decomposition?
    # f(S) = g(S) where g is some simple function of S?

    sample_count = min(1 << n_vars, 200)
    import random
    random.seed(42)
    if sample_count < (1 << n_vars):
        samples = random.sample(range(1 << n_vars), sample_count)
    else:
        samples = list(range(1 << n_vars))

    # Check if f(S)+f(S^c) equals the same value for all S
    universal_pair_count = 0
    total_checked = 0

    for mask in samples:
        vals = {}
        for idx, (a, b) in enumerate(all_vars):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        s = {x: 1 - T(x, I) - T(J, x) for x in others}
        neg_sum_s = -sum(s[x] for x in others)

        # Compute f(S) for each S
        f_vals = {}
        for smask in all_subsets:
            S = subset_from_mask(smask)
            R = subset_from_mask(full_mask ^ smask)

            Si = S + [I]
            Sj = S + [J]
            Rj = [J] + R
            Ri = [I] + R

            he_i = h_end_v(T, Si, I)
            hs_j = h_start_v(T, Rj, J)
            he_j = h_end_v(T, Sj, J)
            hs_i = h_start_v(T, Ri, I)

            f_vals[smask] = he_i * hs_j - he_j * hs_i

        # Check total
        total_f = sum(f_vals.values())
        # This should equal H(T) - H(T')
        ht = sum(1 for p in permutations(range(n))
                 if all(T(p[k], p[k+1]) for k in range(n-1)))

        def Tp(a, b):
            if a == b: return 0
            if a == I and b == J: return 0
            if a == J and b == I: return 1
            return T(a, b)
        htp = sum(1 for p in permutations(range(n))
                  if all(Tp(p[k], p[k+1]) for k in range(n-1)))

        delta_h = ht - htp
        if total_f != delta_h:
            print(f"ERROR: f sum {total_f} != delta_H {delta_h} at mask {mask}")
            return

        # Check pairings
        pair_values = []
        for smask in range(1 << (num_others - 1)):
            # Only iterate over half (|S| <= |S^c|)
            comp = full_mask ^ smask
            pair_val = f_vals[smask] + f_vals[comp]
            pair_values.append(pair_val)

        # Are all pair values the same?
        total_checked += 1
        if len(set(pair_values)) == 1:
            universal_pair_count += 1

        # Record pair values by |S|
        for smask in range(1 << num_others):
            comp = full_mask ^ smask
            if smask > comp:
                continue
            size_s = bin(smask).count('1')
            pair_val = f_vals[smask] + f_vals[comp]
            pair_stats[size_s].append(pair_val)

    print(f"\nAll pairs f(S)+f(S^c) equal for: {universal_pair_count}/{total_checked} cases")

    print(f"\nPair values f(S)+f(S^c) by |S|:")
    for size in sorted(pair_stats.keys()):
        vals = pair_stats[size]
        if not vals:
            continue
        avg = sum(vals)/len(vals)
        minv = min(vals)
        maxv = max(vals)
        print(f"  |S|={size}: avg={avg:.3f}, min={minv}, max={maxv}, "
              f"n={len(vals)}, distinct={len(set(vals))}")


def prove_pairing_n4():
    """Verify the n=4 pairing proof algebraically."""
    print("=== n=4 Pairing Proof (algebraic) ===\n")

    I, J = 0, 1
    others = [2, 3]
    n = 4
    interface = [(x, I) for x in others] + [(J, x) for x in others]
    internal = [(a, b) for a in others for b in others if a < b]
    all_vars = interface + internal
    n_vars = len(all_vars)

    pair_always_equal = True

    for mask in range(1 << n_vars):
        vals = {}
        for idx, (a, b) in enumerate(all_vars):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        s = {x: 1 - T(x, I) - T(J, x) for x in others}
        neg_sum = -sum(s[x] for x in others)

        # f(empty) + f({x,y})
        # f(empty):
        f_empty = h_start_v(T, [J, 2, 3], J) - h_start_v(T, [I, 2, 3], I)
        # f({x,y}) = f({2,3}):
        f_full = h_end_v(T, [2, 3, I], I) - h_end_v(T, [2, 3, J], J)
        # f({2}):
        f_2 = T(2, I) * T(J, 3) - T(2, J) * T(I, 3)
        # f({3}):
        f_3 = T(3, I) * T(J, 2) - T(3, J) * T(I, 2)

        pair1 = f_empty + f_full
        pair2 = f_2 + f_3

        if pair1 != neg_sum or pair2 != neg_sum:
            pair_always_equal = False
            if mask < 5:
                print(f"mask={mask}: pair1={pair1}, pair2={pair2}, -sum(s)={neg_sum}")

        total = f_empty + f_2 + f_3 + f_full
        if total != 2 * neg_sum:
            print(f"ERROR: total {total} != {2*neg_sum}")

    if pair_always_equal:
        print("CONFIRMED: f(S)+f(S^c) = -sum(s_x) for ALL S, ALL tournaments.")
        print("This gives H(T)-H(T') = 2*(-sum(s_x)) = -2*sum(s_x) QED for n=4.\n")


def deeper_analysis_n5():
    """
    At n=5, f(S)+f(S^c) is NOT always the same.
    What is f(S)+f(S^c) equal to?

    Hypothesis: f(S)+f(S^c) depends on the arcs between S and {i,j}
    and between S^c and {i,j}, plus the 5-cycle structure.
    """
    print("=== n=5 Deeper Pairing Analysis ===\n")

    n = 5
    I, J = 0, 1
    others = [2, 3, 4]
    num_others = 3

    interface = [(x, I) for x in others] + [(J, x) for x in others]
    internal = [(a, b) for a in others for b in others if a < b]
    all_vars = interface + internal
    n_vars = len(all_vars)
    full_mask = (1 << num_others) - 1

    def subset_from_mask(smask):
        return [others[bit] for bit in range(num_others) if smask & (1 << bit)]

    # For each tournament, compute f(S)+f(S^c) for each pair
    # and express it in terms of s-values and other quantities
    print("Checking: what determines f(S)+f(S^c)?")

    # Pair types by |S|: |S|=0 pairs with |S|=3, |S|=1 pairs with |S|=2
    # |S|=0: single pair (empty, {a,b,c})
    # |S|=1: three pairs ({a},{b,c}), ({b},{a,c}), ({c},{a,b})

    # For |S|=1, {x} pairs with B_x = others\{x}.
    # Hypothesis: f({x}) + f(B_x) = -sum(s_z) + correction(x)?

    corrections = defaultdict(list)

    for mask in range(1 << n_vars):
        vals = {}
        for idx, (a, b) in enumerate(all_vars):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        s = {x: 1 - T(x, I) - T(J, x) for x in others}
        neg_sum = -sum(s[x] for x in others)

        # Compute all f(S)
        f_vals = {}
        for smask in range(1 << num_others):
            S = subset_from_mask(smask)
            R = subset_from_mask(full_mask ^ smask)
            Si = S + [I]
            Sj = S + [J]
            Ri = [I] + R
            Rj = [J] + R
            f_vals[smask] = (h_end_v(T, Si, I) * h_start_v(T, Rj, J)
                             - h_end_v(T, Sj, J) * h_start_v(T, Ri, I))

        # Pair: f(empty) + f(full)
        pair_0 = f_vals[0] + f_vals[full_mask]

        # Pairs: f({x}) + f(B_x)
        pair_x = {}
        for bit_idx, x in enumerate(others):
            smask_x = 1 << bit_idx
            smask_Bx = full_mask ^ smask_x
            pair_x[x] = f_vals[smask_x] + f_vals[smask_Bx]

        # What is the correction from -sum(s)?
        corr_0 = pair_0 - neg_sum
        for x in others:
            corr_x = pair_x[x] - neg_sum
            corrections[('x', s[x])].append(corr_x)

        corrections[('0', tuple(s[x] for x in others))].append(corr_0)

        # Check: is correction related to 5-cycles through x?
        # A 5-cycle through {i,j} has form (i,j,a,b,c,i) or similar
        # There's only one 5-cycle possible at n=5 (all vertices)

        # Lost 5-cycle: (I,J,...,I) in T
        lost5 = 0
        for perm in permutations(others):
            if T(J, perm[0]) and T(perm[0], perm[1]) and T(perm[1], perm[2]) and T(perm[2], I):
                lost5 += 1
        # Gained 5-cycle: (J,I,...,J) in T'
        gained5 = 0
        for perm in permutations(others):
            if T(I, perm[0]) and T(perm[0], perm[1]) and T(perm[1], perm[2]) and T(perm[2], J):
                gained5 += 1

        net_5 = gained5 - lost5  # should be C5-D5

        # Check if corrections equal -(C5-D5) / 2?
        # Or some multiple of net_5?
        if mask < 5:
            print(f"mask={mask}: s={[s[x] for x in others]}, -sum(s)={neg_sum}, "
                  f"net5={net_5}, pair_0={pair_0}={neg_sum}+{corr_0}, "
                  f"pairs_x={[pair_x[x] for x in others]}")

    print(f"\nCorrection from -sum(s) for boundary pair (|S|=0,3):")
    sig_corrs = defaultdict(list)
    for key, vals in corrections.items():
        if key[0] == '0':
            sig = key[1]
            sig_corrs[sig].extend(vals)

    # Aggregate: is correction constant for given s-signature?
    for sig in sorted(sig_corrs.keys()):
        vals = sig_corrs[sig]
        if len(set(vals)) == 1:
            print(f"  s={sig}: correction={vals[0]} (constant, {len(vals)} cases)")
        else:
            print(f"  s={sig}: corrections={sorted(set(vals))}, NOT constant ({len(vals)} cases)")

    print(f"\nCorrection for vertex pairs (|S|=1,2) by s_x:")
    for sx in [-1, 0, 1]:
        key = ('x', sx)
        if key in corrections:
            vals = corrections[key]
            uniq = sorted(set(vals))
            if len(uniq) <= 5:
                print(f"  s_x={sx:+d}: values={uniq}, counts={[vals.count(v) for v in uniq]}")
            else:
                print(f"  s_x={sx:+d}: {len(uniq)} distinct values, range=[{min(vals)},{max(vals)}]")


def total_f_by_size_n5():
    """
    Decompose sum f(S) by |S| at n=5.

    Is sum_{|S|=k} f(S) a clean function of s-values and cycles?
    """
    print("\n=== Sum of f(S) by |S| at n=5 ===\n")

    n = 5
    I, J = 0, 1
    others = [2, 3, 4]
    num_others = 3
    full_mask = (1 << num_others) - 1

    interface = [(x, I) for x in others] + [(J, x) for x in others]
    internal = [(a, b) for a in others for b in others if a < b]
    all_vars = interface + internal
    n_vars = len(all_vars)

    def subset_from_mask(smask):
        return [others[bit] for bit in range(num_others) if smask & (1 << bit)]

    # Aggregate by s-signature
    sums_by_sig = defaultdict(lambda: defaultdict(list))

    for mask in range(1 << n_vars):
        vals = {}
        for idx, (a, b) in enumerate(all_vars):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        s = {x: 1 - T(x, I) - T(J, x) for x in others}
        sig = tuple(s[x] for x in others)

        # Compute f(S) for each size
        size_sums = defaultdict(int)
        for smask in range(1 << num_others):
            S = subset_from_mask(smask)
            R = subset_from_mask(full_mask ^ smask)
            f_s = (h_end_v(T, S + [I], I) * h_start_v(T, [J] + R, J)
                   - h_end_v(T, S + [J], J) * h_start_v(T, [I] + R, I))
            size_sums[len(S)] += f_s

        for k in range(num_others + 1):
            sums_by_sig[sig][k].append(size_sums[k])

    print("Sum of f(S) by |S|, grouped by s-signature:")
    print(f"{'sig':>12} {'|S|=0':>8} {'|S|=1':>8} {'|S|=2':>8} {'|S|=3':>8} {'total':>8} {'expected':>8}")

    for sig in sorted(sums_by_sig.keys()):
        entry = sums_by_sig[sig]
        # Check if values are constant for this signature
        vals = {k: entry[k] for k in range(num_others + 1)}
        # Get common value (should be constant per signature since it only depends on interface + internal)
        # Actually no, internal arcs vary. Let me check.
        row = []
        for k in range(num_others + 1):
            v = vals[k]
            if len(set(v)) == 1:
                row.append(f"{v[0]:>8}")
            else:
                row.append(f"[{min(v)},{max(v)}]")

        total_v = [sum(vals[k][idx] for k in range(num_others + 1))
                   for idx in range(len(vals[0]))]
        total_str = f"{total_v[0]:>8}" if len(set(total_v)) == 1 else f"[{min(total_v)},{max(total_v)}]"

        neg_sum = -sum(sig)
        # Expected = -2*sum(s) + 2*(C5-D5); but C5-D5 varies with internal arcs
        print(f"{str(sig):>12} {'  '.join(row):>40} {total_str:>8}  (-2*sum(s)={-2*sum(sig)})")


if __name__ == "__main__":
    prove_pairing_n4()
    analyze_pairing(5)
    deeper_analysis_n5()
    total_f_by_size_n5()
