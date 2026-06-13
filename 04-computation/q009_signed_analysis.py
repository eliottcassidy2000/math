#!/usr/bin/env python3
"""
Signed adjacency analysis for the even-odd split identity.

The even-odd split identity sum_S (-1)^|S| Delta(S,R) = 0 is equivalent to:

  A_signed = B_signed

where:
  B_signed = sum over T-paths using i->j of (-1)^{# W-vertices before i}
  A_signed = sum over T'-paths using j->i of (-1)^{# W-vertices before j}

The swap involution pairs T-paths with T'-paths, preserving sign.
Matched paths cancel. The identity reduces to:

  sum (signed unmatched T-paths) = sum (signed unmatched T'-paths)

Key structural finding:
- Unmatched T-paths are blocked by type-C vertices (p_w=1, q_w=1, s_w=-1)
- Unmatched T'-paths are blocked by type-D vertices (p_w=0, q_w=0, s_w=+1)

This script analyzes the signed unmatched structure.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
from collections import defaultdict


def analyze_signed_structure(n, num_trials=None):
    """Analyze signed unmatched path counts."""
    print(f"=== Signed Unmatched Analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    # Enumerate arc variables
    arc_vars = []
    for a in W:
        arc_vars.append((a, I))  # p_a = T[a][i]
    for a in W:
        arc_vars.append((J, a))  # q_a = T[j][a]
    for idx_a, a in enumerate(W):
        for b in W[idx_a+1:]:
            arc_vars.append((a, b))  # internal
    n_vars = len(arc_vars)

    total_configs = 1 << n_vars
    if num_trials is None or num_trials >= total_configs:
        configs = range(total_configs)
        num_trials = total_configs
    else:
        import random
        random.seed(42)
        configs = random.sample(range(total_configs), num_trials)

    print(f"Variables: {n_vars}, Configs: {num_trials}")

    all_pass = 0
    signed_data = []

    for mask in configs:
        vals = {}
        for idx, (a, b) in enumerate(arc_vars):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        # Vertex types
        p = {w: T(w, I) for w in W}
        q = {w: T(J, w) for w in W}
        s = {w: 1 - p[w] - q[w] for w in W}

        # Classify vertices
        types = {}
        for w in W:
            if p[w] == 1 and q[w] == 0:
                types[w] = 'A'  # w->i, w->j
            elif p[w] == 0 and q[w] == 1:
                types[w] = 'B'  # i->w, j->w
            elif p[w] == 1 and q[w] == 1:
                types[w] = 'C'  # w->i, j->w (s=-1)
            else:
                types[w] = 'D'  # i->w, w->j (s=+1)

        # Find T-paths using i->j
        B_signed = 0
        B_unmatched_signed = 0
        B_unmatched_by_pos = defaultdict(int)
        B_unmatched_by_blocker = defaultdict(int)

        for perm in permutations(range(n)):
            # Check if i->j is consecutive
            ij_pos = None
            for k in range(n-1):
                if perm[k] == I and perm[k+1] == J:
                    ij_pos = k
                    break
            if ij_pos is None:
                continue

            # Check all arcs valid in T
            valid = True
            for k in range(n-1):
                if not T(perm[k], perm[k+1]):
                    valid = False
                    break
            if not valid:
                continue

            # Count W-vertices before i
            pos = sum(1 for k in range(ij_pos) if perm[k] in W)
            sign = (-1) ** pos
            B_signed += sign

            # Check if swap works
            # Swap: replace i->j with j->i, i.e., create perm' with j,i instead of i,j
            perm_swap = list(perm)
            perm_swap[ij_pos] = J
            perm_swap[ij_pos + 1] = I

            # Check if perm_swap is valid in T' (T with i<->j flipped)
            def Tp(a, b):
                if a == I and b == J: return 0
                if a == J and b == I: return 1
                return T(a, b)

            swap_valid = True
            for k in range(n-1):
                if not Tp(perm_swap[k], perm_swap[k+1]):
                    swap_valid = False
                    break

            if not swap_valid:
                B_unmatched_signed += sign
                B_unmatched_by_pos[pos] += sign
                # Identify blocker
                if ij_pos > 0:
                    prev_v = perm[ij_pos - 1]
                    if prev_v in W and types[prev_v] == 'C':
                        B_unmatched_by_blocker['prev_C'] += sign
                if ij_pos + 2 < n:
                    next_v = perm[ij_pos + 2]
                    if next_v in W and types[next_v] == 'C':
                        B_unmatched_by_blocker['next_C'] += sign

        # Find T'-paths using j->i
        A_signed = 0
        A_unmatched_signed = 0
        A_unmatched_by_pos = defaultdict(int)

        for perm in permutations(range(n)):
            ji_pos = None
            for k in range(n-1):
                if perm[k] == J and perm[k+1] == I:
                    ji_pos = k
                    break
            if ji_pos is None:
                continue

            def Tp(a, b):
                if a == I and b == J: return 0
                if a == J and b == I: return 1
                return T(a, b)

            valid = True
            for k in range(n-1):
                if not Tp(perm[k], perm[k+1]):
                    valid = False
                    break
            if not valid:
                continue

            pos = sum(1 for k in range(ji_pos) if perm[k] in W)
            sign = (-1) ** pos
            A_signed += sign

            # Check if swap works (j,i -> i,j)
            perm_swap = list(perm)
            perm_swap[ji_pos] = I
            perm_swap[ji_pos + 1] = J

            swap_valid = True
            for k in range(n-1):
                if not T(perm_swap[k], perm_swap[k+1]):
                    swap_valid = False
                    break

            if not swap_valid:
                A_unmatched_signed += sign
                A_unmatched_by_pos[pos] += sign

        match = (A_signed == B_signed)
        unmatched_match = (A_unmatched_signed == B_unmatched_signed)
        if match:
            all_pass += 1

        type_sig = ''.join(types[w] for w in W)
        signed_data.append({
            'A_signed': A_signed,
            'B_signed': B_signed,
            'A_unmatched': A_unmatched_signed,
            'B_unmatched': B_unmatched_signed,
            'type_sig': type_sig,
            's_values': tuple(s[w] for w in W),
        })

    print(f"\nA_signed == B_signed: {all_pass}/{num_trials}")

    # Aggregate statistics
    print(f"\nSigned unmatched analysis:")
    unmatched_match_count = sum(1 for d in signed_data if d['A_unmatched'] == d['B_unmatched'])
    print(f"A_unmatched_signed == B_unmatched_signed: {unmatched_match_count}/{num_trials}")

    # By type signature
    by_sig = defaultdict(list)
    for d in signed_data:
        by_sig[d['type_sig']].append(d)

    print(f"\nBy vertex type signature:")
    for sig in sorted(by_sig.keys()):
        entries = by_sig[sig]
        a_vals = [d['A_signed'] for d in entries]
        b_vals = [d['B_signed'] for d in entries]
        au_vals = [d['A_unmatched'] for d in entries]
        bu_vals = [d['B_unmatched'] for d in entries]
        match_count = sum(1 for d in entries if d['A_signed'] == d['B_signed'])
        print(f"  {sig}: n={len(entries)}, A_signed range=[{min(a_vals)},{max(a_vals)}], "
              f"match={match_count}/{len(entries)}, "
              f"A_um range=[{min(au_vals)},{max(au_vals)}], B_um range=[{min(bu_vals)},{max(bu_vals)}]")

    # KEY TEST: when all vertices are type A or B (s=0 for all), what happens?
    neutral = [d for d in signed_data if all(v == 0 for v in d['s_values'])]
    if neutral:
        print(f"\nAll-neutral (s=0 for all W-vertices): {len(neutral)} configs")
        for d in neutral[:5]:
            print(f"  types={d['type_sig']}: A_signed={d['A_signed']}, B_signed={d['B_signed']}, "
                  f"A_um={d['A_unmatched']}, B_um={d['B_unmatched']}")
        # When s=0 for all, there are NO type C or D vertices, so NO blockers
        # => all paths are matched => unmatched = 0 => A_signed = B_signed trivially!
        all_zero_um = all(d['A_unmatched'] == 0 and d['B_unmatched'] == 0 for d in neutral)
        print(f"  All unmatched = 0: {all_zero_um}")

    # When exactly one C or D vertex
    one_cd = [d for d in signed_data if sum(1 for v in d['s_values'] if v != 0) == 1]
    if one_cd:
        print(f"\nExactly one C/D vertex: {len(one_cd)} configs")
        a_eq_b = sum(1 for d in one_cd if d['A_signed'] == d['B_signed'])
        print(f"  A_signed == B_signed: {a_eq_b}/{len(one_cd)}")

    return signed_data


def test_signed_unmatched_structure(n, num_trials=None):
    """
    More detailed analysis: decompose unmatched paths by which vertex blocks.

    For each unmatched T-path ... v_k -> i -> j -> v_{k+1} ...,
    the blocker is v_k (if type C, blocking the v_k->j connection)
    or v_{k+1} (if type C, blocking the i->v_{k+1} connection).
    """
    print(f"\n=== Blocker Decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    arc_vars = []
    for a in W:
        arc_vars.append((a, I))
    for a in W:
        arc_vars.append((J, a))
    for idx_a, a in enumerate(W):
        for b in W[idx_a+1:]:
            arc_vars.append((a, b))
    n_vars = len(arc_vars)

    total_configs = 1 << n_vars
    if num_trials is None or num_trials >= total_configs:
        configs = range(total_configs)
        num_trials = total_configs
    else:
        import random
        random.seed(42)
        configs = random.sample(range(total_configs), num_trials)

    # Track: for each blocker vertex w, what is the signed contribution?
    # B_unmatched = sum_w B_blocked_by_w
    # where B_blocked_by_w = sum over paths blocked at w, (-1)^pos

    per_vertex_match = 0

    for mask in configs:
        vals = {}
        for idx, (a, b) in enumerate(arc_vars):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        p = {w: T(w, I) for w in W}
        q = {w: T(J, w) for w in W}

        # For T-paths using i->j, decompose unmatched by blocker vertex
        B_by_w = {w: 0 for w in W}
        for perm in permutations(range(n)):
            ij_pos = None
            for k in range(n-1):
                if perm[k] == I and perm[k+1] == J:
                    ij_pos = k
                    break
            if ij_pos is None:
                continue
            valid = all(T(perm[k], perm[k+1]) for k in range(n-1))
            if not valid:
                continue

            pos = sum(1 for k in range(ij_pos) if perm[k] in W)
            sign = (-1) ** pos

            # Check swap
            blocked_by_prev = (ij_pos > 0 and perm[ij_pos-1] in W and q[perm[ij_pos-1]] == 1)
            blocked_by_next = (ij_pos + 2 < n and perm[ij_pos+2] in W and p[perm[ij_pos+2]] == 1)

            if blocked_by_prev and not blocked_by_next:
                B_by_w[perm[ij_pos-1]] += sign
            elif blocked_by_next and not blocked_by_prev:
                B_by_w[perm[ij_pos+2]] += sign
            elif blocked_by_prev and blocked_by_next:
                # Both block — attribute to prev (or we could split)
                B_by_w[perm[ij_pos-1]] += sign  # convention: attribute to left blocker

        # For T'-paths using j->i, decompose unmatched by blocker vertex
        A_by_w = {w: 0 for w in W}
        for perm in permutations(range(n)):
            ji_pos = None
            for k in range(n-1):
                if perm[k] == J and perm[k+1] == I:
                    ji_pos = k
                    break
            if ji_pos is None:
                continue

            def Tp(a, b):
                if a == I and b == J: return 0
                if a == J and b == I: return 1
                return T(a, b)

            valid = all(Tp(perm[k], perm[k+1]) for k in range(n-1))
            if not valid:
                continue

            pos = sum(1 for k in range(ji_pos) if perm[k] in W)
            sign = (-1) ** pos

            # For T'-path, swap creates T-path with i->j
            # Swap blocked when:
            # prev vertex u_k can't connect to i: p[u_k] = 0 (type D)
            # next vertex u_{k+1} can't be reached from j: q[u_{k+1}] = 0 (type D)
            blocked_by_prev = (ji_pos > 0 and perm[ji_pos-1] in W and p[perm[ji_pos-1]] == 0)
            blocked_by_next = (ji_pos + 2 < n and perm[ji_pos+2] in W and q[perm[ji_pos+2]] == 0)

            if blocked_by_prev and not blocked_by_next:
                A_by_w[perm[ji_pos-1]] += sign
            elif blocked_by_next and not blocked_by_prev:
                A_by_w[perm[ji_pos+2]] += sign
            elif blocked_by_prev and blocked_by_next:
                A_by_w[perm[ji_pos-1]] += sign

        # Check per-vertex matching
        B_total = sum(B_by_w.values())
        A_total = sum(A_by_w.values())
        per_vertex = all(B_by_w[w] == A_by_w[w] for w in W)
        if per_vertex:
            per_vertex_match += 1

        if mask < 3 or (not per_vertex and mask < 50):
            print(f"mask={mask}: B_by_w={dict(B_by_w)}, A_by_w={dict(A_by_w)}, "
                  f"B_total={B_total}, A_total={A_total}, per_vertex={per_vertex}")

    print(f"\nPer-vertex matching: {per_vertex_match}/{num_trials}")
    print("(If per-vertex always matches, the identity decomposes vertex-by-vertex)")


def test_position_parity_structure(n, num_trials=None):
    """
    Test whether the signed count has a clean form.

    Hypothesis: A_signed = B_signed = some function of the type signature.
    """
    print(f"\n=== Position-Parity Structure at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    arc_vars = []
    for a in W:
        arc_vars.append((a, I))
    for a in W:
        arc_vars.append((J, a))
    for idx_a, a in enumerate(W):
        for b in W[idx_a+1:]:
            arc_vars.append((a, b))
    n_vars = len(arc_vars)

    total_configs = 1 << n_vars
    if num_trials is None or num_trials >= total_configs:
        configs = range(total_configs)
        num_trials = total_configs
    else:
        import random
        random.seed(42)
        configs = random.sample(range(total_configs), num_trials)

    # Track: A_signed as function of type signature and internal arcs
    by_typesig = defaultdict(list)

    for mask in configs:
        vals = {}
        for idx, (a, b) in enumerate(arc_vars):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        p = {w: T(w, I) for w in W}
        q = {w: T(J, w) for w in W}

        # Signed adj count for T-paths using i->j
        B_signed = 0
        for perm in permutations(range(n)):
            ij_pos = None
            for k in range(n-1):
                if perm[k] == I and perm[k+1] == J:
                    ij_pos = k
                    break
            if ij_pos is None:
                continue
            if all(T(perm[k], perm[k+1]) for k in range(n-1)):
                pos = sum(1 for k in range(ij_pos) if perm[k] in W)
                B_signed += (-1) ** pos

        type_sig = tuple((p[w], q[w]) for w in W)
        by_typesig[type_sig].append(B_signed)

    print(f"Distinct type signatures: {len(by_typesig)}")
    print(f"\nIs B_signed determined by type signature? (same for all internal arcs)")

    determined = 0
    not_determined = 0
    for sig, vals in sorted(by_typesig.items()):
        if len(set(vals)) == 1:
            determined += 1
        else:
            not_determined += 1
            if not_determined <= 5:
                print(f"  sig={sig}: values={sorted(set(vals))}, counts={[vals.count(v) for v in sorted(set(vals))]}")

    print(f"\nDetermined by type sig: {determined}/{len(by_typesig)}")
    print(f"NOT determined: {not_determined}/{len(by_typesig)}")

    if not_determined == 0:
        print("\nB_signed IS determined by type signature alone!")
        print("This means it depends only on (p_w, q_w) and not on internal arcs.")
        print("This would make the identity much easier to prove.")
    else:
        print("\nB_signed depends on internal arcs too.")
        print("The identity must involve global cancellation across internal arc configurations.")


if __name__ == "__main__":
    # n=4: 32 configs (exhaustive)
    analyze_signed_structure(4)

    # n=5: 512 configs (exhaustive)
    analyze_signed_structure(5)

    # Blocker decomposition
    test_signed_unmatched_structure(4)
    test_signed_unmatched_structure(5)

    # Position-parity structure
    test_position_parity_structure(4)
    test_position_parity_structure(5)

    # n=6: sample
    print("\n" + "="*60)
    test_position_parity_structure(6, num_trials=500)
