#!/usr/bin/env python3
"""
Analyze the structural decomposition of U_T' - U_T = delta_I.

Goal: understand WHY the polynomial identity holds, to guide a general proof.

Key decomposition:
- U_T = pred-blocked + succ-blocked - double-blocked
- Each blocked-T path (...,x,i,j,y,...) requires q_x=1 (pred) or p_y=1 (succ)
- Each blocked-T' path (...,x,j,i,y,...) requires p_x=0 (pred) or q_y=0 (succ)

The pred contribution only involves s_x = +/-1 vertices (s=0 cancels).
Same for succ.

This script decomposes the identity into pred/succ/both components
and compares with the THM-013 formula terms.

Instance: kind-pasteur-2026-03-05-S7
"""

import sys
from itertools import permutations, combinations
from collections import defaultdict


def analyze_n5():
    """Detailed structural analysis at n=5."""
    print("=== n=5 Structural Decomposition ===\n")

    n = 5
    I, J = 0, 1
    others = [2, 3, 4]

    # Variables: p_x = T[x][I], q_x = T[J][x] for x in others (6 vars)
    # Plus internal arcs: T[a][b] for a<b in others (3 vars)
    # Total: 9 vars
    interface = [(x, I) for x in others] + [(J, x) for x in others]
    internal = [(a, b) for a in others for b in others if a < b]
    all_vars = interface + internal
    n_vars = len(all_vars)
    assert n_vars == 9

    var_idx = {}
    for idx, (a, b) in enumerate(all_vars):
        var_idx[(a, b)] = idx

    total_eq = 0
    pred_total = 0
    succ_total = 0
    both_total = 0

    # Track stats by s-signature
    stats_by_sig = defaultdict(lambda: {'count': 0, 'delta': 0, 'pred': 0, 'succ': 0, 'both': 0,
                                        'formula': 0, 'cycle': 0})

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

        def Tp(a, b):
            if a == b: return 0
            if a == I and b == J: return 0
            if a == J and b == I: return 1
            return T(a, b)

        # s_x values
        s = {x: 1 - T(x, I) - T(J, x) for x in others}
        p = {x: T(x, I) for x in others}
        q = {x: T(J, x) for x in others}
        sig = tuple(s[x] for x in others)

        # Enumerate all Ham paths
        def enum_paths(Tf, n):
            paths = []
            for perm in permutations(range(n)):
                ok = True
                for k in range(n - 1):
                    if not Tf(perm[k], perm[k+1]):
                        ok = False
                        break
                if ok:
                    paths.append(perm)
            return paths

        T_paths = enum_paths(T, n)
        Tp_paths = enum_paths(Tp, n)

        # Count unmatched paths by component
        def count_unmatched_components(paths, Tf, src, dst, is_T):
            """Count pred-blocked, succ-blocked, both-blocked."""
            pred_count = 0
            succ_count = 0
            both_count = 0
            total_adj = 0
            for path in paths:
                for k in range(n - 1):
                    if path[k] == src and path[k+1] == dst:
                        total_adj += 1
                        pred_blocked = False
                        succ_blocked = False
                        if k > 0:
                            x = path[k - 1]
                            if is_T:
                                # T-path: blocked if T[x][j]=0
                                if T(x, J) == 0:
                                    pred_blocked = True
                            else:
                                # T'-path: blocked if T[x][i]=0
                                if T(x, I) == 0:
                                    pred_blocked = True
                        if k + 1 < n - 1:
                            y = path[k + 2]
                            if is_T:
                                # T-path: blocked if T[i][y]=0
                                if T(I, y) == 0:
                                    succ_blocked = True
                            else:
                                # T'-path: blocked if T[j][y]=0
                                if T(J, y) == 0:
                                    succ_blocked = True
                        if pred_blocked:
                            pred_count += 1
                        if succ_blocked:
                            succ_count += 1
                        if pred_blocked and succ_blocked:
                            both_count += 1
                        break
            unmatched = pred_count + succ_count - both_count
            return total_adj, pred_count, succ_count, both_count, unmatched

        _, pT, sT, bT, uT = count_unmatched_components(T_paths, T, I, J, True)
        _, pTp, sTp, bTp, uTp = count_unmatched_components(Tp_paths, Tp, J, I, False)

        delta = uTp - uT
        delta_pred = pTp - pT
        delta_succ = sTp - sT
        delta_both = bTp - bT

        # THM-013 formula
        formula_sum = 0
        for x in others:
            Bx = [v for v in others if v != x]
            hbx = sum(1 for p in permutations(Bx) if all(T(p[k], p[k+1]) for k in range(len(p)-1)))
            formula_sum += s[x] * hbx

        # 5-cycle counts
        lost5 = 0
        gained5 = 0
        for perm in permutations(others):
            v1, v2, v3 = perm
            if T(J, v1) and T(v1, v2) and T(v2, v3) and T(v3, I):
                lost5 += 1
            if T(I, v1) and T(v1, v2) and T(v2, v3) and T(v3, J):
                gained5 += 1

        expected = 2 * formula_sum + 2 * (gained5 - lost5)

        if delta != expected:
            print(f"MISMATCH at mask={mask}!")
            return

        entry = stats_by_sig[sig]
        entry['count'] += 1
        entry['delta'] += delta
        entry['pred'] += delta_pred
        entry['succ'] += delta_succ
        entry['both'] += delta_both
        entry['formula'] += 2 * formula_sum
        entry['cycle'] += 2 * (gained5 - lost5)
        total_eq += 1

    print(f"All {total_eq}/512 cases match.\n")
    print("Decomposition by s-signature (s_2, s_3, s_4):")
    print(f"{'sig':>12} {'count':>5} {'avg_delta':>10} {'avg_pred':>10} {'avg_succ':>10} "
          f"{'avg_both':>10} {'avg_form':>10} {'avg_cyc':>10}")

    for sig in sorted(stats_by_sig.keys()):
        e = stats_by_sig[sig]
        c = e['count']
        print(f"{str(sig):>12} {c:>5} {e['delta']/c:>10.2f} {e['pred']/c:>10.2f} "
              f"{e['succ']/c:>10.2f} {e['both']/c:>10.2f} "
              f"{e['formula']/c:>10.2f} {e['cycle']/c:>10.2f}")

    print("\n--- Key observation ---")
    print("Does delta_pred = delta_succ? (by symmetry of the identity)")

    sym_count = 0
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

        def Tp(a, b):
            if a == b: return 0
            if a == I and b == J: return 0
            if a == J and b == I: return 1
            return T(a, b)

        T_paths = enum_paths(T, n)
        Tp_paths = enum_paths(Tp, n)
        _, pT, sT, bT, uT = count_unmatched_components(T_paths, T, I, J, True)
        _, pTp, sTp, bTp, uTp = count_unmatched_components(Tp_paths, Tp, J, I, False)

        if pTp - pT == sTp - sT:
            sym_count += 1

    print(f"delta_pred == delta_succ in {sym_count}/512 cases")


def analyze_blocking_by_vertex_n5():
    """Check: which vertices cause blocking, and how does it relate to s_x?"""
    print("\n=== Blocking by vertex at n=5 ===\n")

    n = 5
    I, J = 0, 1
    others = [2, 3, 4]
    interface = [(x, I) for x in others] + [(J, x) for x in others]
    internal = [(a, b) for a in others for b in others if a < b]
    all_vars = interface + internal
    n_vars = len(all_vars)

    # Check: for s_x = 0, does x ever cause blocking?
    s0_blocks = 0
    s0_total = 0
    sm1_blocks = 0
    sp1_blocks = 0

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

        for perm in permutations(range(n)):
            ok = True
            for k in range(n - 1):
                if not T(perm[k], perm[k + 1]):
                    ok = False
                    break
            if not ok:
                continue

            for k in range(n - 1):
                if perm[k] == I and perm[k + 1] == J:
                    # Check pred-blocking
                    if k > 0:
                        x = perm[k - 1]
                        if x in others and T(x, J) == 0:
                            # x pred-blocks this T-path
                            if s[x] == 0:
                                s0_blocks += 1
                            elif s[x] == -1:
                                sm1_blocks += 1
                    # Check succ-blocking
                    if k + 1 < n - 1:
                        y = perm[k + 2]
                        if y in others and T(I, y) == 0:
                            if s[y] == 0:
                                s0_blocks += 1
                            elif s[y] == -1:
                                sm1_blocks += 1
                    break

        # Count T'-paths blocked by s=0, s=+1
        def Tp(a, b):
            if a == b: return 0
            if a == I and b == J: return 0
            if a == J and b == I: return 1
            return T(a, b)

        for perm in permutations(range(n)):
            ok = True
            for k in range(n - 1):
                if not Tp(perm[k], perm[k + 1]):
                    ok = False
                    break
            if not ok:
                continue

            for k in range(n - 1):
                if perm[k] == J and perm[k + 1] == I:
                    if k > 0:
                        x = perm[k - 1]
                        if x in others and T(x, I) == 0:
                            if s[x] == 0:
                                s0_blocks += 1
                            elif s[x] == 1:
                                sp1_blocks += 1
                    if k + 1 < n - 1:
                        y = perm[k + 2]
                        if y in others and T(J, y) == 0:
                            if s[y] == 0:
                                s0_blocks += 1
                            elif s[y] == 1:
                                sp1_blocks += 1
                    break

    print(f"s=0 vertices cause blocking in {s0_blocks} path instances")
    print(f"s=-1 vertices cause T-blocking in {sm1_blocks} instances")
    print(f"s=+1 vertices cause T'-blocking in {sp1_blocks} instances")

    # Check THM-014 claim: blocking iff s=-1 (for T), s=+1 (for T')
    print(f"\nTHM-014 claims s=0 never blocks. Actually: {'CORRECT' if s0_blocks==0 else 'WRONG ('+str(s0_blocks)+')'}")


def analyze_factored_counts_n5():
    """
    Compute PC_x(i), PC_x(j) and compare with H(B_x).

    Key question: is there a simple relationship between
    PC_x(i) + PC_x(j) and H(B_x)?
    """
    print("\n=== Factored adjacency counts at n=5 ===\n")

    n = 5
    I, J = 0, 1
    others = [2, 3, 4]
    interface = [(x, I) for x in others] + [(J, x) for x in others]
    internal = [(a, b) for a in others for b in others if a < b]
    all_vars = interface + internal
    n_vars = len(all_vars)

    # For a subset of variable assignments, compute PC_x(I), PC_x(J), H(B_x)
    import random
    random.seed(42)
    samples = random.sample(range(1 << n_vars), min(100, 1 << n_vars))

    relations = defaultdict(list)

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
        p = {x: T(x, I) for x in others}
        q = {x: T(J, x) for x in others}

        def h_paths(verts, Tf=T):
            """Count Ham paths on vertex set."""
            if len(verts) <= 1:
                return 1
            count = 0
            for perm in permutations(verts):
                if all(Tf(perm[k], perm[k+1]) for k in range(len(perm)-1)):
                    count += 1
            return count

        def h_end(verts, v, Tf=T):
            """Count paths on verts ending at v."""
            if len(verts) == 1:
                return 1 if v in verts else 0
            count = 0
            for perm in permutations(verts):
                if perm[-1] != v:
                    continue
                if all(Tf(perm[k], perm[k+1]) for k in range(len(perm)-1)):
                    count += 1
            return count

        def h_start(verts, v, Tf=T):
            """Count paths on verts starting at v."""
            if len(verts) == 1:
                return 1 if v in verts else 0
            count = 0
            for perm in permutations(verts):
                if perm[0] != v:
                    continue
                if all(Tf(perm[k], perm[k+1]) for k in range(len(perm)-1)):
                    count += 1
            return count

        for x in others:
            Bx = [v for v in others if v != x]
            hBx = h_paths(Bx)

            # PC_x(I) = sum_S h_end(S+{x}, x) * h_start({I}+(Bx\S), I)
            pc_i = 0
            pc_j = 0
            for r in range(len(Bx) + 1):
                for S in combinations(Bx, r):
                    S_set = set(S)
                    rest = [v for v in Bx if v not in S_set]
                    he = h_end(list(S) + [x], x)
                    hs_i = h_start([I] + rest, I)
                    hs_j = h_start([J] + rest, J)
                    pc_i += he * hs_i
                    pc_j += he * hs_j

            relations[s[x]].append((pc_i, pc_j, hBx, p[x], q[x]))

    print("Relationships between PC_x(I), PC_x(J), H(B_x) grouped by s_x:\n")
    for sx in [-1, 0, 1]:
        if not relations[sx]:
            continue
        data = relations[sx]
        print(f"s_x = {sx:+d}: ({len(data)} samples)")
        # Check if PC_x(I) or PC_x(J) = H(B_x)
        i_eq = sum(1 for d in data if d[0] == d[2])
        j_eq = sum(1 for d in data if d[1] == d[2])
        ij_sum_eq = sum(1 for d in data if d[0] + d[1] == d[2])
        print(f"  PC_x(I) == H(B_x): {i_eq}/{len(data)}")
        print(f"  PC_x(J) == H(B_x): {j_eq}/{len(data)}")
        print(f"  PC_x(I) + PC_x(J) == H(B_x): {ij_sum_eq}/{len(data)}")

        # Check ratios
        if data:
            ratios_i = [d[0]/d[2] if d[2] > 0 else None for d in data]
            ratios_j = [d[1]/d[2] if d[2] > 0 else None for d in data]
            valid_i = [r for r in ratios_i if r is not None]
            valid_j = [r for r in ratios_j if r is not None]
            if valid_i:
                print(f"  PC_x(I)/H(B_x): min={min(valid_i):.3f}, max={max(valid_i):.3f}, "
                      f"mean={sum(valid_i)/len(valid_i):.3f}")
            if valid_j:
                print(f"  PC_x(J)/H(B_x): min={min(valid_j):.3f}, max={max(valid_j):.3f}, "
                      f"mean={sum(valid_j)/len(valid_j):.3f}")

        # For s=-1: the term q_x*p_x*PC_x(J) contributes to U_T^pred
        # We need: contribution = PC_x(J) (since q_x=p_x=1)
        # For s=+1: the term (1-p_x)(1-q_x)*PC_x(I) contributes to U_{T'}^pred
        # We need: contribution = PC_x(I) (since p_x=q_x=0)
        if sx == -1:
            print(f"  These contribute PC_x(J) to U_T^pred (need to subtract)")
        elif sx == 1:
            print(f"  These contribute PC_x(I) to U_T'^pred (need to add)")
        print()


def verify_pred_succ_symmetry_n4():
    """
    At n=4, verify the full decomposition of U_T' - U_T.

    n=4: others = {a, b}. B_x = single vertex, H(B_x) = 1.
    Expected: U_T' - U_T = 2*(s_a + s_b).
    """
    print("=== n=4 Full Decomposition ===\n")

    n = 4
    I, J = 0, 1
    others = [2, 3]

    # 5 variables: p_2, p_3, q_2, q_3, T[2][3]
    vars_list = [(2, I), (3, I), (J, 2), (J, 3), (2, 3)]
    n_vars = len(vars_list)

    for mask in range(1 << n_vars):
        vals = {}
        for idx, (a, b) in enumerate(vars_list):
            vals[(a, b)] = (mask >> idx) & 1

        def T(a, b):
            if a == b: return 0
            if a == I and b == J: return 1
            if a == J and b == I: return 0
            if (a, b) in vals: return vals[(a, b)]
            return 1 - vals[(b, a)]

        def Tp(a, b):
            if a == b: return 0
            if a == I and b == J: return 0
            if a == J and b == I: return 1
            return T(a, b)

        p = {x: T(x, I) for x in others}
        q = {x: T(J, x) for x in others}
        s = {x: 1 - p[x] - q[x] for x in others}

        # Enumerate all paths
        def get_paths(Tf):
            paths = []
            for perm in permutations(range(n)):
                if all(Tf(perm[k], perm[k+1]) for k in range(n-1)):
                    paths.append(perm)
            return paths

        T_paths = get_paths(T)
        Tp_paths = get_paths(Tp)

        # Decompose unmatched
        def decompose(paths, src, dst, is_T):
            pred_b = succ_b = both_b = 0
            for path in paths:
                for k in range(n - 1):
                    if path[k] == src and path[k + 1] == dst:
                        pb = sb = False
                        if k > 0:
                            x = path[k - 1]
                            if is_T and T(x, J) == 0:
                                pb = True
                            elif not is_T and T(x, I) == 0:
                                pb = True
                        if k + 1 < n - 1:
                            y = path[k + 2]
                            if is_T and T(I, y) == 0:
                                sb = True
                            elif not is_T and T(J, y) == 0:
                                sb = True
                        if pb: pred_b += 1
                        if sb: succ_b += 1
                        if pb and sb: both_b += 1
                        break
            return pred_b, succ_b, both_b

        pT, sT, bT = decompose(T_paths, I, J, True)
        pTp, sTp, bTp = decompose(Tp_paths, J, I, False)

        uT = pT + sT - bT
        uTp = pTp + sTp - bTp
        delta = uTp - uT

        expected = 2 * sum(s[x] for x in others)

        if delta != expected:
            print(f"FAIL at mask={mask}")
            return

        if mask < 4:
            print(f"mask={mask}: s={[s[x] for x in others]}, "
                  f"pred=({pTp}-{pT}={pTp-pT}), succ=({sTp}-{sT}={sTp-sT}), "
                  f"both=({bTp}-{bT}={bTp-bT}), delta={delta}, expected={expected}")

    print(f"\nAll {1 << n_vars} cases match for n=4.")
    print("Note: pred and succ components individually may not equal expected/2,")
    print("but their sum (minus double-counting) always equals 2*sum(s_x).")


if __name__ == "__main__":
    verify_pred_succ_symmetry_n4()
    analyze_blocking_by_vertex_n5()
    analyze_n5()
    analyze_factored_counts_n5()
