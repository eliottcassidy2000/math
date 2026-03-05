#!/usr/bin/env python3
"""
Symbolic proof of delta_H = delta_I via polynomial identity.

At n=4, I proved by hand that U_{T'} - U_T = 4 - 2(p+q+r+t) = 2*sum(s_x),
which equals delta_I. This script extends the approach to n=5 and n=6.

Method:
1. Enumerate all Ham path "shapes" using arc (i,j) consecutively
2. For each shape, write validity and blocking as polynomial in arc variables
3. Sum to get U_T, apply i<->j symmetry for U_T', compute difference
4. Compare with delta_I formula from THM-013

Instance: kind-pasteur-2026-03-05-S6
"""

from itertools import permutations
from collections import defaultdict


def arc_var(a, b, n):
    """Return a symbolic name for the arc a->b."""
    return f"T{a}{b}"


def prove_n4():
    """Prove the identity at n=4 symbolically."""
    print("=== n=4 Symbolic Proof ===\n")
    # Vertices: i=0, j=1, others={2,3}
    # Arc variables: T[x][0] for x in {2,3}, T[1][x] for x in {2,3}, T[2][3]
    # Using: p_x = T[x][0] (x beats i), q_x = T[1][x] (j beats x)
    # s_x = 1 - p_x - q_x

    # All Ham paths using 0->1 consecutively:
    # Position of (0,1) in path of length 4: pos 0-1, 1-2, 2-3
    # Remaining vertices {2,3} fill the other 2 spots

    i, j = 0, 1
    others = [2, 3]
    n = 4

    # Generate all path shapes
    shapes = []
    for pos in range(n - 1):  # position of i in the path
        # i at position pos, j at position pos+1
        remaining_positions = [k for k in range(n) if k != pos and k != pos + 1]
        for perm in permutations(others):
            path = [None] * n
            path[pos] = i
            path[pos + 1] = j
            for idx, rp in enumerate(remaining_positions):
                path[rp] = perm[idx]
            shapes.append(tuple(path))

    print(f"Total path shapes: {len(shapes)}")

    # For each shape, compute validity polynomial and blocking polynomial
    # Validity: product of T[path[k]][path[k+1]] for k != pos (the i->j arc is guaranteed)
    # Blocking: at least one blocker (pred or succ of i,j)

    # Represent polynomials as dicts: monomial (frozenset of arc vars) -> coefficient
    # An arc T[a][b] is represented; T[b][a] = 1 - T[a][b]

    # For simplicity, represent each arc as a variable or (1-variable)
    # Variables: T20, T30, T12, T13, T23 (and their complements)

    # Actually let me just evaluate numerically over all 2^5 = 32 assignments
    # and verify the identity holds for all.

    vars_list = ['T20', 'T30', 'T12', 'T13', 'T23']

    def T(a, b, vals):
        """Get T[a][b] from variable assignment."""
        if a == 0 and b == 1:
            return 1  # arc i->j
        if a == 1 and b == 0:
            return 0
        key = f"T{a}{b}"
        if key in vals:
            return vals[key]
        # Complement
        comp_key = f"T{b}{a}"
        if comp_key in vals:
            return 1 - vals[comp_key]
        raise ValueError(f"Unknown arc {a}->{b}")

    def Tp(a, b, vals):
        """Get T'[a][b] (flipped arc)."""
        if a == 0 and b == 1:
            return 0  # flipped
        if a == 1 and b == 0:
            return 1  # flipped
        return T(a, b, vals)

    ok = 0
    total = 0

    for mask in range(32):
        vals = {}
        for idx, v in enumerate(vars_list):
            vals[v] = (mask >> idx) & 1

        # Compute U_T
        uT = 0
        for path in shapes:
            # Check validity (all arcs exist in T)
            valid = True
            for k in range(n - 1):
                if not T(path[k], path[k + 1], vals):
                    valid = False
                    break
            if not valid:
                continue

            # Check blocking
            pos = path.index(i)
            assert path[pos + 1] == j
            blocked = False
            if pos > 0:
                pred = path[pos - 1]
                if not T(pred, j, vals):  # pred doesn't beat j
                    blocked = True
            if pos + 1 < n - 1:
                succ = path[pos + 2]
                if not T(i, succ, vals):  # i doesn't beat succ
                    blocked = True
            if blocked:
                uT += 1

        # Compute U_T' (paths of T' using j->i, blocked)
        uTp = 0
        shapes_p = []
        for pos in range(n - 1):
            remaining_positions = [k for k in range(n) if k != pos and k != pos + 1]
            for perm in permutations(others):
                path = [None] * n
                path[pos] = j  # j before i in T'
                path[pos + 1] = i
                for idx, rp in enumerate(remaining_positions):
                    path[rp] = perm[idx]
                shapes_p.append(tuple(path))

        for path in shapes_p:
            valid = True
            for k in range(n - 1):
                if not Tp(path[k], path[k + 1], vals):
                    valid = False
                    break
            if not valid:
                continue

            pos = path.index(j)
            assert path[pos + 1] == i
            blocked = False
            if pos > 0:
                pred = path[pos - 1]
                if not Tp(pred, i, vals):
                    blocked = True
            if pos + 1 < n - 1:
                succ = path[pos + 2]
                if not Tp(j, succ, vals):
                    blocked = True
            if blocked:
                uTp += 1

        # Compute delta_I formula: -2*(s_2 + s_3)
        s2 = 1 - T(2, 0, vals) - T(1, 2, vals)
        s3 = 1 - T(3, 0, vals) - T(1, 3, vals)
        delta_I = -2 * (s2 + s3)  # THM-013 convention: H(T)-H(T')

        # Our convention: delta = U_T' - U_T = H(T') - H(T) = -delta_I(THM-013)
        delta = uTp - uT
        expected = 2 * (s2 + s3)  # = -delta_I(THM-013) = H(T')-H(T)

        total += 1
        if delta == expected:
            ok += 1
        else:
            print(f"  FAIL: vals={vals}, uT={uT}, uTp={uTp}, delta={delta}, expected={expected}")

    print(f"Identity verified: {ok}/{total}")
    if ok == total:
        print("PROVED: U_T' - U_T = 2*sum(s_x) for ALL n=4 tournaments.")
        print("This is equivalent to delta_H = delta_I at n=4.")


def prove_n5():
    """Prove the identity at n=5 symbolically."""
    print("\n=== n=5 Symbolic Proof ===\n")
    i, j = 0, 1
    others = [2, 3, 4]
    n = 5

    # Variables: T[x][0], T[1][x] for x in {2,3,4}, T[2][3], T[2][4], T[3][4]
    # Total: 6 + 3 = 9 variables, 2^9 = 512 assignments
    vars_list = ['T20', 'T30', 'T40', 'T12', 'T13', 'T14', 'T23', 'T24', 'T34']

    def T_val(a, b, vals):
        if a == 0 and b == 1: return 1
        if a == 1 and b == 0: return 0
        key = f"T{a}{b}"
        if key in vals: return vals[key]
        return 1 - vals[f"T{b}{a}"]

    def Tp_val(a, b, vals):
        if a == 0 and b == 1: return 0
        if a == 1 and b == 0: return 1
        return T_val(a, b, vals)

    # Generate T-shapes (using i->j)
    shapes_T = []
    for pos in range(n - 1):
        remaining_positions = [k for k in range(n) if k != pos and k != pos + 1]
        for perm in permutations(others):
            path = [None] * n
            path[pos] = i
            path[pos + 1] = j
            for idx, rp in enumerate(remaining_positions):
                path[rp] = perm[idx]
            shapes_T.append(tuple(path))

    # Generate T'-shapes (using j->i)
    shapes_Tp = []
    for pos in range(n - 1):
        remaining_positions = [k for k in range(n) if k != pos and k != pos + 1]
        for perm in permutations(others):
            path = [None] * n
            path[pos] = j
            path[pos + 1] = i
            for idx, rp in enumerate(remaining_positions):
                path[rp] = perm[idx]
            shapes_Tp.append(tuple(path))

    print(f"Path shapes: {len(shapes_T)} for T, {len(shapes_Tp)} for T'")

    ok = 0
    total = 0

    for mask in range(512):
        vals = {}
        for idx, v in enumerate(vars_list):
            vals[v] = (mask >> idx) & 1

        # U_T
        uT = 0
        for path in shapes_T:
            valid = all(T_val(path[k], path[k+1], vals) for k in range(n-1))
            if not valid:
                continue
            pos = path.index(i)
            blocked = False
            if pos > 0 and not T_val(path[pos-1], j, vals):
                blocked = True
            if pos+1 < n-1 and not T_val(i, path[pos+2], vals):
                blocked = True
            if blocked:
                uT += 1

        # U_T'
        uTp = 0
        for path in shapes_Tp:
            valid = all(Tp_val(path[k], path[k+1], vals) for k in range(n-1))
            if not valid:
                continue
            pos = path.index(j)
            blocked = False
            if pos > 0 and not Tp_val(path[pos-1], i, vals):
                blocked = True
            if pos+1 < n-1 and not Tp_val(j, path[pos+2], vals):
                blocked = True
            if blocked:
                uTp += 1

        # delta_I from THM-013: at n=5, delta_I = 2*sum(DL-CL) for all L
        # = 2*(D3-C3) + 2*(D5-C5) = -2*sum(s_x) + 2*(D5-C5)
        # But simpler: 2*(gained_cycles - lost_cycles)

        # Compute s values
        s_vals = {}
        for x in others:
            s_vals[x] = 1 - T_val(x, 0, vals) - T_val(1, x, vals)

        # D5: 5-cycles in T' using j->i = (1,0,v1,v2,v3,1)
        # C5: 5-cycles in T using i->j = (0,1,v1,v2,v3,0)
        D5 = 0
        C5 = 0
        for perm in permutations(others):
            v1, v2, v3 = perm
            # Lost 5-cycle: (0,1,v1,v2,v3,0)
            if (T_val(1, v1, vals) and T_val(v1, v2, vals) and
                T_val(v2, v3, vals) and T_val(v3, 0, vals)):
                C5 += 1
            # Gained 5-cycle: (1,0,v1,v2,v3,1)
            if (T_val(0, v1, vals) and T_val(v1, v2, vals) and
                T_val(v2, v3, vals) and T_val(v3, 1, vals)):
                D5 += 1

        # Wait: gained means exists in T' but not T.
        # (1,0,v1,v2,v3,1) in T': needs T'[1][0]=1 ✓, T'[0][v1]=T[0][v1],
        # T[v1][v2], T[v2][v3], T'[v3][1]=T[v3][1].
        # This uses arc 1->0 which only exists in T'. So it's a gained cycle.
        # But we need to check if arcs 0->v1, v1->v2, v2->v3, v3->1 all exist in T.
        # Actually T'[0][v1] = T[0][v1] since 0,v1 are not the flipped pair.
        # So: gained 5-cycle needs T[0][v1]=1, T[v1][v2]=1, T[v2][v3]=1, T[v3][1]=1.

        # Lost 5-cycle (0,1,v1,v2,v3,0): uses 0->1 in T.
        # Needs T[1][v1]=1, T[v1][v2]=1, T[v2][v3]=1, T[v3][0]=1.

        # So my computation above looks right. Let me recalculate:
        C5 = 0
        D5 = 0
        for perm in permutations(others):
            v1, v2, v3 = perm
            # Lost: (0,1,v1,v2,v3,0) needs T[1][v1], T[v1][v2], T[v2][v3], T[v3][0]
            if (T_val(1, v1, vals) and T_val(v1, v2, vals) and
                T_val(v2, v3, vals) and T_val(v3, 0, vals)):
                C5 += 1
            # Gained: (1,0,v1,v2,v3,1) needs T[0][v1], T[v1][v2], T[v2][v3], T[v3][1]
            if (T_val(0, v1, vals) and T_val(v1, v2, vals) and
                T_val(v2, v3, vals) and T_val(v3, 1, vals)):
                D5 += 1

        sum_s = sum(s_vals[x] for x in others)
        expected = 2 * sum_s + 2 * (D5 - C5)  # H(T')-H(T) convention

        delta = uTp - uT
        total += 1
        if delta == expected:
            ok += 1
        else:
            print(f"  FAIL at mask={mask}: delta={delta}, expected={expected}, "
                  f"sum_s={sum_s}, D5-C5={D5-C5}")
            if total - ok > 5:
                print("  Too many failures, stopping.")
                break

    print(f"Identity verified: {ok}/{total}")
    if ok == total:
        print("PROVED: U_T' - U_T = 2*sum(s_x) + 2*(D5-C5) for ALL n=5 tournaments.")
        print("This proves delta_H = delta_I at n=5, hence OCF at n=5.")


def prove_n6():
    """Prove the identity at n=6 symbolically."""
    print("\n=== n=6 Symbolic Proof ===\n")
    i, j = 0, 1
    others = [2, 3, 4, 5]
    n = 6

    # Variables: T[x][0], T[1][x] for x in {2,3,4,5}, T[2][3], T[2][4], T[2][5], T[3][4], T[3][5], T[4][5]
    # Total: 8 + 6 = 14 variables, 2^14 = 16384 assignments
    vars_list = ['T20', 'T30', 'T40', 'T50',
                 'T12', 'T13', 'T14', 'T15',
                 'T23', 'T24', 'T25', 'T34', 'T35', 'T45']

    def T_val(a, b, vals):
        if a == 0 and b == 1: return 1
        if a == 1 and b == 0: return 0
        key = f"T{a}{b}"
        if key in vals: return vals[key]
        return 1 - vals[f"T{b}{a}"]

    def Tp_val(a, b, vals):
        if a == 0 and b == 1: return 0
        if a == 1 and b == 0: return 1
        return T_val(a, b, vals)

    # Generate shapes (using i->j for T, j->i for T')
    shapes_T = []
    shapes_Tp = []
    for pos in range(n - 1):
        remaining_positions = [k for k in range(n) if k != pos and k != pos + 1]
        for perm in permutations(others):
            path_t = [None] * n
            path_t[pos] = i
            path_t[pos + 1] = j
            path_tp = [None] * n
            path_tp[pos] = j
            path_tp[pos + 1] = i
            for idx, rp in enumerate(remaining_positions):
                path_t[rp] = perm[idx]
                path_tp[rp] = perm[idx]
            shapes_T.append(tuple(path_t))
            shapes_Tp.append(tuple(path_tp))

    print(f"Path shapes: {len(shapes_T)} each")
    print(f"Variable assignments to check: {2**len(vars_list)}")

    ok = 0
    total = 0
    fail_count = 0

    for mask in range(2**len(vars_list)):
        vals = {}
        for idx, v in enumerate(vars_list):
            vals[v] = (mask >> idx) & 1

        # U_T
        uT = 0
        for path in shapes_T:
            valid = True
            for k in range(n-1):
                if not T_val(path[k], path[k+1], vals):
                    valid = False
                    break
            if not valid:
                continue
            pos = path.index(i)
            blocked = False
            if pos > 0 and not T_val(path[pos-1], j, vals):
                blocked = True
            if pos+1 < n-1 and not T_val(i, path[pos+2], vals):
                blocked = True
            if blocked:
                uT += 1

        # U_T'
        uTp = 0
        for path in shapes_Tp:
            valid = True
            for k in range(n-1):
                if not Tp_val(path[k], path[k+1], vals):
                    valid = False
                    break
            if not valid:
                continue
            pos = path.index(j)
            blocked = False
            if pos > 0 and not Tp_val(path[pos-1], i, vals):
                blocked = True
            if pos+1 < n-1 and not Tp_val(j, path[pos+2], vals):
                blocked = True
            if blocked:
                uTp += 1

        # THM-013 formula for n=6:
        # delta_I = -2*sum(s_x*H(B_x)) + 2*(D5-C5)
        # H(T')-H(T) = 2*sum(s_x*H(B_x)) - 2*(D5-C5)
        # Wait, need to be careful with sign convention.
        # THM-013: DeltaH = H(T) - H(T') = sum 2^k Delta(alpha_k)
        # Our delta: U_T' - U_T = H(T') - H(T)

        # Compute using full cycle counting
        s_vals = {}
        for x in others:
            s_vals[x] = 1 - T_val(x, 0, vals) - T_val(1, x, vals)

        # H(B_x) for 3-vertex sub-tournaments
        def h3(verts, vals):
            """H of 3-vertex tournament."""
            a, b, c = verts
            # Count Ham paths
            count = 0
            for p in permutations(verts):
                if T_val(p[0], p[1], vals) and T_val(p[1], p[2], vals):
                    count += 1
            return count

        # 5-cycle counts
        C5 = 0  # lost 5-cycles (using 0->1)
        D5 = 0  # gained 5-cycles (using 1->0)
        for perm in permutations(others):
            v1, v2, v3, v4 = perm[:4]  # wait, others has 4 elements for 5-cycle
            # A 5-cycle uses 5 vertices. At n=6, a 5-cycle uses {0,1} + 3 of {2,3,4,5}
            pass

        # Actually 5-cycles at n=6: use 5 of 6 vertices.
        # A 5-cycle using i->j: (0, 1, v1, v2, v3, 0) with {v1,v2,v3} subset of others
        from itertools import combinations
        C5 = 0
        D5 = 0
        for subset in combinations(others, 3):
            for perm in permutations(subset):
                v1, v2, v3 = perm
                # Lost: (0,1,v1,v2,v3,0)
                if (T_val(1, v1, vals) and T_val(v1, v2, vals) and
                    T_val(v2, v3, vals) and T_val(v3, 0, vals)):
                    C5 += 1
                # Gained: (1,0,v1,v2,v3,1)
                if (T_val(0, v1, vals) and T_val(v1, v2, vals) and
                    T_val(v2, v3, vals) and T_val(v3, 1, vals)):
                    D5 += 1

        # THM-013 formula (H(T')-H(T) convention):
        # = 2*sum(s_x * H(B_x)) - 2*(D5-C5) ... no wait.
        # THM-013: H(T)-H(T') = -2*sum(s_x*H(B_x)) + 2*(D5-C5)
        # So H(T')-H(T) = 2*sum(s_x*H(B_x)) - 2*(D5-C5)

        formula_sum = 0
        for x in others:
            B_x = [v for v in others if v != x]
            hbx = h3(B_x, vals)
            formula_sum += s_vals[x] * hbx

        # THM-013 convention: D=destroyed(lost), C=created(gained)
        # My code: C5=lost, D5=gained. So THM-013's (D5-C5) = my (C5-D5).
        # H(T')-H(T) = 2*sum(s_x*H(B_x)) + 2*(D5-C5) [my D/C convention]
        expected = 2 * formula_sum + 2 * (D5 - C5)

        delta = uTp - uT
        total += 1
        if delta == expected:
            ok += 1
        else:
            fail_count += 1
            if fail_count <= 3:
                print(f"  FAIL at mask={mask}: delta={delta}, expected={expected}")
            if fail_count > 10:
                print(f"  Too many failures ({fail_count}), aborting n=6.")
                break

        if total % 4000 == 0:
            print(f"  Progress: {total}/{2**len(vars_list)}, ok={ok}")

    print(f"Identity verified: {ok}/{total}")
    if ok == total:
        print("PROVED: delta_H = delta_I for ALL n=6 tournaments (16384 cases).")
        print("Combined with base case (transitive), this proves OCF at n=6.")


if __name__ == "__main__":
    prove_n4()
    prove_n5()
    # n=6 takes longer, run separately if needed
    import sys
    if '--n6' in sys.argv:
        prove_n6()
