#!/usr/bin/env python3
"""
Attempt to prove D(-1) = 0 for all n using the alternating subset convolution.

Identity: B(L_i, R_j) = B(L_j, R_i)
where B(f,g) = sum_S (-1)^|S| f(S) g(W\S) is the alternating subset convolution.

KEY IDEA: Define phi_v(S) = L_v(S) + R_v(S) for v in {i,j}.
  L_v(S) = h_end(S + {v}, v) = paths on S union {v} ending at v
  R_v(S) = h_start({v} + S, v) = paths on {v} union S starting at v

Note: L_v(S) + R_v(S) counts paths ending at v PLUS paths starting at v.
For |S| = 0: L_v(empty) = 1, R_v(empty) = 1, so phi_v = 2.
For |S| = 1, {w}: L_v({w}) = T[w][v], R_v({w}) = T[v][w] = 1 - T[w][v].
  phi_v({w}) = 1 for ALL v, w! This is the "complementary arc" identity.

For |S| >= 2: phi_v(S) counts the total "boundary" Hamiltonian paths
(either starting or ending at v). Is this symmetric in v?

TEST: Is phi_i(S) = phi_j(S) for all S? If so, B(L_i, R_j) = B(L_j, R_i) follows easily.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
import random


def h_end(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def h_start(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[0] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def test_phi_symmetry(n, num_trials=100):
    """Test: is phi_i(S) = phi_j(S) for all S?"""
    print(f"=== phi_v(S) = L_v(S) + R_v(S) Symmetry Test at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)
    all_sym = True

    for trial in range(min(num_trials, 50)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.random()
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        def T(a, b):
            if a == b:
                return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        for smask in range(1 << m):
            S = [W[bit] for bit in range(m) if smask & (1 << bit)]

            phi_i = h_end(T, S + [I], I) + h_start(T, [I] + S, I)
            phi_j = h_end(T, S + [J], J) + h_start(T, [J] + S, J)

            if abs(phi_i - phi_j) > 1e-8:
                all_sym = False
                if trial == 0:
                    print(f"  S={S}: phi_i={phi_i:.4f}, phi_j={phi_j:.4f}, diff={phi_i-phi_j:.4f}")

    if all_sym:
        print(f"  phi_i(S) = phi_j(S) for ALL S, ALL tournaments: TRUE!")
        print(f"  This would prove D(-1) = 0 trivially!")
    else:
        print(f"  phi_i(S) != phi_j(S) in general. Not symmetric.")


def test_H_boundary(n, num_trials=100):
    """
    For a vertex set S, H_boundary(S, v) = h_end(S+{v}, v) + h_start({v}+S, v)
    counts "boundary-weighted" Hamiltonian paths.

    For |S| = 0: H_boundary = 2 (trivially independent of v)
    For |S| = 1, {w}: H_boundary = T[w][v] + T[v][w] = 1 (independent of v!)

    For |S| = 2, {a,b}: h_end({a,b,v}, v) = T[a][b]*T[b][v] + T[b][a]*T[a][v]
    h_start({v,a,b}, v) = T[v][a]*T[a][b] + T[v][b]*T[b][a]

    For v=i: h_end = T[a][b]*p_b + T[b][a]*p_a, h_start = (1-p_a)*T[a][b] + (1-p_b)*T[b][a]
    Sum = T[a][b]*(p_b + 1 - p_a) + T[b][a]*(p_a + 1 - p_b)
        = T[a][b] + T[b][a] + T[a][b]*(p_b - p_a) + T[b][a]*(p_a - p_b)
        = 1 + (p_b - p_a)*(T[a][b] - T[b][a])
        = 1 + (p_b - p_a)*(2*T[a][b] - 1)

    For v=j: similarly = 1 + ((1-q_b)-(1-q_a))*(2*T[a][b]-1)
        = 1 + (q_a - q_b)*(2*T[a][b]-1)

    These are DIFFERENT unless p_b - p_a = q_a - q_b, i.e., p_a + q_a = p_b + q_b,
    i.e., s_a = s_b. Not always true.

    So phi_i != phi_j for |S| >= 2. The direct phi approach fails.
    """
    print(f"\n=== H_boundary detailed analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    # Check |S| = 2 specifically
    if m >= 2:
        asym_count = 0
        total_count = 0
        for trial in range(min(num_trials, 50)):
            arc_values = {(I, J): 1.0, (J, I): 0.0}
            for (a, b) in free_arcs:
                v = random.random()
                arc_values[(a, b)] = v
                arc_values[(b, a)] = 1 - v

            def T(a, b):
                if a == b:
                    return 0
                return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

            for smask in range(1 << m):
                S = [W[bit] for bit in range(m) if smask & (1 << bit)]
                if len(S) != 2:
                    continue
                total_count += 1

                phi_i = h_end(T, S + [I], I) + h_start(T, [I] + S, I)
                phi_j = h_end(T, S + [J], J) + h_start(T, [J] + S, J)

                if abs(phi_i - phi_j) > 1e-8:
                    asym_count += 1

        print(f"  |S|=2: asymmetric in {asym_count}/{total_count} cases")
    print(f"  Conclusion: phi_v(S) depends on v for |S| >= 2.")
    print(f"  The identity B(L_i, R_j) = B(L_j, R_i) does NOT follow from phi_i = phi_j.")


def test_deeper_identity(n, num_trials=100):
    """
    Even though phi_i != phi_j, the alternating convolution B(L_i, R_j) = B(L_j, R_i)
    still holds. This means the asymmetry in phi cancels out in the alternating sum.

    Test: define Delta_phi(S) = phi_i(S) - phi_j(S). Then:
    D(-1) = B(L_i, R_j) - B(L_j, R_i)
           = sum_S (-1)^|S| [L_i(S)*R_j(R) - L_j(S)*R_i(R)]

    Since L_v = (phi_v - R_v)/1... no, phi_v = L_v + R_v, so L_v = phi_v - R_v.

    B(L_i, R_j) = B(phi_i - R_i, R_j) = B(phi_i, R_j) - B(R_i, R_j)
    B(L_j, R_i) = B(phi_j - R_j, R_i) = B(phi_j, R_i) - B(R_j, R_i)

    Using B(f,g) = (-1)^m B(g,f):
    B(R_i, R_j) = (-1)^m B(R_j, R_i)

    So: D(-1) = B(phi_i, R_j) - B(R_i, R_j) - B(phi_j, R_i) + B(R_j, R_i)
              = B(phi_i, R_j) - B(phi_j, R_i) - B(R_i, R_j) + B(R_j, R_i)

    If m is even: B(R_i,R_j) = B(R_j,R_i), so last two terms cancel.
    Then D(-1) = B(phi_i, R_j) - B(phi_j, R_i).

    If m is odd: B(R_i,R_j) = -B(R_j,R_i), so last two terms = -2*B(R_i,R_j).
    Then D(-1) = B(phi_i, R_j) - B(phi_j, R_i) - 2*B(R_i, R_j).

    Hmm, this doesn't simplify cleanly. Let me try yet another decomposition.

    Actually, let me test: is B(R_i, R_j) = B(R_j, R_i) always? (Not just for even m)
    """
    print(f"\n=== Deeper Identity Tests at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 20)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.random()
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        def T(a, b):
            if a == b:
                return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        # Compute B(f, g) for various f, g
        def compute_B(f, g):
            total = 0.0
            for smask in range(1 << m):
                S = [W[bit] for bit in range(m) if smask & (1 << bit)]
                R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
                total += ((-1) ** len(S)) * f(S) * g(R)
            return total

        def Li(S):
            return h_end(T, S + [I], I)
        def Lj(S):
            return h_end(T, S + [J], J)
        def Ri(S):
            return h_start(T, [I] + S, I)
        def Rj(S):
            return h_start(T, [J] + S, J)
        def phi_i(S):
            return Li(S) + Ri(S)
        def phi_j(S):
            return Lj(S) + Rj(S)

        B_Li_Rj = compute_B(Li, Rj)
        B_Lj_Ri = compute_B(Lj, Ri)
        B_Ri_Rj = compute_B(Ri, Rj)
        B_Rj_Ri = compute_B(Rj, Ri)
        B_Li_Lj = compute_B(Li, Lj)
        B_Lj_Li = compute_B(Lj, Li)

        # B(f,g) = (-1)^m * B(g,f)
        check_sym = abs(B_Ri_Rj - ((-1)**m) * B_Rj_Ri)

        # The key identity
        D = B_Li_Rj - B_Lj_Ri  # should be 0

        if trial < 5:
            print(f"  trial {trial}: D={D:.2e}, "
                  f"B(Li,Rj)={B_Li_Rj:.4f}, B(Lj,Ri)={B_Lj_Ri:.4f}, "
                  f"B(Ri,Rj)={B_Ri_Rj:.4f}, B(Rj,Ri)={B_Rj_Ri:.4f}, "
                  f"B_sym_check={check_sym:.2e}")

        # Test: is B(Li, Rj) + B(Ri, Lj) = B(phi_i, Rj)? (YES by linearity)
        # Test: B(Li, Rj) + B(Li, Lj) = B(Li, phi_j)? (YES)

        # CRITICAL TEST: does B(Li, Lj) = B(Ri, Rj)?
        check_LR = abs(B_Li_Lj - B_Ri_Rj)
        if trial < 5:
            print(f"    B(Li,Lj)={B_Li_Lj:.4f}, B(Ri,Rj)={B_Ri_Rj:.4f}, "
                  f"B(Li,Lj)-B(Ri,Rj)={check_LR:.2e}")

        # Test: B(Li, Ri) = B(Lj, Rj)?  [alternating self-convolution]
        B_Li_Ri = compute_B(Li, Ri)
        B_Lj_Rj = compute_B(Lj, Rj)
        check_self = abs(B_Li_Ri - B_Lj_Rj)
        if trial < 5:
            print(f"    B(Li,Ri)={B_Li_Ri:.4f}, B(Lj,Rj)={B_Lj_Rj:.4f}, "
                  f"diff={check_self:.2e}")

    print()
    print("  Legend: B(f,g) = sum_S (-1)^|S| f(S) g(W\\S)")
    print("  Identity to prove: B(Li, Rj) = B(Lj, Ri) [equivalently D = 0]")


if __name__ == "__main__":
    test_phi_symmetry(4)
    test_phi_symmetry(5)
    test_H_boundary(5)
    test_deeper_identity(4)
    test_deeper_identity(5)
