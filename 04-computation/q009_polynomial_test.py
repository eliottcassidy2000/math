#!/usr/bin/env python3
"""
Test whether D(-1) = 0 holds as a POLYNOMIAL identity at n=5.

If it holds for non-binary values of arc variables, it's a polynomial identity.
If it only holds over {0,1}, we need Boolean-specific arguments.

This is the KEY discriminant for choosing a proof strategy.

Also tests n=3 and n=4 as sanity checks.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
import random


def compute_D_neg1(n, arc_values):
    """
    Compute D(-1) = sum_S (-1)^|S| [L_i(S)*R_j(R) - L_j(S)*R_i(R)]
    using real-valued arc variables.

    arc_values: dict (a,b) -> real value for T[a][b]
    """
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    def T(a, b):
        if a == b: return 0
        if (a, b) in arc_values: return arc_values[(a, b)]
        return 1 - arc_values[(b, a)]

    total = 0.0
    for smask in range(1 << m):
        S = [W[bit] for bit in range(m) if smask & (1 << bit)]
        R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
        sign = (-1) ** len(S)

        # L_i(S) = h_end(S + {i}, i)
        L_i = h_end(T, S + [I], I)
        # L_j(S) = h_end(S + {j}, j)
        L_j = h_end(T, S + [J], J)
        # R_i(R) = h_start({i} + R, i)
        R_i = h_start(T, [I] + R, I)
        # R_j(R) = h_start({j} + R, j)
        R_j = h_start(T, [J] + R, J)

        Delta = L_j * R_i - L_i * R_j
        total += sign * Delta

    return total


def h_end(T, verts, v):
    """Sum of path weights on verts ending at v."""
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        weight = 1.0
        for k in range(len(p) - 1):
            weight *= T(p[k], p[k+1])
        total += weight
    return total


def h_start(T, verts, v):
    """Sum of path weights on verts starting at v."""
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[0] != v:
            continue
        weight = 1.0
        for k in range(len(p) - 1):
            weight *= T(p[k], p[k+1])
        total += weight
    return total


def test_polynomial_identity(n, num_trials=100):
    """Test D(-1) = 0 with random real-valued arc variables."""
    print(f"=== Testing D(-1) = 0 as polynomial identity at n={n} ===")

    I, J = 0, 1
    W = list(range(2, n))

    # Build list of independent arc variables
    arc_pairs = [(I, J)]  # fixed to 1
    free_arcs = []
    for a in range(n):
        for b in range(a+1, n):
            if (a, b) == (I, J):
                continue
            free_arcs.append((a, b))

    random.seed(42)
    max_err = 0.0

    for trial in range(num_trials):
        # Random real values in [0, 1]
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.random()  # uniform [0,1]
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        d = compute_D_neg1(n, arc_values)
        max_err = max(max_err, abs(d))

        if abs(d) > 1e-8 and trial < 5:
            print(f"  Trial {trial}: D(-1) = {d:.6e} (NON-ZERO!)")

    if max_err < 1e-8:
        print(f"  All {num_trials} trials: D(-1) < {max_err:.2e} -- POLYNOMIAL IDENTITY")
    else:
        print(f"  Max |D(-1)| = {max_err:.6e} -- NOT a polynomial identity")

    # Also test with values outside [0,1]
    max_err2 = 0.0
    for trial in range(num_trials):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-2, 3)  # outside [0,1]
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        d = compute_D_neg1(n, arc_values)
        max_err2 = max(max_err2, abs(d))

    if max_err2 < 1e-6:
        print(f"  Extended range [-2,3]: max |D(-1)| = {max_err2:.2e} -- STRONG POLYNOMIAL IDENTITY")
    else:
        print(f"  Extended range [-2,3]: max |D(-1)| = {max_err2:.6e} -- FAILS outside [0,1]")

    return max_err < 1e-8


def test_D_at_various_x(n, num_trials=50):
    """Test what values of x make D(x) = 0."""
    print(f"\n=== D(x) = 0 locus at n={n} ===")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    free_arcs = [(a, b) for a in range(n) for b in range(a+1, n) if (a,b) != (I,J)]

    random.seed(42)

    test_x_values = [-2, -1, -0.5, 0, 0.5, 1, 2]

    for x in test_x_values:
        max_D = 0.0
        for trial in range(num_trials):
            arc_values = {(I, J): 1.0, (J, I): 0.0}
            for (a, b) in free_arcs:
                v = random.uniform(-1, 2)
                arc_values[(a, b)] = v
                arc_values[(b, a)] = 1 - v

            def T(a, b):
                if a == b: return 0
                return arc_values.get((a,b), 1 - arc_values.get((b,a), 0))

            # Compute D(x) = sum_k D[k] * x^k where D[k] = F[k] - G[k]
            D_coeffs = [0.0] * (m + 1)
            for perm in permutations(range(n)):
                # Check i->j
                for pos in range(n-1):
                    if perm[pos] == I and perm[pos+1] == J:
                        weight = 1.0
                        for k in range(n-1):
                            weight *= T(perm[k], perm[k+1])
                        nw = sum(1 for r in range(pos) if perm[r] in W)
                        D_coeffs[nw] += weight
                        break
                # Check j->i (in T')
                for pos in range(n-1):
                    if perm[pos] == J and perm[pos+1] == I:
                        # T' weight: same as T except T'[j][i]=1, T'[i][j]=0
                        weight = 1.0
                        for k in range(n-1):
                            a, b = perm[k], perm[k+1]
                            if a == J and b == I:
                                weight *= 1.0  # T'[j][i] = 1
                            elif a == I and b == J:
                                weight *= 0.0  # T'[i][j] = 0
                            else:
                                weight *= T(a, b)
                        nw = sum(1 for r in range(pos) if perm[r] in W)
                        D_coeffs[nw] -= weight
                        break

            D_at_x = sum(D_coeffs[k] * (x ** k) for k in range(m+1))
            max_D = max(max_D, abs(D_at_x))

        status = "= 0" if max_D < 1e-6 else f"!= 0 (max={max_D:.4f})"
        print(f"  x={x:>5.1f}: D(x) {status}")


if __name__ == "__main__":
    test_polynomial_identity(3)
    test_polynomial_identity(4)
    test_polynomial_identity(5)

    test_D_at_various_x(4)
    test_D_at_various_x(5)
