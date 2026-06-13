#!/usr/bin/env python3
"""
Test the inductive structure of the D(-1) = 0 identity.

The identity: sum_{S⊆W} (-1)^|S| [L_j(S)*R_i(R) - L_i(S)*R_j(R)] = 0
where R = W\S, i=0, j=1, W = {2,...,n-1}.

Key reformulation:
  B(L_i, R_j) = B(L_j, R_i)  [alternating subset convolution symmetry]

where B(f, g) = sum_{S⊆W} (-1)^|S| f(S) g(W\S).

APPROACH: Test whether the identity has an inductive structure by:
1. Decomposing by a fixed vertex w (condition on w ∈ S vs w ∈ R)
2. Checking if each piece relates to the (n-1)-vertex identity

Also test a key structural observation:
  L_i(S) = sum_{w∈S} p_w * E_w(S)     [paths on S ending at w, then w→i]
  L_j(S) = sum_{w∈S} (1-q_w) * E_w(S) [paths on S ending at w, then w→j]

This means L_i - L_j = -sum s_w * E_w(S), which factors the "i↔j difference."

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
from collections import defaultdict
import random


def h_end_weighted(T, verts, v):
    """Ham path weight on verts ending at v (real-valued arcs)."""
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


def h_start_weighted(T, verts, v):
    """Ham path weight on verts starting at v (real-valued arcs)."""
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


def test_vertex_decomposition(n, num_trials=100):
    """
    Test: fix vertex w_max = n-1. Split D(-1) into:
      D(-1) = D_L(w on left) + D_R(w on right)

    Does each piece relate to the (n-1)-case identity?
    """
    print(f"=== Vertex Decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    w_max = n - 1
    W_prime = [w for w in W if w != w_max]
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

        p_w = T(w_max, I)
        q_w = T(J, w_max)
        s_w = 1 - p_w - q_w

        # Compute D(-1) split by w_max position
        D_L = 0.0  # w_max ∈ S (left)
        D_R = 0.0  # w_max ∈ R (right)

        for smask in range(1 << m):
            S = [W[bit] for bit in range(m) if smask & (1 << bit)]
            R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
            sign = (-1) ** len(S)

            Li = h_end_weighted(T, S + [I], I)
            Lj = h_end_weighted(T, S + [J], J)
            Ri = h_start_weighted(T, [I] + R, I)
            Rj = h_start_weighted(T, [J] + R, J)

            Delta = Lj * Ri - Li * Rj
            term = sign * Delta

            if w_max in S:
                D_L += term
            else:
                D_R += term

        # Also compute the (n-1)-case D(-1) for the sub-tournament on {0,1} ∪ W'
        D_sub = 0.0
        m_prime = len(W_prime)
        for smask in range(1 << m_prime):
            S = [W_prime[bit] for bit in range(m_prime) if smask & (1 << bit)]
            R = [W_prime[bit] for bit in range(m_prime) if not (smask & (1 << bit))]
            sign = (-1) ** len(S)

            Li = h_end_weighted(T, S + [I], I)
            Lj = h_end_weighted(T, S + [J], J)
            Ri = h_start_weighted(T, [I] + R, I)
            Rj = h_start_weighted(T, [J] + R, J)

            Delta = Lj * Ri - Li * Rj
            D_sub += sign * Delta

        if trial < 10:
            print(f"  trial {trial}: D_L={D_L:.6f}, D_R={D_R:.6f}, "
                  f"D_L+D_R={D_L + D_R:.2e}, D_sub(n-1)={D_sub:.2e}, "
                  f"s_w={s_w:.3f}, D_L/D_R={'N/A' if abs(D_R) < 1e-10 else f'{D_L / D_R:.4f}'}")

    # Check: is D_L = -D_R always? (yes, since D_L + D_R = 0)
    print(f"\n  D_L + D_R = 0: always (this is D(-1) = 0)")
    print(f"  D_sub(n-1) = 0: always (this is the (n-1) identity)")
    print(f"\n  Question: does D_L relate to D_sub? Or to s_w * something?")


def test_factored_structure(n, num_trials=100):
    """
    Test: does L_j(S)*R_i(R) - L_i(S)*R_j(R) factor in a useful way?

    We know:
      L_i(S) = sum_{w∈S} p_w * E_w(S)
      L_j(S) = sum_{w∈S} (1-q_w) * E_w(S)

    So L_i(S) = sum p_w E_w, L_j(S) = sum (1-q_w) E_w.

    Similarly:
      R_i(R) = sum_{w∈R} (1-p_w) * ST_w(R)
      R_j(R) = sum_{w∈R} q_w * ST_w(R)

    where E_w(S) = h_end(S, w) and ST_w(R) = h_start(R, w).

    Delta = L_j*R_i - L_i*R_j
          = [sum (1-q_u) E_u][sum (1-p_v) ST_v] - [sum p_u E_u][sum q_v ST_v]
          = sum_{u,v} [(1-q_u)(1-p_v) - p_u*q_v] E_u(S) * ST_v(R)

    Let c(u,v) = (1-q_u)(1-p_v) - p_u*q_v = 1 - q_u - p_v + q_u*p_v - p_u*q_v
    Then: Delta(S,R) = sum_{u∈S, v∈R} c(u,v) * E_u(S) * ST_v(R)

    And D(-1) = sum_S (-1)^|S| sum_{u∈S, v∈R} c(u,v) * E_u(S) * ST_v(R) = 0

    Can we swap the order of summation?
    """
    print(f"\n=== Factored Structure at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    # For each pair (u,v) with u,v ∈ W, u ≠ v, compute
    # C(u,v) = sum_{S: u∈S, v∉S} (-1)^|S| c(u,v) E_u(S) ST_v(W\S)
    # D(-1) = sum_{u≠v, u,v∈W} C(u,v)

    for trial in range(min(num_trials, 10)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.random()
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        def T(a, b):
            if a == b:
                return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        # Compute c(u,v) for each pair
        c_uv = {}
        for u in W:
            for v in W:
                if u == v:
                    continue
                pu = T(u, I)
                qu = T(J, u)
                pv = T(v, I)
                qv = T(J, v)
                c_uv[(u, v)] = (1 - qu) * (1 - pv) - pu * qv

        # For each (u,v) pair, compute the partial alternating sum
        C_pair = {}
        for u in W:
            for v in W:
                if u == v:
                    continue
                # Sum over S containing u but not v
                W_rest = [w for w in W if w != u and w != v]
                C_val = 0.0
                for rest_mask in range(1 << len(W_rest)):
                    S_rest = [W_rest[bit] for bit in range(len(W_rest)) if rest_mask & (1 << bit)]
                    S = S_rest + [u]
                    R_rest = [W_rest[bit] for bit in range(len(W_rest)) if not (rest_mask & (1 << bit))]
                    R = R_rest + [v]

                    sign = (-1) ** len(S)  # |S| = |S_rest| + 1

                    # E_u(S) = h_end(S, u) [paths on S ending at u, internal arcs only]
                    E_u = h_end_weighted(T, S, u)
                    # ST_v(R) = h_start(R, v) [paths on R starting at v, internal arcs only]
                    ST_v = h_start_weighted(T, R, v)

                    C_val += sign * c_uv[(u, v)] * E_u * ST_v

                C_pair[(u, v)] = C_val

        total = sum(C_pair.values())

        if trial < 5:
            print(f"  trial {trial}: total = {total:.2e}")
            # Print pair contributions
            for u in W:
                for v in W:
                    if u != v:
                        cv = C_pair[(u, v)]
                        if abs(cv) > 1e-10:
                            print(f"    C({u},{v}) = {cv:.6f}, c = {c_uv[(u, v)]:.4f}")

        # Check: for each (u,v), is C(u,v) + C(v,u) = 0?
        antisym_ok = True
        for u in W:
            for v in W:
                if u < v:
                    pair_sum = C_pair[(u, v)] + C_pair[(v, u)]
                    if abs(pair_sum) > 1e-8:
                        antisym_ok = False
                        if trial == 0:
                            print(f"    C({u},{v})+C({v},{u}) = {pair_sum:.6f} != 0")

        if trial < 5:
            print(f"    C(u,v) + C(v,u) = 0 for all pairs: {antisym_ok}")


def test_recursion_via_endpoint(n, num_trials=50):
    """
    Test: for a fixed last vertex w_k of the left path (just before i or j),
    does the remaining sum have (n-1)-case structure?

    L_i(S) = sum_{w∈S} p_w * E_w(S)

    So: sum_S (-1)^|S| L_i(S) R_j(R) = sum_S (-1)^|S| [sum_{w∈S} p_w E_w(S)] R_j(R)
      = sum_w p_w * [sum_{S∋w} (-1)^|S| E_w(S) R_j(W\S)]

    For a fixed w, the inner sum is:
      sum_{S_0 ⊆ W\{w}} (-1)^{|S_0|+1} E_w(S_0∪{w}) R_j(W\(S_0∪{w}))
    = -sum_{S_0 ⊆ W\{w}} (-1)^{|S_0|} E_w(S_0∪{w}) R_j((W\{w})\S_0)

    Similarly for the j-terms. So D(-1) can be written as:
    D(-1) = sum_w [-(1-q_w)*A_w + p_w*B_w] + ...

    where A_w and B_w are alternating sums over the (n-1)-case.

    Test if this recursive structure leads anywhere.
    """
    print(f"\n=== Endpoint Recursion at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 10)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.random()
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        def T(a, b):
            if a == b:
                return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        # For each w ∈ W, compute:
        #   A_w = sum_{S_0 ⊆ W\{w}} (-1)^|S_0| E_w(S_0∪{w}) R_j((W\{w})\S_0)
        #   B_w = sum_{S_0 ⊆ W\{w}} (-1)^|S_0| E_w(S_0∪{w}) R_i((W\{w})\S_0)
        # Then contribution of w as left-endpoint:
        #   sum_S (-1)^|S| [L_j * R_i - L_i * R_j]
        #   where L_j has w as endpoint: (1-q_w) * E_w(S) * R_i(R)
        #   where L_i has w as endpoint: p_w * E_w(S) * R_j(R)
        #   so: (1-q_w)*B_w - p_w*A_w (with appropriate signs from S including w)

        per_w = {}
        for w in W:
            W_rest = [u for u in W if u != w]
            A_w = 0.0
            B_w = 0.0
            for rest_mask in range(1 << len(W_rest)):
                S0 = [W_rest[bit] for bit in range(len(W_rest)) if rest_mask & (1 << bit)]
                R0 = [W_rest[bit] for bit in range(len(W_rest)) if not (rest_mask & (1 << bit))]
                sign = (-1) ** len(S0)

                S = S0 + [w]
                # E_w(S) = h_end(S, w)
                Ew = h_end_weighted(T, S, w)
                Rj = h_start_weighted(T, [J] + R0, J)
                Ri = h_start_weighted(T, [I] + R0, I)

                A_w += sign * Ew * Rj
                B_w += sign * Ew * Ri

            pw = T(w, I)
            qw = T(J, w)
            # Contribution when w is last-before-endpoint on the LEFT side
            # For S containing w: sign = (-1)^|S| = (-1)^{|S_0|+1} = -(-1)^|S_0|
            # So actual contribution from w as left-endpoint is:
            # -[(1-qw)*Bw - pw*Aw] = pw*Aw - (1-qw)*Bw
            contrib = pw * A_w - (1 - qw) * B_w
            per_w[w] = (contrib, A_w, B_w)

        # Similarly need contributions when w is the first vertex on the RIGHT side
        # R_j(R) = sum_{v∈R} q_v * ST_v(R)  [for j]
        # R_i(R) = sum_{v∈R} (1-p_v) * ST_v(R)  [for i]
        per_w_right = {}
        for w in W:
            W_rest = [u for u in W if u != w]
            C_w = 0.0
            D_w = 0.0
            for rest_mask in range(1 << len(W_rest)):
                S0 = [W_rest[bit] for bit in range(len(W_rest)) if rest_mask & (1 << bit)]
                R0 = [W_rest[bit] for bit in range(len(W_rest)) if not (rest_mask & (1 << bit))]
                sign = (-1) ** len(S0)

                R = R0 + [w]
                STw = h_start_weighted(T, R, w)
                Lj = h_end_weighted(T, S0 + [J], J)
                Li = h_end_weighted(T, S0 + [I], I)

                C_w += sign * Lj * STw
                D_w += sign * Li * STw

            qw = T(J, w)
            pw = T(w, I)
            # For R containing w: R is W\S where w∉S, so sign = (-1)^|S| = (-1)^|S_0|
            # Contribution from w as first-after-endpoint on RIGHT:
            # (1-pw)*Cw - qw*Dw
            contrib_r = (1 - pw) * C_w - qw * D_w
            per_w_right[w] = (contrib_r, C_w, D_w)

        # Total D(-1) should be sum of all left and right endpoint contributions
        # But we're double-counting (each term has both a left endpoint and right endpoint)
        # The correct decomposition is by left endpoint ONLY or right endpoint ONLY
        total_left = sum(c for c, _, _ in per_w.values())
        total_right = sum(c for c, _, _ in per_w_right.values())

        if trial < 5:
            print(f"  trial {trial}: left_sum={total_left:.6f}, right_sum={total_right:.6f}")
            for w in W:
                pw = T(w, I)
                qw = T(J, w)
                sw = 1 - pw - qw
                cl, al, bl = per_w[w]
                cr, cr2, dr = per_w_right[w]
                print(f"    w={w}: s={sw:.3f}, left_contrib={cl:.4f} "
                      f"(A={al:.4f}, B={bl:.4f}), right_contrib={cr:.4f}")


if __name__ == "__main__":
    test_vertex_decomposition(4)
    test_vertex_decomposition(5, num_trials=50)
    test_factored_structure(4)
    test_factored_structure(5, num_trials=50)
