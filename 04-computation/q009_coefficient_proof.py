#!/usr/bin/env python3
"""
Prove alpha_w^H = alpha_w^I by structural analysis.

KEY FACTS:
1. delta_H = sum_w alpha_w^H * s_w  (LINEAR, no cross terms)
2. delta_I = sum_w alpha_w^I * s_w  (LINEAR, from cycle formula)
3. alpha_w^H = alpha_w^I (verified n=4,...,7)
4. PROVING THIS FOR ALL n WOULD PROVE OCF.

Structure of alpha_w^H (coefficient of s_w in delta_H, at s=0):

alpha_w^H = -h_start(W,w) - h_end(W,w)                              [boundary]
          - sum_{S:w∈S, |S|=1..m-1} h_end(S,w) * Rj(W\S)            [w in left]
          - sum_{S:w∉S, |S|=1..m-1} Li(S) * h_start(W\S,w)          [w in right]

where Rj(R) = sum_b q_b h_start(R,b) and Li(S) = sum_a p_a h_end(S,a) at s=0.

Structure of alpha_w^I:

alpha_w^I = -2*H(B_w) + (higher-L cycle contributions)

where B_w = W\{w} and the cycle contributions involve the same bracket identity
applied to longer cycles.

Strategy: decompose alpha_w^H by vertex-pair (a,b) contributions and match to
cycle structure of alpha_w^I.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations, combinations
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


def H_total(T, verts):
    if len(verts) <= 1:
        return 1.0
    total = 0.0
    for p in permutations(verts):
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def compute_alpha_H_analytical(n, arc_values, w_target):
    """Compute alpha_w^H analytically from the bracket decomposition."""
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    d = {w: T(w, I) - T(J, w) for w in W}

    # Boundary terms
    boundary = -h_start(T, list(W), w_target) - h_end(T, list(W), w_target)

    # Intermediate: w in S (left factor)
    # Contribution: -sum_b (1-d_b)/2 * h_end(S,w) * h_start(R,b)
    # where b ranges over R = W\S, and S ranges over subsets containing w
    left_contrib = 0.0
    W_minus_w = [x for x in W if x != w_target]
    for smask in range(1 << len(W_minus_w)):
        S_rest = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        if len(S_rest) == len(W_minus_w):
            continue  # S = W, skip (handled by boundary R=∅)
        S = S_rest + [w_target]
        R = [x for x in W if x not in S]
        h_S_w = h_end(T, list(S), w_target)
        for b in R:
            left_contrib += -(1 - d[b]) / 2 * h_S_w * h_start(T, list(R), b)

    # Intermediate: w in R (right factor)
    # Contribution: -sum_a (1+d_a)/2 * h_end(S,a) * h_start(R,w)
    right_contrib = 0.0
    for smask in range(1, 1 << len(W_minus_w)):  # exclude S=∅ (handled by boundary S=∅)
        S = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        R_rest = [x for x in W_minus_w if x not in S]
        R = R_rest + [w_target]
        h_R_w = h_start(T, list(R), w_target)
        for a in S:
            right_contrib += -(1 + d[a]) / 2 * h_end(T, list(S), a) * h_R_w

    alpha = boundary + left_contrib + right_contrib
    return alpha, boundary, left_contrib, right_contrib


def compute_alpha_I_analytical(n, arc_values, w_target):
    """Compute alpha_w^I from the cycle formula.

    alpha_w^I = d(delta_I)/d(s_w)

    delta_I = -2*sum_x s_x*H(B_x) + 2*sum_{L>=5}(D_L-C_L) + ...

    The 3-cycle contribution: -2*H(B_w) [coefficient of s_w in -2*sum_x s_x*H(B_x)]

    The 5-cycle contribution: d/ds_w of [2*(D5-C5)].
    5-cycle (i,j,w1,w2,w3): weight = q_{w1}*T[w1][w2]*T[w2][w3]*p_{w3}
    d/ds_w changes q_w or p_w (decreasing by 1/2 each).
    """
    I, J = 0, 1
    W = list(range(2, n))

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    # 3-cycle part
    B_w = [x for x in W if x != w_target]
    H_Bw = H_total(T, B_w)
    alpha_3 = -2 * H_Bw

    # 5-cycle part: d(D5-C5)/ds_w
    # For each 5-cycle through {i,j} involving w_target:
    # Gained: q_{w1}*t12*t23*p_{w3}, Lost: (1-p_{w1})*t12*t23*(1-q_{w3})
    # If w_target = w1: d(gained)/ds_w = d(q_w)/ds_w * rest = (-1/2)*rest
    #                    d(lost)/ds_w = d(1-p_w)/ds_w * rest = (1/2)*rest
    # If w_target = w3: d(gained)/ds_w = rest * d(p_w)/ds_w = rest*(-1/2)
    #                    d(lost)/ds_w = rest * d(1-q_w)/ds_w = rest*(1/2)
    # If w_target = w2 (middle): no interface arcs involved → derivative = 0.
    # Actually, w2 in a 5-cycle uses only internal arcs. The interface arcs are only
    # at the first (leaving j) and last (entering i) positions.
    alpha_5 = 0.0
    if n >= 5 and len(W) >= 3:
        for combo in combinations(W, 3):
            for perm in permutations(combo):
                w1, w2, w3 = perm
                internal = T(w1, w2) * T(w2, w3)

                # Gained: q_{w1} * internal * p_{w3}
                # Lost: (1-p_{w1}) * internal * (1-q_{w3})
                # Net: q_{w1}*p_{w3} - (1-p_{w1})*(1-q_{w3}) times internal
                # = bracket(w3, w1) * internal

                # d(net)/ds_w for the 5-cycle contribution (factor 2):
                if w_target == w1:
                    # d/ds_w of [q_w1*p_w3 - (1-p_w1)*(1-q_w3)]
                    # = d/ds_w1 of bracket(w3,w1)
                    # bracket = -(s_{w3}(1-d_{w1}) + s_{w1}(1+d_{w3}))/2
                    # d/ds_{w1} = -(1+d_{w3})/2 = -p_{w3} (at s=0)
                    d_bracket = -(1 + (T(w3, I) - T(J, w3))) / 2
                    alpha_5 += 2 * d_bracket * internal
                elif w_target == w3:
                    # d/ds_{w3} of bracket(w3,w1)
                    # bracket = -(s_{w3}(1-d_{w1}) + s_{w1}(1+d_{w3}))/2
                    # d/ds_{w3} = -(1-d_{w1})/2 = -q_{w1} (at s=0)
                    d_bracket = -(1 - (T(w1, I) - T(J, w1))) / 2
                    alpha_5 += 2 * d_bracket * internal

    # 7-cycle part: d(D7-C7)/ds_w (n >= 7, uses 5 vertices from W)
    alpha_7 = 0.0
    if n >= 7 and len(W) >= 5:
        for combo in combinations(W, 5):
            for perm in permutations(combo):
                ws = perm  # w1,...,w5
                internal = 1.0
                for k in range(4):
                    internal *= T(ws[k], ws[k + 1])

                if w_target == ws[0]:
                    # d(bracket(w5,w1))/ds_{w1} = -(1+d_{w5})/2
                    d_bracket = -(1 + (T(ws[4], I) - T(J, ws[4]))) / 2
                    alpha_7 += 2 * d_bracket * internal
                elif w_target == ws[4]:
                    # d(bracket(w5,w1))/ds_{w5} = -(1-d_{w1})/2
                    d_bracket = -(1 - (T(ws[0], I) - T(J, ws[0]))) / 2
                    alpha_7 += 2 * d_bracket * internal

    return alpha_3 + alpha_5 + alpha_7, alpha_3, alpha_5 + alpha_7


def verify_analytical(n, num_trials=20):
    """Verify analytical formulas for alpha_H and alpha_I match."""
    print(f"=== Analytical verification at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)
    max_err_H = 0.0
    max_err_I = 0.0
    max_err_match = 0.0

    for trial in range(num_trials):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(0.1, 0.9)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        for w in W[:min(2, len(W))]:
            aH, bnd, left, right = compute_alpha_H_analytical(n, arc_values, w)
            aI, a3, a5 = compute_alpha_I_analytical(n, arc_values, w)
            err = abs(aH - aI)
            max_err_match = max(max_err_match, err)

            if trial < 2:
                print(f"  trial {trial}, w={w}:")
                print(f"    alpha_H = {aH:.6f} (boundary={bnd:.4f}, left={left:.4f}, right={right:.4f})")
                print(f"    alpha_I = {aI:.6f} (3-cycle={a3:.4f}, 5-cycle={a5:.4f})")
                print(f"    err = {err:.2e}")

    print(f"\n  Max |alpha_H - alpha_I|: {max_err_match:.2e}")
    if max_err_match < 1e-8:
        print(f"  CONFIRMED at n={n}\n")
    else:
        print(f"  FAILED at n={n} (may need higher-L cycle terms)\n")
    return max_err_match < 1e-8


if __name__ == "__main__":
    for n in range(4, 8):
        verify_analytical(n, num_trials=20 if n <= 6 else 5)
