#!/usr/bin/env python3
"""
Attempt to prove alpha_w^H = alpha_w^I by induction.

alpha_w^H (the s_w coefficient of delta_H at s=0) can be decomposed as:

  alpha_w^H = boundary + left + right

where:
  boundary = -h_start(W,w) - h_end(W,w)
  left = -sum_{S:w∈S} h_end(S,w) * Rj(R)  [w enters i from left side]
  right = -sum_{S:w∉S} Li(S) * h_start(R,w)  [w leaves j on right side]

At s=0: Rj(R) = sum_b q_b*h_start(R,b) and Li(S) = sum_a p_a*h_end(S,a)
where p_a = (1+d_a)/2, q_b = (1-d_b)/2.

alpha_w^I = -2*H(B_w) + (cycle derivatives).

KEY INSIGHT: Both the alternating sum (THM-016/017) and the unsigned sum (OCF)
involve the SAME structural ingredients. The difference is only in the signs.

For the alternating sum: (-1)^|S| → cancellation gives 0.
For the unsigned sum: (+1)^|S| → gives alpha_w^I.

Can we COMBINE the two results?

From THM-016/017: sum_S (-1)^|S| [stuff(S,w)] = 0
From OCF: sum_S (+1)^|S| [stuff(S,w)] = alpha_w^I

These give: 2*sum_{|S| even} [stuff] = alpha_w^I and sum_{|S| even} = sum_{|S| odd}.

Can we compute sum_{|S| even} [stuff] directly?

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


def compute_alpha_H_unsigned(n, arc_values, w_target):
    """
    Compute alpha_w^H: coefficient of s_w in delta_H (UNSIGNED sum) at s=0.

    Decompose into contributions from each pair (a,b) where
    a is the last W-vertex before i, and b is the first W-vertex after j.
    """
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    p = {w: T(w, I) for w in W}
    q = {w: T(J, w) for w in W}

    total = 0.0

    # Boundary S=∅: coefficient of s_w is -h_start(W, w)
    total += -h_start(T, list(W), w_target)

    # Boundary R=∅ (S=W): coefficient of s_w is -h_end(W, w)
    total += -h_end(T, list(W), w_target)

    # Intermediate: w ∈ S (w is the "last before i" vertex or contributes through h_end)
    # For each S containing w (1 ≤ |S| ≤ m-1), for each b ∈ R:
    # contribution = -(1-d_b)/2 * h_end(S,w) * h_start(R,b) = -q_b * h_end(S,w) * h_start(R,b)
    W_minus_w = [x for x in W if x != w_target]
    for smask in range(1 << len(W_minus_w)):
        S_rest = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        if len(S_rest) == len(W_minus_w):
            continue  # S would be W (already handled by boundary R=∅)
        S = S_rest + [w_target]
        R = [x for x in W if x not in S]
        h_S_w = h_end(T, list(S), w_target)
        for b in R:
            total += -q[b] * h_S_w * h_start(T, list(R), b)

    # Intermediate: w ∈ R (w is the "first after j" vertex or contributes through h_start)
    for smask in range(1, 1 << len(W_minus_w)):  # exclude empty (already handled by boundary S=∅)
        S = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        R_rest = [x for x in W_minus_w if x not in S]
        R = R_rest + [w_target]
        h_R_w = h_start(T, list(R), w_target)
        for a in S:
            total += -p[a] * h_end(T, list(S), a) * h_R_w

    return total


def compute_alpha_H_alternating(n, arc_values, w_target):
    """
    Compute alpha_w: coefficient of s_w in B(Li,Rj)-B(Lj,Ri) (ALTERNATING sum) at s=0.
    Should be 0 by THM-016/017.
    """
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    p = {w: T(w, I) for w in W}
    q = {w: T(J, w) for w in W}

    total = 0.0

    # Same structure but with (-1)^|S| signs
    # Boundary S=∅: (-1)^0 = 1. Contribution: -h_start(W,w)
    total += -h_start(T, list(W), w_target)

    # Boundary S=W: (-1)^m. Contribution: (-1)^m * (-h_end(W,w))
    total += ((-1) ** m) * (-h_end(T, list(W), w_target))

    W_minus_w = [x for x in W if x != w_target]
    for smask in range(1 << len(W_minus_w)):
        S_rest = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        if len(S_rest) == len(W_minus_w):
            continue
        S = S_rest + [w_target]
        R = [x for x in W if x not in S]
        sign = (-1) ** len(S)
        h_S_w = h_end(T, list(S), w_target)
        for b in R:
            total += sign * (-q[b]) * h_S_w * h_start(T, list(R), b)

    for smask in range(1, 1 << len(W_minus_w)):
        S = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        R_rest = [x for x in W_minus_w if x not in S]
        R = R_rest + [w_target]
        sign = (-1) ** len(S)
        h_R_w = h_start(T, list(R), w_target)
        for a in S:
            total += sign * (-p[a]) * h_end(T, list(S), a) * h_R_w

    return total


def test_sum_and_difference(n, num_trials=30):
    """
    We have two identities:
      sum_S (-1)^|S| f(S,w) = 0  (alternating, THM-016/017)
      sum_S (+1)   * f(S,w) = alpha_w^I  (unsigned, OCF)

    Adding: 2 * sum_{|S| even} f(S,w) = alpha_w^I
    Subtracting: 2 * sum_{|S| odd} f(S,w) = alpha_w^I

    So: alpha_w^I = 2 * sum_{|S| even} f(S,w) = 2 * sum_{|S| odd} f(S,w)

    What is sum_{|S| even} f(S,w)? Can we express it in terms of cycles?
    """
    print(f"=== Sum decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 3)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(0.1, 0.9)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        for w in W[:2]:
            alpha_unsigned = compute_alpha_H_unsigned(n, arc_values, w)
            alpha_alternating = compute_alpha_H_alternating(n, arc_values, w)

            print(f"  trial {trial}, w={w}:")
            print(f"    alpha_unsigned = {alpha_unsigned:.6f}")
            print(f"    alpha_alternating = {alpha_alternating:.2e}")
            print(f"    half_unsigned = {alpha_unsigned/2:.6f}")

            # Decompose by even/odd |S|
            def T(a, b):
                if a == b: return 0
                return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))
            p = {x: T(x, I) for x in W}
            q = {x: T(J, x) for x in W}

            even_sum = 0.0
            odd_sum = 0.0
            # ... for boundary and intermediate terms
            # This is complex. Let me compute from the total:
            # even + odd = alpha_unsigned, even - odd = alpha_alternating = 0
            # So even = odd = alpha_unsigned / 2
            even_sum = alpha_unsigned / 2
            odd_sum = alpha_unsigned / 2

            # What is -H(B_w)?
            B_w = [x for x in W if x != w]
            H_Bw = H_total(T, B_w)
            print(f"    -H(B_w) = {-H_Bw:.6f}")
            print(f"    alpha/2 = {alpha_unsigned/2:.6f}")
            print(f"    alpha/2 + H(B_w) = {alpha_unsigned/2 + H_Bw:.6f}")
            print()


def test_recursive_structure(n, num_trials=20):
    """
    Is alpha_w^H related to sub-tournament quantities in a way that
    naturally matches alpha_w^I = -2*H(B_w) + cycle terms?

    At n=4: alpha_w^H = -2 = -2*1 = -2*H(B_w) since B_w has 1 vertex.
    At n=5: alpha_w^H = -2*H(B_w) + 5-cycle derivative.

    Hypothesis: alpha_w^H/(-2) = H(B_w) - [5-cycle derivative]/2
    Can we express alpha_w^H as a formula involving H of sub-tournaments?
    """
    print(f"\n=== Recursive structure at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 5)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(0.1, 0.9)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        def T(a, b):
            if a == b: return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        for w in W[:2]:
            alpha = compute_alpha_H_unsigned(n, arc_values, w)
            B_w = [x for x in W if x != w]
            H_Bw = H_total(T, B_w)

            # Compute sub-tournament delta_H for each sub-tournament W\{w}
            # with the same interface (i,j) but on n-1 vertices
            # delta_H_{n-1} at the same d,internal values but without vertex w

            # The "alpha_w^H" should relate to delta_H on the reduced tournament
            # through the bracket decomposition

            # What is alpha_w^H + 2*H(B_w)? This is the "residual" beyond the
            # 3-cycle contribution
            residual = alpha + 2 * H_Bw
            print(f"  trial {trial}, w={w}: alpha={alpha:.4f}, -2*H(B_w)={-2*H_Bw:.4f}, "
                  f"residual={residual:.4f}")

    print()


if __name__ == "__main__":
    for n in [4, 5, 6]:
        test_sum_and_difference(n, 3)

    for n in [4, 5, 6]:
        test_recursive_structure(n, 5)
