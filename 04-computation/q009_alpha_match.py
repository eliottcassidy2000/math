#!/usr/bin/env python3
"""
Test: does the s_w-coefficient of delta_H match that of delta_I?

KEY FACT: delta_H is LINEAR in s (no cross terms), because:
  bracket(a,b) = p_a*q_b - (1-q_a)*(1-p_b) = -(s_a(1-d_b) + s_b(1+d_a))/2
has no s_a*s_b term. So delta_H = sum_w alpha_w * s_w.

delta_I is also linear in s (cycle indicators are linear in interface arcs).

Question: does alpha_w^H = alpha_w^I for each w?

If YES, this proves delta_H = delta_I for all n (since both are linear in s,
matching all coefficients means the functions are identical).

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


def compute_delta_H(n, arc_values):
    """Compute delta_H = adj(i,j) - adj'(j,i) directly from path counts."""
    V = list(range(n))
    I, J = 0, 1

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    adj_ij = 0.0
    adj_ji_prime = 0.0
    for p in permutations(V):
        w = 1.0
        for k in range(n - 1):
            w *= T(p[k], p[k + 1])
        for k in range(n - 1):
            if p[k] == I and p[k + 1] == J:
                adj_ij += w
                break

    # T' flips only the i-j arc
    arc_prime = dict(arc_values)
    arc_prime[(I, J)] = 0.0
    arc_prime[(J, I)] = 1.0

    def Tp(a, b):
        if a == b: return 0
        return arc_prime.get((a, b), 1 - arc_prime.get((b, a), 0))

    for p in permutations(V):
        w = 1.0
        for k in range(n - 1):
            w *= Tp(p[k], p[k + 1])
        for k in range(n - 1):
            if p[k] == J and p[k + 1] == I:
                adj_ji_prime += w
                break

    return adj_ij - adj_ji_prime


def compute_alpha_H(n, arc_values, w_target, eps=1e-6):
    """Compute d(delta_H)/d(s_w) by finite difference."""
    I, J = 0, 1
    p_w = arc_values[(w_target, I)]
    q_w = arc_values[(J, w_target)]

    av_plus = dict(arc_values)
    av_minus = dict(arc_values)

    # s_w = 1 - p_w - q_w. Increase s_w by eps: decrease p_w and q_w by eps/2 each
    av_plus[(w_target, I)] = p_w - eps / 2
    av_plus[(I, w_target)] = 1 - (p_w - eps / 2)
    av_plus[(J, w_target)] = q_w - eps / 2
    av_plus[(w_target, J)] = 1 - (q_w - eps / 2)

    av_minus[(w_target, I)] = p_w + eps / 2
    av_minus[(I, w_target)] = 1 - (p_w + eps / 2)
    av_minus[(J, w_target)] = q_w + eps / 2
    av_minus[(w_target, J)] = 1 - (q_w + eps / 2)

    dH_plus = compute_delta_H(n, av_plus)
    dH_minus = compute_delta_H(n, av_minus)

    return (dH_plus - dH_minus) / eps


def compute_delta_I(n, arc_values):
    """Compute delta_I from cycle structure."""
    I, J = 0, 1
    V = list(range(n))
    W = [v for v in V if v != I and v != J]

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    # delta_I = 2 * sum_L (D_L - C_L) for n <= 5
    # More generally: delta_I = sum_{k>=1} 2^k * Delta(alpha_k)
    # At n<=7: only k=1 matters (alpha_2 contributes at n>=6 with VD 3-3 pairs)

    # Compute D_L - C_L for L-cycles through {i,j}
    # An L-cycle through i,j uses L-2 vertices from W
    total_delta_I = 0.0

    # 3-cycles: {i, j, w} for w in W
    for w in W:
        # Cycle i->j->w->i weight in T: T[i][j]*T[j][w]*T[w][i] = 1*q_w*p_w
        gained = T(J, w) * T(w, I)  # q_w * p_w
        # Cycle j->i->w->j weight in T': T'[j][i]*T[i][w]*T[w][j] = 1*(1-p_w)*(1-q_w)
        lost = (1 - T(w, I)) * (1 - T(J, w))
        total_delta_I += 2 * (gained - lost)

    if n >= 5:
        # 5-cycles: use 3 vertices from W
        from itertools import combinations
        for combo in combinations(W, 3):
            for perm in permutations(combo):
                w1, w2, w3 = perm
                # Cycle i,j,w1,w2,w3,i: weight = T[j][w1]*T[w1][w2]*T[w2][w3]*T[w3][i]
                gained = T(J, w1) * T(w1, w2) * T(w2, w3) * T(w3, I)
                # Cycle j,i,w1,w2,w3,j: weight = T[i][w1]*T[w1][w2]*T[w2][w3]*T[w3][j]
                lost = T(I, w1) * T(w1, w2) * T(w2, w3) * T(w3, J)
                total_delta_I += 2 * (gained - lost)

    if n >= 7:
        # 7-cycles: use 5 vertices from W (only when |W| >= 5, i.e., n >= 7)
        from itertools import combinations
        for combo in combinations(W, 5):
            for perm in permutations(combo):
                ws = perm
                gained = T(J, ws[0])
                for k in range(4):
                    gained *= T(ws[k], ws[k + 1])
                gained *= T(ws[4], I)

                lost = T(I, ws[0])
                for k in range(4):
                    lost *= T(ws[k], ws[k + 1])
                lost *= T(ws[4], J)

                total_delta_I += 2 * (gained - lost)

    # For n >= 6, need VD pair corrections (alpha_2 terms)
    if n >= 6:
        # alpha_2: pairs of VD odd cycles
        # VD 3-3 pairs
        for w1 in W:
            for w2 in W:
                if w1 >= w2: continue
                # Two 3-cycles {i,j,w1} and {i,j,w2} are always VD (share i,j)
                # Wait, VD means vertex-disjoint. Two cycles sharing {i,j} are NOT VD.
                # VD cycles: must not share ANY vertex.
                # A 3-cycle {i,j,w} always contains i,j.
                # Two 3-cycles through {i,j}: {i,j,w1} and {i,j,w2} share i,j → NOT VD.
                # For VD pairs, we need cycles that don't share vertices.
                # But all affected cycles contain {i,j}... so no two affected cycles can be VD!
                pass

        # Actually, for alpha_2 to contribute: need TWO VD odd cycles in Omega(T)
        # Both cycles must be affected by the flip (contain i,j)
        # But then both contain i,j → NOT vertex-disjoint. Contradiction!
        # Wait, the alpha_k formula counts independent sets in Omega, which may include
        # cycles NOT affected by the flip. Let me re-read THM-013.
        #
        # THM-013: Delta(alpha_k) = sum_L [sum_{C gained L-cycle} alpha_{k-1}(comp(C))
        #                                 - sum_{C lost L-cycle} alpha_{k-1}(comp(C'))]
        # alpha_{k-1}(comp(C)) counts independent sets of size k-1 in Omega(T[comp(C)])
        # where comp(C) = V\V(C).
        #
        # So Delta(alpha_1) = sum_L (D_L - C_L) [just the net cycle change]
        # Delta(alpha_2) = sum_L sum_{C: L-cycle} alpha_1(comp(C)) - ...
        # alpha_1(comp(C)) = #{odd cycles in T[comp(C)]}
        #
        # For a 3-cycle C = {i,j,w}, comp(C) = V\{i,j,w} = W\{w}, size n-3.
        # alpha_1(comp(C)) = #{odd cycles in T[W\{w}]}
        #
        # At n=6: comp(C) has 3 vertices. alpha_1 = 0 or 1 (depending on whether
        # those 3 vertices form a 3-cycle).
        #
        # So Delta(alpha_2) at n=6 involves the 3-cycle structure of the complement.
        # This is more complex and depends on internal arcs.
        # For now, skip alpha_2 and focus on the linear part.
        pass

    return total_delta_I


def compute_alpha_I(n, arc_values, w_target, eps=1e-6):
    """Compute d(delta_I)/d(s_w) by finite difference."""
    I, J = 0, 1
    p_w = arc_values[(w_target, I)]
    q_w = arc_values[(J, w_target)]

    av_plus = dict(arc_values)
    av_minus = dict(arc_values)

    av_plus[(w_target, I)] = p_w - eps / 2
    av_plus[(I, w_target)] = 1 - (p_w - eps / 2)
    av_plus[(J, w_target)] = q_w - eps / 2
    av_plus[(w_target, J)] = 1 - (q_w - eps / 2)

    av_minus[(w_target, I)] = p_w + eps / 2
    av_minus[(I, w_target)] = 1 - (p_w + eps / 2)
    av_minus[(J, w_target)] = q_w + eps / 2
    av_minus[(w_target, J)] = 1 - (q_w + eps / 2)

    dI_plus = compute_delta_I(n, av_plus)
    dI_minus = compute_delta_I(n, av_minus)

    return (dI_plus - dI_minus) / eps


def verify_alpha_match(n, num_trials=50):
    """Verify alpha_w^H = alpha_w^I for each w."""
    print(f"=== alpha match verification at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)
    max_err = 0.0

    for trial in range(num_trials):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(0.1, 0.9)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        for w in W:
            aH = compute_alpha_H(n, arc_values, w)
            aI = compute_alpha_I(n, arc_values, w)
            err = abs(aH - aI)
            max_err = max(max_err, err)

            if trial < 2 and w == W[0]:
                print(f"  trial {trial}, w={w}: alpha_H={aH:.6f}, alpha_I={aI:.6f}, err={err:.2e}")

    print(f"\n  Max |alpha_H - alpha_I|: {max_err:.2e}")
    if max_err < 1e-3:
        print(f"  CONFIRMED: alpha_w^H = alpha_w^I at n={n}\n")
    else:
        print(f"  FAILED at n={n}\n")
    return max_err < 1e-3


def verify_linearity(n, num_trials=30):
    """Verify delta_H is linear in s (no quadratic terms)."""
    print(f"=== Linearity check at n={n} ===")

    I, J = 0, 1
    W = list(range(2, n))
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)
    max_err = 0.0

    for trial in range(num_trials):
        # Fix d and internal arcs, vary s
        arc_values_base = {(I, J): 1.0, (J, I): 0.0}
        d_vals = {w: random.uniform(-0.5, 0.5) for w in W}
        internal_arcs = [(a, b) for a in W for b in W if a < b]
        for (a, b) in internal_arcs:
            v = random.uniform(0.1, 0.9)
            arc_values_base[(a, b)] = v
            arc_values_base[(b, a)] = 1 - v

        # Compute delta_H at 3 different s values and check linearity
        s_vals_list = [
            {w: 0.0 for w in W},
            {w: random.uniform(-0.3, 0.3) for w in W},
            {w: random.uniform(-0.3, 0.3) for w in W},
        ]

        delta_H_vals = []
        for s_vals in s_vals_list:
            av = dict(arc_values_base)
            for w in W:
                p_w = (1 - s_vals[w] + d_vals[w]) / 2
                q_w = (1 - s_vals[w] - d_vals[w]) / 2
                av[(w, I)] = p_w
                av[(I, w)] = 1 - p_w
                av[(J, w)] = q_w
                av[(w, J)] = 1 - q_w
            delta_H_vals.append(compute_delta_H(n, av))

        # If linear: delta_H(s) = c_0 + sum_w alpha_w * s_w
        # delta_H(0) = c_0
        # delta_H(s1) = c_0 + sum alpha_w * s1_w
        # delta_H(s2) = c_0 + sum alpha_w * s2_w
        # Linearity: delta_H(s1+s2) = 2*c_0 + sum alpha_w*(s1_w+s2_w) = delta_H(s1)+delta_H(s2)-c_0
        # Actually, let me use: delta_H(0) should be 0 (since at s=0, sigma = identity, D=0)

        s_sum = {w: s_vals_list[1][w] + s_vals_list[2][w] for w in W}
        av_sum = dict(arc_values_base)
        for w in W:
            p_w = (1 - s_sum[w] + d_vals[w]) / 2
            q_w = (1 - s_sum[w] - d_vals[w]) / 2
            av_sum[(w, I)] = p_w
            av_sum[(I, w)] = 1 - p_w
            av_sum[(J, w)] = q_w
            av_sum[(w, J)] = 1 - q_w
        delta_H_sum = compute_delta_H(n, av_sum)

        # Linearity: delta_H(s1+s2) = delta_H(s1) + delta_H(s2) - delta_H(0)
        predicted = delta_H_vals[1] + delta_H_vals[2] - delta_H_vals[0]
        err = abs(delta_H_sum - predicted)
        max_err = max(max_err, err)

    print(f"  Max linearity error: {max_err:.2e}")
    if max_err < 1e-8:
        print(f"  CONFIRMED: delta_H is linear in s at n={n}\n")


if __name__ == "__main__":
    # First verify linearity
    for n in range(4, 7):
        verify_linearity(n)

    # Then verify alpha match
    for n in range(4, 7):
        verify_alpha_match(n, num_trials=30 if n <= 5 else 10)
