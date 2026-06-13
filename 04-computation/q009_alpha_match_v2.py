#!/usr/bin/env python3
"""
Verify alpha_w^H = alpha_w^I using the FULL delta_I formula.

delta_H is LINEAR in s (proved). delta_I is also linear in s.
Both = sum_w alpha_w * s_w. We verify alpha_w^H = alpha_w^I.

delta_I (full formula) = sum_{k>=1} 2^k * Delta(alpha_k)
At n<=7: delta_I = -2*sum_w s_w*H(B_w) + 2*sum_{L>=5}(D_L-C_L)
where B_w = V\{i,j,w} (complement of the 3-cycle), H(B_w) = Ham path count.

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


def compute_delta_H(n, arc_values):
    """Compute delta_H via subset decomposition."""
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    delta_H = 0.0
    for smask in range(1 << m):
        S = [W[bit] for bit in range(m) if smask & (1 << bit)]
        R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
        Li = h_end(T, S + [I], I)
        Rj = h_start(T, [J] + R, J)
        Lj = h_end(T, S + [J], J)
        Ri = h_start(T, [I] + R, I)
        delta_H += Li * Rj - Lj * Ri
    return delta_H


def compute_delta_I_full(n, arc_values):
    """Compute delta_I using the FULL formula from THM-013.

    At n<=7: delta_I = -2*sum_w s_w*H(B_w) + 2*sum_{L>=5}(D_L-C_L)
    where B_w = W\{w} (the complement of 3-cycle {i,j,w}).
    """
    I, J = 0, 1
    W = list(range(2, n))

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    # s_w and H(B_w)
    total = 0.0
    for w in W:
        s_w = 1 - T(w, I) - T(J, w)
        B_w = [x for x in W if x != w]
        H_Bw = H_total(T, B_w)
        total += -2 * s_w * H_Bw

    # 5-cycle terms: D5-C5
    if n >= 5:
        from itertools import combinations
        for combo in combinations(W, 3):
            for perm in permutations(combo):
                w1, w2, w3 = perm
                # Gained 5-cycle (i,j,w1,w2,w3): T[j][w1]*T[w1][w2]*T[w2][w3]*T[w3][i]
                gained = T(J, w1) * T(w1, w2) * T(w2, w3) * T(w3, I)
                # Lost 5-cycle: T[i][w1]*T[w1][w2]*T[w2][w3]*T[w3][j]
                lost = T(I, w1) * T(w1, w2) * T(w2, w3) * T(w3, J)
                total += 2 * (gained - lost)

    # 7-cycle terms (n >= 7)
    if n >= 7:
        from itertools import combinations
        for combo in combinations(W, 5):
            for perm in permutations(combo):
                ws = perm
                gained = T(J, ws[0])
                lost = T(I, ws[0])
                for k in range(4):
                    gained *= T(ws[k], ws[k + 1])
                    lost *= T(ws[k], ws[k + 1])
                gained *= T(ws[4], I)
                lost *= T(ws[4], J)
                total += 2 * (gained - lost)

    return total


def verify_delta_match(n, num_trials=50):
    """Verify delta_H = delta_I using the full formula."""
    print(f"=== delta_H = delta_I verification at n={n} ===")

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

        dH = compute_delta_H(n, arc_values)
        dI = compute_delta_I_full(n, arc_values)
        err = abs(dH - dI)
        max_err = max(max_err, err)

        if trial < 3:
            print(f"  trial {trial}: delta_H={dH:.6f}, delta_I={dI:.6f}, err={err:.2e}")

    print(f"  Max error: {max_err:.2e}")
    if max_err < 1e-8:
        print(f"  CONFIRMED at n={n}\n")
    else:
        print(f"  FAILED at n={n}\n")
    return max_err < 1e-8


def compute_alpha_by_diff(n, arc_values, w_target, compute_fn, eps=1e-6):
    """Compute d(fn)/d(s_w) by finite difference."""
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

    return (compute_fn(n, av_plus) - compute_fn(n, av_minus)) / eps


def verify_alpha_match(n, num_trials=30):
    """Verify alpha_w^H = alpha_w^I for each w, using FULL delta_I."""
    print(f"=== alpha match at n={n} ===")

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

        for w in W[:min(3, len(W))]:
            aH = compute_alpha_by_diff(n, arc_values, w, compute_delta_H)
            aI = compute_alpha_by_diff(n, arc_values, w, compute_delta_I_full)
            err = abs(aH - aI)
            max_err = max(max_err, err)

            if trial < 2 and w == W[0]:
                print(f"  trial {trial}, w={w}: alpha_H={aH:.4f}, alpha_I={aI:.4f}, err={err:.2e}")

    print(f"  Max |alpha_H - alpha_I|: {max_err:.2e}")
    if max_err < 1e-3:
        print(f"  CONFIRMED at n={n}\n")
    else:
        print(f"  FAILED at n={n}\n")
    return max_err < 1e-3


def check_alpha_formula(n, num_trials=30):
    """
    The alpha_w coefficient of delta_I at n<=7 is:
      alpha_w^I = -2*H(B_w) + 2*d(D5-C5)/ds_w + ...

    The alpha_w coefficient of delta_H should equal this.

    Can we express alpha_w^H analytically?

    alpha_w^H = boundary + sum over (a,b) pairs involving w

    From the bracket decomposition:
    alpha_w^H = -h_start(W,w) - h_end(W,w)
              + sum_{S:w∈S} h_end(S,w) * (-Rj(W\S))  [w enters i, last step -(1-d_b)/2*h_start]
              + sum_{S:w∉S} (-Li(S)) * h_start(W\S,w) [w leaves j, first step -(1+d_a)/2*h_end]

    At s=0: -Rj(R) = -sum_b q_b h_start(R,b), -Li(S) = -sum_a p_a h_end(S,a).
    """
    print(f"=== Alpha formula check at n={n} ===")

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

        def T(a, b):
            if a == b: return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        for w in W[:2]:
            # Numerical alpha
            aH = compute_alpha_by_diff(n, arc_values, w, compute_delta_H)
            aI = compute_alpha_by_diff(n, arc_values, w, compute_delta_I_full)

            # Analytical alpha_I from the n<=7 formula
            s_w = 1 - T(w, I) - T(J, w)
            B_w = [x for x in W if x != w]
            H_Bw = H_total(T, B_w)
            # alpha_w^I ≈ -2*H(B_w) + 5-cycle derivative contribution
            alpha_I_3cycle = -2 * H_Bw

            print(f"  trial {trial}, w={w}: alpha_H={aH:.4f}, alpha_I={aI:.4f}, "
                  f"-2*H(B_w)={alpha_I_3cycle:.4f}, "
                  f"5-cycle contrib={aI-alpha_I_3cycle:.4f}")

    print()


if __name__ == "__main__":
    # Verify delta match with full formula
    for n in range(4, 8):
        ok = verify_delta_match(n, num_trials=30 if n <= 6 else 10)
        if not ok:
            break

    # Verify alpha match with full formula
    for n in range(4, 8):
        ok = verify_alpha_match(n, num_trials=20 if n <= 6 else 5)
        if not ok:
            break

    # Check alpha formula
    for n in [4, 5, 6]:
        check_alpha_formula(n, 3)
