#!/usr/bin/env python3
"""
Direct proof of B(Li,Rj) = B(Lj,Ri) via the alternating sum structure.

KEY INSIGHT: B(Li,Rj) - B(Lj,Ri) is LINEAR in the s-variables.
Each coefficient alpha_v = 0 for all v in W.

From the bracket decomposition: for u in S, w in R (u != w always):
  p_u*q_w - (1-q_u)*(1-p_w) = -[s_u*(1-d_w) + s_w*(1+d_u)] / 2

The s_v-linear coefficient alpha_v gets contributions:
1. Boundary S=empty: -h_start(W,v)
2. Boundary S=W: (-1)^m * (-h_end(W,v))
3. Intermediate (v in R): -(1+d_u)/2 weighted by h_end/h_start
4. Intermediate (v in S): -(1-d_w)/2 weighted by h_end/h_start

The CRUCIAL CANCELLATION: for each vertex u != v, the d_u-dependent parts
from contributions 3 and 4 cancel because (1+d_u)/2 + (1-d_u)/2 = 1.
This requires a PAIRING between contributions where u is in S (before v)
and where u is in R (after v).

Test: verify this cancellation numerically and look for the proof structure.

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


def compute_alpha_v_detailed(n, v_target, arc_values):
    """
    Compute the s_v coefficient alpha_v in B(Li,Rj) - B(Lj,Ri) in detail.
    """
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    # Interface parameters
    d = {}
    for w in W:
        pw = arc_values[(w, I)]
        qw = arc_values[(J, w)]
        d[w] = pw - qw

    # Boundary term 1: S=empty
    # Coefficient of s_v in B_j(W) - B_i(W) = -sum_w s_w * h_start(W, w)
    boundary1 = -h_start(T, list(W), v_target) if v_target in W else 0

    # Boundary term 2: S=W
    # Coefficient of s_v in (-1)^m * (Li(W) - Lj(W)) = (-1)^m * (-sum_u s_u h_end(W,u))
    boundary2 = ((-1) ** m) * (-h_end(T, list(W), v_target)) if v_target in W else 0

    # Intermediate terms
    # For each S with 1 <= |S| <= m-1:
    intermediate_d_dependent = 0.0  # terms with d_u factor
    intermediate_d_independent = 0.0  # terms without d factor (from the 1/2 part)

    for smask in range(1, (1 << m) - 1):  # exclude empty and full
        S = [W[bit] for bit in range(m) if smask & (1 << bit)]
        R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
        sign = (-1) ** len(S)

        if v_target in R:
            # v is in R (right side). Contribution from bracket s_w term with w=v.
            # Coeff: sign * sum_{u in S} -(1+d_u)/2 * h_end(S,u) * h_start(R, v)
            h_R_v = h_start(T, list(R), v_target)
            for u in S:
                h_S_u = h_end(T, list(S), u)
                coeff = sign * (-(1 + d[u]) / 2) * h_S_u * h_R_v
                intermediate_d_dependent += sign * (-d[u] / 2) * h_S_u * h_R_v
                intermediate_d_independent += sign * (-1 / 2) * h_S_u * h_R_v

        elif v_target in S:
            # v is in S (left side). Contribution from bracket s_u term with u=v.
            # Coeff: sign * sum_{w in R} -(1-d_w)/2 * h_end(S, v) * h_start(R, w)
            h_S_v = h_end(T, list(S), v_target)
            for w in R:
                h_R_w = h_start(T, list(R), w)
                coeff = sign * (-(1 - d[w]) / 2) * h_S_v * h_R_w
                intermediate_d_dependent += sign * (d[w] / 2) * h_S_v * h_R_w
                intermediate_d_independent += sign * (-1 / 2) * h_S_v * h_R_w

    alpha = boundary1 + boundary2 + intermediate_d_dependent + intermediate_d_independent

    return {
        'alpha': alpha,
        'boundary1': boundary1,
        'boundary2': boundary2,
        'd_dep': intermediate_d_dependent,
        'd_indep': intermediate_d_independent,
        'boundary_total': boundary1 + boundary2,
        'intermediate_total': intermediate_d_dependent + intermediate_d_independent,
    }


def test_alpha_decomposition(n, num_trials=50):
    """Test the decomposition of alpha_v and verify cancellation."""
    print(f"=== alpha_v decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 10)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-1, 2)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        for v in W[:2]:  # Just show first 2 vertices
            res = compute_alpha_v_detailed(n, v, arc_values)
            if trial < 3:
                print(f"  trial {trial}, v={v}:")
                print(f"    boundary = {res['boundary_total']:.6f} "
                      f"(start={res['boundary1']:.4f}, end={res['boundary2']:.4f})")
                print(f"    intermediate = {res['intermediate_total']:.6f} "
                      f"(d_dep={res['d_dep']:.4f}, d_indep={res['d_indep']:.4f})")
                print(f"    alpha = {res['alpha']:.2e}")
                print(f"    boundary + d_indep = {res['boundary_total']+res['d_indep']:.6f}")
                print(f"    (should equal -d_dep = {-res['d_dep']:.6f})")

    # Check: is alpha = 0 for all v, all trials?
    max_alpha = 0.0
    max_boundary_plus_indep = 0.0
    for trial in range(num_trials):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-2, 3)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        for v in W:
            res = compute_alpha_v_detailed(n, v, arc_values)
            max_alpha = max(max_alpha, abs(res['alpha']))
            max_boundary_plus_indep = max(max_boundary_plus_indep,
                                          abs(res['boundary_total'] + res['d_indep']))

    print(f"\n  Max |alpha_v|: {max_alpha:.2e}")
    print(f"  Max |boundary + d_indep|: {max_boundary_plus_indep:.2e}")


def test_d_independent_part(n, num_trials=50):
    """
    The d-independent part of alpha_v:
    alpha_v^{d-indep} = boundary + sum_S (-1)^|S| (-1/2) * (stuff)

    This should equal MINUS the d-dependent part.
    Check: what is the d-independent part by itself?
    """
    print(f"\n=== d-independent structure at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 5)):
        # Set all d = 0 (p_w = q_w = 1/2) to isolate d-independent behavior
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for w in W:
            arc_values[(w, I)] = 0.5  # p_w = 1/2
            arc_values[(I, w)] = 0.5
            arc_values[(J, w)] = 0.5  # q_w = 1/2
            arc_values[(w, J)] = 0.5
        # Random internal arcs
        for a in W:
            for b in W:
                if a < b:
                    v = random.uniform(-1, 2)
                    arc_values[(a, b)] = v
                    arc_values[(b, a)] = 1 - v

        for v in W[:2]:
            res = compute_alpha_v_detailed(n, v, arc_values)
            print(f"  trial {trial}, v={v} (all d=0):")
            print(f"    boundary = {res['boundary_total']:.6f}")
            print(f"    d_dep = {res['d_dep']:.6f} (should be 0 since all d=0)")
            print(f"    d_indep = {res['d_indep']:.6f}")
            print(f"    alpha = {res['alpha']:.2e}")
            print(f"    boundary + d_indep = {res['boundary_total'] + res['d_indep']:.2e}")


def test_pairing_structure(n, v_target, num_trials=10):
    """
    For each vertex u != v, trace the paired contributions from:
    - u in S (before v): gives (1-d_w)/2 factor for some w
    - u in R (after v): gives (1+d_u)/2 factor

    Check: do these pair up to give u-independent contributions?
    """
    print(f"\n=== Pairing structure at n={n}, v={v_target} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 3)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(0, 1)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        def T(a, b):
            if a == b: return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        d = {w: arc_values[(w, I)] - arc_values[(J, w)] for w in W}

        print(f"  trial {trial}:")

        # For each u != v, compute the contribution from subsets where u is the boundary
        for u in W:
            if u == v_target:
                continue

            # Contribution where u is LEFT boundary (u in S, v in R):
            # sum over S with u in S, v not in S, sign * -(1+d_u)/2 * h_end(S,u) * h_start(R,v)
            left_contrib = 0.0
            for smask in range(1, (1 << m) - 1):
                S = [W[bit] for bit in range(m) if smask & (1 << bit)]
                R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
                if u not in S or v_target not in R:
                    continue
                sign = (-1) ** len(S)
                left_contrib += sign * h_end(T, list(S), u) * h_start(T, list(R), v_target)

            # Contribution where u is RIGHT boundary (u in R, v in S):
            # sum over S with v in S, u not in S, sign * -(1-d_u)/2 * h_end(S,v) * h_start(R,u)
            right_contrib = 0.0
            for smask in range(1, (1 << m) - 1):
                S = [W[bit] for bit in range(m) if smask & (1 << bit)]
                R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
                if v_target not in S or u not in R:
                    continue
                sign = (-1) ** len(S)
                right_contrib += sign * h_end(T, list(S), v_target) * h_start(T, list(R), u)

            # The paired contribution to alpha_v from vertex u:
            # left: -(1+d_u)/2 * left_contrib
            # right: -(1-d_u)/2 * right_contrib
            # Total: -(1+d_u)/2 * left_contrib - (1-d_u)/2 * right_contrib
            # = -1/2 * (left_contrib + right_contrib) - d_u/2 * (left_contrib - right_contrib)

            total = -(1 + d[u]) / 2 * left_contrib - (1 - d[u]) / 2 * right_contrib
            d_indep_part = -1 / 2 * (left_contrib + right_contrib)
            d_dep_part = -d[u] / 2 * (left_contrib - right_contrib)

            print(f"    u={u}: left={left_contrib:.4f}, right={right_contrib:.4f}, "
                  f"L+R={left_contrib + right_contrib:.4f}, L-R={left_contrib - right_contrib:.4f}, "
                  f"total={total:.4f}")

        # Also compute boundary contributions
        boundary = -h_start(T, list(W), v_target) + ((-1) ** m) * (-h_end(T, list(W), v_target))
        print(f"    boundary = {boundary:.4f}")
        print()


if __name__ == "__main__":
    test_alpha_decomposition(4)
    test_alpha_decomposition(5)
    test_alpha_decomposition(6)
    test_d_independent_part(5)
    test_pairing_structure(5, 2)
