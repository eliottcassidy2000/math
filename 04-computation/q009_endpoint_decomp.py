#!/usr/bin/env python3
"""
Decompose alpha_w^H into boundary + interior parts and compare with alpha_w^I.

KEY OBSERVATION: Define f(v) = T(w,v) + q_v for v in B_w.
Then g(v) = T(v,w) + p_v = 2 - f(v).

alpha_w^H = -[start_sum + interior_sum + end_sum] where:
  start_sum = sum_v f(v) * h_start(v, B_w)
  end_sum = sum_v (2-f(v)) * h_end(v, B_w) = 2*H(B_w) - sum_v f(v)*h_end(v)

So: alpha_w^H = -2*H(B_w) - sum_v f(v)*E(v) - interior_sum

where E(v) = h_start(v) - h_end(v) is the "endpoint excess".

The identity alpha_w^H = alpha_w^I becomes:
  interior_sum + sum_v f(v)*E(v) + chain_terms = 0

This script verifies this decomposition and explores the structure of each piece.

Instance: opus-2026-03-05-S6
"""

from itertools import permutations, combinations
import random


def make_T(n, arc_values):
    def T(a, b):
        if a == b: return 0.0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))
    return T


def random_s0_tournament(n, seed=None):
    if seed is not None:
        random.seed(seed)
    I, J = 0, 1
    W = list(range(2, n))
    arc_values = {(I, J): 1.0, (J, I): 0.0}
    for w in W:
        q_w = random.uniform(0.1, 0.9)
        arc_values[(J, w)] = q_w
        arc_values[(w, J)] = 1 - q_w
        arc_values[(w, I)] = 1 - q_w
        arc_values[(I, w)] = q_w
    for a in W:
        for b in W:
            if a < b:
                v = random.uniform(0.1, 0.9)
                arc_values[(a, b)] = v
                arc_values[(b, a)] = 1 - v
    return arc_values


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


def h_start(T, verts, v):
    """Sum of Ham paths starting at v."""
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


def h_end(T, verts, v):
    """Sum of Ham paths ending at v."""
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


def decompose_alpha_H(T, W, w_target):
    """Decompose alpha_w^H into start_sum, interior_sum, end_sum.

    Returns (start_sum, interior_sum, end_sum, alpha_H).
    """
    I, J = 0, 1
    p = {v: T(v, I) for v in W}
    q = {v: T(J, v) for v in W}
    B_w = [x for x in W if x != w_target]

    start_total = 0.0
    interior_total = 0.0
    end_total = 0.0

    for pi in permutations(B_w):
        B_wt = 1.0
        for k in range(len(pi) - 1):
            B_wt *= T(pi[k], pi[k + 1])

        m1 = len(pi)
        for pos in range(m1 + 1):
            if pos == 0:
                val = (T(w_target, pi[0]) + q[pi[0]]) * B_wt
                start_total += val
            elif pos == m1:
                val = (T(pi[-1], w_target) + p[pi[-1]]) * B_wt
                end_total += val
            else:
                u_k = pi[pos - 1]
                u_k1 = pi[pos]
                numer = T(u_k, w_target) * q[u_k1] + p[u_k] * T(w_target, u_k1)
                val = numer * B_wt / T(u_k, u_k1)
                interior_total += val

    alpha_H = -(start_total + interior_total + end_total)
    return start_total, interior_total, end_total, alpha_H


def compute_chain_terms(T, W, w_target):
    """Compute the chain terms from alpha_w^I (excluding -2*H(B_w))."""
    I, J = 0, 1
    p = {v: T(v, I) for v in W}
    q = {v: T(J, v) for v in W}

    chain_sum = 0.0
    for chain_len in range(3, len(W) + 1, 2):
        for combo in combinations(W, chain_len):
            if w_target not in combo:
                continue
            comp = [x for x in W if x not in combo]
            H_comp = H_total(T, comp)
            for perm in permutations(combo):
                internal = 1.0
                for k in range(chain_len - 1):
                    internal *= T(perm[k], perm[k + 1])
                if perm[0] == w_target:
                    chain_sum += 2 * (-p[perm[-1]]) * internal * H_comp
                elif perm[-1] == w_target:
                    chain_sum += 2 * (-q[perm[0]]) * internal * H_comp

    return chain_sum


def test_decomposition(n, num_trials=5):
    """Test the decomposition: interior + f*E + chains = 0."""
    print(f"\n=== Endpoint decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    max_err = 0.0

    for trial in range(num_trials):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p = {v: T(v, I) for v in W}
        q = {v: T(J, v) for v in W}

        for w in W[:2]:
            B_w = [x for x in W if x != w]

            # Decompose alpha_H
            start_sum, interior_sum, end_sum, alpha_H = decompose_alpha_H(T, W, w)

            # f(v) and endpoint excess
            f = {v: T(w, v) + q[v] for v in B_w}
            E = {v: h_start(T, B_w, v) - h_end(T, B_w, v) for v in B_w}
            fE_sum = sum(f[v] * E[v] for v in B_w)

            # Chain terms
            chains = compute_chain_terms(T, W, w)

            # Verify decomposition
            H_Bw = H_total(T, B_w)

            # alpha_H should = -2*H(B_w) - fE_sum - interior_sum
            alpha_H_check = -2 * H_Bw - fE_sum - interior_sum
            err1 = abs(alpha_H - alpha_H_check)

            # Identity: interior + fE + chains = 0
            identity_val = interior_sum + fE_sum + chains
            err2 = abs(identity_val)

            max_err = max(max_err, err2)

            if trial < 2:
                print(f"  trial {trial}, w={w}:")
                print(f"    H(B_w) = {H_Bw:.4f}")
                print(f"    start = {start_sum:.4f}, interior = {interior_sum:.4f}, "
                      f"end = {end_sum:.4f}")
                print(f"    fE = {fE_sum:.4f}, chains = {chains:.4f}")
                print(f"    interior + fE + chains = {identity_val:.2e}")
                print(f"    decomp check: {err1:.2e}")

                # Show individual f and E values
                for v in B_w[:3]:
                    print(f"      v={v}: f={f[v]:.3f}, E={E[v]:.4f}, "
                          f"h_s={h_start(T, B_w, v):.4f}, h_e={h_end(T, B_w, v):.4f}")

    print(f"\n  Max |interior + fE + chains|: {max_err:.2e}")
    return max_err < 1e-10


def analyze_interior_by_edge(T, W, w_target):
    """Break down interior_sum by which edge (u,v) is being split.

    For each directed edge u->v in B_w, the interior contribution from
    inserting w between u and v is:
      [T(u,w)*q_v + p_u*T(w,v)] * (product of other path edges)

    Sum over all paths containing u->v at an interior position.
    """
    I, J = 0, 1
    p = {v: T(v, I) for v in W}
    q = {v: T(J, v) for v in W}
    B_w = [x for x in W if x != w_target]

    edge_contrib = {}
    for u in B_w:
        for v in B_w:
            if u == v:
                continue
            edge_contrib[(u, v)] = 0.0

    for pi in permutations(B_w):
        B_wt = 1.0
        for k in range(len(pi) - 1):
            B_wt *= T(pi[k], pi[k + 1])

        for pos in range(1, len(pi)):  # interior positions
            u_k = pi[pos - 1]
            u_k1 = pi[pos]
            numer = T(u_k, w_target) * q[u_k1] + p[u_k] * T(w_target, u_k1)
            val = numer * B_wt / T(u_k, u_k1)
            edge_contrib[(u_k, u_k1)] += val

    return edge_contrib


def analyze_edge_structure(n, num_trials=3):
    """Analyze interior contributions by edge and relate to chain terms."""
    print(f"\n=== Edge-level analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    for trial in range(min(num_trials, 2)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p = {v: T(v, I) for v in W}
        q = {v: T(J, v) for v in W}

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            edge_contrib = analyze_interior_by_edge(T, W, w)

            print(f"  trial {trial}, w={w}:")

            # For each edge, compute the "insertion weight" = perm / T(u,v)
            # perm(u,v) = T(u,w)*q_v + p_u*T(w,v)
            for u in B_w:
                for v in B_w:
                    if u == v:
                        continue
                    perm_uv = T(u, w) * q[v] + p[u] * T(w, v)
                    ratio = perm_uv / T(u, v) if T(u, v) > 0.01 else float('inf')

                    # Paths through u->v (weight of all other edges)
                    paths_through = 0.0
                    for pi in permutations(B_w):
                        for k in range(len(pi) - 1):
                            if pi[k] == u and pi[k + 1] == v:
                                other_wt = 1.0
                                for j in range(len(pi) - 1):
                                    if j != k:
                                        other_wt *= T(pi[j], pi[j + 1])
                                paths_through += other_wt

                    if len(B_w) <= 4:
                        print(f"    edge ({u},{v}): perm/T = {ratio:.3f}, "
                              f"paths_through = {paths_through:.4f}, "
                              f"contrib = {edge_contrib[(u,v)]:.4f}")


def test_interior_rewrite(n, num_trials=5):
    """Test if interior_sum can be rewritten in terms of h-functions.

    Interior_sum = sum_{u,v in B_w} perm(u,v) * G(u,v)

    where perm(u,v) = T(u,w)*q_v + p_u*T(w,v)
    and G(u,v) = sum of paths through u->v / T(u,v) = "split paths weight"

    G(u,v) = sum_{S+{u}, T+{v} partition of B_w} h_end(u, S+{u}) * h_start(v, T+{v})

    Can we express G(u,v) in terms of h functions on smaller sets?
    """
    print(f"\n=== Interior rewrite at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    for trial in range(min(num_trials, 2)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p = {v: T(v, I) for v in W}
        q = {v: T(J, v) for v in W}

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            m = len(B_w)

            # Compute G(u,v) = sum of products of path edges / T(u,v)
            # for all paths containing edge u->v
            G = {}
            for u in B_w:
                for v in B_w:
                    if u == v:
                        continue
                    total = 0.0
                    rest = [x for x in B_w if x != u and x != v]
                    # Split rest into prefix (ending at u) and suffix (starting at v)
                    for r in range(len(rest) + 1):
                        for prefix_set in combinations(rest, r):
                            suffix_set = tuple(x for x in rest if x not in prefix_set)
                            prefix_verts = list(prefix_set) + [u]
                            suffix_verts = [v] + list(suffix_set)
                            he = h_end(T, prefix_verts, u)
                            hs = h_start(T, suffix_verts, v)
                            total += he * hs
                    G[(u, v)] = total

            # Verify: interior_sum = sum_{u,v} perm(u,v) * G(u,v)
            interior_from_G = 0.0
            for u in B_w:
                for v in B_w:
                    if u == v:
                        continue
                    perm_uv = T(u, w) * q[v] + p[u] * T(w, v)
                    interior_from_G += perm_uv * G[(u, v)]

            _, interior_sum, _, _ = decompose_alpha_H(T, W, w)
            err = abs(interior_from_G - interior_sum)
            print(f"  trial {trial}, w={w}: interior={interior_sum:.6f}, "
                  f"from_G={interior_from_G:.6f}, err={err:.2e}")

            # Now test: does G(u,v) have a nice form?
            # G(u,v) = sum_S h_end(u, S+{u}) * h_start(v, complement+{v})
            # Note: G(u,v) + G(v,u) involves SAME splits (just reversed roles)
            if m <= 4:
                print(f"    G values:")
                for u in B_w:
                    for v in B_w:
                        if u == v:
                            continue
                        print(f"      G({u},{v}) = {G[(u,v)]:.4f}")


def test_perm_decomposition(n, num_trials=5):
    """Test if perm(u,v) = T(u,w)*q_v + p_u*T(w,v) has a nice decomposition.

    Note: perm(u,v) = f(v) - T(w,u)*q_v - q_u*T(w,v)
    where f(v) = T(w,v) + q_v.

    Also: perm(u,v) = T(u,w)*q_v + p_u*T(w,v)
    = (1-T(w,u))*(1-p_v) + (1-q_u)*T(w,v)    [at s=0]
    = 1 - p_v - T(w,u) + T(w,u)*p_v + T(w,v) - q_u*T(w,v)
    = 1 - p_v + T(w,v) - T(w,u) + T(w,u)*p_v - q_u*T(w,v)

    Another form: perm(u,v) = the permanent of [[T(u,w), T(w,v)], [p_u, q_v]]
    = T(u,w)*q_v + p_u*T(w,v)
    """
    print(f"\n=== perm(u,v) decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    for trial in range(1):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p = {v: T(v, I) for v in W}
        q = {v: T(J, v) for v in W}

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            f = {v: T(w, v) + q[v] for v in B_w}

            print(f"  w={w}:")
            for u in B_w:
                for v in B_w:
                    if u == v:
                        continue
                    perm_uv = T(u, w) * q[v] + p[u] * T(w, v)
                    # Alternative: f(v) - cross(u,v)
                    cross = T(w, u) * q[v] + q[u] * T(w, v)
                    alt = f[v] - cross
                    # perm should equal alt? Let's check:
                    # perm = T(u,w)*q_v + p_u*T(w,v) = (1-T(w,u))*q_v + (1-q_u)*T(w,v)
                    # = q_v - T(w,u)*q_v + T(w,v) - q_u*T(w,v)
                    # = f(v) - (T(w,u)*q_v + q_u*T(w,v))
                    # = f(v) - cross(u,v) ✓
                    print(f"    perm({u},{v}) = {perm_uv:.4f}, "
                          f"f(v)-cross = {alt:.4f}, "
                          f"cross = {cross:.4f}, "
                          f"f(v) = {f[v]:.4f}")


if __name__ == "__main__":
    print("=" * 60)
    print("Testing endpoint decomposition identity:")
    print("  interior_sum + sum_v f(v)*E(v) + chain_terms = 0")
    print("=" * 60)

    for n in [5, 6, 7, 8]:
        ok = test_decomposition(n, num_trials=5 if n <= 7 else 3)
        if not ok:
            print(f"FAILED at n={n}")
            break

    # Deeper analysis at small n
    for n in [5, 6]:
        test_perm_decomposition(n, 1)
        test_interior_rewrite(n, 2)

    for n in [5, 6]:
        analyze_edge_structure(n, 1)
