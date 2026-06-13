#!/usr/bin/env python3
"""
Investigate the identity G(u,v) + G(v,u) = 2*H(comp) and its corrections.

G(u,v) = sum_{L+R = comp} h_end(u, L+{u}) * h_start(v, R+{v})
       = sum of paths through edge u->v, with weight of u->v removed.

At n=6: G(u,v) + G(v,u) = 2*H(comp) exactly (comp = 1 vertex).

At n>=7: G(u,v) + G(v,u) != 2*H(comp) because some paths don't have u,v adjacent.
The CORRECTION should relate to longer chains in the OCF formula!

HYPOTHESIS: The correction comes from paths where u and v are separated,
and these correspond to chain terms of length >= 5.

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
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[0] != v: continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def h_end(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[-1] != v: continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def compute_G(T, B_w, u, v):
    """G(u,v) = sum over partitions of comp into (L,R),
    h_end(u, L+{u}) * h_start(v, R+{v}).

    Equivalently: paths of B_w through u->v, with T(u,v) weight removed.
    """
    comp = [x for x in B_w if x != u and x != v]
    total = 0.0
    m = len(comp)
    for r in range(m + 1):
        for L in combinations(comp, r):
            R = tuple(x for x in comp if x not in L)
            prefix = list(L) + [u]
            suffix = [v] + list(R)
            he = h_end(T, prefix, u)
            hs = h_start(T, suffix, v)
            total += he * hs
    return total


def compute_G_direct(T, B_w, u, v):
    """Compute G(u,v) directly from path enumeration."""
    total = 0.0
    for pi in permutations(B_w):
        for k in range(len(pi) - 1):
            if pi[k] == u and pi[k+1] == v:
                # Weight = B_wt / T(u,v)
                wt = 1.0
                for j in range(len(pi) - 1):
                    if j != k:
                        wt *= T(pi[j], pi[j+1])
                total += wt
    return total


def test_G_sum_identity(n, num_trials=5):
    """Test: G(u,v) + G(v,u) = 2*H(comp) ?"""
    print(f"\n=== G(u,v) + G(v,u) vs 2*H(comp) at n={n} ===\n")

    W = list(range(2, n))

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            print(f"  trial {trial}, w={w}, B_w={B_w}:")

            for u in B_w:
                for v in B_w:
                    if u >= v: continue
                    comp = [x for x in B_w if x != u and x != v]
                    H_comp = H_total(T, comp)

                    G_uv = compute_G(T, B_w, u, v)
                    G_vu = compute_G(T, B_w, v, u)
                    G_sum = G_uv + G_vu

                    # Also compute directly
                    G_uv_d = compute_G_direct(T, B_w, u, v)
                    G_vu_d = compute_G_direct(T, B_w, v, u)

                    err_direct = abs(G_uv - G_uv_d) + abs(G_vu - G_vu_d)
                    diff = G_sum - 2 * H_comp

                    print(f"    ({u},{v}): G+G'={G_sum:.6f}, 2H={2*H_comp:.6f}, "
                          f"diff={diff:.6f}, direct_err={err_direct:.2e}")


def compute_non_adjacent_paths(T, B_w, u, v):
    """Sum of paths where u and v are NOT adjacent, with weight.

    For each path pi where u and v are not consecutive,
    compute: B_wt(pi) / ???

    Actually, the non-adjacent paths contribute the correction to G(u,v)+G(v,u).
    """
    total = 0.0
    count = 0
    for pi in permutations(B_w):
        # Check if u and v are adjacent
        adjacent = False
        for k in range(len(pi) - 1):
            if (pi[k] == u and pi[k+1] == v) or (pi[k] == v and pi[k+1] == u):
                adjacent = True
                break
        if not adjacent:
            wt = 1.0
            for k in range(len(pi) - 1):
                wt *= T(pi[k], pi[k+1])
            total += wt
            count += 1
    return total, count


def analyze_correction(n, num_trials=3):
    """Analyze the correction G(u,v)+G(v,u) - 2*H(comp) at n>=7."""
    print(f"\n=== Correction analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    for trial in range(min(num_trials, 2)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p = {v: T(v, I) for v in W}
        q = {v: T(J, v) for v in W}

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            print(f"  trial {trial}, w={w}:")

            for u_idx, u in enumerate(B_w):
                for v in B_w[u_idx+1:]:
                    comp = [x for x in B_w if x != u and x != v]
                    H_comp = H_total(T, comp)
                    G_uv = compute_G_direct(T, B_w, u, v)
                    G_vu = compute_G_direct(T, B_w, v, u)
                    correction = G_uv + G_vu - 2 * H_comp

                    # Non-adjacent paths
                    non_adj, n_count = compute_non_adjacent_paths(T, B_w, u, v)

                    # What are the non-adjacent paths?
                    # Path has u at position i, v at position j, |i-j| > 1
                    # The vertices between u and v form a "separator"
                    # This separator is a chain from u to v (or v to u) through interior vertices

                    # For length-5 chains: chain (w, a, b, c, w) uses w and 3 others
                    # The "separator" between u and v in the path corresponds to
                    # a chain in the OCF formula

                    if abs(correction) > 1e-10:
                        print(f"    ({u},{v}): G+G'={G_uv+G_vu:.4f}, 2H={2*H_comp:.4f}, "
                              f"corr={correction:.4f}, non_adj={non_adj:.4f}")

            # Now test: can the full identity be decomposed by pairs?
            # interior = sum_{u,v} perm(u,v) * G(u,v)
            # = sum_{u<v} [perm(u,v)*G(u,v) + perm(v,u)*G(v,u)]
            #
            # Using perm(u,v) = f(v) - cross(u,v):
            # perm(u,v)*G(u,v) + perm(v,u)*G(v,u)
            # = [f(v)-cross(u,v)]*G(u,v) + [f(u)-cross(v,u)]*G(v,u)
            # = f(v)*G(u,v) + f(u)*G(v,u) - cross(u,v)*G(u,v) - cross(v,u)*G(v,u)

            # Key quantities per pair
            print(f"\n    Pair-level analysis:")
            for u_idx, u in enumerate(B_w[:3]):
                for v in B_w[u_idx+1:]:
                    comp = [x for x in B_w if x != u and x != v]
                    H_comp = H_total(T, comp)
                    G_uv = compute_G_direct(T, B_w, u, v)
                    G_vu = compute_G_direct(T, B_w, v, u)

                    f_u = T(w, u) + q[u]
                    f_v = T(w, v) + q[v]
                    cross_uv = T(w, u) * q[v] + q[u] * T(w, v)
                    cross_vu = T(w, v) * q[u] + q[v] * T(w, u)  # same as cross_uv!

                    perm_uv = f_v - cross_uv
                    perm_vu = f_u - cross_vu

                    pair_interior = perm_uv * G_uv + perm_vu * G_vu

                    # Chain terms for this pair (length-3 chains using {u,v} with w at endpoint)
                    # (w, u, v): 2*(-p_v)*T(w,u)*T(u,v)*H(comp)
                    # (w, v, u): 2*(-p_u)*T(w,v)*T(v,u)*H(comp)
                    # (u, v, w): 2*(-q_u)*T(u,v)*T(v,w)*H(comp)
                    # (v, u, w): 2*(-q_v)*T(v,u)*T(u,w)*H(comp)
                    pair_chains_3 = (
                        2*(-p[v])*T(w,u)*T(u,v)*H_comp +
                        2*(-p[u])*T(w,v)*T(v,u)*H_comp +
                        2*(-q[u])*T(u,v)*T(v,w)*H_comp +
                        2*(-q[v])*T(v,u)*T(u,w)*H_comp
                    )

                    # What about longer chains through this pair?
                    pair_chains_5 = 0.0
                    for chain_len in range(5, len(W) + 1, 2):
                        for combo in combinations(W, chain_len):
                            if w not in combo or u not in combo or v not in combo:
                                continue
                            comp5 = [x for x in W if x not in combo]
                            H_comp5 = H_total(T, comp5)
                            for perm in permutations(combo):
                                internal = 1.0
                                for k in range(chain_len - 1):
                                    internal *= T(perm[k], perm[k+1])
                                if perm[0] == w:
                                    pair_chains_5 += 2*(-p[perm[-1]])*internal*H_comp5
                                elif perm[-1] == w:
                                    pair_chains_5 += 2*(-q[perm[0]])*internal*H_comp5

                    print(f"    ({u},{v}): pair_int={pair_interior:.4f}, "
                          f"ch3={pair_chains_3:.4f}, ch5={pair_chains_5:.4f}, "
                          f"sum={pair_interior+pair_chains_3+pair_chains_5:.4f}")


def test_symmetric_decomposition(n, num_trials=3):
    """Test: perm(u,v)*G(u,v) + perm(v,u)*G(v,u) = ?

    Note: cross(u,v) = cross(v,u) = T(w,u)*q_v + q_u*T(w,v).

    So: perm(u,v)*G(u,v) + perm(v,u)*G(v,u)
      = (f_v - cross)*G(u,v) + (f_u - cross)*G(v,u)
      = f_v*G(u,v) + f_u*G(v,u) - cross*[G(u,v)+G(v,u)]

    If G(u,v)+G(v,u) = 2*H(comp), this simplifies to:
      = f_v*G(u,v) + f_u*G(v,u) - 2*cross*H(comp)

    And: f_v*G(u,v) + f_u*G(v,u) = f_v*G(u,v) + f_u*(2*H(comp)-G(u,v))
      = (f_v-f_u)*G(u,v) + 2*f_u*H(comp)

    So pair_interior = (f_v-f_u)*G(u,v) + 2*(f_u - cross)*H(comp)
                     = (f_v-f_u)*G(u,v) + 2*perm(v,u)*H(comp)

    Hmm, this is interesting but only valid when G+G' = 2*H.
    """
    print(f"\n=== Symmetric decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    for trial in range(min(num_trials, 2)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p = {v: T(v, I) for v in W}
        q = {v: T(J, v) for v in W}

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            print(f"  trial {trial}, w={w}:")

            for u_idx, u in enumerate(B_w[:3]):
                for v in B_w[u_idx+1:]:
                    comp = [x for x in B_w if x != u and x != v]
                    H_comp = H_total(T, comp)
                    G_uv = compute_G_direct(T, B_w, u, v)
                    G_vu = compute_G_direct(T, B_w, v, u)

                    f_u = T(w, u) + q[u]
                    f_v = T(w, v) + q[v]
                    cross = T(w, u) * q[v] + q[u] * T(w, v)

                    pair_int = (f_v - cross)*G_uv + (f_u - cross)*G_vu

                    # Alternative using G+G' = 2*H (check if valid)
                    G_sum = G_uv + G_vu
                    alt = (f_v - f_u)*G_uv + 2*(f_u - cross)*H_comp
                    err_alt = abs(pair_int - alt) if abs(G_sum - 2*H_comp) < 1e-10 else float('inf')

                    # The chain terms: -2*[p_v*T(w,u)*T(u,v) + p_u*T(w,v)*T(v,u)
                    #                      + q_u*T(u,v)*T(v,w) + q_v*T(v,u)*T(u,w)]*H(comp)
                    chain_3 = -2*H_comp*(
                        p[v]*T(w,u)*T(u,v) + p[u]*T(w,v)*T(v,u) +
                        q[u]*T(u,v)*T(v,w) + q[v]*T(v,u)*T(u,w)
                    )

                    # Check: chain_3 = -2*H_comp*[T(u,v)*(...) + T(v,u)*(...)]
                    coeff_uv = p[v]*T(w,u) + q[u]*(1-T(w,v))
                    coeff_vu = p[u]*T(w,v) + q[v]*(1-T(w,u))
                    chain_3_alt = -2*H_comp*(T(u,v)*coeff_uv + T(v,u)*coeff_vu)

                    # Note: coeff_uv = p_v*(1-T(u,w)) + q_u*T(v,w)... hmm
                    # Actually: p_v*T(w,u) + q_u*(1-T(w,v))
                    # = (1-q_v)*T(w,u) + q_u - q_u*T(w,v)
                    # = T(w,u) - q_v*T(w,u) + q_u - q_u*T(w,v)
                    # = cross? No, cross = T(w,u)*q_v + q_u*T(w,v)
                    # So coeff_uv = T(w,u) + q_u - cross = ... hmm.

                    # Let's verify: perm(v,u) = T(v,w)*q_u + p_v*T(w,u)
                    # = (1-T(w,v))*q_u + (1-q_v)*T(w,u)
                    # = q_u - q_u*T(w,v) + T(w,u) - q_v*T(w,u)
                    # = T(w,u) + q_u - cross
                    # So coeff_uv = perm(v,u)! Let me check:
                    perm_vu = T(v,w)*q[u] + p[v]*T(w,u)
                    err_coeff = abs(coeff_uv - perm_vu)

                    print(f"    ({u},{v}): cross={cross:.4f}, perm_vu={perm_vu:.4f}, "
                          f"coeff_uv={coeff_uv:.4f}, err={err_coeff:.2e}")

                    # So chain_3 = -2*H_comp*(T(u,v)*perm(v,u) + T(v,u)*perm(u,v))!
                    chain_3_check = -2*H_comp*(T(u,v)*perm_vu + T(v,u)*(f_v-cross))
                    err_chain = abs(chain_3 - chain_3_check)
                    print(f"           chain_3={chain_3:.4f}, check={chain_3_check:.4f}, "
                          f"err={err_chain:.2e}")


if __name__ == "__main__":
    # First check basic G identity
    for n in [5, 6, 7, 8]:
        test_G_sum_identity(n, 3 if n <= 7 else 2)

    # Analyze corrections at n>=7
    for n in [7, 8]:
        analyze_correction(n, 2)

    # Test symmetric decomposition
    for n in [6, 7]:
        test_symmetric_decomposition(n, 2)
