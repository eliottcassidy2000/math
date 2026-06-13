#!/usr/bin/env python3
"""
Test whether G_in(v) = sum_u G(u,v) has a simple formula.

At m=2: G_in(v) = 1
At m=3: G_in(v) = d_out(v) + 1 (verified algebraically)
At m>=4: ???

Also test: Can the full identity interior + fE + chains = 0
be reorganized into a per-vertex identity?

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


def compute_G_direct(T, B_w, u, v):
    """G(u,v) = paths through u->v / T(u,v)."""
    total = 0.0
    for pi in permutations(B_w):
        for k in range(len(pi) - 1):
            if pi[k] == u and pi[k + 1] == v:
                wt = 1.0
                for j in range(len(pi) - 1):
                    if j != k:
                        wt *= T(pi[j], pi[j + 1])
                total += wt
    return total


def test_G_in(n, num_trials=5):
    """Test whether G_in(v) = d_out(v) + C for some constant C."""
    print(f"\n=== G_in(v) analysis at n={n}, |B_w|={n-3} ===\n")

    W = list(range(2, n))

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            m = len(B_w)

            print(f"  trial {trial}, w={w}, m={m}:")
            for v in B_w:
                G_in = sum(compute_G_direct(T, B_w, u, v)
                           for u in B_w if u != v)
                d_out = sum(T(v, u) for u in B_w if u != v)
                d_in = sum(T(u, v) for u in B_w if u != v)
                # d_out + d_in = m - 1

                # Test G_in = d_out + 1? (works at m=3)
                residual_1 = G_in - d_out - 1
                # Test G_in = d_out + (m-2)?
                residual_2 = G_in - d_out - (m - 2)
                # Test G_in = some function of d_out?
                # At m=2: G_in = 1, d_out = T(v,u). Is G_in = d_out + (1-d_out)? No, = 1.
                # At m=3: G_in = d_out + 1
                # At m=4: G_in = ???

                print(f"    v={v}: G_in={G_in:.4f}, d_out={d_out:.4f}, "
                      f"d_in={d_in:.4f}, G_in-d_out={G_in-d_out:.4f}, "
                      f"residual(d+1)={residual_1:.4f}")

            # Also compute G_out
            print(f"  G_out:")
            for v in B_w:
                G_out = sum(compute_G_direct(T, B_w, v, u)
                            for u in B_w if u != v)
                d_out = sum(T(v, u) for u in B_w if u != v)
                print(f"    v={v}: G_out={G_out:.4f}, d_out={d_out:.4f}, "
                      f"G_out-d_out={G_out-d_out:.4f}")

            # Check: sum_v G_in(v) = sum_v G_out(v) = (m-1)*H
            H = H_total(T, B_w)
            sum_G_in = sum(
                sum(compute_G_direct(T, B_w, u, v) for u in B_w if u != v)
                for v in B_w
            )
            print(f"  sum G_in = {sum_G_in:.4f}, (m-1)*H = {(m-1)*H:.4f}")
            print()


def test_E_and_G(n, num_trials=3):
    """Test relationship between E(v) and G-functions.

    E(v) = h_start(v) - h_end(v)
    h_start(v) = sum_u T(v,u)*h_start(u, B_w\{v})  (first-step decomp)
    h_end(v) = sum_u T(u,v)*h_end(u, B_w\{v})  (last-step decomp)

    Also: sum_v E(v) = 0.

    At m=2: E(v) = 2*T(v,u) - 1 = 2*d_out - 1
    At m=3: E(v) = d_out - 1 (verified)
    General: ???
    """
    print(f"\n=== E(v) analysis at n={n}, |B_w|={n-3} ===\n")

    W = list(range(2, n))

    for trial in range(min(num_trials, 2)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            m = len(B_w)

            print(f"  trial {trial}, w={w}, m={m}:")
            sum_E = 0
            for v in B_w:
                hs = 0
                he = 0
                for pi in permutations(B_w):
                    wt = 1.0
                    for k in range(len(pi) - 1):
                        wt *= T(pi[k], pi[k + 1])
                    if pi[0] == v:
                        hs += wt
                    if pi[-1] == v:
                        he += wt
                E_v = hs - he
                sum_E += E_v
                d_out = sum(T(v, u) for u in B_w if u != v)

                # Test: E(v) = d_out - (m-1)/2?
                res_half = E_v - d_out + (m - 1) / 2

                # Test: E(v) depends only on d_out?
                print(f"    v={v}: E={E_v:.4f}, d_out={d_out:.4f}, "
                      f"E - (d_out - (m-1)/2) = {res_half:.4f}")

            print(f"    sum E = {sum_E:.6f}")


def test_telescoping_decomp(n, num_trials=3):
    """Test a telescoping decomposition of S(pi) along the path.

    S(pi) = f(v_1) + sum_k perm(v_k,v_{k+1})/T(v_k,v_{k+1}) + (2-f(v_m))

    Can we write perm(u,v)/T(u,v) = [something that telescopes in f] + [correction]?

    Try: perm(u,v)/T(u,v) = a(u,v) - f(u)/T(u,v) + f(v)/T(u,v) + ...?
    Or: perm(u,v)/T(u,v) = 1 + [f(v)-f(u)-cross(u,v)+T(u,v)]/T(u,v)?

    Actually: perm(u,v)/T(u,v) = [f(v) - cross(u,v)] / T(u,v)
    and f(v) = T(w,v) + q_v.

    So perm/T = f(v)/T(u,v) - cross(u,v)/T(u,v).

    In S(pi):
    sum_k f(v_{k+1})/T(v_k,v_{k+1}) - sum_k cross(v_k,v_{k+1})/T(v_k,v_{k+1}) + f(v_1) + 2 - f(v_m)

    The f terms are: f(v_1) + sum_k f(v_{k+1})/T(v_k,v_{k+1}) + (2-f(v_m))
    = f(v_1) + [f(v_2)/T(v_1,v_2) + ... + f(v_m)/T(v_{m-1},v_m)] + 2 - f(v_m)

    Hmm, f(v_k) appears in both the boundary and interior but with different weights.
    Not obviously telescoping.

    What if we define g(v_k) = f(v_k) * [cumulative something]?

    Actually, let me try a different approach. Define:
    R(k) = product T(v_1,v_2)*...*T(v_{k-1},v_k) (partial path weight)

    Then B_wt * perm(v_k,v_{k+1})/T(v_k,v_{k+1}) = R(k) * perm(v_k,v_{k+1}) * [product after k+1]

    This doesn't factor nicely either.

    Let me just compute S(pi) for several random paths and look for patterns.
    """
    print(f"\n=== Telescoping at n={n} ===\n")
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
            m = len(B_w)

            perms_list = list(permutations(B_w))
            for pi in perms_list[:min(6, len(perms_list))]:
                B_wt = 1.0
                for k in range(m - 1):
                    B_wt *= T(pi[k], pi[k + 1])

                # S(pi) terms
                S_start = f[pi[0]]
                S_end = 2 - f[pi[-1]]
                S_interior = []
                for k in range(m - 1):
                    u, v = pi[k], pi[k + 1]
                    perm_uv = T(u, w) * q[v] + p[u] * T(w, v)
                    cross_uv = T(w, u) * q[v] + q[u] * T(w, v)
                    S_k = perm_uv / T(u, v)
                    f_ratio = f[v] / T(u, v)
                    cross_ratio = cross_uv / T(u, v)
                    S_interior.append((S_k, f_ratio, cross_ratio, perm_uv, T(u, v)))

                S_total = S_start + sum(s[0] for s in S_interior) + S_end
                defect = S_total - 2

                print(f"  pi={pi}: S={S_total:.4f}, defect={defect:.4f}, B_wt={B_wt:.4f}")
                print(f"    start={S_start:.3f}, end={S_end:.3f}")
                for k, (sk, fr, cr, pm, t) in enumerate(S_interior):
                    print(f"    edge {k}: S_k={sk:.3f}, f/T={fr:.3f}, cross/T={cr:.3f}, "
                          f"perm={pm:.3f}, T={t:.3f}")


if __name__ == "__main__":
    for n in [5, 6, 7, 8]:
        test_G_in(n, 3 if n <= 7 else 2)

    for n in [5, 6, 7, 8]:
        test_E_and_G(n, 2)
