#!/usr/bin/env python3
"""
Attempt a general proof of alpha_w^H = alpha_w^I for all n.

PROOF STRATEGY: Strong induction on n.

Base case: n=4 (trivial, S(pi)=2 always).

Inductive step: Assume OCF for all m < n. Then:
  alpha_w^I = -2*H(B_w) + sum_{chains with w at endpoint} 2*(-p/-q)*internal*H(comp)

We need: alpha_w^H = alpha_w^I where alpha_w^H = -sum_pi B_wt*S(pi).

KEY OBSERVATION: Both sides can be decomposed by "removing a vertex u from B_w".

For alpha_w^H: group perms of B_w by position of u, then the remaining vertices
form a perm of B_w \ {u}. The contribution involves H of sub-parts.

For alpha_w^I: H(B_w) decomposes by position of u, and the chain terms either
include u or not.

The INDUCTIVE STRUCTURE suggests:
  alpha_w^H(n) - alpha_w^H(n-1) involves u-dependent terms
  alpha_w^I(n) - alpha_w^I(n-1) involves u-dependent terms
And by induction at n-1, these differences match.

Actually, let me think about this differently. At n, we have W with m vertices.
The identity is a polynomial in the arcs of the complete graph on W plus the
interface arcs p_v, q_v.

When we add vertex u to W (going from n-1 to n):
- New arcs: T(u,v) for all v in W, plus p_u, q_u
- New perms: those including u
- New chains: those including u

For alpha_w^H: the perms of B_w' = B_w ∪ {u} include u at each position.
For alpha_w^I: H(B_w') and the new chains that include u.

Can we express the DIFFERENCE alpha_w^H(n) - alpha_w^H(n-1) in terms of
alpha_w^H at smaller sizes? And similarly for alpha_w^I?

Actually, there's a cleaner approach. Let me think about what happens when
we FIX the position of u in the B_w permutation.

For a perm pi' = (v_1, ..., v_j, u, v_{j+1}, ..., v_{m-2}) of B_w':
the S(pi') value involves w inserted at each of m positions.

The contributions split into:
1. w inserted BEFORE u: involves only (v_1,...,v_j) and u, v_{j+1},...
2. w inserted AFTER u: involves v_1,...,v_j, u and (v_{j+1},...,v_{m-2})
3. w inserted at u's position (before or after): adjacent to u

In cases 1 and 2, the w-insertion only interacts with ONE side of u.
This is where the inductive structure might help: the "left" and "right"
parts of pi' are perms of subsets of B_w, and by induction their S-values
are related to H and cycle derivatives of those subsets.

Let me test this numerically.

Instance: opus-2026-03-05-S5b
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


def compute_alpha_H_by_u_position(T, W, w_target, u):
    """Decompose alpha_w^H by position of u in B_w permutations.

    For each perm pi of B_w = B_w_inner ∪ {u} where B_w_inner = B_w\{u}:
    u is at some position j. The other vertices form a perm of B_w_inner.

    Then w is inserted at each position in pi. Group the contributions by
    whether w is to the LEFT of u, at u's position, or to the RIGHT.
    """
    I, J = 0, 1
    p = {v: T(v, I) for v in W}
    q = {v: T(J, v) for v in W}

    B_w = [x for x in W if x != w_target]
    B_w_inner = [x for x in B_w if x != u]

    # For each perm of B_w with u at position j:
    contributions = {'left_of_u': 0.0, 'right_of_u': 0.0, 'at_u': 0.0}

    for pi_inner in permutations(B_w_inner):
        for j in range(len(pi_inner) + 1):
            # Insert u at position j to get pi
            pi = list(pi_inner[:j]) + [u] + list(pi_inner[j:])
            m1 = len(pi)

            # Path weight
            B_wt = 1.0
            for k in range(m1 - 1):
                B_wt *= T(pi[k], pi[k + 1])

            # Insert w at each position
            for pos in range(m1 + 1):
                if pos == 0:
                    val = -(T(w_target, pi[0]) + q[pi[0]]) * B_wt
                elif pos == m1:
                    val = -(T(pi[-1], w_target) + p[pi[-1]]) * B_wt
                else:
                    u_k = pi[pos - 1]
                    u_k1 = pi[pos]
                    numer = T(u_k, w_target) * q[u_k1] + p[u_k] * T(w_target, u_k1)
                    val = -numer * B_wt / T(u_k, u_k1)

                # Classify: is w to the left of u, right of u, or at u's position?
                if pos <= j:
                    contributions['left_of_u'] += val
                elif pos == j + 1:
                    # w is right after u (but this is "at u's right boundary")
                    contributions['at_u'] += val
                else:
                    contributions['right_of_u'] += val

    return contributions


def test_u_decomposition(n, num_trials=5):
    """Test how alpha_w^H decomposes by u's position."""
    print(f"\n=== u-decomposition at n={n} ===\n")

    W = list(range(2, n))

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            for u in B_w[:1]:
                contribs = compute_alpha_H_by_u_position(T, W, w, u)
                total = sum(contribs.values())

                # Compare with full alpha_H
                alpha_total = 0.0
                for pi in permutations(B_w):
                    B_wt = 1.0
                    for k in range(len(pi) - 1):
                        B_wt *= T(pi[k], pi[k + 1])
                    m1 = len(pi)
                    for pos in range(m1 + 1):
                        if pos == 0:
                            val = -(T(w, pi[0]) + T(J, pi[0])) * B_wt
                        elif pos == m1:
                            val = -(T(pi[-1], w) + T(pi[-1], I)) * B_wt
                        else:
                            u_k = pi[pos - 1]
                            u_k1 = pi[pos]
                            numer = T(u_k, w) * T(J, u_k1) + T(u_k, I) * T(w, u_k1)
                            val = -numer * B_wt / T(u_k, u_k1)
                        alpha_total += val

                print(f"  trial {trial}, w={w}, u={u}:")
                print(f"    left={contribs['left_of_u']:.4f}, at_u={contribs['at_u']:.4f}, "
                      f"right={contribs['right_of_u']:.4f}")
                print(f"    total={total:.4f}, alpha_H={alpha_total:.4f}, "
                      f"err={abs(total-alpha_total):.2e}")
    print()


def analyze_first_step(n, num_trials=5):
    """Analyze alpha_w^H using a FIRST-STEP decomposition of B_w perms.

    For perm pi of B_w, the first vertex pi[0] = u determines:
    - The first arc T(u, pi[1])
    - The rest is a perm of B_w\{u} starting from pi[1]

    This decomposes H(B_w) as: H(B_w) = sum_u sum_v T(u,v) * h_start(B_w\{u}, v)

    Similarly for alpha_w^H, the first-step decomposition groups by the first vertex.

    For the chain terms in alpha_w^I:
    - 3-cycle part: -2*H(B_w)
    - Chains starting with w then u: 2*(-p_{last})*T(w,u)*... * H(comp)
    - Chains ending with u then w: 2*(-q_u)*T(u,...)*... * H(comp)

    If the first-step decomposition matches at each level, we get an inductive proof!
    """
    print(f"\n=== First-step decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    p = lambda v, T=None: T(v, I)
    q = lambda v, T=None: T(J, v)

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p_dict = {v: T(v, I) for v in W}
        q_dict = {v: T(J, v) for v in W}

        for w in W[:1]:
            B_w = [x for x in W if x != w]

            # Decompose alpha_w^H by first vertex of B_w perm
            first_vertex_contrib = {}
            for u in B_w:
                total_u = 0.0
                rest = [x for x in B_w if x != u]
                for pi_rest in permutations(rest):
                    pi = [u] + list(pi_rest)
                    B_wt = 1.0
                    for k in range(len(pi) - 1):
                        B_wt *= T(pi[k], pi[k + 1])
                    m1 = len(pi)
                    for pos in range(m1 + 1):
                        if pos == 0:
                            val = -(T(w, pi[0]) + q_dict[pi[0]]) * B_wt
                        elif pos == m1:
                            val = -(T(pi[-1], w) + p_dict[pi[-1]]) * B_wt
                        else:
                            u_k = pi[pos - 1]
                            u_k1 = pi[pos]
                            numer = T(u_k, w) * q_dict[u_k1] + p_dict[u_k] * T(w, u_k1)
                            val = -numer * B_wt / T(u_k, u_k1)
                        total_u += val
                first_vertex_contrib[u] = total_u

            # Decompose alpha_w^I by first vertex
            # -2*H(B_w) decomposes as -2*sum_u sum_v T(u,v)*H_start(rest, v)
            # Chains: group by first W-vertex after w (or last before w)

            # For now, just compare totals
            alpha_H = sum(first_vertex_contrib.values())

            # Full alpha_I
            alpha_I = -2 * H_total(T, B_w)
            for chain_len in range(3, len(W) + 1, 2):
                for combo in combinations(W, chain_len):
                    if w not in combo:
                        continue
                    comp = [x for x in W if x not in combo]
                    H_comp = H_total(T, comp)
                    for perm in permutations(combo):
                        internal = 1.0
                        for k in range(chain_len - 1):
                            internal *= T(perm[k], perm[k + 1])
                        if perm[0] == w:
                            alpha_I += 2 * (-p_dict[perm[-1]]) * internal * H_comp
                        elif perm[-1] == w:
                            alpha_I += 2 * (-q_dict[perm[0]]) * internal * H_comp

            print(f"  trial {trial}, w={w}:")
            print(f"    alpha_H = {alpha_H:.6f}, alpha_I = {alpha_I:.6f}, "
                  f"err = {abs(alpha_H - alpha_I):.2e}")
            for u in B_w[:3]:
                print(f"    contrib[first={u}] = {first_vertex_contrib[u]:.4f}")
    print()


def test_inductive_relation(n, num_trials=5):
    """Test the inductive relation between n and n-1.

    At n: alpha_w^H(n) = insertion decomposition using B_w with m-1 vertices.
    At n-1: alpha_w^H(n-1) = insertion decomposition using B_w\{u} with m-2 vertices.

    Is alpha_w^H(n) = alpha_w^H(n-1) * (something involving u)
                    + (correction from chains involving u)?
    """
    print(f"\n=== Inductive relation n={n} vs n-1={n-1} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    W_prev = list(range(2, n - 1))

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p_dict = {v: T(v, I) for v in W}
        q_dict = {v: T(J, v) for v in W}

        for w in W_prev[:1]:
            # alpha at n
            B_w_n = [x for x in W if x != w]
            alpha_n = 0.0
            for pi in permutations(B_w_n):
                B_wt = 1.0
                for k in range(len(pi) - 1):
                    B_wt *= T(pi[k], pi[k + 1])
                m1 = len(pi)
                for pos in range(m1 + 1):
                    if pos == 0:
                        val = -(T(w, pi[0]) + q_dict[pi[0]]) * B_wt
                    elif pos == m1:
                        val = -(T(pi[-1], w) + p_dict[pi[-1]]) * B_wt
                    else:
                        u_k = pi[pos - 1]
                        u_k1 = pi[pos]
                        numer = T(u_k, w) * q_dict[u_k1] + p_dict[u_k] * T(w, u_k1)
                        val = -numer * B_wt / T(u_k, u_k1)
                    alpha_n += val

            # alpha at n-1 (without the last vertex)
            u_new = n - 1  # the vertex added going from n-1 to n
            B_w_prev = [x for x in W_prev if x != w]
            alpha_prev = 0.0
            for pi in permutations(B_w_prev):
                B_wt = 1.0
                for k in range(len(pi) - 1):
                    B_wt *= T(pi[k], pi[k + 1])
                m1 = len(pi)
                for pos in range(m1 + 1):
                    if pos == 0:
                        val = -(T(w, pi[0]) + q_dict[pi[0]]) * B_wt
                    elif pos == m1:
                        val = -(T(pi[-1], w) + p_dict[pi[-1]]) * B_wt
                    else:
                        u_k = pi[pos - 1]
                        u_k1 = pi[pos]
                        numer = T(u_k, w) * q_dict[u_k1] + p_dict[u_k] * T(w, u_k1)
                        val = -numer * B_wt / T(u_k, u_k1)
                    alpha_prev += val

            # alpha_I at n
            alpha_I_n = -2 * H_total(T, B_w_n)
            for chain_len in range(3, len(W) + 1, 2):
                for combo in combinations(W, chain_len):
                    if w not in combo:
                        continue
                    comp = [x for x in W if x not in combo]
                    H_comp = H_total(T, comp)
                    for perm in permutations(combo):
                        internal = 1.0
                        for k in range(chain_len - 1):
                            internal *= T(perm[k], perm[k + 1])
                        if perm[0] == w:
                            alpha_I_n += 2 * (-p_dict[perm[-1]]) * internal * H_comp
                        elif perm[-1] == w:
                            alpha_I_n += 2 * (-q_dict[perm[0]]) * internal * H_comp

            # alpha_I at n-1
            alpha_I_prev = -2 * H_total(T, B_w_prev)
            for chain_len in range(3, len(W_prev) + 1, 2):
                for combo in combinations(W_prev, chain_len):
                    if w not in combo:
                        continue
                    comp = [x for x in W_prev if x not in combo]
                    H_comp = H_total(T, comp)
                    for perm in permutations(combo):
                        internal = 1.0
                        for k in range(chain_len - 1):
                            internal *= T(perm[k], perm[k + 1])
                        if perm[0] == w:
                            alpha_I_prev += 2 * (-p_dict[perm[-1]]) * internal * H_comp
                        elif perm[-1] == w:
                            alpha_I_prev += 2 * (-q_dict[perm[0]]) * internal * H_comp

            delta_H = alpha_n - alpha_prev
            delta_I = alpha_I_n - alpha_I_prev

            print(f"  trial {trial}, w={w}, u_new={u_new}:")
            print(f"    alpha_H(n) = {alpha_n:.4f}, alpha_H(n-1) = {alpha_prev:.4f}, "
                  f"delta_H = {delta_H:.4f}")
            print(f"    alpha_I(n) = {alpha_I_n:.4f}, alpha_I(n-1) = {alpha_I_prev:.4f}, "
                  f"delta_I = {delta_I:.4f}")
            print(f"    delta_H - delta_I = {delta_H - delta_I:.2e}")
    print()


if __name__ == "__main__":
    for n in [5, 6, 7]:
        test_inductive_relation(n, 3)

    for n in [5, 6]:
        analyze_first_step(n, 3)
