#!/usr/bin/env python3
"""
Attempt an INDUCTIVE PROOF of alpha_w^H = alpha_w^I for all n.

The full OCF formula (THM-013) gives:
  delta_I = sum_{k>=1} 2^k * Delta(alpha_k)

where Delta(alpha_k) = sum_C [alpha_{k-1}(comp(C_gained)) - alpha_{k-1}(comp(C_lost))]
and comp(C) = V \ V(C) (the complement of the cycle).

By STRONG INDUCTION: assume OCF for all m < n. Then for sub-tournaments on comp(C):
  alpha_{k-1}(comp(C)) is determined by H values of sub-sub-tournaments.

At the k=1 level (individual cycles):
  Delta(alpha_1) = sum_L (D_L - C_L)
  = sum_w (-s_w) + sum_{L>=5} (D_L - C_L)   [3-cycle + longer cycles]

At the k=2 level (VD cycle pairs):
  Delta(alpha_2) = sum_C [alpha_1(comp(C_gained)) - alpha_1(comp(C_lost))]
  = sum_{3-cycles C} (-s_w) * alpha_1(B_w)
    + sum_{5-cycles C'} bracket(C') * alpha_1(comp(C'))
    + ...

By inductive OCF on comp(C): alpha_1(comp(C)) = (H(comp(C)) - 1) / 2.

So 2*alpha_1(comp(C)) = H(comp(C)) - 1.

Combining: delta_I = 2*Delta(alpha_1) + 4*Delta(alpha_2) + 8*Delta(alpha_3) + ...

For the 3-cycle contribution to all levels:
  2*(-s_w) + 4*(-s_w)*alpha_1(B_w) + 8*(-s_w)*alpha_2(B_w) + ...
  = (-s_w) * 2 * sum_{k>=0} 2^k * alpha_k(B_w)
  = (-s_w) * 2 * I(Omega(B_w), 2)
  = -2*s_w * H(B_w)   [by inductive OCF on B_w]

This is EXACTLY the -2*sum s_w * H(B_w) term! And it holds for ALL n by induction!

For the 5-cycle contribution:
  2*(D5-C5) + 4*sum_{5-cycles} bracket * alpha_1(comp) + ...
  = 2 * sum_{chains of length 3} bracket(w_3, w_1) * sum_{k>=0} 2^k * alpha_k(comp)
  = 2 * sum_chains bracket * H(comp)   [by inductive OCF]

The coefficient of s_w in this:
  2 * sum_chains d(bracket)/ds_w * H(comp)
  = 2 * sum_chains (-p or -q) * H(comp)

where comp has n-5 vertices (W \ {w1,w2,w3} if cycle is {i,j,w1,w2,w3}).

This is a GENERALIZED cycle derivative that includes H(comp) factors!

So the FULL alpha_w^I formula for all n is:

  alpha_w^I = -2*H(B_w)
            + sum_{L=5,7,...} 2 * sum_{chains of length L-2 with w at endpoint}
              (-p or -q) * internal * H(comp(chain))

where comp(chain) = B_w \ {vertices in chain other than w} for chains starting with w,
or = W \ {vertices in chain} for chains ending with w.

Let me verify this formula numerically at n=8.

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


def compute_alpha_H(T, W, w_target):
    """Compute alpha_w^H via insertion decomposition (verified at all n)."""
    I, J = 0, 1
    p = {w: T(w, I) for w in W}
    q = {w: T(J, w) for w in W}
    B_w = [x for x in W if x != w_target]

    total = 0.0
    for pi in permutations(B_w):
        B_wt = 1.0
        for k in range(len(pi) - 1):
            B_wt *= T(pi[k], pi[k + 1])

        m1 = len(pi)
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
            total += val
    return total


def compute_alpha_I_full(T, W, w_target):
    """Compute alpha_w^I using the FULL inductive formula.

    alpha_w^I = -2*H(B_w)
              + sum_{chains of odd length >= 3 with w at endpoint}
                2 * (-p or -q) * internal * H(comp)

    where comp = W \ {chain vertices} (complement in W, not V).
    """
    I, J = 0, 1
    p = {w: T(w, I) for w in W}
    q = {w: T(J, w) for w in W}
    m = len(W)

    B_w = [x for x in W if x != w_target]
    alpha = -2 * H_total(T, B_w)

    # Chains of odd length >= 3 with w at endpoint
    for chain_len in range(3, m + 1, 2):  # 3, 5, 7, ...
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
                    # d/ds_w of bracket at w = first: -p_{last}
                    alpha += 2 * (-p[perm[-1]]) * internal * H_comp
                elif perm[-1] == w_target:
                    # d/ds_w of bracket at w = last: -q_{first}
                    alpha += 2 * (-q[perm[0]]) * internal * H_comp

    return alpha


def verify_full_formula(n, num_trials=20):
    """Verify the full inductive alpha_w^I formula matches alpha_w^H."""
    print(f"\n=== Full inductive formula verification at n={n} ===\n")

    W = list(range(2, n))
    max_err = 0.0

    for trial in range(num_trials):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:min(3, len(W))]:
            aH = compute_alpha_H(T, W, w)
            aI = compute_alpha_I_full(T, W, w)
            err = abs(aH - aI)
            max_err = max(max_err, err)

            if trial < 2 and w == W[0]:
                print(f"  trial {trial}, w={w}: alpha_H={aH:.6f}, alpha_I_full={aI:.6f}, "
                      f"err={err:.2e}")

    print(f"  Max |alpha_H - alpha_I|: {max_err:.2e}")
    if max_err < 1e-8:
        print(f"  CONFIRMED at n={n}\n")
    else:
        print(f"  FAILED at n={n}\n")
    return max_err < 1e-8


if __name__ == "__main__":
    for n in range(4, 10):
        ok = verify_full_formula(n, num_trials=20 if n <= 7 else 10)
        if not ok:
            break
