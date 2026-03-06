#!/usr/bin/env python3
"""
Attack the unsigned sum: prove delta_H = delta_I for all n.

We know: delta_H = sum_S Delta(S,R) where Delta(S,R) = Li(S)*Rj(R) - Lj(S)*Ri(R).

Key idea: Delta(S,R) can be decomposed using the bracket identity:
  B(u,w) := p_u*q_w - (1-q_u)*(1-p_w) = (p_u+q_w-1) - (q_u+p_w-1)*...
Actually: B(u,w) = p_u*q_w - (1-q_u)(1-p_w)
         = p_u*q_w - 1 + q_u + p_w - q_u*p_w
         = (p_u*q_w - q_u*p_w) + (q_u + p_w - 1)

In terms of s,d: B(u,w) = -[s_u(1-d_w) + s_w(1+d_u)]/2

Approach: express delta_H = sum_S Li(S)*Rj(R) - Lj(S)*Ri(R) explicitly
using the bracket decomposition, then match to delta_I = sum_k 2^k Delta(alpha_k).

The relationship Li(S)*Rj(R) counts "Li paths on S ending at i, then i→j, then Rj paths
on R starting at j". The key transition is i→j, contributing factor 1 (since T[i][j]=1).

Let me try a completely different approach: use the TRANSFER MATRIX.

For a tournament T on V with arc i→j, define the transfer matrix M where M[a][b] = T(a,b).
Then H(T) = sum of all Hamiltonian path weights = permanent-like quantity.

adj(i,j) = sum over permutations with i immediately before j.

Can we express delta_H = adj(i,j) - adj(j,i,T') using determinantal identities?

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
import random
from fractions import Fraction


def h_end(T, verts, v):
    if len(verts) == 1:
        return 1 if v == verts[0] else 0
    total = 0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        w = 1
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def h_start(T, verts, v):
    if len(verts) == 1:
        return 1 if v == verts[0] else 0
    total = 0
    for p in permutations(verts):
        if p[0] != v:
            continue
        w = 1
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def test_bracket_expansion(n, num_trials=20):
    """
    Expand Delta(S,R) using bracket decomposition.

    Delta(S,R) = Li(S)*Rj(R) - Lj(S)*Ri(R)

    For |S|=1, S={u}, R=W\{u}:
    Delta({u}, R) = p_u * Rj(R) - (1-q_u) * Ri(R)

    For |S|=0: Delta(∅, W) = 1*Rj(W) - 1*Ri(W) = sum_w q_w*h(W\{w}→w) - sum_w (1-p_w)*h(W\{w}→w)
    Wait, Rj(W) = h_start(T, {j}∪W, j) = sum_{w∈W} T(j,w)*h_start(T, W\{w}∪..., ...)
    Hmm, this is getting complicated. Let me just look at specific cases.
    """
    print(f"=== Bracket expansion at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a+1, n) if (a,b) != (I,J)]

    random.seed(42)

    for trial in range(min(num_trials, 3)):
        arc_values = {(I,J): 1, (J,I): 0}
        for (a,b) in free_arcs:
            v = random.randint(0, 1)
            arc_values[(a,b)] = v
            arc_values[(b,a)] = 1 - v

        def T(a, b):
            if a == b: return 0
            return arc_values.get((a,b), 1 - arc_values.get((b,a), 0))

        # Compute s and d values
        s = {w: 1 - T(w, I) - T(J, w) for w in W}
        d = {w: T(w, I) - T(J, w) for w in W}

        # Compute delta_H
        delta_H = 0
        for smask in range(1 << m):
            S = [W[bit] for bit in range(m) if smask & (1 << bit)]
            R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]

            Li_S = h_end(T, S + [I], I)
            Rj_R = h_start(T, [J] + R, J)
            Lj_S = h_end(T, S + [J], J)
            Ri_R = h_start(T, [I] + R, I)

            delta_H += Li_S * Rj_R - Lj_S * Ri_R

        # Compute delta_I = 2 * sum_x (-s_x) for n=4,5 (simplified formula)
        delta_I_simple = -2 * sum(s[w] for w in W)

        # For n >= 6, also need 5-cycle terms
        # For now, just compare with delta_H

        print(f"  trial {trial}:")
        print(f"    s = {s}")
        print(f"    d = {d}")
        print(f"    delta_H = {delta_H}")
        if n <= 5:
            print(f"    -2*sum(s) = {delta_I_simple}")
            print(f"    Match: {delta_H == delta_I_simple}")


def test_inductive_approach(n, num_trials=30):
    """
    Try to prove delta_H = delta_I by induction on n.

    For the inductive step, remove a vertex w from W.
    Does delta_H for (n) decompose in terms of delta_H for (n-1)?

    More precisely: delta_H = sum_S Delta(S,R).
    Split by whether w ∈ S or w ∈ R:

    delta_H = sum_{S not containing w} Delta(S,R) + sum_{S containing w} Delta(S,R)

    In the first case, w ∈ R. In the second, w ∈ S.
    """
    print(f"\n=== Inductive decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a+1, n) if (a,b) != (I,J)]

    random.seed(42)

    for trial in range(min(num_trials, 5)):
        arc_values = {(I,J): 1, (J,I): 0}
        for (a,b) in free_arcs:
            v = random.randint(0, 1)
            arc_values[(a,b)] = v
            arc_values[(b,a)] = 1 - v

        def T(a, b):
            if a == b: return 0
            return arc_values.get((a,b), 1 - arc_values.get((b,a), 0))

        s = {w: 1 - T(w, I) - T(J, w) for w in W}

        # Full delta_H
        delta_H = 0
        for smask in range(1 << m):
            S = [W[bit] for bit in range(m) if smask & (1 << bit)]
            R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
            Li = h_end(T, S + [I], I)
            Rj = h_start(T, [J] + R, J)
            Lj = h_end(T, S + [J], J)
            Ri = h_start(T, [I] + R, I)
            delta_H += Li * Rj - Lj * Ri

        # Split by vertex w = W[-1]
        w_target = W[-1]
        W_minus_w = [x for x in W if x != w_target]
        m2 = len(W_minus_w)

        # Part 1: w in R (not in S)
        part_wR = 0
        for smask in range(1 << m2):
            S = [W_minus_w[bit] for bit in range(m2) if smask & (1 << bit)]
            R_minus_w = [x for x in W_minus_w if x not in S]
            R = R_minus_w + [w_target]
            Li = h_end(T, S + [I], I)
            Rj = h_start(T, [J] + R, J)
            Lj = h_end(T, S + [J], J)
            Ri = h_start(T, [I] + R, I)
            part_wR += Li * Rj - Lj * Ri

        # Part 2: w in S (in S)
        part_wS = 0
        for smask in range(1 << m2):
            S_minus_w = [W_minus_w[bit] for bit in range(m2) if smask & (1 << bit)]
            S = S_minus_w + [w_target]
            R = [x for x in W_minus_w if x not in S_minus_w]
            Li = h_end(T, S + [I], I)
            Rj = h_start(T, [J] + R, J)
            Lj = h_end(T, S + [J], J)
            Ri = h_start(T, [I] + R, I)
            part_wS += Li * Rj - Lj * Ri

        # delta_H for the reduced tournament (W minus w)
        delta_H_reduced = 0
        for smask in range(1 << m2):
            S = [W_minus_w[bit] for bit in range(m2) if smask & (1 << bit)]
            R = [x for x in W_minus_w if x not in S]
            Li = h_end(T, S + [I], I)
            Rj = h_start(T, [J] + R, J)
            Lj = h_end(T, S + [J], J)
            Ri = h_start(T, [I] + R, I)
            delta_H_reduced += Li * Rj - Lj * Ri

        print(f"  trial {trial}: s_w={s[w_target]}")
        print(f"    delta_H = {delta_H}")
        print(f"    part(w in R) = {part_wR}")
        print(f"    part(w in S) = {part_wS}")
        print(f"    delta_H_reduced (W\\{{w}}) = {delta_H_reduced}")

        # Can we express part_wR and part_wS in terms of delta_H_reduced and s_w?
        # If delta_H = delta_H_reduced + something depending on s_w and d_w...
        # That would give an inductive formula.

        if delta_H_reduced != 0:
            ratio_R = part_wR / delta_H_reduced if delta_H_reduced != 0 else "N/A"
            ratio_S = part_wS / delta_H_reduced if delta_H_reduced != 0 else "N/A"
        else:
            ratio_R = "div0"
            ratio_S = "div0"
        print(f"    part_wR / delta_H_red = {ratio_R}")
        print(f"    part_wS / delta_H_red = {ratio_S}")
        print()


def test_last_step_decomposition(n, num_trials=20):
    """
    Use first/last step decomposition to relate delta_H to sub-tournament quantities.

    adj(i,j) = sum over Ham paths (..., a, i, j, b, ...) where a is the vertex before i
    and b is the vertex after j.

    Decompose by the pair (a, b):
    adj(i,j) = sum_{a,b ∈ W, a≠b} T(a,i) * T(j,b) * [Ham paths on W\{a,b} connecting a's predecessor to b's successor]
    + boundary terms (a=none when i is first, b=none when j is last)

    More precisely:
    adj(i,j) = sum_S h_end(S∪{i}, i) * h_start({j}∪R, j)

    Using last-step: h_end(S∪{i}, i) = sum_{a∈S} T(a,i) * h_end(S, a) + [S=∅]
    (i.e., the vertex just before i in the left part is a, or S is empty)

    Using first-step: h_start({j}∪R, j) = sum_{b∈R} T(j,b) * h_start(R, b) + [R=∅]

    So: Li(S)*Rj(R) = (sum_a T(a,i)*h_end(S,a) + [S=∅]) * (sum_b T(j,b)*h_start(R,b) + [R=∅])
    = sum_{a∈S,b∈R} p_a * q_b * h_end(S,a) * h_start(R,b)
    + [S=∅] * sum_b q_b * h_start(R,b)
    + [R=∅] * sum_a p_a * h_end(S,a)
    + [S=∅,R=∅]  (impossible since m≥1 typically)

    Similarly for Lj and Ri, with p_a→(1-q_a) and q_b→(1-p_b):
    Lj(S)*Ri(R) = sum_{a,b} (1-q_a)(1-p_b) h_end(S,a) h_start(R,b) + boundaries

    So Delta(S,R) = sum_{a∈S,b∈R} [p_a*q_b - (1-q_a)(1-p_b)] h_end(S,a) h_start(R,b) + boundary Delta

    Using bracket: p_a*q_b - (1-q_a)(1-p_b) = -[s_a(1-d_b)+s_b(1+d_a)]/2

    This means: delta_H = sum_S sum_{a∈S,b∈R} bracket(a,b) * h_end(S,a)*h_start(R,b) + boundary

    Exchanging summation: for each pair (a,b):
    delta_H contribution from (a,b) = bracket(a,b) * sum_{S: a∈S, b∉S} h_end(S,a)*h_start(R,b)

    THIS is what we need to evaluate and match to delta_I!
    """
    print(f"\n=== Last-step decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a+1, n) if (a,b) != (I,J)]

    random.seed(42)

    for trial in range(min(num_trials, 3)):
        arc_values = {(I,J): 1, (J,I): 0}
        for (a,b) in free_arcs:
            v = random.randint(0, 1)
            arc_values[(a,b)] = v
            arc_values[(b,a)] = 1 - v

        def T(a, b):
            if a == b: return 0
            return arc_values.get((a,b), 1 - arc_values.get((b,a), 0))

        s = {w: 1 - T(w, I) - T(J, w) for w in W}
        d = {w: T(w, I) - T(J, w) for w in W}

        # Full delta_H
        delta_H = 0
        for smask in range(1 << m):
            S = [W[bit] for bit in range(m) if smask & (1 << bit)]
            R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
            Li = h_end(T, S + [I], I)
            Rj = h_start(T, [J] + R, J)
            Lj = h_end(T, S + [J], J)
            Ri = h_start(T, [I] + R, I)
            delta_H += Li * Rj - Lj * Ri

        # Decompose by (a,b) pair
        pair_contrib = {}
        for a in W:
            for b in W:
                if a == b:
                    continue
                # bracket(a,b) = -(s_a(1-d_b) + s_b(1+d_a))/2
                bracket = -(s[a]*(1-d[b]) + s[b]*(1+d[a])) / 2

                # Contribution: bracket * sum_{S: a∈S, b∉S} h_end(S,a)*h_start(R,b)
                # where R = W\S, and S contains a but not b.
                # S ranges over subsets of W\{b} containing a.
                W_minus_ab = [w for w in W if w != a and w != b]
                contrib_sum = 0
                for smask in range(1 << len(W_minus_ab)):
                    S_rest = [W_minus_ab[bit] for bit in range(len(W_minus_ab)) if smask & (1 << bit)]
                    S = S_rest + [a]
                    R = [w for w in W if w not in S]  # contains b
                    contrib_sum += h_end(T, list(S), a) * h_start(T, list(R), b)

                pair_contrib[(a,b)] = bracket * contrib_sum

        # Boundary terms: S=∅ or R=∅
        # S=∅: Delta = Rj(W) - Ri(W) = sum_b q_b h_start(W\{b},b) - sum_b (1-p_b) h_start(W\{b},b)
        #     = sum_b (q_b - 1 + p_b) h_start(W\{b},b) = sum_b (-s_b) h_start(W\{b},b)
        boundary_S_empty = 0
        for b in W:
            hR = h_start(T, [w for w in W if w != b] + [b], b) if m > 1 else 1
            # Actually h_start(W, b) since W = R when S=∅
            hR = h_start(T, list(W), b)
            boundary_S_empty += (-s[b]) * hR

        # R=∅: Delta = Li(W) - Lj(W) = sum_a p_a h_end(W\{a},a) - sum_a (1-q_a) h_end(W\{a},a)
        #     = sum_a (p_a - 1 + q_a) h_end(W\{a},a) = sum_a (-s_a) h_end(W\{a},a)
        boundary_R_empty = 0
        for a in W:
            hL = h_end(T, list(W), a)
            boundary_R_empty += (-s[a]) * hL

        total_from_pairs = sum(pair_contrib.values())
        total_decomp = total_from_pairs + boundary_S_empty + boundary_R_empty

        print(f"  trial {trial}:")
        print(f"    delta_H = {delta_H}")
        print(f"    pair contributions = {total_from_pairs}")
        print(f"    boundary (S=∅) = {boundary_S_empty}")
        print(f"    boundary (R=∅) = {boundary_R_empty}")
        print(f"    total decomp = {total_decomp}")
        print(f"    err = {abs(delta_H - total_decomp):.2e}")

        # Now, can we group this by s_w?
        # delta_H = sum_w s_w * [...] + sum_{w1,w2} s_w1 * s_w2 * [...]
        # From bracket: bracket(a,b) involves s_a and s_b linearly.
        # So the s-structure of delta_H is:
        # delta_H = sum_w (-s_w) * [boundary_w + sum_u cross_wu] + higher order

        # Let's compute the s_w-coefficient of delta_H
        print(f"    s-coefficients of delta_H:")
        for w in W:
            # Coefficient of s_w: comes from bracket(a,b) where a=w or b=w
            # plus boundary terms
            coeff = 0

            # From bracket(w,b) for b != w: bracket = -(s_w(1-d_b) + s_b(1+d_w))/2
            # The s_w-linear part: -(1-d_b)/2 * contrib_sum(w,b)
            for b in W:
                if b == w: continue
                W_minus_wb = [x for x in W if x != w and x != b]
                cs = 0
                for smask in range(1 << len(W_minus_wb)):
                    S_rest = [W_minus_wb[bit] for bit in range(len(W_minus_wb)) if smask & (1 << bit)]
                    S = S_rest + [w]
                    R = [x for x in W if x not in S]
                    cs += h_end(T, list(S), w) * h_start(T, list(R), b)
                coeff += -(1-d[b])/2 * cs

            # From bracket(a,w) for a != w: bracket = -(s_a(1-d_w) + s_w(1+d_a))/2
            # The s_w-linear part: -(1+d_a)/2 * contrib_sum(a,w)
            for a in W:
                if a == w: continue
                W_minus_aw = [x for x in W if x != a and x != w]
                cs = 0
                for smask in range(1 << len(W_minus_aw)):
                    S_rest = [W_minus_aw[bit] for bit in range(len(W_minus_aw)) if smask & (1 << bit)]
                    S = S_rest + [a]
                    R = [x for x in W if x not in S]
                    cs += h_end(T, list(S), a) * h_start(T, list(R), w)
                coeff += -(1+d[a])/2 * cs

            # Boundary: -s_w * h_start(W, w) (from S=∅)
            coeff += -h_start(T, list(W), w)
            # Boundary: -s_w * h_end(W, w) (from R=∅)
            coeff += -h_end(T, list(W), w)

            print(f"      s_{w}: coeff = {coeff:.4f}, s_{w} = {s[w]}")

        # delta_I (simplified, n<=5)
        if n <= 5:
            delta_I = -2 * sum(s[w] for w in W)
            print(f"    delta_I = {delta_I} = -2*sum(s)")
            print(f"    → s_w coefficient in delta_I = -2 for all w")

        print()


if __name__ == "__main__":
    for n in [4, 5]:
        test_bracket_expansion(n, num_trials=5)

    for n in [4, 5]:
        test_last_step_decomposition(n, num_trials=3)
