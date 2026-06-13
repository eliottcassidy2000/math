#!/usr/bin/env python3
"""
Prove alpha_w^H = alpha_w^I for all n by algebraic decomposition.

STRATEGY: At s=0, express alpha_w^H as a sum over permutations of W,
then reorganize to match the cycle structure of alpha_w^I.

alpha_w^H = -h_start(W,w) - h_end(W,w)
           - sum_{S:w in S, R!=empty} h_end(S,w) * sum_b q_b * h_start(R,b)
           - sum_{S:w not in S, S!=empty} sum_a p_a * h_end(S,a) * h_start(R,w)

At s=0: p_a + q_a = 1, so p_a = T(a,i), q_a = T(j,a) = 1 - T(a,j).

alpha_w^I = -2*H(B_w) + sum_{L>=5, odd} cycle_derivative_L(w)

KEY QUESTION: Can we decompose alpha_w^H into -2*H(B_w) + (something that matches cycle derivatives)?

APPROACH 1: Rewrite alpha_w^H by combining the boundary and intermediate terms
into a sum over "paths through W that pass through w" times connection factors.

APPROACH 2: Use the identity h_start(V,v) + h_end(V,v) = H(V) - H(V\{v}) + ???
(This doesn't hold in general.)

APPROACH 3: Try induction on |W|. At n=4 (|W|=2): alpha = -2. At n=5: add vertex u
and track how alpha changes.

APPROACH 4: Expand everything as sums over permutations and identify cancellations.

Instance: opus-2026-03-05-S5
"""

from itertools import permutations, combinations
import random
from collections import defaultdict


def make_T(n, arc_values):
    """Create tournament function with s=0 constraint."""
    def T(a, b):
        if a == b: return 0.0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))
    return T


def random_s0_tournament(n, seed=None):
    """Generate a random tournament with s_w = 0 for all w in W.

    At s=0: p_w + q_w = 1, i.e., T(w,i) + T(j,w) = 1 for all w.
    This means T(w,0) + T(1,w) = 1, i.e., T(w,0) = 1 - T(1,w).
    """
    if seed is not None:
        random.seed(seed)
    I, J = 0, 1
    W = list(range(2, n))

    arc_values = {(I, J): 1.0, (J, I): 0.0}

    # For each w in W: T(w,I) + T(J,w) = 1
    # Free parameter: T(J,w) = q_w in (0,1), then T(w,I) = 1 - q_w
    for w in W:
        q_w = random.uniform(0.1, 0.9)
        arc_values[(J, w)] = q_w
        arc_values[(w, J)] = 1 - q_w
        arc_values[(w, I)] = 1 - q_w  # p_w = 1 - q_w
        arc_values[(I, w)] = q_w      # T(I,w) = 1 - T(w,I)

    # Internal arcs: free
    for a in W:
        for b in W:
            if a < b:
                v = random.uniform(0.1, 0.9)
                arc_values[(a, b)] = v
                arc_values[(b, a)] = 1 - v

    return arc_values


def h_paths(T, verts):
    """Return all Hamiltonian paths as (weight, start, end) tuples."""
    if len(verts) == 0:
        return []
    if len(verts) == 1:
        return [(1.0, verts[0], verts[0])]
    results = []
    for p in permutations(verts):
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        results.append((w, p[0], p[-1]))
    return results


def h_end(T, verts, v):
    return sum(w for w, s, e in h_paths(T, verts) if e == v)


def h_start(T, verts, v):
    return sum(w for w, s, e in h_paths(T, verts) if s == v)


def H_total(T, verts):
    return sum(w for w, s, e in h_paths(T, verts))


def compute_alpha_H(T, W, w_target):
    """Compute alpha_w^H at s=0.

    Returns (total, boundary, left, right) and also a detailed decomposition.
    """
    I, J = 0, 1
    p = {w: T(w, I) for w in W}
    q = {w: T(J, w) for w in W}

    # Verify s=0
    for w in W:
        assert abs(p[w] + q[w] - 1) < 1e-10, f"s_{w} = {1-p[w]-q[w]} != 0"

    boundary = -h_start(T, list(W), w_target) - h_end(T, list(W), w_target)

    W_minus_w = [x for x in W if x != w_target]

    left = 0.0
    for smask in range(1 << len(W_minus_w)):
        S_rest = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        if len(S_rest) == len(W_minus_w):
            continue
        S = S_rest + [w_target]
        R = [x for x in W if x not in S]
        h_S_w = h_end(T, list(S), w_target)
        for b in R:
            left += -q[b] * h_S_w * h_start(T, list(R), b)

    right = 0.0
    for smask in range(1, 1 << len(W_minus_w)):
        S = [W_minus_w[bit] for bit in range(len(W_minus_w)) if smask & (1 << bit)]
        R_rest = [x for x in W_minus_w if x not in S]
        R = R_rest + [w_target]
        h_R_w = h_start(T, list(R), w_target)
        for a in S:
            right += -p[a] * h_end(T, list(S), a) * h_R_w

    return boundary + left + right, boundary, left, right


def compute_alpha_I(T, W, w_target):
    """Compute alpha_w^I at s=0 from the cycle formula.

    alpha_w^I = -2*H(B_w) + sum of cycle derivatives for L>=5.

    For an L-cycle (i,j,w1,...,w_{L-2}), the derivative with respect to s_w
    hits only when w is at position w1 (first, adjacent to j) or w_{L-2} (last,
    adjacent to i).

    d(bracket(w_{L-2}, w_1))/ds_{w1} = -p_{w_{L-2}} at s=0
    d(bracket(w_{L-2}, w_1))/ds_{w_{L-2}} = -q_{w_1} at s=0
    """
    I, J = 0, 1
    p = {w: T(w, I) for w in W}
    q = {w: T(J, w) for w in W}

    B_w = [x for x in W if x != w_target]
    alpha_3 = -2 * H_total(T, B_w)

    alpha_cycles = 0.0
    m = len(W)

    # For each odd L = 2k+1 where k >= 2 (L >= 5):
    # L-2 vertices from W form the internal chain
    for chain_len in range(3, m + 1, 2):  # 3, 5, 7, ... (for L=5,7,9,...)
        for combo in combinations(W, chain_len):
            for perm in permutations(combo):
                # perm = (w1, w2, ..., w_{chain_len})
                # Internal product: T(w1,w2)*T(w2,w3)*...*T(w_{k-1},w_k)
                internal = 1.0
                for k in range(chain_len - 1):
                    internal *= T(perm[k], perm[k + 1])

                # Bracket derivative with respect to s_w
                if w_target == perm[0]:
                    # d/ds_{w1} of bracket(w_last, w1) = -p_{w_last}
                    alpha_cycles += 2 * (-p[perm[-1]]) * internal
                elif w_target == perm[-1]:
                    # d/ds_{w_last} of bracket(w_last, w1) = -q_{w1}
                    alpha_cycles += 2 * (-q[perm[0]]) * internal

    return alpha_3 + alpha_cycles, alpha_3, alpha_cycles


def test_decomposition_by_permutation(n, num_trials=10):
    """Decompose alpha_w^H as a sum over PERMUTATIONS of W.

    Each term in alpha_w^H is a product of tournament arcs along a path in S
    times a product along a path in R times a connection factor (p or q).

    The concatenation (path in S)·(connection)·(path in R) forms a
    "decorated permutation" of W with w at a specific position and a
    connection at the S|R boundary.

    Can we rewrite this as a sum over permutations of W?
    """
    print(f"\n=== Permutation decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:2]:
            aH, bnd, left, right = compute_alpha_H(T, W, w)
            aI, a3, acyc = compute_alpha_I(T, W, w)

            B_w = [x for x in W if x != w]
            H_Bw = H_total(T, B_w)

            # Decompose alpha_w^H by permutation of W
            # A permutation sigma of W has w at position k (0-indexed).
            # Group contributions by k.
            perm_contrib = defaultdict(float)

            for perm in permutations(W):
                pos_w = list(perm).index(w)

                # Path weight in W
                path_wt = 1.0
                for k in range(len(perm) - 1):
                    path_wt *= T(perm[k], perm[k + 1])

                # What's the contribution of this permutation to alpha_w^H?
                # boundary: -h_start(W,w) contributes -path_wt when pos_w=0
                #           -h_end(W,w) contributes -path_wt when pos_w=m-1
                if pos_w == 0:
                    perm_contrib['start'] += -path_wt
                if pos_w == m - 1:
                    perm_contrib['end'] += -path_wt

            # Intermediate terms: not directly a single permutation of W
            # because they involve splitting W into S and R with a p/q connection

            print(f"  trial {trial}, w={w}:")
            print(f"    alpha_H = {aH:.6f} = {bnd:.4f} + {left:.4f} + {right:.4f}")
            print(f"    alpha_I = {aI:.6f} = {a3:.4f} + {acyc:.4f}")
            print(f"    -2*H(B_w) = {-2*H_Bw:.4f}")
            print(f"    err = {abs(aH-aI):.2e}")
            print(f"    boundary = -h_start - h_end = {bnd:.4f}")
            for key in sorted(perm_contrib.keys()):
                print(f"    perm[{key}] = {perm_contrib[key]:.4f}")
    print()


def test_alpha_by_cut(n, num_trials=10):
    """Decompose alpha_w^H by the CUT POSITION.

    For each subset S containing w (and R = W\S):
    - The path in S ends at w: path (...→a→w) in S
    - Connection: q_b = T(j,b) for some b in R
    - The path in R starts at b: path (b→...→c) in R

    The full "trajectory" is: (start_S→...→a→w) | q_b | (b→...→end_R)

    This is like a permutation of W with a "cut" after position |S|-1,
    where w is at the end of the left part, and the cut introduces a
    factor of q_b instead of T(w, b).

    Similarly for w in R: (start_S→...→a) | p_a | (w→...→end_R)
    The cut introduces p_a instead of T(a, w).

    KEY INSIGHT: The "cut" replaces T(w,b) with q_b = T(j,b) or T(a,w) with p_a = T(a,i).
    So we can write:

    alpha_w^H = -sum over perms of W where w is at position k *
                  [path_weight * SOMETHING involving the cut at position k]
    """
    print(f"\n=== Cut decomposition at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:2]:
            aH, bnd, left, right = compute_alpha_H(T, W, w)
            aI, a3, acyc = compute_alpha_I(T, W, w)

            # Rewrite intermediate terms as sums over permutations with cuts
            # For w in S: S-path ends at w, R-path starts at b
            #   Contribution: -q_b * weight(S-path ending at w) * weight(R-path starting at b)
            #   = -T(j,b) * [product of T along S-path] * [product of T along R-path]
            #   Compare with permutation (S-path)(R-path) where arc from w to b would be T(w,b):
            #   That perm weight = [S-path weight] * T(w,b) * [R-path weight]
            #   Ratio: contribution / perm_weight = -T(j,b) / T(w,b) ... messy

            # Alternative: express as sum over perm of W with w at position k,
            #   and for each cut position c after element k, multiply by correction factor.

            # Actually, let's think about this differently.
            # For w in S with |S| = s_size, the S-path is a permutation of S ending at w,
            # and the R-path is a permutation of R starting at some b.
            # These together form a permutation of W where:
            #   - w is at position s_size - 1 (0-indexed)
            #   - The "cut" is after position s_size - 1
            #   - The arc from w to b is REPLACED by -q_b

            # So alpha_w^H(left) = sum over perms sigma of W with w at position k (0<=k<=m-2)
            #   * [product of T(sigma(i), sigma(i+1)) for i != k] * (-q_{sigma(k+1)})

            # And alpha_w^H(right) = sum over perms sigma of W with w at position k (1<=k<=m-1)
            #   * [product of T(sigma(i), sigma(i+1)) for i != k-1] * (-p_{sigma(k-1)})

            # And boundary:
            #   -h_start(W,w): w at position 0, full path weight, factor -1
            #   -h_end(W,w): w at position m-1, full path weight, factor -1

            # Let's verify this reformulation
            rewritten = 0.0
            for perm in permutations(W):
                pos_w = list(perm).index(w)

                # Full path weight (product of all arcs)
                full_wt = 1.0
                for k in range(m - 1):
                    full_wt *= T(perm[k], perm[k + 1])

                # Boundary contributions
                if pos_w == 0:
                    rewritten += -full_wt  # -h_start(W,w)
                if pos_w == m - 1:
                    rewritten += -full_wt  # -h_end(W,w)

                # Left contribution: cut after w (w at position pos_w, next is perm[pos_w+1])
                # Replace T(w, perm[pos_w+1]) with -q_{perm[pos_w+1]}
                if pos_w < m - 1:
                    b = perm[pos_w + 1]
                    # Path weight with T(w,b) removed and -q_b inserted
                    cut_wt = 1.0
                    for k in range(m - 1):
                        if k == pos_w:
                            continue  # skip the w→b arc
                        cut_wt *= T(perm[k], perm[k + 1])
                    rewritten += cut_wt * (-T(J, b))

                # Right contribution: cut before w (w at position pos_w, prev is perm[pos_w-1])
                # Replace T(perm[pos_w-1], w) with -p_{perm[pos_w-1]}
                if pos_w > 0:
                    a = perm[pos_w - 1]
                    cut_wt = 1.0
                    for k in range(m - 1):
                        if k == pos_w - 1:
                            continue  # skip the a→w arc
                        cut_wt *= T(perm[k], perm[k + 1])
                    rewritten += cut_wt * (-T(a, I))

            err_rewrite = abs(rewritten - aH)
            print(f"  trial {trial}, w={w}: alpha_H={aH:.6f}, rewritten={rewritten:.6f}, "
                  f"err={err_rewrite:.2e}")
    print()


def test_combined_formula(n, num_trials=10):
    """After the cut-based rewriting, alpha_w^H becomes:

    alpha_w^H = sum over perms sigma of W, with w at position k:
      CASE k=0:       -full_weight - full_weight_minus_arc(0) * q_{sigma(1)}
      CASE 0<k<m-1:   -cut_left(k) * q_{next} - cut_right(k) * p_{prev}
      CASE k=m-1:     -full_weight - full_weight_minus_arc(m-2) * p_{sigma(m-2)}

    Combine: for each perm, the contribution involves the full path weight times
    replacement factors at the arcs adjacent to w.

    For a perm (..., a, w, b, ...):
      -q_b * (path weight without T(w,b)) - p_a * (path weight without T(a,w))
      = -(path weight / T(w,b)) * q_b - (path weight / T(a,w)) * p_a
      = -path_except_w * [T(a,w)*q_b + T(w,b)*p_a]  ... not quite

    Actually:
      path_weight = L * T(a,w) * T(w,b) * R   (L=left of a, R=right of b)
      Left cut: L * T(a,w) * (-q_b) * R  (remove T(w,b), insert -q_b)
      Right cut: L * (-p_a) * T(w,b) * R  (remove T(a,w), insert -p_a)
      Sum: -L*R * [T(a,w)*q_b + p_a*T(w,b)]

    Plus boundary (w at start or end):
      w at start: perm = (w, b, ...). -full - (remove T(w,b), insert -q_b)
        = -T(w,b)*R - (-q_b)*R = -R*(T(w,b) - q_b) = -R*(T(w,b)-T(j,b))
        Hmm, -full + cut = -T(w,b)*R + q_b*R = R*(q_b - T(w,b))
        No wait: boundary gives -full_wt. Left cut gives -q_b * R (where R is the rest).
        Total: -T(w,b)*R + (-q_b)*R = -(T(w,b) + q_b)*R = -(T(w,b) + T(j,b))*R

    But q_b = T(j,b). And T(w,b) + T(j,b) doesn't simplify further in general.

    Hmm. Let me try: for w at start, the combined factor per perm is:
      -(T(w,b) + T(j,b)) * [rest of path from b]

    For w in interior (...,a,w,b,...):
      -[T(a,w)*T(j,b) + T(a,i)*T(w,b)] * [left of a] * [right of b]

    For w at end:
      -(T(a,w) + T(a,i)) * [rest of path up to a]

    Note: T(a,i) = p_a = T(a,I) = T(a,0).

    OBSERVATION: T(a,w) + T(a,i) = T(a,w) + T(a,0).
    And T(w,b) + T(j,b) = T(w,b) + T(1,b).

    These are sums of tournament values from a common source a (or to a common target b).
    They don't simplify to anything universal.

    UNLESS we use the specific structure of tournaments more carefully.

    KEY: At s=0, T(w,i) + T(j,w) = 1, so T(j,w) = 1 - T(w,i).
    But T(j,b) is not directly related to T(w,b) unless w=i or b has special properties.
    """
    print(f"\n=== Combined formula at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            # Compute using the combined formula
            combined = 0.0

            for perm in permutations(W):
                pos_w = list(perm).index(w)

                if m == 1:
                    # Only vertex, boundary contributions: -1 -1 = -2
                    combined += -2.0
                    continue

                if pos_w == 0:
                    # w at start: next vertex is b = perm[1]
                    b = perm[1]
                    R_wt = 1.0
                    for k in range(1, m - 1):
                        R_wt *= T(perm[k], perm[k + 1])
                    factor = -(T(w, b) + T(J, b))
                    combined += factor * R_wt
                elif pos_w == m - 1:
                    # w at end: prev vertex is a = perm[m-2]
                    a = perm[m - 2]
                    L_wt = 1.0
                    for k in range(0, m - 2):
                        L_wt *= T(perm[k], perm[k + 1])
                    factor = -(T(a, w) + T(a, I))
                    combined += factor * L_wt
                else:
                    # w in interior: prev=a, next=b
                    a = perm[pos_w - 1]
                    b = perm[pos_w + 1]
                    L_wt = 1.0
                    for k in range(0, pos_w - 1):
                        L_wt *= T(perm[k], perm[k + 1])
                    R_wt = 1.0
                    for k in range(pos_w + 1, m - 1):
                        R_wt *= T(perm[k], perm[k + 1])
                    factor = -(T(a, w) * T(J, b) + T(a, I) * T(w, b))
                    combined += factor * L_wt * R_wt

            aH, _, _, _ = compute_alpha_H(T, W, w)
            print(f"  trial {trial}, w={w}: alpha_H={aH:.6f}, combined={combined:.6f}, "
                  f"err={abs(aH-combined):.2e}")
    print()


def test_factor_analysis(n, num_trials=10):
    """Analyze the factors in the combined formula.

    For w in interior (...,a,w,b,...):
      factor = -(T(a,w)*T(j,b) + T(a,i)*T(w,b))

    Note: T(a,w)*T(j,b) + T(a,i)*T(w,b)
      = T(a,w)*T(j,b) + (1-T(w,a))*(1-T(b,w))  ... no, T(a,i) != 1-T(w,a)

    At s=0: p_a = T(a,i), q_b = T(j,b), and p_a + q_a = 1 for all a in W.

    So: T(a,w)*q_b + p_a*T(w,b)

    Compare with the cycle formula's contribution at this position:
    In a (2k+1)-cycle through i,j with w at the boundary:
      derivative contribution = -p_{w_last} or -q_{w_first}

    The cycle derivative hits when w is at the START or END of the internal chain.
    The alpha_H formula involves w at ANY position.

    So the internal positions (w not at boundary of the internal chain) must somehow
    cancel or combine to give -2*H(B_w).

    HYPOTHESIS: The sum over permutations where w is at internal positions gives -2*H(B_w),
    and the sum where w is at boundary positions gives the cycle derivatives.

    Actually for boundary: w at position 0 gives -(T(w,b) + T(j,b))*(rest)
    and position m-1 gives -(T(a,w) + T(a,i))*(rest).

    Split T(w,b) + T(j,b) = T(w,b) + (1-T(b,j)):
      If b = some W vertex, T(j,b) is just another arc value.

    Let me try a DIFFERENT decomposition:

    factor at interior = -(T(a,w)*q_b + p_a*T(w,b))
      = -T(a,w)*q_b - p_a*T(w,b)
      = -q_b*T(a,w) - p_a*T(w,b)

    Compare with:
      -2*T(a,w)*T(w,b) = -T(a,w)*T(w,b) - T(a,w)*T(w,b)

    If we had q_b = T(w,b) and p_a = T(a,w), the factor would be -2*T(a,w)*T(w,b) =
    twice the full path weight through w. This would give -2*H(W) restricted somehow.

    But q_b != T(w,b) in general. The DIFFERENCE:
      factor + 2*T(a,w)*T(w,b) = -q_b*T(a,w) - p_a*T(w,b) + 2*T(a,w)*T(w,b)
        = T(a,w)*(T(w,b) - q_b) + T(w,b)*(T(a,w) - p_a)
        = T(a,w)*(T(w,b) - T(j,b)) + T(w,b)*(T(a,w) - T(a,i))

    This is a "correction" that depends on the difference between internal arcs and
    interface arcs. If a,b are not i or j, these differences are just arc differences.

    Let me try numerically to see if this correction gives the cycle derivatives.
    """
    print(f"\n=== Factor analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            # Compute "base" contribution: -2 * sum over perms with factor T(a,w)*T(w,b) at w
            base_contribution = 0.0
            correction = 0.0

            for perm in permutations(W):
                pos_w = list(perm).index(w)

                # Compute the path weight excluding the two arcs adjacent to w
                outer_wt = 1.0
                for k in range(m - 1):
                    if k == pos_w - 1 or k == pos_w:
                        continue
                    outer_wt *= T(perm[k], perm[k + 1])

                if pos_w == 0:
                    b = perm[1]
                    # factor = -(T(w,b) + T(j,b))
                    # base: -2*T(w,b)
                    # correction: -(T(j,b) - T(w,b)) = T(w,b) - T(j,b)
                    base_contribution += -2 * T(w, b) * outer_wt
                    correction += (T(w, b) - T(J, b)) * outer_wt
                elif pos_w == m - 1:
                    a = perm[m - 2]
                    # factor = -(T(a,w) + T(a,i))
                    # base: -2*T(a,w)
                    # correction: -(T(a,i) - T(a,w)) = T(a,w) - T(a,i)
                    base_contribution += -2 * T(a, w) * outer_wt
                    correction += (T(a, w) - T(a, I)) * outer_wt
                else:
                    a = perm[pos_w - 1]
                    b = perm[pos_w + 1]
                    # factor = -(T(a,w)*T(j,b) + T(a,i)*T(w,b))
                    # base: -2*T(a,w)*T(w,b)
                    # correction: factor - base
                    #   = -(T(a,w)*T(j,b) + T(a,i)*T(w,b)) + 2*T(a,w)*T(w,b)
                    #   = T(a,w)*(T(w,b) - T(j,b)) + T(w,b)*(T(a,w) - T(a,i))
                    base_contribution += -2 * T(a, w) * T(w, b) * outer_wt
                    corr = (T(a, w) * (T(w, b) - T(J, b)) +
                            T(w, b) * (T(a, w) - T(a, I))) * outer_wt
                    correction += corr

            aH, _, _, _ = compute_alpha_H(T, W, w)
            aI, a3, acyc = compute_alpha_I(T, W, w)
            B_w = [x for x in W if x != w]

            print(f"  trial {trial}, w={w}:")
            print(f"    alpha_H = {aH:.6f}")
            print(f"    base (-2*path-through-w sum) = {base_contribution:.6f}")
            print(f"    correction = {correction:.6f}")
            print(f"    base + correction = {base_contribution + correction:.6f}")
            print(f"    -2*H(B_w) = {-2*H_total(T, B_w):.6f}")
            print(f"    cycle derivs = {acyc:.6f}")
            print(f"    base vs -2*H(W) = {base_contribution:.6f} vs {-2*H_total(T, W):.6f}")

            # Is base = -2*H(W)?
            # No: base = -2 * sum_perm T(a,w)*T(w,b)*(rest) = -2*H(W)
            # Wait, the sum over all perms of the full path weight is H(W).
            # And T(a,w)*T(w,b)*(rest) IS the full path weight for that perm.
            # So base = -2*H(W)? Let me check.
            print(f"    -2*H(W) = {-2*H_total(T, W):.6f}")

            # For boundary cases (pos_w=0 or m-1), the "path weight" is
            # T(w,b)*outer_wt (with only one arc from w) or T(a,w)*outer_wt.
            # So base = -2 * sum_perm [full_path_weight] = -2*H(W).
            # The outer_wt * T(a,w) * T(w,b) or outer_wt * T(w,b) or outer_wt * T(a,w)
            # IS the full path weight. So base = -2*H(W). Let me verify.
            H_W = H_total(T, W)
            err_base = abs(base_contribution - (-2 * H_W))
            print(f"    err(base vs -2*H(W)) = {err_base:.2e}")
            print(f"    => correction = alpha_H + 2*H(W) = {aH + 2*H_W:.6f}")
            print(f"    => cycle derivs = {acyc:.6f}")
            print(f"    => correction vs cycle_derivs = {abs(aH + 2*H_W - acyc):.2e}")

            # ALSO: alpha_I = -2*H(B_w) + cycle_derivs
            # So: alpha_H - alpha_I = (-2*H(W) + correction) - (-2*H(B_w) + cycle_derivs)
            #   = -2*(H(W) - H(B_w)) + (correction - cycle_derivs)
            print(f"    H(W) - H(B_w) = {H_W - H_total(T, B_w):.6f}")
            print()
    print()


if __name__ == "__main__":
    # First verify the cut-based rewriting
    for n in [4, 5, 6]:
        test_alpha_by_cut(n, 5)

    # Then the combined formula
    for n in [4, 5, 6]:
        test_combined_formula(n, 5)

    # Then the factor analysis
    for n in [4, 5, 6]:
        test_factor_analysis(n, 5)
