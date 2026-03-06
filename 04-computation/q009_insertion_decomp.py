#!/usr/bin/env python3
"""
Decompose alpha_w^H by INSERTING w into permutations of B_w = W\{w}.

For each permutation pi = (u_1, ..., u_{m-1}) of B_w with weight B_wt,
inserting w at each position gives a contribution:

  alpha_w^H(pi) = -B_wt * S(pi)

where S(pi) = boundary terms + sum of intermediate terms involving
T(u_k,w), T(w,u_{k+1}), T(u_k,i), T(j,u_{k+1}), T(u_k,u_{k+1}).

If S(pi) = 2 for all pi, then alpha_w^H = -2*H(B_w) = alpha_w^I (at n=4).
For n>=5, S(pi) deviates from 2 and the deviation gives cycle derivatives.

Instance: opus-2026-03-05-S5
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


def compute_S_pi(T, pi, w, I, J):
    """Compute S(pi) for a permutation pi of B_w.

    S(pi) = T(w,u_1) + T(j,u_1)                    [pos 0]
          + T(u_{m-1},w) + T(u_{m-1},i)              [pos m-1]
          + sum_{k=1}^{m-2} [T(u_k,w)*T(j,u_{k+1}) + T(u_k,i)*T(w,u_{k+1})] / T(u_k,u_{k+1})
                                                       [interior positions]

    For m-1 = 1 (n=4, B_w has 1 vertex):
      pi = (u,). S = T(w,u) + T(j,u) + T(u,w) + T(u,i). No interior sum.
    """
    m_minus_1 = len(pi)
    if m_minus_1 == 0:
        return 2.0  # Only w in W, no pi vertices

    S = 0.0

    # Position 0: w inserted before u_1
    S += T(w, pi[0]) + T(J, pi[0])

    # Position m-1: w inserted after u_{m-1}
    S += T(pi[-1], w) + T(pi[-1], I)

    # Interior positions: w between u_k and u_{k+1} (k=0,...,m-3 in 0-indexed pi)
    for k in range(m_minus_1 - 1):
        u_k = pi[k]
        u_k1 = pi[k + 1]
        t_uk_uk1 = T(u_k, u_k1)
        if abs(t_uk_uk1) < 1e-15:
            # Edge case: if T(u_k, u_{k+1}) = 0, the B_wt is 0, so the contribution is 0
            # S(pi) is undefined but it doesn't matter
            return float('inf')  # Will be multiplied by 0
        numer = T(u_k, w) * T(J, u_k1) + T(u_k, I) * T(w, u_k1)
        S += numer / t_uk_uk1

    return S


def test_S_pi_values(n, num_trials=5):
    """Compute S(pi) for each permutation pi of B_w and check if it equals 2."""
    print(f"\n=== S(pi) analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            total_check = 0.0
            max_S_dev = 0.0
            S_values = []

            for pi in permutations(B_w):
                B_wt = 1.0
                for k in range(len(pi) - 1):
                    B_wt *= T(pi[k], pi[k + 1])

                S = compute_S_pi(T, pi, w, I, J)
                S_values.append((S, B_wt, pi))
                total_check += -B_wt * S
                if abs(B_wt) > 1e-15:
                    max_S_dev = max(max_S_dev, abs(S - 2))

            H_Bw = H_total(T, B_w)
            # alpha_w^H should equal total_check
            # alpha_w^I = -2*H(B_w) + cycle_derivs

            print(f"  trial {trial}, w={w}:")
            print(f"    sum(-B_wt*S) = {total_check:.6f}")
            print(f"    -2*H(B_w) = {-2*H_Bw:.6f}")
            print(f"    defect = sum(-B_wt*S) - (-2*H(B_w)) = {total_check + 2*H_Bw:.6f}")
            print(f"    max |S-2| = {max_S_dev:.6f}")
            if m <= 4:
                for S, B_wt, pi in sorted(S_values, key=lambda x: -abs(x[1])):
                    if abs(B_wt) > 1e-10:
                        print(f"      pi={pi}, B_wt={B_wt:.4f}, S={S:.4f}, S-2={S-2:.4f}, "
                              f"-B_wt*(S-2)={-B_wt*(S-2):.4f}")
    print()


def test_defect_vs_cycles(n, num_trials=5):
    """Check if the defect sum(-B_wt*(S-2)) equals the cycle derivative contribution."""
    print(f"\n=== Defect vs cycles at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(min(num_trials, 5)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)
        p = {w: T(w, I) for w in W}
        q = {w: T(J, w) for w in W}

        for w in W[:2]:
            B_w = [x for x in W if x != w]
            H_Bw = H_total(T, B_w)

            # Compute defect
            defect = 0.0
            for pi in permutations(B_w):
                B_wt = 1.0
                for k in range(len(pi) - 1):
                    B_wt *= T(pi[k], pi[k + 1])
                S = compute_S_pi(T, pi, w, I, J)
                defect += -B_wt * (S - 2)

            # Compute cycle derivatives
            cycle_derivs = 0.0
            for chain_len in range(3, m + 1, 2):
                for combo in combinations(W, chain_len):
                    for perm in permutations(combo):
                        internal = 1.0
                        for k in range(chain_len - 1):
                            internal *= T(perm[k], perm[k + 1])
                        if w == perm[0]:
                            cycle_derivs += 2 * (-p[perm[-1]]) * internal
                        elif w == perm[-1]:
                            cycle_derivs += 2 * (-q[perm[0]]) * internal

            err = abs(defect - cycle_derivs)
            print(f"  trial {trial}, w={w}: defect={defect:.6f}, cycles={cycle_derivs:.6f}, "
                  f"err={err:.2e}")
    print()


def analyze_S_interior(n, num_trials=3):
    """Analyze the interior terms of S(pi) more carefully.

    Interior term for position k (w between u_k and u_{k+1}):
      [T(u_k,w)*T(j,u_{k+1}) + T(u_k,i)*T(w,u_{k+1})] / T(u_k,u_{k+1})

    At s=0: T(u_k,i) = p_{u_k} and T(j,u_{k+1}) = q_{u_{k+1}}.

    So the term is: [T(u_k,w)*q_{u_{k+1}} + p_{u_k}*T(w,u_{k+1})] / T(u_k,u_{k+1})

    This is a "weighted average" of q and p terms.

    Can we split this as 1 + correction?

    If T(u_k,w)*q_{u_{k+1}} + p_{u_k}*T(w,u_{k+1}) = T(u_k,u_{k+1}) + something,
    then the interior term = 1 + something/T(u_k,u_{k+1}).

    T(u_k,w)*q_{u_{k+1}} + p_{u_k}*T(w,u_{k+1})
    vs T(u_k,u_{k+1}):

    No obvious relation since these involve different arcs.

    But what if we decompose using T(u_k,w) = T(u_k,w), T(w,u_{k+1}) = T(w,u_{k+1})?
    These are internal arcs, not directly related to interface arcs.

    KEY: At s=0, for any vertex v in W:
      p_v = T(v,i) = T(v,0)
      q_v = T(j,v) = T(1,v)
      p_v + q_v = 1

    So q_v = 1 - p_v = 1 - T(v,0) = T(0,v) = T(i,v).

    Therefore: q_{u_{k+1}} = T(i, u_{k+1}).

    Interior term = [T(u_k,w)*T(i,u_{k+1}) + T(u_k,i)*T(w,u_{k+1})] / T(u_k,u_{k+1})

    This looks like a "3-path through i": u_k → w → i → u_{k+1} gives T(u_k,w)*T(w,i)*T(i,u_{k+1})?
    No, that's different.

    Actually: T(u_k,w)*T(i,u_{k+1}) + T(u_k,i)*T(w,u_{k+1})

    This is the sum of two "detour" paths that go through either (u_k→w, then i→u_{k+1})
    or (u_k→i, then w→u_{k+1}).

    This looks like a permanent/determinant of a 2x2 matrix!

    [T(u_k,w)  T(u_k,i)] [T(w,u_{k+1})  ] + cross term
    But it's T(u_k,w)*T(i,u_{k+1}) + T(u_k,i)*T(w,u_{k+1})
    = det or perm of [[T(u_k,w), T(u_k,i)], [T(w,u_{k+1}), T(i,u_{k+1})]]???

    No. The permanent of [[a,b],[c,d]] = ac + bd.
    The determinant of [[a,b],[c,d]] = ad - bc.

    What we have: T(u_k,w)*T(i,u_{k+1}) + T(u_k,i)*T(w,u_{k+1})
    = perm([[T(u_k,w), T(u_k,i)], [T(i,u_{k+1}), T(w,u_{k+1})]])

    YES! This is the permanent of the 2x2 matrix M where:
    M[0][0] = T(u_k, w),   M[0][1] = T(u_k, i)
    M[1][0] = T(i, u_{k+1}), M[1][1] = T(w, u_{k+1})

    Rows indexed by source {u_k}, columns by target {u_{k+1}}.
    The "routing" vertices are {w, i}.

    perm(M) = sum over matchings of {w,i} to {sources, targets}:
    - w routes from u_k (T(u_k,w)), i routes to u_{k+1} (T(i,u_{k+1}))
    - i routes from u_k (T(u_k,i)), w routes to u_{k+1} (T(w,u_{k+1}))

    This is actually the total weight of all Hamiltonian paths from u_k to u_{k+1}
    through the set {w, i}!

    H({w,i}, start=u_k, end=u_{k+1})?
    Well, paths through {w,i} from u_k: u_k→w→i→u_{k+1} or u_k→i→w→u_{k+1}.
    Weight: T(u_k,w)*T(w,i)*T(i,u_{k+1}) + T(u_k,i)*T(i,w)*T(w,u_{k+1}).

    That's NOT the same as the permanent. The permanent has T(u_k,w)*T(i,u_{k+1})
    without T(w,i). So it's a "two-step" path, not a three-step path.

    Hmm. Let me just check numerically if this permanent structure helps.
    """
    print(f"\n=== Interior term analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(num_trials):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            B_w = [x for x in W if x != w]

            # For each pair (u_k, u_{k+1}) in B_w, compute the interior term
            for u_k in B_w:
                for u_k1 in B_w:
                    if u_k == u_k1:
                        continue
                    t = T(u_k, u_k1)
                    numer = T(u_k, w) * T(I, u_k1) + T(u_k, I) * T(w, u_k1)
                    if abs(t) > 1e-10:
                        ratio = numer / t
                    else:
                        ratio = float('inf')

                    # What about the boundary-like version?
                    # For u_k → w at end: T(u_k, w) + T(u_k, I)
                    # For w → u_k1 at start: T(w, u_k1) + T(J, u_k1) = T(w, u_k1) + T(I, u_k1)
                    # Wait, q = T(J,v) = T(1,v). And at s=0, q_v = T(I,v) since q_v = 1-p_v = T(I,v).
                    # Actually no: q_v = T(J,v) = T(1,v). And T(I,v) = T(0,v) = 1-T(v,0) = 1-p_v = q_v.
                    # Yes! q_v = T(I,v) since T(I,v) = 1 - T(v,I) = 1 - p_v = q_v.

                    # So the interior term is [T(u_k,w)*T(I,u_{k+1}) + T(u_k,I)*T(w,u_{k+1})] / T(u_k,u_{k+1})

                    if trial == 0:
                        print(f"    u_k={u_k}, u_{u_k1}={u_k1}: t={t:.4f}, numer={numer:.4f}, "
                              f"ratio={ratio:.4f}")
    print()


def test_telescoping(n, num_trials=5):
    """Test if S(pi) - 2 can be written as a telescoping sum.

    S(pi) = [T(w,u_1) + T(j,u_1)]
          + sum_{k=0}^{m-3} [T(u_k,w)*T(i,u_{k+1}) + T(u_k,i)*T(w,u_{k+1})] / T(u_k,u_{k+1})
          + [T(u_{m-1},w) + T(u_{m-1},i)]

    Note: at s=0, T(j,v) = T(i,v) = q_v for all v in W.
    Wait: q_v = T(J,v) = T(1,v). And T(I,v) = T(0,v) = 1 - T(v,0) = 1 - p_v = q_v.
    So T(J,v) = T(I,v) = q_v? NO!

    T(J,v) = T(1,v), T(I,v) = T(0,v). These are generally different!

    At s=0: p_v + q_v = 1, i.e., T(v,0) + T(1,v) = 1, i.e., T(v,I) + T(J,v) = 1.
    So T(v,I) = 1 - T(J,v). This means T(I,v) = 1 - T(v,I) = 1 - (1-T(J,v)) = T(J,v).
    So T(I,v) = T(J,v) for all v in W at s=0!

    WAIT. That's a KEY insight. At s=0:
      T(v,I) + T(J,v) = 1  =>  T(I,v) = 1 - T(v,I) = T(J,v).

    So: T(I,v) = T(J,v) = q_v for all v in W.
    And: T(v,I) = T(v,J) = p_v? No: T(v,J) = 1 - T(J,v) = 1 - q_v = p_v. Yes!

    So at s=0: T(I,v) = T(J,v) and T(v,I) = T(v,J) for all v in W.

    This means in the S(pi) formula:
    - T(j,u_1) = T(J,u_1) = T(I,u_1)
    - T(u_{m-1},i) = T(u_{m-1},I) = T(u_{m-1},J)
    - The interior numerator: T(u_k,w)*T(I,u_{k+1}) + T(u_k,I)*T(w,u_{k+1})

    But also: T(u_k,I) = T(u_k,J) = p_{u_k}, and T(I,u_{k+1}) = T(J,u_{k+1}) = q_{u_{k+1}}.

    So i and j are INTERCHANGEABLE in these formulas (at s=0)!

    This means S(pi) is SYMMETRIC under swapping i and j (at s=0).

    Hmm, but that doesn't immediately help with proving S(pi) = 2 or computing the defect.

    Let me try another angle. Define for each vertex v:
      f(v) = T(v,w) / T(v,I)  (ratio of how much v "prefers" w over i)

    Then T(v,w) = f(v) * T(v,I) = f(v) * p_v.

    Interior numerator = p_{u_k} * f(u_k) * q_{u_{k+1}} + p_{u_k} * T(w, u_{k+1})
    = p_{u_k} * [f(u_k) * q_{u_{k+1}} + T(w, u_{k+1})]

    Hmm, doesn't simplify.

    Let me try: what if we use T(u_k,w) + T(u_k,i) = T(u_k,w) + p_{u_k}?

    At the boundary:
      pos 0: T(w,u_1) + q_{u_1}
      pos m-1: T(u_{m-1},w) + p_{u_{m-1}}

    Sum of boundary: T(w,u_1) + q_{u_1} + T(u_{m-1},w) + p_{u_{m-1}}

    For m-1=1: this is T(w,u) + q_u + T(u,w) + p_u = 1 + 1 = 2. So S=2. ✓

    For m-1=2: pi = (u,v). S = T(w,u)+q_u + [T(u,w)*q_v + p_u*T(w,v)]/T(u,v) + T(v,w)+p_v
    = 1 + (q_u + p_v) + [T(u,w)*q_v + p_u*T(w,v)]/T(u,v)
    Hmm wait: T(w,u)+T(u,w) = 1 always. So boundary = 1 + (q_u + p_v).

    At s=0: q_u + p_v = q_u + p_v. No simplification unless u=v.

    For m-1=2: S = 1 + q_u + p_v + [T(u,w)*q_v + p_u*T(w,v)]/T(u,v)

    We need S = 2 + defect. So defect = q_u + p_v - 1 + [T(u,w)*q_v + p_u*T(w,v)]/T(u,v)
    = (q_u - p_u) + (p_v + p_u - 1) + ...

    This is getting messy. Let me just verify the key equality numerically for many trials.
    """
    print(f"\n=== S(pi) = 2 test at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        # Verify T(I,v) = T(J,v) at s=0
        for v in W:
            assert abs(T(I, v) - T(J, v)) < 1e-10, f"T(I,{v})={T(I,v)} != T(J,{v})={T(J,v)}"

        for w in W[:1]:
            B_w = [x for x in W if x != w]
            defect_total = 0.0
            for pi in permutations(B_w):
                B_wt = 1.0
                for k in range(len(pi) - 1):
                    B_wt *= T(pi[k], pi[k + 1])
                S = compute_S_pi(T, pi, w, I, J)
                defect_total += -B_wt * (S - 2)

            # Compare with cycle derivatives
            p = {v: T(v, I) for v in W}
            q = {v: T(J, v) for v in W}
            cycle_derivs = 0.0
            for chain_len in range(3, m + 1, 2):
                for combo in combinations(W, chain_len):
                    for perm in permutations(combo):
                        internal = 1.0
                        for k in range(chain_len - 1):
                            internal *= T(perm[k], perm[k + 1])
                        if w == perm[0]:
                            cycle_derivs += 2 * (-p[perm[-1]]) * internal
                        elif w == perm[-1]:
                            cycle_derivs += 2 * (-q[perm[0]]) * internal

            print(f"  trial {trial}, w={w}: defect_sum={defect_total:.6f}, "
                  f"cycle_derivs={cycle_derivs:.6f}, err={abs(defect_total - cycle_derivs):.2e}")
    print()


def test_key_relation(n, num_trials=10):
    """At s=0, T(I,v) = T(J,v) = q_v for all v in W.

    This means the interior term simplifies:
    [T(u_k,w)*q_{u_{k+1}} + p_{u_k}*T(w,u_{k+1})] / T(u_k,u_{k+1})

    With p+q=1, this is:
    [T(u_k,w)*(1-p_{u_{k+1}}) + (1-q_{u_k})*T(w,u_{k+1})] / T(u_k,u_{k+1})
    = [T(u_k,w) - T(u_k,w)*p_{u_{k+1}} + T(w,u_{k+1}) - q_{u_k}*T(w,u_{k+1})] / T(u_k,u_{k+1})
    = [1 - T(u_k,w)*p_{u_{k+1}} - q_{u_k}*T(w,u_{k+1})] / T(u_k,u_{k+1})

    since T(u_k,w) + T(w,u_{k+1}) = T(u_k,w) + 1 - T(u_{k+1},w).

    Actually that's only true if u_k and u_{k+1} are the same vertex, which they're not.
    T(u_k,w) + T(w,u_{k+1}) != 1 in general.

    ALTERNATIVE: use the bracket identity.
    bracket(a,b) = p_a*q_b - (1-q_a)*(1-p_b) = p_a*q_b - p_a_bar*q_b_bar
    where p_a_bar = 1-q_a = T(a,J) = T(a,I) at s=0, and q_b_bar = 1-p_b.

    Wait: p_a = T(a,I), q_a = T(J,a) = T(I,a).
    1-q_a = T(a,I) = p_a. And 1-p_b = T(I,b) = q_b.

    So bracket(a,b) = p_a*q_b - p_a*q_b = 0 at s=0! Of course, that's the whole point.

    The DERIVATIVE d(bracket)/ds_w = -(1+d_w)/2 or -(1-d_w)/2 depending on which variable.

    OK let me try yet another approach. The key equation to prove is:

    sum_pi B_wt * (S(pi) - 2) = -cycle_derivs

    The LHS sums over permutations of B_w. The RHS involves permutations of subsets of W
    containing w. Can we match them?

    The cycle derivative for a chain (w1,...,w_k) with w at endpoint:
    - If w = w1: contribution = 2*(-p_{w_k})*T(w1,w2)*...*T(w_{k-1},w_k)
    - If w = w_k: contribution = 2*(-q_{w_1})*T(w1,w2)*...*T(w_{k-1},w_k)

    The key insight: each such chain uses k-1 vertices from B_w (the ones that aren't w).
    The remaining m-1-(k-1) = m-k vertices of B_w form the "rest" of the permutation.

    HYPOTHESIS: The defect sum over pi can be decomposed by which vertices of B_w
    participate in the "cycle chain" and which are just passing through.

    Let me verify: for n=5 (m=3), B_w has 2 vertices. The 5-cycle chain uses 2 more W-vertices
    beyond w, so all of B_w. So the chain IS a permutation of B_w.
    """
    print(f"\n=== Key relation test at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    # Just verify T(I,v) = T(J,v)
    arc_values = random_s0_tournament(n, seed=42)
    T = make_T(n, arc_values)
    for v in W:
        print(f"  T(I,{v}) = {T(I,v):.4f}, T(J,{v}) = {T(J,v):.4f}, "
              f"diff = {T(I,v)-T(J,v):.2e}")
    print()


if __name__ == "__main__":
    # Test S(pi) values
    for n in [4, 5, 6]:
        test_S_pi_values(n, 3)

    # Test defect = cycle derivatives
    for n in [4, 5, 6, 7]:
        test_defect_vs_cycles(n, 5)

    # Verify T(I,v)=T(J,v)
    for n in [4, 5]:
        test_key_relation(n, 1)

    # Telescoping analysis
    for n in [4, 5]:
        test_telescoping(n, 3)
