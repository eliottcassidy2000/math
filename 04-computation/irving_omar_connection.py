#!/usr/bin/env python3
"""
Irving-Omar connection to transfer matrix M[a,b].

KEY FORMULA (Irving-Omar Proposition 2):
  ham(D) = sum_S det(Ā[S]) · per(A[S^c])

Our transfer matrix:
  M[a,b] = sum_S (-1)^|S| E_a(S+a) · B_b(R+b)

QUESTION: Can we express M[a,b] as a det×per-type formula that
makes the even-r property transparent?

TOURNAMENT IDENTITY (Corollary 19):
  U_D = sum over σ with all odd cycles of 2^{ψ(σ)} p^{cyc}(σ)

Even-length cycles DON'T contribute to U_D for tournaments!
This is the symmetric-function version of even r-powers.

CONNECTION TO TEST:
1. Express M[a,b] in terms of cycle covers of a modified digraph
2. Show that even-cycle contributions cancel via the same mechanism
3. The "all odd cycles" constraint forces even r-powers

Instance: opus-2026-03-06-S24
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Symbol, Matrix, det, Poly
from collections import defaultdict
import sys

def make_symbols(n):
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    return r, s

def edge_weight(i, j, r, s):
    if i == j:
        return 0  # no self-loops
    return r + s[(i,j)]

def ham_paths_ending_at(vertex_set, target, r, s):
    vs = list(vertex_set)
    if len(vs) == 1: return 1
    total = 0
    for perm in permutations([v for v in vs if v != target]):
        path = list(perm) + [target]
        w = 1
        for k in range(len(path)-1): w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def ham_paths_beginning_at(vertex_set, source, r, s):
    vs = list(vertex_set)
    if len(vs) == 1: return 1
    total = 0
    for perm in permutations([v for v in vs if v != source]):
        path = [source] + list(perm)
        w = 1
        for k in range(len(path)-1): w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def compute_M(a, b, n, r, s):
    V = list(range(n))
    U = [v for v in V if v != a and v != b]
    total = 0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])
        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)
        total += sign * Ea * Bb
    return expand(total)


def build_weight_matrix(vertices, r, s):
    """Build the weight matrix A[i,j] = t_{ij} = r + s_{ij}."""
    vlist = sorted(vertices)
    n = len(vlist)
    A = Matrix.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = edge_weight(vlist[i], vlist[j], r, s)
    return A, vlist


def permanent(M):
    """Compute permanent of matrix M."""
    n = M.shape[0]
    total = 0
    for perm in permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= M[i, perm[i]]
        total += prod
    return expand(total)


def test_det_per_decomposition(n):
    """
    Test if M[a,b] can be expressed via det×per of submatrices.

    Irving-Omar: ham(D) = sum_S det(Ā[S]) * per(A[S^c])

    For c-tournaments: Ā_{ij} = t_{ji} = r + s_{ji} = r - s_{ij}
    So Ā = A^T (as a matrix).

    Can we adapt this to M[a,b]?
    """
    print(f"\n{'='*60}")
    print(f"DET × PER DECOMPOSITION: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    a, b = 0, 1
    M_ab = compute_M(a, b, n, r, s)

    # Build A on all vertices
    A, vlist = build_weight_matrix(range(n), r, s)
    AT = A.T  # = Ā for tournaments (complement = transpose)

    # Irving-Omar: ham = sum_S det(A^T[S]) * per(A[S^c])
    # This sums over ALL subsets S ⊆ [n].
    # For our M[a,b], we need endpoint-conditioned paths.

    # IDEA: Modify the matrix to condition on endpoints.
    # Define A' where we force paths to end at specific vertices.
    #
    # One approach: use a BORDERED matrix.
    # Add extra rows/columns for the endpoints a, b.
    # Or: use the "source-sink" formulation of LGV.

    # Actually, E_a(S+a) = sum of ham paths ending at a through S+a
    # = per of a certain submatrix conditioned on last vertex = a.
    # = sum over permutations σ of S where σ starts from S and last step → a.
    # This is: sum_{last step u→a} per(A[S, S\{u}]) * t(u,a)
    # Or more directly: E_a(W) = permanent of the "path matrix" for W ending at a.

    # KEY REALIZATION: E_a(W) is NOT a permanent of a square submatrix of A.
    # It's a permanent of a RECTANGULAR arrangement:
    # rows = positions 1,...,|W|-1, cols = vertices of W\{a},
    # with the last entry fixed as "→a".
    # This is actually the (a,*)-row-conditioned permanent.

    # Let me try a different angle: think of M as endpoint-conditioned
    # version of the Irving-Omar formula.

    # The I-O formula computes ham(D) = total Hamiltonian paths.
    # We want M[a,b] = "signed Hamiltonian 2-path-cover count with endpoints a,b".

    # APPROACH: Define the "endpoint matrix"
    # Phi[a,b] = sum over Ham paths from * to a, * to b through V
    #          = some function of A
    # Then M[a,b] might be a cofactor of Phi.

    # Actually, let's directly test:
    # M[a,b] = cofactor of the PERMANENT matrix?
    # Define P where P[i,j] = permanent of A with row i, col j deleted.
    # Then per(A) = sum_j A[i,j] * P[i,j] (Laplace expansion of permanent).
    # M[a,b] ≠ P[a,b] (already tested). But maybe a signed version?

    # Let me instead look at the CYCLE COVER interpretation.
    # A permutation σ on V corresponds to a set of directed cycles covering V.
    # The weight of σ is prod_{i} A[i,σ(i)] = prod_i t(i,σ(i)).

    # per(A) = sum_σ prod_i t(i,σ(i)) = sum over all cycle covers.
    # det(A) = sum_σ sgn(σ) prod_i t(i,σ(i)).

    # For tournaments: A^T_{ij} = t(j,i) = 2r - t(i,j) at c=2r.
    # So per(A^T) = per(A) and det(A^T) = det(A).

    # The Irving-Omar Corollary 19 says: only ODD-length cycles contribute
    # to U_D for tournaments. This means: for each permutation σ with an
    # even-length cycle, the contribution cancels.

    # WHY? Because for a k-cycle C = (i_1, i_2, ..., i_k) in σ:
    # Reversing C gives C^rev = (i_1, i_k, i_{k-1}, ..., i_2).
    # The weight changes: prod t(i_j, i_{j+1}) → prod t(i_{j+1}, i_j).
    # For c-tournaments: t(j,i) = c - t(i,j), so each factor changes.
    # For a k-cycle: prod t(i_{j+1}, i_j) = prod (c - t(i_j, i_{j+1})).

    # Now: the SIGN of the permutation changes by (-1)^{k-1} when we reverse
    # a k-cycle (reversal = conjugation by flip, changes parity by k-1).

    # For EVEN k: (-1)^{k-1} = -1, so the sign flips.
    # The weight: prod (c - t) vs prod t. In the det, the sign flip means
    # the two contributions have OPPOSITE signs but DIFFERENT weights.
    # They don't cancel directly... unless we expand in powers of c.

    # AH! This is the key: in the p-basis expansion of U_D,
    # the cycle-reversal involution PAIRS permutations with even cycles.
    # For each pair, the p-monomial is the same (same cycle type),
    # but the sign changes. So they cancel in the p-expansion.

    # This is EXACTLY the "even cycle vanishing theorem" (THM-148).

    # Now: how does this connect to M[a,b]?
    # M[a,b] is not a cycle cover — it's a 2-PATH-cover.
    # A 2-path-cover has no cycles — it's two vertex-disjoint paths.

    # But the INCLUSION-EXCLUSION in M converts path products into
    # something that CAN be expressed via cycle covers!

    # Specifically: E_a(S+a) * B_b(R+b) involves two separate permanents.
    # The (-1)^|S| sign creates a signed sum.
    # The TOTAL sum over all S is like a "broken" determinant.

    # DEEPER INSIGHT: M[a,b] is the (a,b) entry of a matrix that
    # can be expressed as:
    # M = E^T · diag((-1)^|S|) · B
    # where E[S,v] = E_v(S+v), B[S,v] = B_v(R+v).

    # For even r-powers, we need E^T Λ B = (E^T Λ B)^T,
    # i.e., E^T Λ B = B^T Λ E.

    # Under path reversal: E_v(W; r, s) at -r gives B_v(W; r, -s)
    # times (-1)^|W|-1.

    # So E[S,v; -r] = (-1)^|S| B[S,v; r, -s], meaning the E and B
    # matrices are related by r → -r and s → -s.

    # If M(r) = E^T(r) Λ B(r), then
    # M(-r) = E^T(-r) Λ B(-r) = ((-1)^|·| B(r,-s))^T Λ ((-1)^|·| E(r,-s))
    # The signs from (-1)^|S| get squared (= 1), so:
    # M(-r) = B^T(r,-s) Λ E(r,-s) = M^T(r,-s)
    # = M(r,-s) if M is symmetric.

    # So M(-r) = M(r,-s), which is the tautology again.
    # We need M(-r) = M(r), which requires M(r,-s) = M(r,s) for even n
    # and M(r,-s) = -M(r,s) for odd n. These are the s-parity claims.

    # NEW ANGLE: Use Irving-Omar's Theorem 10(b):
    # U_D = L_n exp( sum_k (1/k) p_k [tr(XA)^k + (-1)^{k-1} tr(XĀ)^k] )
    # For tournaments: tr(XĀ)^k = tr(XA)^k for k > 1.
    # So: tr(XA)^k + (-1)^{k-1} tr(XA)^k = (1 + (-1)^{k-1}) tr(XA)^k
    # = 0 for k even, 2·tr(XA)^k for k odd!

    # This is WHY only odd cycles contribute!

    # FOR M[a,b]: We need the endpoint-conditioned version.
    # M[a,b] extracts the coefficient of x_a * x_b^{n-1} (or similar)
    # from U_D. Wait, no: M[a,b] is not directly a coefficient of U_D.

    # Actually, M[a,b] relates to [s_{(1^n)}] conditioned on endpoints.
    # The Schur function s_{(1^n)} extracts the Hamiltonian path count.
    # Conditioning on endpoints a and b...

    # Let me just compute and compare.

    # Test: compute per(A) and det(A) and relate to M
    per_A = permanent(A)
    det_A = expand(det(A))

    print(f"\n  per(A) = {per_A}")
    print(f"  det(A) = {det_A}")

    # Total Hamiltonian path count H(T) = sum_{a,b, a≠b} (something involving M)
    # Actually at r=1/2 (c=1, tournament), M[a,b] at binary weights gives the
    # signed path cover. sum_a M[a,a] = H(T) for odd n.

    # Let me compute the actual cycle-cover expansion.
    # For each permutation σ on V, compute:
    #   weight(σ) = prod_i t(i, σ(i))
    #   sign(σ) = (-1)^{n - #cycles(σ)}
    # and check which σ contribute to M[a,b] via some formula.

    print(f"\n  Cycle cover analysis:")

    V = list(range(n))
    total_det = 0
    total_per = 0
    cycle_type_weights = defaultdict(lambda: 0)

    for perm in permutations(V):
        weight = 1
        for i in range(n):
            weight *= edge_weight(i, perm[i], r, s)
        weight = expand(weight)

        # Compute cycle type
        visited = [False] * n
        cycles = []
        for i in range(n):
            if not visited[i]:
                cycle = []
                j = i
                while not visited[j]:
                    visited[j] = True
                    cycle.append(j)
                    j = perm[j]
                cycles.append(tuple(cycle))

        cycle_lengths = tuple(sorted([len(c) for c in cycles], reverse=True))
        sign = (-1)**(n - len(cycles))

        cycle_type_weights[cycle_lengths] = expand(
            cycle_type_weights[cycle_lengths] + weight)

        total_det = expand(total_det + sign * weight)
        total_per = expand(total_per + weight)

    print(f"  Cycle type decomposition of per(A):")
    for ct in sorted(cycle_type_weights.keys()):
        w = cycle_type_weights[ct]
        has_even = any(l % 2 == 0 for l in ct)
        label = " [HAS EVEN]" if has_even else " [all odd]"
        print(f"    {ct}: {w}{label}")

    # Check: for cycle types with even parts, what is the weight?
    even_total = 0
    odd_total = 0
    for ct, w in cycle_type_weights.items():
        if any(l % 2 == 0 for l in ct):
            even_total = expand(even_total + w)
        else:
            odd_total = expand(odd_total + w)

    print(f"\n  Total weight (all odd cycles): {odd_total}")
    print(f"  Total weight (has even cycle): {even_total}")
    print(f"  per(A) = {expand(odd_total + even_total)}")
    print(f"  per(A) matches? {expand(per_A - odd_total - even_total) == 0}")

    # KEY TEST: Is even_total always odd-in-r?
    # For tournaments: the even-cycle contributions should have opposite
    # r-parity from the odd-cycle contributions.
    p_even = Poly(even_total, r) if even_total != 0 else None
    p_odd = Poly(odd_total, r) if odd_total != 0 else None

    if p_even:
        even_coeffs = {m[0]: c for m, c in p_even.as_dict().items()}
        print(f"\n  Even-cycle r-powers present: {sorted(even_coeffs.keys())}")
    if p_odd:
        odd_coeffs = {m[0]: c for m, c in p_odd.as_dict().items()}
        print(f"  Odd-cycle r-powers present: {sorted(odd_coeffs.keys())}")

    # Now: what about M[a,b]?
    print(f"\n  M[{a},{b}] = {M_ab}")

    # Is M[a,b] related to the "odd-cycle-only" permanent?
    # Or to a specific Schur coefficient of U_D?

    return M_ab


def test_endpoint_conditioned_cycles(n):
    """
    Express M[a,b] in terms of cycle covers of a MODIFIED digraph.

    Idea: Add edges (b→a) and (a→b) to create cycles from path pairs.
    A 2-path-cover (P_a ending at a, P_b starting at b) can be completed
    to a permutation by adding the edge b→a (connecting end of P_b to start
    of P_a) and a→b (connecting end of P_a to start of P_b).
    Wait, that's not right. Let me think...

    P_a: v_1 → v_2 → ... → a  (Ham path through S+a ending at a)
    P_b: b → w_1 → w_2 → ... → w_k  (Ham path through R+b starting at b)

    To make a permutation: connect end of P_b (w_k) to start of P_a (v_1),
    and end of P_a (a) to start of P_b (b).

    This creates a SINGLE CYCLE: v_1 → v_2 → ... → a → b → w_1 → ... → w_k → v_1

    So EACH 2-path-cover corresponds to a SINGLE HAMILTONIAN CYCLE
    (with the extra edges a→b and w_k→v_1 added).

    Wait, but M[a,b] sums over DIFFERENT partitions (S, R).
    And the inclusion-exclusion sign (-1)^|S| is crucial.

    Hmm, the single-cycle interpretation doesn't directly work because
    the partition (S, R) is part of the structure.

    But maybe: the sum over all S of (-1)^|S| * (single Ham cycle through
    fixed edge a→b) gives M[a,b].

    NO — the 2-path-cover has TWO separate paths, not one cycle with a→b.
    The vertex partition S∪{a}, R∪{b} is essential.

    Let me think again: M[a,b] involves E_a (paths ending at a) and B_b
    (paths starting at b). If we add the "virtual edge" a→b, each 2-path-
    cover becomes a single Hamiltonian cycle containing the edge a→b.

    Specifically: compose P_a (ending at a) + edge(a→b) + P_b (starting at b)
    + edge(w_k → v_1) to get a Ham cycle through edge a→b.

    But the last edge w_k→v_1 is NOT part of the original tournament.
    It's an artifact of cycle completion.

    ALTERNATIVE: M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b).
    Fix the edge a→b and consider a permutation σ where σ(a) = b.
    Then σ restricted to S+a gives the E_a path structure,
    and σ restricted to R+b gives the B_b path structure.

    Actually: if σ(a) = b, then in the cycle decomposition of σ,
    a and b are in the same cycle, with a→b as an edge.
    The restriction of σ to V\{a→b edge} splits into:
    - a path from ?(→...→a) in part of V
    - a path from b(→...→?) in the rest of V

    This IS a 2-path-cover! The vertices in the "a-side" are those
    that appear before a in the cycle, and the "b-side" are those after b.

    So: M[a,b] relates to cycle covers where a→b is a fixed edge!
    """
    print(f"\n{'='*60}")
    print(f"ENDPOINT-CONDITIONED CYCLE ANALYSIS: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    a, b = 0, 1
    V = list(range(n))

    M_ab = compute_M(a, b, n, r, s)

    # Consider permutations σ with σ(a) = b.
    # Such σ has a cycle containing a→b.
    # The cycle is: a → b → σ(b) → σ²(b) → ... → a.

    # For each such σ, the weight is prod_i t(i, σ(i)).
    # But we want the SIGNED weight from the inclusion-exclusion.

    # The key: the (-1)^|S| sign in M relates to... what?
    # |S| = number of vertices in the a-side path (excluding a).
    # In the cycle, |S| = number of vertices between b's predecessor
    # and a (going backwards from a).
    # Actually: the cycle is ..., u, a, b, w, ...
    # So in the 2-path-cover: S+a = {the vertices from cycle start to a},
    # R+b = {b, the vertices after b in cycle, ..., back to start}.

    # For a FULL cycle of length n (Hamiltonian cycle through a→b):
    # a → b → w_1 → w_2 → ... → w_{n-2} → a
    # Here S = {w_1, ..., w_{n-2}} and R = ∅ (or vice versa).
    # Wait no: S+a contains a and all predecessors of a,
    # R+b contains b and all successors of b.

    # In the cycle a → b → w_1 → ... → w_{n-2} → a:
    # Predecessors of a = {w_{n-2}} (the vertex before a in the cycle).
    # But in our formula, S+a is the vertex set of the E_a path.
    # E_a path ends at a. If the Ham cycle is a→b→w_1→...→w_{n-2}→a,
    # then removing edge a→b gives two paths:
    #   P_a: w_{n-2} → ... → w_1 → b... wait that doesn't end at a.

    # Let me reconsider. Remove edge a→b from the cycle.
    # We get a Hamiltonian PATH from b to a:
    # b → w_1 → w_2 → ... → w_{n-2} → a.
    # This is a single path, not two paths!

    # So a Ham cycle through a→b corresponds to a single Ham PATH from b to a
    # (after removing the a→b edge). This counts TOTAL Ham paths from b to a,
    # not the 2-path-cover that M computes.

    # So M[a,b] ≠ Ham cycle count through a→b.
    # M[a,b] is MORE structured: it's the inclusion-exclusion sum.

    # OK let me try yet another connection.
    # Consider permutations σ on V. Each σ is a cycle cover.
    # For σ with σ(a) = b (a fixed), this constrains the cycle structure.

    # The "signed permanent with σ(a)=b fixed" is:
    # sum_{σ: σ(a)=b} sgn(σ) * prod_i t(i, σ(i))
    # = (-1)^{a+b} * cofactor(a,b) of A.

    # We already know M ≠ cofactor (degree mismatch).

    # BUT: what if we restrict to permutations where a→b is in a cycle
    # of SPECIFIC length? The length-2 cycle (a→b→a) gives σ(a)=b, σ(b)=a.
    # And the rest of σ acts on U = V\{a,b}.

    # For length-2 cycle: σ(a)=b, σ(b)=a, and σ|_U is an arbitrary perm of U.
    # Weight = t(a,b) * t(b,a) * prod_{u∈U} t(u, σ(u))
    #        = (r+s_{ab})(r-s_{ab}) * prod = (r²-s_{ab}²) * prod.

    # For length-k cycle through a→b: σ = (a, b, u_1, ..., u_{k-2}).
    # Weight = t(a,b) * t(b,u_1) * ... * t(u_{k-2}, a) * rest.

    # Summing over all σ with σ(a)=b, weighted by sgn:
    # = cofactor(a,b) of A.
    # cofactor(a,b) has degree n-1 (includes the a→b edge weight).
    # M has degree n-2 (doesn't include a→b edge weight).

    # INSIGHT: M[a,b] is like cofactor(a,b) but WITHOUT the edge a→b.
    # It's an "amputated" cofactor!

    # Compute: cofactor(a,b) / t(a,b) =? M[a,b]
    A_mat, _ = build_weight_matrix(V, r, s)
    cof_ab = expand((-1)**(a+b) * det(A_mat.minor_submatrix(a, b)))
    t_ab = edge_weight(a, b, r, s)

    # This division won't work symbolically unless cof_ab is divisible by t_ab.
    # Let me check.
    print(f"\n  cofactor({a},{b}) = {cof_ab}")
    print(f"  t({a},{b}) = {t_ab}")
    print(f"  M[{a},{b}] = {M_ab}")

    # Check if cofactor = t_ab * M + something
    diff = expand(cof_ab - t_ab * M_ab)
    print(f"  cofactor - t_ab * M = {diff}")

    # Maybe: cofactor(a,b) = t_ab * M[a,b] + correction?
    # Or: M[a,b] = (cofactor(a,b) - something) / t_ab?


def test_cycle_type_in_M(n):
    """
    Decompose M[a,b] by the cycle type of the underlying permutation.

    Each 2-path-cover (P_a, P_b) can be completed to a permutation σ
    by adding edges a→b and end(P_b)→start(P_a), but the resulting
    cycle structure depends on the partition.

    Instead: think of the PRODUCT E_a(S+a) * B_b(R+b) as summing
    over pairs of paths. Each pair is a permutation of V\{a,b} restricted
    to S and R respectively, composed with the endpoint constraints.

    The full permutation σ on V is obtained by:
    σ(v_i) = v_{i+1} for consecutive vertices in P_a (ending at a)
    σ(v_j) = v_{j+1} for consecutive vertices in P_b (starting at b)
    σ(a) = ??? and σ(last(P_b)) = ???

    These aren't permutations of V — they're PATHS. To make them into
    permutations, we'd need to close the cycles.

    ALTERNATIVE APPROACH: Think of M[a,b] via the EXPONENTIAL formula.

    Irving-Omar Theorem 10(b):
    U_D = L_n exp(sum_k (1/k) p_k [tr(XA)^k + (-1)^{k-1} tr(XĀ)^k])

    For tournaments: the exponent simplifies to
    sum_k_odd (2/k) p_k tr(XA)^k  (only odd k survive)

    Now: M[a,b] is related to the coefficient [x_a x_b prod_{u∈U} x_u]
    in s_{(1^n)} coefficient of U_D. But s_{(1^n)} extracts the permanent.

    Wait: [s_{(1^n)}] U_D = ham(D) = total number of Ham paths.
    And M[a,b] is a SIGNED count, not the total.

    Actually: sum_{a≠b} M[a,b] = 0 (odd n) or 2H (even n).
    And M[a,b] is the (a,b) entry of the transfer matrix.

    The connection to U_D is through the SCHUR expansion:
    [s_{(2,1^{n-2})}] U_D might relate to the transfer matrix.
    The hook partition (2,1^{n-2}) has a nice character theory.
    """
    print(f"\n{'='*60}")
    print(f"SCHUR COEFFICIENT CONNECTION: n={n}")
    print(f"{'='*60}")

    # For now, let's just verify the Irving-Omar even-cycle cancellation
    # at the level of INDIVIDUAL cycle types, using our weighted adjacency.

    r, s = make_symbols(n)
    V = list(range(n))

    # Weight of each cycle type
    ct_weights = defaultdict(lambda: 0)

    for perm in permutations(V):
        weight = 1
        for i in range(n):
            weight *= edge_weight(i, perm[i], r, s)

        # Cycle decomposition
        visited = [False] * n
        cycles = []
        for i in range(n):
            if not visited[i]:
                cycle = []
                j = i
                while not visited[j]:
                    visited[j] = True
                    cycle.append(j)
                    j = perm[j]
                cycles.append(len(cycle))

        ct = tuple(sorted(cycles, reverse=True))
        ct_weights[ct] = expand(ct_weights[ct] + weight)

    # For each cycle type with even parts, check r-parity
    print(f"\n  Cycle type weights and r-parity:")
    for ct in sorted(ct_weights.keys()):
        w = ct_weights[ct]
        has_even = any(l % 2 == 0 for l in ct)

        # Check if w is even in r
        if w != 0:
            p = Poly(w, r)
            r_degrees = sorted(p.as_dict().keys())
            odd_r = [d for d in r_degrees if d[0] % 2 == 1]
            even_r = [d for d in r_degrees if d[0] % 2 == 0]

            parity_label = "ODD r-powers" if odd_r else "even r-powers only"
            even_label = "HAS_EVEN_CYCLE" if has_even else "all_odd_cycles"

            print(f"  {ct} [{even_label}]: {parity_label}, r-degrees: {[d[0] for d in r_degrees]}")

            if has_even:
                # Check if this weight is ODD in r (only odd r-powers)
                even_r_terms = sum(c * r**d[0] for d, c in p.as_dict().items() if d[0] % 2 == 0)
                print(f"    Even-r-power part: {expand(even_r_terms)}")

    # KEY PREDICTION: For tournament-weighted A:
    # - Cycle types with ALL ODD parts → weight has EVEN r-parity (only r^0, r^2, ...)
    # - Cycle types with ANY EVEN part → weight has ??? r-parity
    # The even-cycle cancellation in Irving-Omar says the SUM of all even-cycle-type
    # weights vanishes in the p-expansion. But individually they might not.


def main():
    for n in [3, 4]:
        test_det_per_decomposition(n)

    for n in [3, 4, 5]:
        test_cycle_type_in_M(n)

    test_endpoint_conditioned_cycles(4)
    test_endpoint_conditioned_cycles(5)

    print(f"\n\n{'='*60}")
    print("SYNTHESIS")
    print(f"{'='*60}")
    print("""
IRVING-OMAR KEY INSIGHTS FOR EVEN R-POWERS:

1. EVEN CYCLE CANCELLATION (Corollary 19):
   For tournaments, tr(XĀ)^k = tr(XA)^k for k > 1.
   So: tr(XA)^k + (-1)^{k-1} tr(XĀ)^k = (1+(-1)^{k-1}) tr(XA)^k
   = 0 for k even, 2·tr(XA)^k for k odd.
   This makes U_D = sum over odd-cycle-only σ of 2^{ψ(σ)} p^{cyc}(σ).

2. CONNECTION TO EVEN R-POWERS:
   The even-cycle cancellation in the symmetric function expansion
   is the EXACT SAME mechanism as the even r-powers in M[a,b].
   Both arise from the T↔T^op symmetry:
   - For cycles: reversing an even cycle flips the sign
   - For r-powers: r → -r combined with T → T^op preserves M

3. THE MISSING LINK:
   M[a,b] is an ENDPOINT-CONDITIONED version of the Hamiltonian
   path generating function. The Irving-Omar formula gives the
   TOTAL (unconditioned) count. We need the endpoint-conditioned
   analogue that preserves the even-cycle cancellation mechanism.

4. AMPUTATED COFACTOR:
   cofactor(a,b) of A has degree n-1 (includes t(a,b)).
   M[a,b] has degree n-2 (excludes t(a,b)).
   If cofactor(a,b) = t(a,b) · M[a,b] + correction,
   this would express M as a "residue" of the cofactor.
""")


if __name__ == '__main__':
    main()
