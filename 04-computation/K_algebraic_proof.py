#!/usr/bin/env python3
"""
K_algebraic_proof.py -- Algebraic proof that K is constant for regular tournaments.

KEY RESULT FROM PRIOR ANALYSIS:
  K = c5 + butterfly - n*L*(L-1) - 3*c3
  where butterfly = 2*ov2 + 3*c3, L = (n^2-1)/8 = c3(v) per vertex.
  Simplifies to: K constant iff c5 + 2*ov2 = constant.

ALGEBRAIC APPROACH:
  For regular tournament A on n=2m+1 vertices:
    A + A^T = J - I
    A^T = J - I - A
    (A^T)^2 = I + 2A + A^2 - J

  c5 = tr(A^5)/5
  ov2 = sum_{i->j} C(lambda_{ij}, 2) where lambda_{ij} = (A^2)_{ji}

  Express c5 + 2*ov2 in terms of traces and use regularity.

Author: kind-pasteur-2026-03-12-S60
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
import random


def random_regular_tournament(n):
    assert n % 2 == 1
    m = (n - 1) // 2
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    for _ in range(10000):
        scores = A.sum(axis=1)
        if scores.max() == m and scores.min() == m:
            return A
        high = int(scores.argmax())
        low = int(scores.argmin())
        if high != low and A[high][low]:
            A[high][low] = 0
            A[low][high] = 1
        else:
            for j in range(n):
                if j != high and A[high][j]:
                    if scores[j] < m:
                        A[high][j] = 0
                        A[j][high] = 1
                        break
    scores = A.sum(axis=1)
    if all(s == m for s in scores):
        return A
    return None


def main():
    print("=" * 70)
    print("ALGEBRAIC PROOF: K CONSTANT FOR REGULAR TOURNAMENTS")
    print("=" * 70)

    # ====== PART 1: Verify the regular tournament matrix identity ======
    print("\n" + "=" * 60)
    print("PART 1: Matrix identities for regular tournaments")
    print("=" * 60)

    n = 7
    m = 3
    random.seed(42)

    A = random_regular_tournament(n)
    I_n = np.eye(n, dtype=int)
    J = np.ones((n, n), dtype=int)

    AT = A.T
    A2 = A @ A
    AT2 = AT @ AT  # = (A^2)^T

    # Verify: A + A^T = J - I
    print(f"  A + A^T = J - I: {np.array_equal(A + AT, J - I_n)}")

    # Verify: (A^T)^2 = I + 2A + A^2 - J
    formula = I_n + 2*A + A2 - J
    print(f"  (A^T)^2 = I + 2A + A^2 - J: {np.array_equal(AT2, formula)}")
    print(f"  (A^T)^2 = I + 2A + A^2 (WITHOUT -J): {np.array_equal(AT2, I_n + 2*A + A2)}")

    # Additional identity: A*J = J*A = m*J
    print(f"  A*J = m*J: {np.array_equal(A @ J, m * J)}")
    print(f"  J*A = m*J: {np.array_equal(J @ A, m * J)}")

    # ====== PART 2: Butterfly in terms of matrix products ======
    print("\n" + "=" * 60)
    print("PART 2: Butterfly sum decomposition")
    print("=" * 60)

    # butterfly = sum_{i,j} A_{ij} * [(A^2)_{ji}]^2
    #           = sum_{i,j} A_{ij} * [(A^T)^2_{ij}]^2
    # Using (A^T)^2 = I + 2A + A^2 - J, let B = I + 2A + A^2 - J:
    # butterfly = sum_{ij} A_{ij} * B_{ij}^2
    # = sum_{ij} A_{ij} * (I_{ij} + 2A_{ij} + (A^2)_{ij} - J_{ij})^2

    # Since A has 0 diagonal and J_{ij}=1 always, I_{ij}=delta_{ij}:
    # For i != j: B_{ij} = 2*A_{ij} + (A^2)_{ij} - 1
    # For i = j: B_{ii} = 1 + 0 + (A^2)_{ii} - 1 = (A^2)_{ii}
    # But A_{ii}=0, so diagonal terms don't contribute to butterfly.

    # For i != j with A_{ij}=1 (i->j):
    # B_{ij} = 2*1 + (A^2)_{ij} - 1 = 1 + (A^2)_{ij}
    # lambda_{ij} = (A^2)_{ji} = B_{ij} = 1 + (A^2)_{ij}

    # Hmm wait, this is interesting: for edge i->j,
    # lambda_{ij} = (A^2)_{ji} = B_{ij} = 2*A_{ij} + (A^2)_{ij} - 1 = 2 + (A^2)_{ij} - 1 = 1 + (A^2)_{ij}

    # Verify this: lambda_{ij} = # {k: j->k->i} = (A^2)_{ji}
    # And (A^2)_{ij} = # {k: i->k->j} = # 2-paths from i to j
    # For edge i->j: lambda_{ij} = 1 + (A^2)_{ij}?

    # Check: for Paley T_7 with edge 0->1:
    # (A^2)_{10} = # 2-paths from 1 to 0 = # {k: 1->k->0}
    # (A^2)_{01} = # 2-paths from 0 to 1 = # {k: 0->k->1}
    # Are these related by lambda = 1 + (A^2)_{ij}? Let's check.

    # Direct computation on actual edges and non-edges:
    edges_shown = 0
    non_edges_shown = 0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            lam = int(A2[j][i])
            a2_ij = int(A2[i][j])
            b_ij = 2*int(A[i][j]) + a2_ij - 1
            if A[i][j] == 1 and edges_shown < 3:
                print(f"  Edge {i}->{j}: lambda={lam}, (A^2)_{{ij}}={a2_ij}, "
                      f"B_{{ij}}=2A+A2-1={b_ij}, match={lam==b_ij}")
                edges_shown += 1
            elif A[i][j] == 0 and non_edges_shown < 3:
                print(f"  Non-edge {i}->{j}: (A^2)_{{ji}}={lam}, (A^2)_{{ij}}={a2_ij}, "
                      f"B_{{ij}}=A2-1={b_ij}, match={lam==b_ij}")
                non_edges_shown += 1

    # So for ANY i != j in a regular tournament:
    # (A^2)_{ji} = 2*A_{ij} + (A^2)_{ij} - 1
    # i.e., (A^2)_{ji} + 1 = 2*A_{ij} + (A^2)_{ij}

    # This is equivalent to: (A^2)_{ji} - (A^2)_{ij} = 2*A_{ij} - 1
    # When i->j (A_{ij}=1): (A^2)_{ji} = (A^2)_{ij} + 1
    # When j->i (A_{ij}=0): (A^2)_{ji} = (A^2)_{ij} - 1

    # Verify this identity for ALL edges:
    all_match = True
    for i in range(n):
        for j in range(n):
            if i != j:
                diff = int(A2[j][i]) - int(A2[i][j])
                expected = 2*int(A[i][j]) - 1
                if diff != expected:
                    all_match = False
                    print(f"  FAIL: ({i},{j}): diff={diff}, expected={expected}")
    print(f"\n  (A^2)_{{ji}} - (A^2)_{{ij}} = 2*A_{{ij}} - 1 for all i!=j: {all_match}")

    # This is a KNOWN IDENTITY for regular tournaments!
    # It follows from: A^2 - (A^T)^2 = A^2 - (I + 2A + A^2 - J) = -(I + 2A - J)
    #                                = J - I - 2A = (J-I) - 2A = (A+A^T) - 2A = A^T - A
    # So: (A^2 - A^{T2})_{ij} = (A^T-A)_{ij} = A_{ji} - A_{ij} = -(2A_{ij}-1) for i!=j
    # I.e., (A^2)_{ij} - (A^{T2})_{ij} = -(2A_{ij}-1) = 1-2A_{ij}
    # So (A^{T2})_{ij} = (A^2)_{ij} + 2A_{ij} - 1
    # And (A^{T2})_{ij} = (A^2)_{ji} (since A^{T2} = (A^2)^T)
    # Therefore: (A^2)_{ji} = (A^2)_{ij} + 2A_{ij} - 1. Confirmed!

    # ====== PART 3: Expand butterfly using the identity ======
    print("\n" + "=" * 60)
    print("PART 3: Expand butterfly algebraically")
    print("=" * 60)

    # For edge i->j (A_{ij}=1):
    # lambda_{ij} = (A^2)_{ji} = (A^2)_{ij} + 1
    # Let mu_{ij} = (A^2)_{ij} (number of 2-paths i->k->j)
    # Then lambda_{ij} = mu_{ij} + 1

    # butterfly = sum_{i->j} lambda^2 = sum_{i->j} (mu+1)^2
    #           = sum_{i->j} (mu^2 + 2*mu + 1)
    #           = sum_{i->j} mu^2 + 2*sum_{i->j} mu + n*m

    # Now, sum_{i->j} mu_{ij} = sum_{i->j} (A^2)_{ij}
    # = sum_{ij} A_{ij} * (A^2)_{ij}
    # = sum_{ijk} A_{ij} * A_{ik} * A_{kj}   (expanding (A^2)_{ij} = sum_k A_{ik}*A_{kj})

    # This counts: oriented paths i->j where there exists k with i->k->j and i->j.
    # I.e., directed triangles where we mark one "direct" edge.
    # Each directed triangle {i->j->k->i} contributes when (the direct edge) is
    # one of: i->j (with 2-path i->k->j), j->k (with 2-path j->i->k), k->i (with 2-path k->j->i).
    # So each directed triangle contributes 3 to this sum.
    # But wait, each UNDIRECTED triangle contributes from BOTH orientations...
    # Actually there are two directed 3-cycles per triangle vertex set: i->j->k->i and i->k->j->i.
    # For each, the 3 edges contribute. So 6 per triangle vertex set.
    # sum_{i->j} mu_{ij} = 6 * c3_undirected = 6c3? No...

    # Actually wait. A directed 3-cycle is one of TWO orientations. Each triangle (vertex set)
    # is exactly ONE directed 3-cycle (since tournament). So we have c3 directed 3-cycles.
    # Each has 3 directed edges. For each edge e in the cycle, the other two edges form a
    # 2-path connecting the endpoints of e in the right direction. So mu_e >= 1 for each edge in a cycle.

    # sum_{i->j} mu_{ij} = sum_{ij} A_{ij} * (A^2)_{ij}
    # = trace of A * Hadamard * A^2... no, this is element-wise, not trace.
    # = sum of entries of A .* A^2 (Hadamard product)

    sum_A_hadamard_A2 = int((A * A2).sum())
    print(f"  sum A.*A^2 = {sum_A_hadamard_A2}")
    print(f"  3*c3 = {3 * int(np.trace(A2 @ A)) // 3}")

    # Hmm, let me check: sum_{ij} A_{ij}*A_{ik}*A_{kj} summed over i,j,k
    # = sum_i (sum_j A_{ij})*(sum_k A_{ik}) - ... no that's not right
    # sum_{ij} A_{ij} * (A^2)_{ij} = sum_{ijk} A_{ij}*A_{ik}*A_{kj}

    # For k=j: A_{ij}*A_{ij}*A_{jj} = 0 (diagonal of A is 0)
    # For k != j, k != i: this is i->j AND i->k AND k->j.
    # So it counts ordered triples (i,j,k) with i->j, i->k, k->j, all distinct.
    # This is the number of "transitive triples" where i dominates both j and k, and k dominates j.
    # Each such triple is a transitive triangle i->j, i->k, k->j.
    # Each transitive triple contributes once.
    # Number of transitive triples = C(n,3) - c3 (all triples are either transitive or cyclic).
    # But: each transitive triple has a UNIQUE vertex that dominates both others.
    # The count sum_{ijk} counts ORDERED (i,j,k) so each transitive triple with top vertex i
    # contributes 2 (for the two orderings of j,k among the bottom two).

    # Wait, let me be more careful.
    # sum_{ijk, all distinct} A_{ij}*A_{ik}*A_{kj} counts:
    # Ordered triples (i,j,k) where i->j, i->k, k->j.
    # This means i is at the top, k is in the middle, j is at the bottom.
    # The triple {i,j,k} is TRANSITIVE with i winning against j and k, and k winning against j.
    # For each transitive triple with dominance order i > k > j, this is counted once (i,j,k).
    # Are there other orderings? Let's see:
    # We need i->j, i->k, k->j. In the transitive triple i>k>j:
    # - (i,j,k): i->j=yes, i->k=yes, k->j=yes. COUNTED.
    # - (i,k,j): i->k=yes, i->j=yes, j->k=? No, we need k->j. And j->k? Since k>j, k->j. So A_{jk}=0, not counted in this formulation.
    # Actually the sum is over distinct i,j,k. The condition is A_{ij}*A_{ik}*A_{kj}.
    # (i,j,k): i->j, i->k, k->j. Top=i, mid=k, bot=j.
    # (i,k,j): i->k, i->j, j->k. But j->k means j>k which contradicts k>j. So NOT counted.
    # Only ONE ordered triple per transitive triple.
    # So: sum_{ijk distinct} A_{ij}*A_{ik}*A_{kj} = # transitive triples = C(n,3) - c3.

    # Including k=i or k=j: k=i gives A_{ij}*A_{ii}*A_{ij}=0. k=j gives A_{ij}*A_{ij}*A_{jj}=0.
    # So the sum over all i,j,k (with i!=j but k free) = same as distinct.

    # Therefore: sum_{ij} A_{ij}*(A^2)_{ij} = C(n,3) - c3

    c3_val = int(np.trace(A2 @ A)) // 3
    transitive = n*(n-1)*(n-2)//6 - c3_val
    print(f"  c3 = {c3_val}, C(n,3) = {n*(n-1)*(n-2)//6}, transitive = {transitive}")
    print(f"  sum A.*A^2 = {sum_A_hadamard_A2}, C(n,3)-c3 = {transitive}")
    print(f"  Match: {sum_A_hadamard_A2 == transitive}")

    # GREAT! So for ANY tournament:
    # sum_{i->j} mu_{ij} = C(n,3) - c3

    # For REGULAR tournament: c3 = n(n^2-1)/24
    # C(n,3) - c3 = n(n-1)(n-2)/6 - n(n^2-1)/24
    # = n(n-1)[4(n-2) - (n+1)] / 24
    # = n(n-1)(3n-9) / 24
    # = n(n-1)(n-3) / 8
    # = n * C(n-1,2) * ... hmm
    # = n(n-1)(n-3)/8

    print(f"  For regular: C(n,3)-c3 = n(n-1)(n-3)/8 = {n*(n-1)*(n-3)//8}")

    # ====== PART 4: The second moment of mu ======
    print("\n" + "=" * 60)
    print("PART 4: sum_{i->j} mu^2 (second moment)")
    print("=" * 60)

    # We need: sum_{i->j} mu_{ij}^2 = sum_{ij} A_{ij} * [(A^2)_{ij}]^2
    # = sum_{ij} A_{ij} * [sum_k A_{ik}*A_{kj}]^2
    # = sum_{ij} A_{ij} * sum_{kl} A_{ik}*A_{kj}*A_{il}*A_{lj}
    # = sum_{ijkl} A_{ij}*A_{ik}*A_{kj}*A_{il}*A_{lj}

    # This counts: i->j, i->k->j, i->l->j (with k,l possibly equal).
    # Geometrically: vertex i beats j directly, and has two 2-paths through k and l.

    sum_mu_sq = int((A * A2 * A2).sum())
    print(f"  sum mu^2 = {sum_mu_sq}")

    # Now: butterfly = sum (mu+1)^2 = sum mu^2 + 2*sum mu + n*m
    butterfly_formula = sum_mu_sq + 2*(n*(n-1)*(n-2)//6 - c3_val) + n*m
    butterfly_direct = int((A * (A2.T)**2).sum())

    # Wait, butterfly = sum_{i->j} lambda^2 where lambda = (A^2)_{ji}
    # And we showed lambda = mu + 1.
    # So butterfly = sum (mu+1)^2. But let me verify.
    butterfly_check = 0
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                lam = int(A2[j][i])
                mu = int(A2[i][j])
                assert lam == mu + 1, f"lambda={lam}, mu+1={mu+1}"
                butterfly_check += lam**2

    print(f"  butterfly (direct) = {butterfly_check}")
    print(f"  butterfly (formula sum(mu+1)^2) = {butterfly_formula}")
    print(f"  Match: {butterfly_check == butterfly_formula}")

    # ====== PART 5: c5 + 2*ov2 expansion ======
    print("\n" + "=" * 60)
    print("PART 5: c5 + 2*ov2 = c5 + butterfly - 3*c3")
    print("=" * 60)

    # From butterfly = 2*ov2 + 3*c3 (verified in previous script):
    # 2*ov2 = butterfly - 3*c3

    # c5 + 2*ov2 = c5 + butterfly - 3*c3
    # = tr(A^5)/5 + sum(mu+1)^2 - 3*c3
    # = tr(A^5)/5 + sum mu^2 + 2*sum mu + n*m - 3*c3
    # = tr(A^5)/5 + sum mu^2 + 2*(C(n,3)-c3) + n*m - 3*c3
    # = tr(A^5)/5 + sum mu^2 + 2*C(n,3) - 2*c3 + n*m - 3*c3
    # = tr(A^5)/5 + sum mu^2 + 2*C(n,3) + n*m - 5*c3

    # For regular: c3 = n(n^2-1)/24, C(n,3) = n(n-1)(n-2)/6, n*m = n(n-1)/2
    # 2*C(n,3) + n*m - 5*c3 = n(n-1)(n-2)/3 + n(n-1)/2 - 5n(n^2-1)/24
    # = n(n-1)[8(n-2) + 12 - 5(n+1)] / 24
    # = n(n-1)[8n-16+12-5n-5] / 24
    # = n(n-1)(3n-9) / 24
    # = n(n-1)(n-3) / 8
    # Hmm, same as C(n,3) - c3.

    reg_const = n*(n-1)*(n-3)//8
    alg_const = 2*n*(n-1)*(n-2)//6 + n*m - 5*c3_val
    print(f"  2*C(n,3) + n*m - 5*c3 = {alg_const}")
    print(f"  n(n-1)(n-3)/8 = {reg_const}")

    # Hmm that's not equal. Let me recompute:
    # 2*C(n,3) = 2*7*6*5/6 = 70
    # n*m = 7*3 = 21
    # 5*c3 = 5*14 = 70
    # Sum: 70 + 21 - 70 = 21
    # n(n-1)(n-3)/8 = 7*6*4/8 = 21. YES they match!

    print(f"  Match: {alg_const == reg_const}")

    # So: c5 + 2*ov2 = tr(A^5)/5 + sum_{i->j} mu^2 + n(n-1)(n-3)/8
    # And this should be constant. So: tr(A^5)/5 + sum mu^2 = constant.

    # Let's verify this!
    target_const = -3*n*(n**2-1)*(n**2-9)//320 + n*(n**2-1)//8 * ((n**2-1)//8 - 1)
    target_c5_2ov2 = target_const  # K + n*L*(L-1)

    # So: tr(A^5)/5 + sum mu^2 = target_c5_2ov2 - n(n-1)(n-3)/8
    target_trA5_sum_mu2 = target_c5_2ov2 - reg_const

    print(f"\n  Target for c5 + 2*ov2 = {target_c5_2ov2}")
    print(f"  Target for tr(A^5)/5 + sum mu^2 = {target_trA5_sum_mu2}")

    # Verify across tournaments
    random.seed(42)
    for trial in range(15):
        A = random_regular_tournament(n)
        if A is None:
            continue
        A2 = A @ A
        A5 = A @ A2 @ A2
        c5 = int(np.trace(A5)) // 5
        sum_mu2 = int((A * A2 * A2).sum())
        val = c5 + sum_mu2
        print(f"    Trial {trial}: c5={c5}, sum_mu2={sum_mu2}, "
              f"c5+sum_mu2={val}, target={target_trA5_sum_mu2}, "
              f"{'OK' if val == target_trA5_sum_mu2 else 'FAIL'}")

    # ====== PART 6: Express sum mu^2 in trace terms ======
    print("\n" + "=" * 60)
    print("PART 6: sum mu^2 as a trace expression")
    print("=" * 60)

    # sum_{ij} A_{ij} * (A^2)_{ij}^2 = sum_{ij} A_{ij} * (A^2 .* A^2)_{ij}
    # This is NOT a standard trace but a Hadamard-weighted sum.
    #
    # However, we can try to expand:
    # (A^2)_{ij}^2 = [sum_k A_{ik}A_{kj}]^2 = sum_{kl} A_{ik}A_{kj}A_{il}A_{lj}
    # So: sum_{ij} A_{ij}*(A^2)_{ij}^2 = sum_{ijkl} A_{ij}A_{ik}A_{kj}A_{il}A_{lj}
    #
    # This can be split: k=l and k!=l.
    # k=l: sum_{ijk} A_{ij}A_{ik}^2*A_{kj}^2 = sum_{ijk} A_{ij}A_{ik}A_{kj}
    # (since A is 0-1, A^2=A)
    # = sum_{ij} A_{ij}*(A^2)_{ij} = C(n,3) - c3

    # k!=l: sum_{ijkl, k!=l} A_{ij}A_{ik}A_{kj}A_{il}A_{lj}
    # This counts: i->j, i->k->j, i->l->j with k!=l.
    # I.e., two DISTINCT 2-paths from i to j, plus the direct edge.

    # So: sum mu^2 = (C(n,3)-c3) + sum_{i->j} mu*(mu-1)
    # But sum mu*(mu-1) = sum mu^2 - sum mu. So this is circular!
    # That's expected: sum mu^2 = k=l part + k!=l part = sum mu + (sum mu^2 - sum mu) = sum mu^2.

    # Let me try a different decomposition.
    # sum_{ijkl} A_{ij}A_{ik}A_{kj}A_{il}A_{lj}
    # = sum_i [sum_j A_{ij} * prod over two paths from i through j]
    # = sum_i sum_{j,k,l} A_{ij}*A_{ik}*A_{kj}*A_{il}*A_{lj}
    # = sum_i sum_{k,l} A_{ik}*A_{il} * [sum_j A_{ij}*A_{kj}*A_{lj}]

    # For fixed i,k,l: sum_j A_{ij}*A_{kj}*A_{lj} = # vertices j beaten by all of i,k,l
    # This is the "common out-neighbor" count.

    # Alternative: think of this as a trace of something.
    # sum_{ijkl} A_{ij}*A_{ik}*A_{kj}*A_{il}*A_{lj}
    # Let C = A .* A^2 (Hadamard product). Then C_{ij} = A_{ij}*(A^2)_{ij} = mu_{ij}*A_{ij}.
    # sum mu^2 = sum_{ij} A_{ij}*(A^2)_{ij}^2 = sum_{ij} C_{ij}*(A^2)_{ij}
    # = sum_{ijk} C_{ij}*A_{ik}*A_{kj}
    # = tr(C*A^2)... wait, sum_{ij} C_{ij}*(A^2)_{ij} is Hadamard inner product, not matrix product.

    # Actually, sum_{ij} X_{ij}*Y_{ij} = tr(X^T * Y) when X,Y are matrices.
    # Wait: tr(X^T Y) = sum_j (X^T Y)_{jj} = sum_{jk} X_{kj}*Y_{kj} = sum_{kj} X_{kj}*Y_{kj}
    # = sum_{ij} X_{ij}*Y_{ij}. Yes!

    # So: sum mu^2 = sum_{ij} A_{ij}*(A^2)_{ij}^2 = sum_{ij} (A.*A^2)_{ij}*(A^2)_{ij}
    # = tr((A.*A^2)^T * A^2) = tr(C^T * A^2)

    C = A * A2  # Hadamard product
    val_trace = int(np.trace(C.T @ A2))
    print(f"  sum mu^2 = {sum_mu_sq}")
    print(f"  tr(C^T * A^2) = {val_trace}")
    print(f"  Match: {sum_mu_sq == val_trace}")

    # Good! But tr(C^T * A^2) = tr((A.*A^2)^T * A^2) is still a Hadamard expression.

    # Let's try to express sum_mu^2 as a polynomial in standard traces.
    # Use the expansion:
    # sum_{ijkl} A_{ij}A_{ik}A_{kj}A_{il}A_{lj}
    # = sum_{i} sum_{kl} A_{ik}A_{il} sum_j A_{ij}A_{kj}A_{lj}

    # For fixed i: let out(i) be the out-neighborhood.
    # sum_{kl} A_{ik}A_{il} sum_j A_{ij}A_{kj}A_{lj}
    # = sum_{k in out(i)} sum_{l in out(i)} |out(i) cap out(k) cap out(l)|

    # This is getting complex. Let me instead try to find the value
    # by fitting a polynomial in known traces.

    # We need: c5 + sum_mu2 = constant for regular tournaments.
    # c5 = tr(A^5)/5.
    # Let's check: is tr(A^5)/5 + tr(C^T A^2) = const?

    # Already verified above. Now let's try different trace combos:
    random.seed(42)
    data = []
    for _ in range(500):
        A = random_regular_tournament(n)
        if A is None:
            continue
        A2 = A @ A
        A3 = A @ A2
        AT = A.T
        AT2 = AT @ AT
        AAT = A @ AT

        trA5 = int(np.trace(A3 @ A2))
        trA3AT = int(np.trace(A3 @ AT))
        trA2AT2 = int(np.trace(A2 @ AT2))
        trAAT2 = int(np.trace(AAT @ AAT))
        trA3 = int(np.trace(A3))
        sum_mu2 = int((A * A2 * A2).sum())

        data.append({
            'trA5': trA5, 'trA3AT': trA3AT, 'trA2AT2': trA2AT2,
            'trAAT2': trAAT2, 'sum_mu2': sum_mu2, 'trA3': trA3
        })
        if len(data) >= 20:
            break

    # Check if sum_mu2 is a linear combination of standard degree-5 traces
    print(f"\n  Checking if sum_mu^2 = a*trA5 + b*trA3AT + c*trA2AT2 + d*trAAT^2 + const")
    print(f"  Data points:")
    for d in data[:10]:
        print(f"    trA5={d['trA5']}, trA3AT={d['trA3AT']}, trA2AT2={d['trA2AT2']}, "
              f"trAAT2={d['trAAT2']}, sum_mu2={d['sum_mu2']}")

    # Use numpy to solve: sum_mu2 = a*trA5 + b*trA3AT + c*trA2AT2 + d*trAAT^2 + e
    if len(data) >= 5:
        X = np.array([[d['trA5'], d['trA3AT'], d['trA2AT2'], d['trAAT2'], 1] for d in data], dtype=float)
        y = np.array([d['sum_mu2'] for d in data], dtype=float)
        coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
        print(f"\n  Fit: sum_mu2 = {coeffs[0]:.6f}*trA5 + {coeffs[1]:.6f}*trA3AT "
              f"+ {coeffs[2]:.6f}*trA2AT2 + {coeffs[3]:.6f}*trAAT2 + {coeffs[4]:.6f}")
        if len(residuals) > 0:
            print(f"  Residuals: {residuals}")
        else:
            print(f"  Residuals: (perfect fit or underdetermined)")

        # Check fit
        for d in data[:5]:
            pred = (coeffs[0]*d['trA5'] + coeffs[1]*d['trA3AT'] +
                    coeffs[2]*d['trA2AT2'] + coeffs[3]*d['trAAT2'] + coeffs[4])
            print(f"    Predicted={pred:.2f}, Actual={d['sum_mu2']}, "
                  f"Error={abs(pred-d['sum_mu2']):.6f}")

    # ====== PART 7: Does the identity extend to non-prime odd n? ======
    print("\n" + "=" * 60)
    print("PART 7: Does K constant hold for non-prime odd n?")
    print("=" * 60)

    for nn in [5, 7, 9, 11, 13, 15]:
        mm = (nn - 1) // 2
        K_formula = -3 * nn * (nn**2 - 1) * (nn**2 - 9) // 320
        K_vals = set()
        random.seed(42)
        count = 0
        for _ in range(3000):
            A = random_regular_tournament(nn)
            if A is None:
                continue
            count += 1
            A2 = A @ A
            c5 = int(np.trace(A @ A2 @ A2)) // 5

            # Compute ov2 directly
            ov2 = 0
            for i in range(nn):
                for j in range(nn):
                    if A[i][j]:
                        lam = int(A2[j][i])
                        ov2 += lam * (lam - 1) // 2

            # Compute ov1 via c3(v) and ov2
            c3_val = int(np.trace(A2 @ A)) // 3
            c3v = np.zeros(nn, dtype=int)
            for v in range(nn):
                for a in range(nn):
                    if a == v: continue
                    for b in range(a+1, nn):
                        if b == v: continue
                        if (A[v][a] and A[a][b] and A[b][v]) or (A[v][b] and A[b][a] and A[a][v]):
                            c3v[v] += 1
            sum_c3v_c2 = int(sum(cv*(cv-1)//2 for cv in c3v))
            ov1 = sum_c3v_c2 - 2*ov2

            K = c5 - 2*ov1 - 2*ov2
            K_vals.add(K)

            if count >= 50:
                break

        status = "CONSTANT" if len(K_vals) == 1 else f"VARIES ({len(K_vals)} values)"
        K_list = sorted(K_vals)
        match_str = ""
        if len(K_vals) == 1:
            match_str = f", formula={K_formula}, match={K_list[0]==K_formula}"

        print(f"  n={nn}: {count} tournaments, K={K_list[:5]}{'...' if len(K_list)>5 else ''} -- {status}{match_str}")

    # ====== PART 8: The lambda -> mu identity more carefully ======
    print("\n" + "=" * 60)
    print("PART 8: The key identity A^2 - A^{T2} = A^T - A")
    print("=" * 60)

    # For ANY tournament: A + A^T = J - I (complete tournament property)
    # So A^T = J - I - A.
    # A^2 - (A^T)^2 = A^2 - (J-I-A)^2
    # Let's expand (J-I-A)^2:
    # = J^2 - J - JA - J + I + A - AJ + A + A^2
    # For REGULAR: JA = AJ = mJ.
    # = nJ - J - mJ - J + I + A - mJ + A + A^2
    # = (n-2-2m)J + I + 2A + A^2
    # With n=2m+1: (2m+1-2-2m)J = -J.
    # So (J-I-A)^2 = -J + I + 2A + A^2
    #
    # A^2 - (A^T)^2 = A^2 - (-J + I + 2A + A^2) = J - I - 2A = (A+A^T) - 2A = A^T - A
    #
    # So: A^2 - A^{T2} = A^T - A = (J-I) - 2A    ... (*)
    #
    # This means: for regular tournaments, (A^2)_{ij} - (A^2)_{ji} = A_{ji} - A_{ij}
    # For i->j: (A^2)_{ij} - (A^2)_{ji} = 0 - 1 = -1
    # I.e., mu - lambda = -1, so lambda = mu + 1. Confirmed!

    print("  For regular tournaments:")
    print("  A^2 - (A^T)^2 = A^T - A = (J-I) - 2A")
    print("  Entry-wise for i!=j: (A^2)_{ij} - (A^2)_{ji} = A_{ji} - A_{ij}")
    print("  For edge i->j: lambda_{ij} = mu_{ij} + 1")
    print()

    # ====== PART 9: Why c5/5 + sum mu^2 = const ======
    print("=" * 60)
    print("PART 9: The central identity")
    print("=" * 60)

    # We've reduced to: for regular tournaments,
    # tr(A^5)/5 + sum_{ij} A_{ij}*(A^2)_{ij}^2 = constant
    #
    # Both terms involve degree-5 monomials in A.
    # tr(A^5)/5 = sum_{i0,...,i4} A_{i0,i1}*A_{i1,i2}*A_{i2,i3}*A_{i3,i4}*A_{i4,i0} / 5
    # sum mu^2 = sum_{i,j,k,l} A_{ij}*A_{ik}*A_{kj}*A_{il}*A_{lj}
    #
    # The first counts 5-cycles (up to starting point).
    # The second counts "double-path + direct edge" configurations.
    #
    # For the proof, we need to show that the sum of these is determined by n alone
    # when A is regular.
    #
    # Key approach: use A^T = J - I - A to express everything in terms of A and J.
    # For regular A, any polynomial in A and J can be evaluated using:
    # - A*J = J*A = m*J
    # - J^2 = n*J
    # - tr(J) = n
    # - A + A^T = J - I
    #
    # The trace tr(A^k) for a regular tournament depends on the specific tournament
    # for k >= 5 (it's constant for k=1,2,3,4 but varies for k=5+).
    #
    # However, sum mu^2 also varies, and the SUM c5 + sum_mu2 is constant.
    # This suggests a CANCELLATION that follows from the regularity identities.

    # Let me try to prove this using the A^T substitution.
    # sum_mu2 = sum_{ijkl} A_{ij}A_{ik}A_{kj}A_{il}A_{lj}
    # Fix an index pattern. The vertices are i,j,k,l (not necessarily distinct).
    #
    # When all 5 indices in tr(A^5) are distinct, we get 5-cycles.
    # When some coincide, we get shorter structures.
    # Similarly for sum_mu2.

    # Actually, maybe a more productive approach: compute both as polynomials in
    # the eigenvalues of A.
    #
    # For a regular tournament on n=2m+1 vertices:
    # A has eigenvalue m (eigenvector: all-ones) and n-1 other eigenvalues
    # whose real parts are all -1/2 (well-known).
    #
    # Let lambda_1,...,lambda_{n-1} be the non-trivial eigenvalues.
    # tr(A^5) = m^5 + sum lambda_k^5
    # But sum_mu^2 is NOT a simple trace, so it doesn't decompose by eigenvalues alone.

    # The fact that sum_mu2 + c5 is constant despite NEITHER being individually constant
    # is truly remarkable. Let me compute the constant value for general n.

    print(f"\n  Values of c5 + sum_mu^2 for various n:")
    for nn in [5, 7, 9]:
        mm = (nn - 1) // 2
        random.seed(42)
        vals = set()
        for _ in range(2000):
            A = random_regular_tournament(nn)
            if A is None:
                continue
            A2 = A @ A
            c5 = int(np.trace(A @ A2 @ A2)) // 5
            sum_mu2 = int((A * A2 * A2).sum())
            vals.add(c5 + sum_mu2)
            if len(vals) > 1:
                break
        print(f"  n={nn}: c5+sum_mu2 = {sorted(vals)}")

    # The constant should be: target_c5_2ov2 - n(n-1)(n-3)/8
    # where target_c5_2ov2 = K + n*L*(L-1)
    # = -3n(n^2-1)(n^2-9)/320 + n*(n^2-1)/8*((n^2-1)/8 - 1)

    print(f"\n  Predicted constants:")
    for nn in [5, 7, 9, 11, 13]:
        mm = (nn - 1) // 2
        LL = (nn**2 - 1) // 8
        K = -3*nn*(nn**2-1)*(nn**2-9)//320
        target = K + nn*LL*(LL-1) - nn*(nn-1)*(nn-3)//8
        print(f"  n={nn}: c5+sum_mu2 = {target}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
