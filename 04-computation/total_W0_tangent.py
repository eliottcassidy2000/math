#!/usr/bin/env python3
"""
Total sum of W(0) over all tournaments — tangent number connection.

At n=5 (odd): W(0) = constant term of W(r) = F_4(0) + 2*F_2(0)*t3 + 2*F_0(0)*t5
= 1 - t3 + 2*t5.

Sum over all 2^m tournaments of W(0) = sum (1 - t3 + 2*t5)
= 2^m - sum(t3) + 2*sum(t5)

Average t3 at n=5 with backbone: need to compute.
Average t5 at n=5 with backbone: need to compute.

But also: sum W(0) = sum_T prod_{perm} prod_step s_i
Can we compute this from the transfer matrix?

W(0) = (1/2)^{n-1} * sum_perm (-1)^{descents(P,T)}
where descent means backward edge.

Sum over T: sum_T W_T(0) = (1/2)^{n-1} * sum_T sum_perm (-1)^{desc}
= (1/2)^{n-1} * sum_perm sum_T (-1)^{desc(P,T)}

For a fixed perm P, each non-backbone edge (i,j) in the path contributes
independently: A[P(i),P(i+1)] is 1 or 0.
If (P(i),P(i+1)) is a backbone edge: always forward, contributes s_i = +1/2.
If not: contributes +1/2 or -1/2 with equal probability over tournaments.

So for fixed perm P:
sum_T prod s_i = prod_step [if backbone: (1/2)] * [if not backbone: (1/2 + (-1/2))/2 = 0]

Wait! If any step is a non-backbone edge, the product averages to 0!

But that would make sum_T W_T(0) = 0 (since most perms have non-backbone edges).
Hmm, but only perms using ONLY backbone edges would contribute.
The only perm using only backbone edges is the identity: 0,1,...,n-1.
For this perm, all edges are forward, so W(0) contribution = (1/2)^{n-1}.

But wait, 2^m tournaments * (1/2)^{n-1} per identity perm = 2^m / 2^{n-1}.
At n=5: 2^6 / 2^4 = 4. But our computed sum was 8.

Let me reconsider. The issue is that sum_T prod s_i for fixed perm P depends
on which edges are backbone.

Actually, for a fixed pair (P(i), P(i+1)):
- If it's a backbone edge (consecutive in 0,1,...,n-1): s_i = +1/2 for ALL T
- If it's a non-backbone edge: s_i = +1/2 for half the tournaments,
  s_i = -1/2 for the other half.
  Average s_i over T = 0.
- BUT we want average of PRODUCT, not product of averages!

Actually for INDEPENDENT non-backbone edges: average of product = product of averages = 0.
So any perm with at least one non-backbone edge has average product = 0.

Only the identity perm (and maybe its reverse) use only backbone edges.

Identity: 0->1->2->...->n-1, all edges are backbone (consecutive), all forward.
Product = (1/2)^{n-1}.

Reverse: n-1->n-2->...->0, all edges are backbone (consecutive) but BACKWARD.
Product = (-1/2)^{n-1}.

These are the only two perms using only backbone edges.

So: sum_T W_T(0) = 2^m * [(1/2)^{n-1} + (-1/2)^{n-1}]
At odd n: = 2^m * [(1/2)^{n-1} + (-1/2)^{n-1}] = 2^m * 2 * (1/2)^{n-1} = 2^{m-n+2}
At even n: = 2^m * [(1/2)^{n-1} - (1/2)^{n-1}] = 0

Wait, let me be more careful.
At odd n: (-1/2)^{n-1} = (1/2)^{n-1} (since n-1 is even)
At even n: (-1/2)^{n-1} = -(1/2)^{n-1} (since n-1 is odd)

So sum_T W_T(0):
- Odd n: 2^m * 2 * (1/2)^{n-1} = 2^{m-n+2}
- Even n: 0

Let's verify.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations
from fractions import Fraction

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def compute_W0(A, n):
    half = Fraction(1, 2)
    total = Fraction(0)
    for perm in permutations(range(n)):
        prod = Fraction(1)
        for i in range(n-1):
            prod *= (Fraction(A[perm[i]][perm[i+1]]) - half)
        total += prod
    return total

for n in [3, 4, 5]:
    m = num_tiling_bits(n)
    total_W0 = Fraction(0)
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        total_W0 += compute_W0(A, n)

    predicted = 2**(m-n+2) if n % 2 == 1 else 0
    print(f"n={n}: m={m}, sum W(0) = {total_W0} = {float(total_W0):.4f}, predicted = {predicted}")

# Also check sum W(1/2) = sum H(T) = total Ham paths
for n in [3, 4, 5]:
    m = num_tiling_bits(n)
    total_H = 0
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        # Count Ham paths directly
        count = 0
        for perm in permutations(range(n)):
            if all(A[perm[i]][perm[i+1]] for i in range(n-1)):
                count += 1
        total_H += count

    # Expected: each edge is forward with prob 1/2 (for non-backbone)
    # Only backbone-only perms contribute to the average
    # Average H = sum over perms of prod (prob forward)
    # For identity perm: all backbone, all forward: prob 1
    # For other perms: each non-backbone edge has prob 1/2
    avg_H = Fraction(total_H, 2**m)
    # Alternative: avg H = sum_perm prod_step p_i where p_i = 1 for backbone, 1/2 for non-backbone
    # Hmm, this is a more complex calculation

    print(f"n={n}: m={m}, sum H(T) = {total_H}, avg H = {avg_H}")

print(f"\n{'='*60}")
print("INTERPRETATION")
print("="*60)
print("""
sum_T W_T(0) counts, up to (1/2)^{n-1} scaling, the signed Hamiltonian
weight: each backbone-only path contributes ±1.

At odd n: identity and reverse both contribute +1, total = 2.
At even n: they cancel (one + one - = 0).

This is a very natural "boundary" computation — only boundary perms
(those using only the fixed backbone) survive averaging over all
tournament orientations.
""")

# Now: what about the VARIANCE of W(0) over tournaments?
print(f"\n{'='*60}")
print("VARIANCE of W(0) and W(1/2)")
print("="*60)
for n in [3, 5]:
    m = num_tiling_bits(n)
    W0_vals = []
    H_vals = []
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        W0 = compute_W0(A, n)
        W0_vals.append(W0)

        count = sum(1 for perm in permutations(range(n))
                    if all(A[perm[i]][perm[i+1]] for i in range(n-1)))
        H_vals.append(count)

    mean_W0 = sum(W0_vals) / len(W0_vals)
    var_W0 = sum((w - mean_W0)**2 for w in W0_vals) / len(W0_vals)

    mean_H = Fraction(sum(H_vals), len(H_vals))
    var_H = Fraction(sum((h - mean_H)**2 for h in H_vals), len(H_vals))

    print(f"  n={n}: E[W(0)] = {mean_W0}, Var[W(0)] = {var_W0} ({float(var_W0):.4f})")
    print(f"         E[H] = {mean_H}, Var[H] = {var_H} ({float(var_H):.4f})")
    print(f"         Var[W(0)]/Var[H] = {float(var_W0)/float(var_H):.6f}")
