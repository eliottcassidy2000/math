#!/usr/bin/env python3
"""
W(0) = 0 when t3 is odd — does this generalize beyond n=5?

At n=5: W(0) = 1 - t3 + 2*t5 = 0 for ALL tournaments with odd t3.
This is because odd t3 forces t3 = 1 + 2*t5.

At n=7: W(0) = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc
Does this vanish (or have a pattern) for odd vs even t3?

Actually, let me think about this differently.
W(0) = sum_perm prod_{i=0}^{n-2} s_i where s_i in {+1/2, -1/2}.
Each s_i is +1/2 if forward edge, -1/2 if backward.

W(0) = (1/2)^{n-1} * sum_perm (-1)^{backward_edges}
     = (1/2)^{n-1} * sum_perm (-1)^{n-1-k} where k = forward edges
     = (1/2)^{n-1} * (-1)^{n-1} * sum_perm (-1)^{-k}
     = ((-1)/2)^{n-1} * sum_perm (-1)^k

Wait: W(0) = sum_perm prod(s_i) where s_i = +1/2 or -1/2.
Let k = #forward edges in perm. Then prod = (1/2)^k * (-1/2)^{n-1-k}
= (1/2)^{n-1} * (-1)^{n-1-k} = (1/2)^{n-1} * (-1)^{n-1} * (-1)^{-k}
= (-1/2)^{n-1} * (-1)^k ... hmm, let me be more careful.

prod_{i=0}^{n-2} s_i = prod (1/2 if fwd, -1/2 if bwd)
= (1/2)^{n-1} * (-1)^{(n-1-k)} where k = forward edges
= (1/2)^{n-1} * (-1)^{n-1-k}

So W(0) = (1/2)^{n-1} * sum_perm (-1)^{n-1-k(P)}
        = (1/2)^{n-1} * (-1)^{n-1} * sum_perm (-1)^{k(P)}
        ... wait, that doesn't factor nicely since k(P) varies.

Let me just write:
W(0) = (1/2)^{n-1} * sum_perm (-1)^{bwd(P)}
where bwd(P) = n-1-k(P) = number of backward edges.

This is: W(0) = (1/2)^{n-1} * sum_k a_k * (-1)^{n-1-k}
where a_k = #{perms with k forward edges}.

By palindromy a_k = a_{n-1-k}:
W(0) = (1/2)^{n-1} * sum_k a_k * (-1)^{n-1-k}

Let me compute this for n=5:
n-1=4. sum_k a_k * (-1)^{4-k} = a_0 - a_1 + a_2 - a_3 + a_4
With palindromy: = a_0 - a_1 + a_2 - a_1 + a_0 = 2a_0 - 2a_1 + a_2

For transitive n=5: a = [1, 26, 66, 26, 1]
W(0) = (1/16) * (1 - 26 + 66 - 26 + 1) = (1/16)*16 = 1. ✓ (t3=0)

For t3=3 tournament: a = [9, 30, 42, 30, 9]
W(0) = (1/16) * (9 - 30 + 42 - 30 + 9) = (1/16)*0 = 0. ✓ (t3=3, odd)

Interesting! The alternating sum of Eulerian numbers is 0 when t3 is odd.

Let me check at n=6 and n=7.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict

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

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def compute_forward_dist(A, n):
    """Distribution of forward edge counts over all permutations."""
    dist = defaultdict(int)
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        dist[fwd] += 1
    return dict(dist)

import random
random.seed(42)

for n in [5, 6]:
    m = num_tiling_bits(n)
    total = 2**m
    test_set = range(total) if total <= 1024 else random.sample(range(total), 200)

    print(f"\n{'='*60}")
    print(f"n={n}: W(0) vs t3 parity")
    print(f"{'='*60}")

    # Compute W(0) for each tournament
    t3_parity_w0 = defaultdict(list)
    for bits in test_set:
        A = tournament_from_tiling(n, bits)
        t3 = count_3cycles(A, n)

        # W(0) = (1/2)^{n-1} * sum_k a_k * (-1)^{n-1-k}
        dist = compute_forward_dist(A, n)
        alt_sum = sum(dist.get(k, 0) * ((-1)**(n-1-k)) for k in range(n))
        W0 = Fraction(alt_sum, 2**(n-1))

        t3_parity_w0[t3 % 2].append((t3, W0))

    for parity in [0, 1]:
        entries = t3_parity_w0[parity]
        W0_vals = sorted(set(float(w0) for _, w0 in entries))
        zero_count = sum(1 for _, w0 in entries if w0 == 0)
        print(f"  t3 {'odd' if parity == 1 else 'even'}: {len(entries)} tournaments")
        print(f"    W(0) values: {W0_vals[:10]}{'...' if len(W0_vals)>10 else ''}")
        print(f"    W(0)=0 count: {zero_count}/{len(entries)}")

# n=6: odd powers only, so W(0) = 0 for ALL tournaments!
# Because W has only odd powers of r, so W(0) = 0 trivially.
print(f"\n{'='*60}")
print("KEY INSIGHT: At even n, W has only ODD r-powers, so W(0) = 0 ALWAYS.")
print("At odd n, W has only EVEN r-powers, so W(0) is the constant term.")
print("="*60)

# At n=7: W(0) = constant term (non-zero in general)
# Let me check a sample at n=7
n = 7
m = num_tiling_bits(n)
print(f"\nn={n}: W(0) vs t3 parity (sample)")

for bits in [0, 1, 42, 100, 1000, 2**m-1]:
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    dist = compute_forward_dist(A, n)
    alt_sum = sum(dist.get(k, 0) * ((-1)**(n-1-k)) for k in range(n))
    W0 = Fraction(alt_sum, 2**(n-1))
    print(f"  bits={bits}: t3={t3}, t3%2={t3%2}, W(0)={W0} ({float(W0):.4f})")
    print(f"    Forward dist: {dict(sorted(dist.items()))}")
    print(f"    Alt sum = {alt_sum}")

# At n=7 (odd): W(0) = alternating sum / 64
# Does the alternating sum have a specific parity related to t3?
print(f"\nn=7: Checking alt_sum parity vs t3 parity (sample of 200)")
alt_sum_parities = defaultdict(list)
for _ in range(200):
    bits = random.randint(0, 2**m - 1)
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    dist = compute_forward_dist(A, n)
    alt_sum = sum(dist.get(k, 0) * ((-1)**(n-1-k)) for k in range(n))
    alt_sum_parities[t3 % 2].append(alt_sum)

for parity in [0, 1]:
    entries = alt_sum_parities[parity]
    vals = sorted(set(entries))
    print(f"  t3 {'odd' if parity == 1 else 'even'}: {len(entries)} samples")
    print(f"    Alt sum values: {vals[:15]}...")
    parities = set(v % 2 for v in entries)
    print(f"    Alt sum parities: {parities}")
    # Check mod 4
    mod4s = set(v % 4 for v in entries)
    print(f"    Alt sum mod 4: {mod4s}")
