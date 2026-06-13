#!/usr/bin/env python3
"""
c3_constant_and_qr_parity.py -- kind-pasteur-2026-03-13-S60

Two questions:
1. WHY is c3=55 constant for all orientations at p=11?
2. What distinguishes the two H-classes within #QR ∈ {2,3}?

For Q1: A circulant tournament on Z_p has c3 = sum over 3-element
orbits of n(orbit). Since all 3-sets form orbits of size p under
Z_p translation, there are C(p,3)/p orbits (when p prime).
For each orbit, n(orbit) = # directed 3-cycles on the canonical
representative. Is this always 2 for tournaments?

Actually: every 3-set in a tournament supports exactly 0 or 2
directed 3-cycles (transitive vs cyclic triple). A REGULAR
tournament has c3 = p(p-1)(p-3)/24 (classical formula).
But our tournaments ARE regular (all circulant T on Z_p are regular
since S has size (p-1)/2 for each vertex).

For Q2: The split of {2,3}-QR into two H-classes must come from
some FINER invariant of the orientation, beyond just counting
how many QR elements are in S.
"""

from itertools import combinations
from collections import defaultdict

p = 11
m = (p - 1) // 2
S_qr_set = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

print(f"p={p}, QR={sorted(S_qr_set)}")
print(f"Pairs (QR, NQR):")
pairs = []
for j in range(1, m+1):
    qr_j = j if j in S_qr_set else p-j
    nqr_j = p-j if j in S_qr_set else j
    pairs.append((qr_j, nqr_j, j))
    print(f"  j={j}: pair ({j}, {p-j}), QR={qr_j}, NQR={nqr_j}")

# Q1: c3 for regular tournaments
# Classical: For ANY regular tournament on n=2k+1 vertices,
# c3 = n(n-1)(n-3)/24 ... wait, let me derive.
# Each vertex has out-degree k = (n-1)/2.
# A 3-cycle (a->b->c->a) exists iff a->b, b->c, c->a.
# Total directed 3-cycles = (1/3) * sum_v d+(v) * "how many of out-neighbors
#   of v have edges closing cycles"
# Actually: c3(directed) = C(n,3) * 2 * P(cyclic) where P depends on the tournament.
# For regular: c3 = C(n,3) * (n^2 - 4n + 3) / (n choose 3) ... no.

# Simpler: For regular tournament on n vertices with d+ = d- = (n-1)/2 = k:
# Number of 3-sets that are CYCLIC (support 2 directed cycles) =
#   C(n,3) - number of transitive triples
# Transitive triples = sum_v C(d+(v), 2) = n * C(k, 2)
# So cyclic triples = C(n,3) - n*C(k,2)
# Directed 3-cycles = 2 * cyclic triples

n = p
k = (n-1) // 2
from math import comb
transitive = n * comb(k, 2)
cyclic = comb(n, 3) - transitive
c3_directed = 2 * cyclic
c3_vsets = cyclic

print(f"\nQ1: WHY c3=55 IS CONSTANT")
print(f"  n={n}, k={k}")
print(f"  C(n,3) = {comb(n,3)}")
print(f"  Transitive triples = n*C(k,2) = {n}*{comb(k,2)} = {transitive}")
print(f"  Cyclic triples = {cyclic}")
print(f"  Directed 3-cycles = 2 * {cyclic} = {c3_directed}")
print(f"  This matches c3=55 vertex sets (c3_directed = 2*55 = 110)")

# VERIFIED: c3 depends ONLY on the degree sequence, which is constant
# (regular with d+=5) for ALL circulant tournaments on Z_11.
# So c3 = C(11,3) - 11*C(5,2) = 165 - 110 = 55 cyclic 3-sets,
# giving 110 directed 3-cycles.

# General formula: For any regular tournament on n = 2k+1 vertices:
# c3 = C(n,3) - n*C(k,2) cyclic triples
# = n(n-1)(n-2)/6 - n*k(k-1)/2
# = n * [(n-1)(n-2)/6 - k(k-1)/2]
# = n * [(n-1)(n-2) - 3k(k-1)] / 6
# With k = (n-1)/2: 3k(k-1) = 3(n-1)(n-3)/4
# So: [(n-1)(n-2) - 3(n-1)(n-3)/4] / 6
# = (n-1)/6 * [(n-2) - 3(n-3)/4]
# = (n-1)/6 * [(4n-8-3n+9)/4]
# = (n-1)(n+1) / 24
# So c3 (cyclic triples) = n(n-1)(n+1)/24 for regular tournaments!
c3_formula = n * (n-1) * (n+1) // 24
print(f"\n  Formula: c3_cyclic = n(n-1)(n+1)/24 = {c3_formula}")
print(f"  Match: {c3_formula == cyclic}")

# At p=7: c3 = 7*6*8/24 = 14. ✓ (known)
# At p=11: c3 = 11*10*12/24 = 55. ✓
# At p=13: c3 = 13*12*14/24 = 91
# At p=17: c3 = 17*16*18/24 = 204

for pp in [3, 5, 7, 11, 13, 17, 19, 23]:
    print(f"  p={pp}: c3 = {pp*(pp-1)*(pp+1)//24}")

print(f"\n{'='*70}")
print(f"Q2: WHAT DISTINGUISHES THE TWO #QR=2,3 H-CLASSES?")
print(f"{'='*70}")

# The H-classes at p=11:
# H=93467: bits = [4,5,7,11,15,16,20,24,26,27]  (#QR ∈ {2,3})
# H=92411: bits = [1,8,9,12,14,17,19,22,23,30]  (#QR ∈ {2,3})

# Let's look at which SPECIFIC pairs are flipped
# Pairs: j=1: (1,10), j=2: (2,9), j=3: (3,8), j=4: (4,7), j=5: (5,6)
# QR element in pair j: 1,9,3,4,5

class1 = [4,5,7,11,15,16,20,24,26,27]  # H=93467
class2 = [1,8,9,12,14,17,19,22,23,30]  # H=92411

print(f"\n  H=93467 class (10 members):")
for bits in class1:
    choices = []
    for j in range(m):
        if bits & (1 << j):
            choices.append(j+1)
        else:
            choices.append(p-(j+1))
    n_qr = sum(1 for c in choices if c in S_qr_set)
    which_qr = [j+1 for j in range(m) if (choices[j] in S_qr_set)]
    print(f"    bits={bits:>2} (0b{bits:05b}): S={choices}, #QR={n_qr}, QR_pairs={which_qr}")

print(f"\n  H=92411 class (10 members):")
for bits in class2:
    choices = []
    for j in range(m):
        if bits & (1 << j):
            choices.append(j+1)
        else:
            choices.append(p-(j+1))
    n_qr = sum(1 for c in choices if c in S_qr_set)
    which_qr = [j+1 for j in range(m) if (choices[j] in S_qr_set)]
    print(f"    bits={bits:>2} (0b{bits:05b}): S={choices}, #QR={n_qr}, QR_pairs={which_qr}")

# Let's look at the QR_pairs pattern more carefully
# Which SUBSET of the 5 pairs has the QR element chosen?
# Pair j=1 has QR=1, j=2 has QR=9, j=3 has QR=3, j=4 has QR=4, j=5 has QR=5

print(f"\n  Reformulation: which pairs pick QR element?")
print(f"  Pair labels: j=1(1/10), j=2(9/2), j=3(3/8), j=4(4/7), j=5(5/6)")

def qr_pattern(bits):
    """Returns which pairs (1-indexed) pick the QR element."""
    pattern = set()
    for j in range(m):
        if bits & (1 << j):
            elem = j+1
        else:
            elem = p-(j+1)
        if elem in S_qr_set:
            pattern.add(j+1)
    return frozenset(pattern)

print(f"\n  H=93467 QR patterns:")
for bits in class1:
    pat = qr_pattern(bits)
    n_qr = len(pat)
    print(f"    bits={bits:>2}: QR@pairs={sorted(pat)} (#={n_qr})")

print(f"\n  H=92411 QR patterns:")
for bits in class2:
    pat = qr_pattern(bits)
    n_qr = len(pat)
    print(f"    bits={bits:>2}: QR@pairs={sorted(pat)} (#={n_qr})")

# Let's look at it from the PRODUCT perspective
# At p=11: QR = {1,3,4,5,9}. The products s*s' mod 11 matter.
# Maybe what distinguishes the classes is the product structure?

print(f"\n  Product analysis: sum of s_i * s_j for pairs in S")
for H_val, bits_list in [(95095, [29]), (93467, [4]), (93027, [0]), (92411, [1])]:
    bits = bits_list[0]
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j+1)
        else:
            S.append(p-(j+1))

    prod_sum = sum(S[i]*S[j] % p for i in range(m) for j in range(i+1, m))
    sq_sum = sum(s*s % p for s in S)
    lin_sum = sum(S)
    prod_set = set(S[i]*S[j] % p for i in range(m) for j in range(i+1, m))
    sq_set = set(s*s % p for s in S)

    print(f"  H={H_val}: S={sorted(S)}, sum={lin_sum}, sum_sq={sq_sum}, "
          f"prod_sum={prod_sum}, sq_set={sorted(sq_set)}, prod_set={sorted(prod_set)}")

# Let's try: the distinguishing invariant might be the MULTISET of
# pairwise differences d = s_i - s_j mod p in S
print(f"\n  Difference multiset analysis:")
for H_val, bits_list in [(95095, [29]), (93467, [4]), (93027, [0]), (92411, [1])]:
    bits = bits_list[0]
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j+1)
        else:
            S.append(p-(j+1))

    diffs = sorted((S[i] - S[j]) % p for i in range(m) for j in range(m) if i != j)
    diff_counter = defaultdict(int)
    for d in diffs:
        diff_counter[d] += 1

    # How many differences land in QR vs NQR?
    n_qr_diff = sum(c for d, c in diff_counter.items() if d in S_qr_set)
    n_nqr_diff = sum(c for d, c in diff_counter.items() if d not in S_qr_set and d != 0)

    print(f"  H={H_val}: #QR_diffs={n_qr_diff}, #NQR_diffs={n_nqr_diff}, "
          f"ratio={n_qr_diff/n_nqr_diff if n_nqr_diff else 'inf':.3f}")

# Check if the CLOSURE property matters
# S is "closed under QR-multiplication" if s*r ∈ S for all s in S, r in QR
print(f"\n  QR-closure analysis:")
for H_val, bits_list in [(95095, [29]), (93467, [4]), (93027, [0]), (92411, [1])]:
    bits = bits_list[0]
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j+1)
        else:
            S.append(p-(j+1))
    S_set = set(S)

    # How many s*r mod p remain in S?
    qr_closure = sum(1 for s in S for r in S_qr_set if (s*r) % p in S_set)
    nqr_closure = sum(1 for s in S for r in (set(range(1,p))-S_qr_set) if (s*r) % p in S_set)

    print(f"  H={H_val}: S={sorted(S)}, QR_closure={qr_closure}/{m*m}, "
          f"NQR_closure={nqr_closure}/{m*m}")

print("\nDONE.")
