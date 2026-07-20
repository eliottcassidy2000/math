#!/usr/bin/env python3
"""ocf_low_n_exact_kps_S128c114.py -- kind-pasteur-2026-07-20-S128c114

FOLLOW-UP: the positive that fell out of the failed odd/even duality.

hp = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ... where alpha_k counts
INDEPENDENT SETS of size k in the odd-cycle conflict graph -- i.e. sets of k PAIRWISE
VERTEX-DISJOINT directed odd cycles.

Two vertex-disjoint odd cycles need at least 3 + 3 = 6 vertices.  So

    alpha_k = 0 for all k >= 2 whenever n <= 5,

not as an observation but as a counting fact, and therefore

    hp(T) = 1 + 2 * c_odd(T)   EXACTLY, for every tournament on n <= 5 vertices,

where c_odd is the total number of directed odd cycles.  At n = 6 the bound is met
(3+3) and the identity must break.  Checked exhaustively at n <= 5 and refuted at
n = 6 with an explicit witness, which is the honest way round: the theorem is proved
by the vertex count and the computation only confirms the break lands where predicted.
"""
import sys
from itertools import permutations

NEX = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def from_bits(bits, n):
    adj = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            idx += 1
    return adj


def hp_count(adj, n):
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            m = adj[v] & ~mask
            while m:
                b = m & -m
                dp[mask | b][b.bit_length() - 1] += c
                m ^= b
    return sum(dp[size - 1])


def odd_cycles(adj, n):
    """List of vertex-sets of directed ODD cycles (each cycle once)."""
    out = []
    for k in range(3, n + 1, 2):
        for sub in range(1 << n):
            if bin(sub).count('1') != k:
                continue
            vs = [i for i in range(n) if sub >> i & 1]
            base = vs[0]
            for perm in permutations(vs[1:]):
                seq = [base] + list(perm)
                if all((adj[seq[i]] >> seq[(i + 1) % k]) & 1 for i in range(k)):
                    out.append(sub)
    return out


print("=" * 78)
print("hp = 1 + 2*c_odd  EXACTLY for n <= 5  (because two disjoint odd cycles need 6)")
print("=" * 78)
for n in range(3, NEX + 1):
    m = n * (n - 1) // 2
    bad = 0
    tot = 0
    for bits in range(1 << m):
        adj = from_bits(bits, n)
        tot += 1
        if hp_count(adj, n) != 1 + 2 * len(odd_cycles(adj, n)):
            bad += 1
    print("  n = %d : %d tournaments, violations of hp = 1 + 2*c_odd : %d  -> %s"
          % (n, tot, bad, "HOLDS" if bad == 0 else "FAILS"))

print()
print("=" * 78)
print("AND IT MUST BREAK AT n = 6 -- explicit witness")
print("=" * 78)
n = 6
m = n * (n - 1) // 2
found = None
import random
random.seed(2)
for _ in range(4000):
    bits = random.getrandbits(m)
    adj = from_bits(bits, n)
    cyc = odd_cycles(adj, n)
    hpv = hp_count(adj, n)
    if hpv != 1 + 2 * len(cyc):
        # locate a disjoint pair, the reason for the break
        pair = None
        for i in range(len(cyc)):
            for j in range(i + 1, len(cyc)):
                if cyc[i] & cyc[j] == 0:
                    pair = (cyc[i], cyc[j])
                    break
            if pair:
                break
        found = (bits, hpv, len(cyc), pair)
        break
if found:
    bits, hpv, c, pair = found
    print("  witness bits = %d" % bits)
    print("     hp = %d,  c_odd = %d,  1 + 2*c_odd = %d,  residual = %d"
          % (hpv, c, 1 + 2 * c, hpv - 1 - 2 * c))
    print("     residual divisible by 4 : %s" % ((hpv - 1 - 2 * c) % 4 == 0))
    if pair:
        a = [i for i in range(n) if pair[0] >> i & 1]
        b = [i for i in range(n) if pair[1] >> i & 1]
        print("     a vertex-DISJOINT pair of odd cycles: %s and %s" % (a, b))
        print("     -> alpha_2 >= 1, which is exactly what n = 6 first permits (3+3).")
else:
    print("  no break found in the sample -- unexpected; the identity should fail at n=6")
print()
print("  So the OCF's higher terms switch on at precisely n = 6, and the reason is a")
print("  vertex count, not an arithmetic accident: 3 + 3 = 6 is the first room for two")
print("  vertex-disjoint odd cycles.")
