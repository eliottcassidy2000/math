#!/usr/bin/env python3
"""dsat_is_a003141_halfcube_kps_S128c108.py -- kind-pasteur-2026-07-20-S128c108

TWO CHECKS, both of which correct something.

(A) IS d_sat A NEW INVARIANT?  mac-mini-S126's THM-1390 presents
    d_sat(n) = 2,3,4,7 -- the least d at which the waggly truncation G^(<=d)
    saturates to the complete graph -- as a NEW metagraph invariant, and hands off
    "anyone extending it should compute n=8 before conjecturing".  But
    07-reflections/diameter-is-feedback-arc-set.md (opus-S306, four months earlier)
    already records  diam(G_n) = max_T min-FAS(T) = A003141(n), with growth ~n^2/4,
    and OPEN-QUESTIONS lists it as RESOLVED.  If those are the same object then
    d_sat is not new, n=8 needs no computation, and the asymptotics are known.
    Checked here by computing max_T min-FAS(T) directly.

(B) IS MY OWN HALF-CUBE CLAIM RIGHT?  I told the owner that the half-cube (even-
    weight vectors, Hamming distance 2) "is literally the d=2 waggly layer".  Two
    corrections, the second load-bearing:
      1. distance-2 moves PRESERVE weight parity, so the d=2 layer on all 2^m
         tilings is DISCONNECTED into two components, each a copy of the halved
         cube -- TWO copies, not one;
      2. the even/odd split is not S_n-invariant, so it does NOT descend to the
         metagraph.  Tiling weight = number of chords of the corresponding even
         subgraph, which depends on the choice of base path, so weight parity is
         not a class function.  Verified here by exhibiting two tilings in the SAME
         iso class with OPPOSITE weight parity.
    If (2) holds, the halved-cube structure is destroyed by exactly the quotient the
    project is built on, and my suggestion was worth less than I said.
"""
import sys
from itertools import combinations, permutations

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def fas(adj, n):
    """Minimum feedback arc set size, by subset DP over linear orders."""
    adj_in = [0] * n
    for v in range(n):
        for u in range(n):
            if u != v and (adj[u] >> v) & 1:
                adj_in[v] |= 1 << u
    size = 1 << n
    dp = [0] * size
    for S in range(1, size):
        best = -1
        m = S
        while m:
            b = m & -m
            v = b.bit_length() - 1
            m ^= b
            val = dp[S ^ b] + bin(adj_in[v] & (S ^ b)).count('1')
            if val > best:
                best = val
        dp[S] = best
    return n * (n - 1) // 2 - dp[size - 1]


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


print("=" * 78)
print("(A) max_T min-FAS(T)  -- is this mac-mini's 'new invariant' d_sat?")
print("=" * 78)
seq = []
for n in range(3, NMAX + 1):
    m = n * (n - 1) // 2
    best = 0
    arg = None
    for bits in range(1 << m):
        f = fas(from_bits(bits, n), n)
        if f > best:
            best = f
            arg = bits
    seq.append(best)
    print("  n = %d : max over all %d labelled tournaments of min-FAS = %d"
          % (n, 1 << m, best))
print()
print("  sequence n=3..%d : %s" % (NMAX, seq))
print("  opus-S306 (07-reflections/diameter-is-feedback-arc-set.md) records")
print("     diam(G_n) = max_T min-FAS(T) = A003141(n), growth ~n^2/4,")
print("     and OPEN-QUESTIONS lists it RESOLVED.")
print("  mac-mini-S126 THM-1390 reports d_sat(n) = 2,3,4,7 for n=4..7 as NEW.")
print("  merged-metagraph diameters already in canon (waggly_completeness_s301.out)")
print("     are 1,3,4 at n=4,5,6.")
print("  -> compare the n=4,5,6 entries of the computed sequence above with 1,3,4.")

print()
print("=" * 78)
print("(B1) does the d=2 layer split by weight parity?  (two halved cubes, not one)")
print("=" * 78)
print("  A distance-2 move flips exactly 2 bits, so it PRESERVES Hamming weight")
print("  parity.  Hence the d=2 graph on 2^m tilings has (at least) two components.")
for m in (4, 6, 10):
    ev = sum(1 for b in range(1 << m) if bin(b).count('1') % 2 == 0)
    od = (1 << m) - ev
    deg = m * (m - 1) // 2
    print("  m = %-3d : 2^m = %-6d  even %-5d  odd %-5d  each C(m,2)-regular, deg = %d"
          % (m, 1 << m, ev, od, deg))
    print("            d=2 edges total = 2 * (%d * %d / 2) = %d"
          % (ev, deg, ev * deg))
print("  -> so 'the half-cube IS the d=2 layer' is off by a factor of two:")
print("     the d=2 layer is a DISJOINT UNION of two halved cubes.")

print()
print("=" * 78)
print("(B2) is weight parity a CLASS function?  (does the split survive the quotient)")
print("=" * 78)
# tiling model: base path n-1 -> n-2 -> ... -> 0; tiles = arcs (x,y) with x-y >= 2.
for n in (4, 5):
    tiles = [(x, y) for y in range(n) for x in range(n) if x - y >= 2]
    m = len(tiles)

    def build(bits):
        adj = [0] * n
        for k in range(n - 1, 0, -1):
            adj[k] |= 1 << (k - 1)          # base path arcs k -> k-1
        for i, (x, y) in enumerate(tiles):
            if bits >> i & 1:
                adj[x] |= 1 << y
            else:
                adj[y] |= 1 << x
        return adj

    def canon(adj):
        best = None
        for p in permutations(range(n)):
            e = frozenset((p[a], p[b]) for a in range(n) for b in range(n)
                          if a != b and (adj[a] >> b) & 1)
            key = tuple(sorted(e))
            if best is None or key < best:
                best = key
        return best

    cls = {}
    for bits in range(1 << m):
        c = canon(build(bits))
        cls.setdefault(c, []).append(bits)
    mixed = 0
    witness = None
    for c, bl in cls.items():
        par = {bin(b).count('1') % 2 for b in bl}
        if len(par) > 1:
            mixed += 1
            if witness is None:
                e = [b for b in bl if bin(b).count('1') % 2 == 0][0]
                o = [b for b in bl if bin(b).count('1') % 2 == 1][0]
                witness = (e, o)
    print("  n = %d (m = %d tiles, %d iso classes):" % (n, m, len(cls)))
    print("     classes containing BOTH parities : %d of %d" % (mixed, len(cls)))
    if witness:
        e, o = witness
        print("     witness: tilings %s (weight %d, EVEN) and %s (weight %d, ODD)"
              % (bin(e), bin(e).count('1'), bin(o), bin(o).count('1')))
        print("              are ISOMORPHIC as tournaments -- same class.")
print()
print("  -> if any class is mixed, weight parity is NOT a class function, the")
print("     even/odd bipartition does NOT descend to the metagraph, and the halved-")
print("     cube structure is destroyed by exactly the S_n quotient the project is")
print("     built on.  My 'half-cube = d=2 layer' handoff needs both corrections.")
