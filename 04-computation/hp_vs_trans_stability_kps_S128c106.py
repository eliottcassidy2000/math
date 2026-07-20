#!/usr/bin/env python3
"""hp_vs_trans_stability_kps_S128c106.py -- kind-pasteur-2026-07-20-S128c106

QUESTION.  The repo's central invariant is hp(T) = the number of Hamiltonian paths
(Redei: always odd; hp = 1 iff T is transitive).  The Erdos-Hajnal quantity for
tournaments is tr(T) = the size of the largest TRANSITIVE subtournament.  Both
measure "distance from transitivity", from opposite ends, and as far as I can tell
nobody in this repo has related them.

Two extremes agree: hp = 1 <=> transitive <=> tr = n, and the hp-MAXIMISER is
quasi-random (Paley/circulant), which is exactly the regime where tr is SMALL.
That suggests a monotone tension.  The sharp form worth having is a STABILITY
statement:

    (S)   hp(T) small  =>  tr(T) large,

i.e. a lower bound tr >= phi(n, hp) decreasing in hp.  Stability statements of this
shape are how Erdos-Hajnal-type results are usually approached, so a clean phi
would be worth something, and a clean COUNTEREXAMPLE to monotonicity would be
worth knowing before anyone builds on the intuition.

WHAT IS COMPUTED
  (A) n = 3..6: EXHAUSTIVE over all 2^C(n,2) labelled tournaments.  For every
      attained value of hp, the minimum and maximum of tr.  This is the exact
      profile, no sampling.
  (B) n = 7,8,9: the near-transitive regime exactly, by flipping every set of
      k <= 3 arcs of the transitive tournament (C(36,3) = 7140 at n = 9).  This is
      where (S) lives, and it is exactly enumerable even though the full space is
      not.
  (C) A direct monotonicity test: is min{tr : hp(T) = h} non-increasing in h?

hp is counted by subset DP; tr by "largest 3-cycle-free subset", using the fact
that an induced subtournament is transitive iff it contains no 3-cycle.
"""
import sys
from itertools import combinations

NMAX_EXH = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def hp_count(adj, n):
    """Number of Hamiltonian paths, by DP over subsets."""
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
            outs = adj[v] & ~mask
            m = outs
            while m:
                b = m & -m
                u = b.bit_length() - 1
                dp[mask | b][u] += c
                m ^= b
    full = size - 1
    return sum(dp[full])


def cycle_masks(adj, n):
    """Masks of the 3-element subsets that induce a 3-cycle."""
    out = []
    for i, j, k in combinations(range(n), 3):
        # a 3-set is a cycle iff each vertex has out-degree 1 inside it
        d = 0
        for a, b in ((i, j), (j, k), (i, k)):
            if adj[a] >> b & 1:
                d += 1
        # out-degrees inside are (d_i,d_j,d_k); cyclic iff sum==3 and none is 2
        di = ((adj[i] >> j) & 1) + ((adj[i] >> k) & 1)
        dj = ((adj[j] >> i) & 1) + ((adj[j] >> k) & 1)
        dk = ((adj[k] >> i) & 1) + ((adj[k] >> j) & 1)
        if di == dj == dk == 1:
            out.append((1 << i) | (1 << j) | (1 << k))
    return out


def trans_size(adj, n):
    """Largest transitive subtournament = largest 3-cycle-free vertex subset."""
    cyc = cycle_masks(adj, n)
    if not cyc:
        return n
    best = 0
    for mask in range(1 << n):
        pc = bin(mask).count('1')
        if pc <= best:
            continue
        ok = True
        for c in cyc:
            if mask & c == c:
                ok = False
                break
        if ok:
            best = pc
    return best


def from_bits(bits, n):
    """Build adj from an upper-triangular orientation bitstring."""
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
print("(A) EXHAUSTIVE PROFILE  min/max tr  for each attained hp,  n = 3..%d" % NMAX_EXH)
print("=" * 78)
profiles = {}
for n in range(3, NMAX_EXH + 1):
    m = n * (n - 1) // 2
    prof = {}
    for bits in range(1 << m):
        adj = from_bits(bits, n)
        h = hp_count(adj, n)
        t = trans_size(adj, n)
        if h not in prof:
            prof[h] = [t, t]
        else:
            if t < prof[h][0]:
                prof[h][0] = t
            if t > prof[h][1]:
                prof[h][1] = t
    profiles[n] = prof
    print("  n = %d   (%d labelled tournaments)" % (n, 1 << m))
    print("     %-8s %-8s %-8s" % ("hp", "min tr", "max tr"))
    for h in sorted(prof):
        print("     %-8d %-8d %-8d" % (h, prof[h][0], prof[h][1]))
    mins = [prof[h][0] for h in sorted(prof)]
    nonincr = all(mins[i] >= mins[i + 1] for i in range(len(mins) - 1))
    print("     min-tr non-increasing in hp : %s" % nonincr)
    print()

print("=" * 78)
print("(B) THE NEAR-TRANSITIVE REGIME, exactly: flip k <= 3 arcs of the")
print("    transitive tournament,  n = 7,8,9")
print("=" * 78)
for n in (7, 8, 9):
    pairs = list(combinations(range(n), 2))
    base = [0] * n
    for i, j in pairs:
        base[i] |= 1 << j            # transitive: i -> j for i < j
    print("  n = %d   transitive baseline: hp = %d, tr = %d"
          % (n, hp_count(base, n), trans_size(base, n)))
    for k in (1, 2, 3):
        prof = {}
        for flip in combinations(range(len(pairs)), k):
            adj = list(base)
            for f in flip:
                i, j = pairs[f]
                adj[i] &= ~(1 << j)
                adj[j] |= 1 << i
            h = hp_count(adj, n)
            t = trans_size(adj, n)
            if h not in prof:
                prof[h] = [t, t]
            else:
                prof[h][0] = min(prof[h][0], t)
                prof[h][1] = max(prof[h][1], t)
        print("     k = %d flips  (%d configurations)"
              % (k, len(list(combinations(range(len(pairs)), k)))))
        for h in sorted(prof):
            print("        hp = %-6d  tr in [%d, %d]   (n - min tr = %d)"
                  % (h, prof[h][0], prof[h][1], n - prof[h][0]))
    print()

print("=" * 78)
print("(C) MONOTONICITY VERDICT")
print("=" * 78)
for n in sorted(profiles):
    prof = profiles[n]
    hs = sorted(prof)
    mins = [prof[h][0] for h in hs]
    nonincr = all(mins[i] >= mins[i + 1] for i in range(len(mins) - 1))
    viol = [(hs[i], mins[i], hs[i + 1], mins[i + 1])
            for i in range(len(mins) - 1) if mins[i] < mins[i + 1]]
    print("  n = %d : min-tr non-increasing in hp : %-5s %s"
          % (n, nonincr, ("" if nonincr else "  violations: %s" % viol[:4])))
print()
print("  If min-tr is NOT monotone in hp, the naive stability intuition")
print("  'more Hamiltonian paths => less transitivity' is FALSE as stated, and any")
print("  bound must be phrased with hp as a one-sided constraint only.")
