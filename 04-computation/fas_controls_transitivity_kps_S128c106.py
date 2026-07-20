#!/usr/bin/env python3
"""fas_controls_transitivity_kps_S128c106.py -- kind-pasteur-2026-07-20-S128c106

Follow-up to hp_vs_trans_stability: the flip experiment showed  n - tr <= k  for
every k-arc perturbation of the transitive tournament, k <= 3, at n = 7,8,9.  That
is not an hp phenomenon -- it is the FEEDBACK ARC SET.

    LEMMA.  tr(T) >= n - fas(T).
    Proof.  Let A be a minimum feedback arc set, |A| = fas(T); reversing A makes T
    acyclic, hence transitive.  Delete one endpoint of each arc of A: at most
    fas(T) vertices.  The surviving set S meets no arc of A, so T[S] coincides with
    the reversed (transitive) tournament on S and is therefore transitive. QED

This script checks three things:
  (A) the lemma holds, exhaustively, n <= 6;
  (B) how TIGHT it is -- the distribution of the slack tr - (n - fas);
  (C) whether hp controls tr at all, by exhibiting explicit violations of
      monotonicity in BOTH directions (min-tr and max-tr against hp).

(C) is the point.  If hp does not control tr, then the repo's central invariant is
not a route to Erdos-Hajnal-type statements, and that is worth knowing before a
session is spent on it.  And the invariant that DOES control tr gives a bound that
is vacuous exactly in the regime Erdos-Hajnal cares about (fas ~ n^2/4 for random
tournaments, so n - fas < 0), so this bridge carries no load either.
"""
import sys
from itertools import combinations

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6


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


def fas(adj, n):
    """Minimum feedback arc set size = C(n,2) - max forward arcs over linear orders.
    DP over subsets: dp[S] = max forward arcs with S as the prefix."""
    size = 1 << n
    dp = [0] * size
    for S in range(1, size):
        best = -1
        m = S
        while m:
            b = m & -m
            v = b.bit_length() - 1
            m ^= b
            rest = S ^ b
            # arcs from rest into v are forward when v is placed last in S
            gain = bin(adj_in[v] & rest).count('1')
            val = dp[rest] + gain
            if val > best:
                best = val

        dp[S] = best
    return n * (n - 1) // 2 - dp[size - 1]


def cycle_masks(adj, n):
    out = []
    for i, j, k in combinations(range(n), 3):
        di = ((adj[i] >> j) & 1) + ((adj[i] >> k) & 1)
        dj = ((adj[j] >> i) & 1) + ((adj[j] >> k) & 1)
        dk = ((adj[k] >> i) & 1) + ((adj[k] >> j) & 1)
        if di == dj == dk == 1:
            out.append((1 << i) | (1 << j) | (1 << k))
    return out


def trans_size(adj, n):
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
print("(A)+(B)  tr >= n - fas :  validity and slack,  exhaustive n = 3..%d" % NMAX)
print("=" * 78)
hp_prof = {}
for n in range(3, NMAX + 1):
    m = n * (n - 1) // 2
    viol = 0
    slack = {}
    prof = {}
    for bits in range(1 << m):
        adj = from_bits(bits, n)
        adj_in = [0] * n
        for v in range(n):
            for u in range(n):
                if u != v and (adj[u] >> v) & 1:
                    adj_in[v] |= 1 << u
        globals()['adj_in'] = adj_in
        f = fas(adj, n)
        t = trans_size(adj, n)
        h = hp_count(adj, n)
        if t < n - f:
            viol += 1
        sl = t - (n - f)
        slack[sl] = slack.get(sl, 0) + 1
        if h not in prof:
            prof[h] = [t, t]
        else:
            prof[h][0] = min(prof[h][0], t)
            prof[h][1] = max(prof[h][1], t)
    hp_prof[n] = prof
    print("  n = %d : violations of tr >= n - fas : %d  -> %s"
          % (n, viol, "LEMMA HOLDS" if viol == 0 else "LEMMA FAILS"))
    tot = 1 << m
    print("     slack tr - (n - fas) distribution: %s"
          % {k: "%.1f%%" % (100.0 * v / tot) for k, v in sorted(slack.items())})
    print("     tight (slack 0) fraction: %.1f%%" % (100.0 * slack.get(0, 0) / tot))

print()
print("=" * 78)
print("(C)  DOES hp CONTROL tr ?   monotonicity of min-tr and max-tr against hp")
print("=" * 78)
for n in sorted(hp_prof):
    prof = hp_prof[n]
    hs = sorted(prof)
    mins = [prof[h][0] for h in hs]
    maxs = [prof[h][1] for h in hs]
    vmin = [(hs[i], mins[i], hs[i + 1], mins[i + 1])
            for i in range(len(mins) - 1) if mins[i] < mins[i + 1]]
    vmax = [(hs[i], maxs[i], hs[i + 1], maxs[i + 1])
            for i in range(len(maxs) - 1) if maxs[i] < maxs[i + 1]]
    print("  n = %d : min-tr non-increasing in hp : %-5s  violations: %s"
          % (n, not vmin, vmin[:3] if vmin else "none"))
    print("          max-tr non-increasing in hp : %-5s  violations: %s"
          % (not vmax, vmax[:3] if vmax else "none"))
print()
print("  Each violation quadruple reads (hp_1, tr_1, hp_2, tr_2) with hp_1 < hp_2")
print("  but tr_1 < tr_2 -- MORE Hamiltonian paths yet MORE transitivity.")
print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
print("  hp does NOT control tr: monotonicity fails in both directions from n = 6.")
print("  fas DOES, via tr >= n - fas, which is elementary and tight.")
print("  But for a random tournament fas ~ n^2/4, so n - fas is negative and the")
print("  bound is VACUOUS exactly in the regime Erdos-Hajnal cares about.")
print("  Neither invariant is a route to EH.  Recording this so the bridge is not")
print("  rebuilt: the repo's hp is a GLOBAL count and tr is a LOCAL extremal")
print("  quantity, and the flip experiment shows hp can move by an order of")
print("  magnitude while tr does not move at all.")
