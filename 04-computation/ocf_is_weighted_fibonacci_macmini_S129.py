#!/usr/bin/env python3
"""
The OCF IS the 2-weighted shifted-Pascal identity, and Redei is its mod-2 shadow
                                                        (mac-mini-2026-07-20-S129)
================================================================================
Owner: "think how fibonacci arises from summing pascals triangle shifted" and
"1001 = three sixties relates to the final digit of Fibonacci repeating every 60."

The bridge, tested here:
  Fibonacci from shifted Pascal   sum_k C(m-k, k) * 1^k = F(m+1)
                                    = # independent sets in the path P_m
  The repo's OCF                  H(T) = sum over sets of vertex-disjoint odd cycles
                                         of 2^(#cycles)
                                    = the SAME independent-set sum, weight 2 not 1,
                                      on the odd-cycle intersection graph.
  So when that graph is a PATH,   H = sum_k C(m-k,k) 2^k = JACOBSTHAL J(m+1)
                                    = (2^m - (-1)^m)/3.

Consequences tested:
  * Redei's theorem is the mod-2 shadow: every nonempty collection contributes an
    EVEN 2^|S|, so H = 1 (mod 2).  Oddness of H is one line from the OCF.
  * Periodicity: Fibonacci's last digit has period 60; the repo's 2-weighted analogue
    has last-digit period 4.  Both computed, not asserted.
"""
from itertools import permutations, combinations
from math import comb

# ------------------------------------------------------------------ tournaments
def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(pairs)
    for mask in range(1 << E):
        w = [[0] * n for _ in range(n)]
        for e, (i, j) in enumerate(pairs):
            if mask >> e & 1: w[j][i] = 1
            else:             w[i][j] = 1
        yield w

def ham_paths(w, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c or not (S >> v & 1): continue
            for u in range(n):
                if S >> u & 1 or not w[v][u]: continue
                dp[S | 1 << u][u] += c
    return sum(dp[(1 << n) - 1])

def directed_cycles(w, n, odd_only=True):
    """all directed cycles as (frozenset of vertices); returns list of vertex sets
       WITH multiplicity = number of distinct directed cycles on that vertex set."""
    out = []
    for size in range(3, n + 1):
        if odd_only and size % 2 == 0: continue
        for verts in combinations(range(n), size):
            # count distinct directed Hamiltonian cycles on these vertices
            v0 = verts[0]; rest = verts[1:]
            cnt = 0
            for perm in permutations(rest):
                seq = (v0,) + perm
                if all(w[seq[i]][seq[(i + 1) % size]] for i in range(size)): cnt += 1
            for _ in range(cnt): out.append(frozenset(verts))
    return out

def ocf(w, n):
    """sum over sets of VERTEX-DISJOINT directed odd cycles of 2^(#cycles)"""
    cyc = directed_cycles(w, n)
    total = 0
    def rec(i, used, k):
        nonlocal total
        if i == len(cyc):
            total += 1 << k; return
        rec(i + 1, used, k)                      # skip cycle i
        if not (cyc[i] & used): rec(i + 1, used | cyc[i], k + 1)
    rec(0, frozenset(), 0)
    return total

print("=" * 74)
print("PART 1 -- is H(T) the 2-weighted disjoint-odd-cycle sum?  (all labelled T)")
print("=" * 74)
for n in range(3, 7):
    bad = 0; tot = 0
    for w in all_tournaments(n):
        tot += 1
        if ham_paths(w, n) != ocf(w, n): bad += 1
    print(f"  n={n}: {tot:>6} tournaments, mismatches = {bad}   "
          f"{'OCF CONFIRMED' if bad == 0 else '*** OCF FAILS ***'}")

print()
print("=" * 74)
print("PART 2 -- the shifted-Pascal identity, weight 1 (Fibonacci) vs weight 2")
print("=" * 74)
def shifted_pascal(m, wt): return sum(comb(m - k, k) * wt**k for k in range(m // 2 + 1))
F = [0, 1]
while len(F) < 30: F.append(F[-1] + F[-2])
J = [0, 1]
while len(J) < 30: J.append(J[-1] + 2 * J[-2])
print(f"{'m':>3} {'sum C(m-k,k)':>13} {'F(m+1)':>8} {'ok':>4}   "
      f"{'sum C(m-k,k)2^k':>16} {'J(m+1)':>8} {'ok':>4} {'(2^m-(-1)^m)/3':>15}")
for m in range(12):
    a, b = shifted_pascal(m, 1), shifted_pascal(m, 2)
    cf = (2**m - (-1)**m) // 3
    print(f"{m:>3} {a:>13} {F[m+1]:>8} {str(a==F[m+1]):>4}   "
          f"{b:>16} {J[m+1]:>8} {str(b==J[m+1]):>4} {cf:>15}")

print()
print("  So: independent sets in the path P_m, weight 1 -> FIBONACCI;")
print("      the same sum with weight 2 -> JACOBSTHAL = (2^m - (-1)^m)/3.")
print("      The OCF is exactly this sum on the odd-cycle intersection graph.")

print()
print("=" * 74)
print("PART 3 -- Redei as the mod-2 shadow of the OCF")
print("=" * 74)
print("  Every NONEMPTY collection contributes 2^|S| == 0 (mod 2); only S = {} gives 1.")
print("  Hence H(T) == 1 (mod 2) -- Redei's oddness is one line from the OCF.")
for n in range(3, 7):
    odd = all(ham_paths(w, n) % 2 == 1 for w in all_tournaments(n))
    print(f"  n={n}: every H(T) odd?  {odd}")

print()
print("=" * 74)
print("PART 4 -- last-digit periods: Fibonacci 60 vs the repo's 2-weighted analogue")
print("=" * 74)
def eventual_period(seq_fn, mod, guard=6000):
    """smallest (pre-period, period): s[i]==s[i+p] for all i >= pre."""
    s = [seq_fn(i) % mod for i in range(guard)]
    for p in range(1, guard // 3):
        for pre in range(0, 40):
            if all(s[i] == s[i + p] for i in range(pre, guard - p - 1)):
                return pre, p
    return None, None
def fib(i):
    a, b = 0, 1
    for _ in range(i): a, b = b, a + b
    return a
def jac(i):
    a, b = 0, 1
    for _ in range(i): a, b = b, b + 2 * a
    return a
for mod in (2, 5, 10, 1001):
    (af, pf), (aj, pj) = eventual_period(fib, mod), eventual_period(jac, mod)
    print(f"  mod {mod:>4}:  Fibonacci period = {str(pf):>5} (pre {af})   "
          f"Jacobsthal period = {str(pj):>5} (pre {aj})")
print()
print()
print("  READ THIS OFF THE TABLE, not from intuition:")
print("    Fibonacci  (weight 1): period  60 mod   10,  period 560 mod 1001")
print("    Jacobsthal (weight 2): period   4 mod   10,  period  60 mod 1001")
print("  So the 60 does NOT vanish under the tournament weighting -- it MOVES,")
print("  from modulus 10 to modulus 1001 = 7*11*13.  Reason: 60 = ord_10(F-recursion)")
print("  on one side and 60 = ord_1001(2) on the other, and the weight-2 recursion has")
print("  J(k) = (2^k - (-1)^k)/3, so its period mod q is governed by ord_q(2).")
print("  That is the honest '1001 <-> three sixties' link: same constant, different modulus,")
print("  swapped by the weight-1 -> weight-2 substitution that turns Fibonacci into the OCF.")
