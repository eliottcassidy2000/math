#!/usr/bin/env python3
"""odd_paths_even_cycles_kps_S128c114.py -- kind-pasteur-2026-07-20-S128c114

THE ODD/EVEN DUALITY THE OWNER IS POINTING AT.

Rédei: hp(T), the number of Hamiltonian PATHS, is always ODD.  The repo's OCF gives
the reason in one line: hp = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + ... , so
hp is odd because the empty independent set contributes 1 and everything else is
even.  THM-165 sharpens the even part: each directed Hamiltonian CYCLE conflicts
with every other odd cycle, so changing hc by delta changes hp by exactly 2*delta.

So paths carry the odd part and cycles carry an even part.  But there is a second,
sharper sense in which cycles are the "even" object, and it is the one that connects
to the repo's even graphs E_n:

    a Hamiltonian PATH has two vertices of degree 1  -> NOT an even subgraph
    a Hamiltonian CYCLE has every vertex of degree 2 -> IS an even subgraph

Hamiltonian cycles are literally elements of the cycle space of K_n, i.e. points of
E_n.  Hamiltonian paths are not.  That makes the natural question:

    Rédei says hp is always ODD.  Is hc always EVEN, or does it obey some other
    parity law?

If it does, that is the even-graph-side dual of Rédei, and it would explain the
coefficient 2 of THM-165 structurally rather than combinatorially.

COMPUTED HERE
  (A) hp and hc for every tournament, n = 3..6 exhaustive (control: hp always odd).
  (B) the parity of hc, and its distribution -- overall and split by strong/non-strong.
  (C) hc mod 4, in case the law is finer than mod 2.
  (D) the same at n = 7, 8 by sampling, to see whether any law is an n <= 6 artifact.
  (E) the directed-odd-cycle census (the Omega vertices), to locate where the odd/even
      split actually sits.
"""
import sys
import random
from itertools import permutations

random.seed(17)
NEX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
NSAMP = int(sys.argv[2]) if len(sys.argv) > 2 else 4000


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


def hc_count(adj, n):
    """Directed Hamiltonian cycles, each counted ONCE (rotation fixed at vertex 0)."""
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    dp[1][0] = 1
    for mask in range(size):
        if not (mask & 1):
            continue
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
    full = size - 1
    return sum(dp[full][v] for v in range(1, n) if (adj[v] >> 0) & 1)


def is_strong(adj, n):
    seen = 1
    stack = [0]
    while stack:
        v = stack.pop()
        m = adj[v] & ~seen
        while m:
            b = m & -m
            seen |= b
            stack.append(b.bit_length() - 1)
            m ^= b
    if seen != (1 << n) - 1:
        return False
    radj = [0] * n
    for a in range(n):
        for b2 in range(n):
            if a != b2 and (adj[a] >> b2) & 1:
                radj[b2] |= 1 << a
    seen = 1
    stack = [0]
    while stack:
        v = stack.pop()
        m = radj[v] & ~seen
        while m:
            b = m & -m
            seen |= b
            stack.append(b.bit_length() - 1)
            m ^= b
    return seen == (1 << n) - 1


def odd_cycle_census(adj, n):
    """Number of directed cycles of each ODD length (the vertices of Omega(T))."""
    out = {}
    for k in range(3, n + 1, 2):
        tot = 0
        for sub in range(1 << n):
            if bin(sub).count('1') != k:
                continue
            vs = [i for i in range(n) if sub >> i & 1]
            base = vs[0]
            for perm in permutations(vs[1:]):
                seq = [base] + list(perm)
                if all((adj[seq[i]] >> seq[(i + 1) % k]) & 1 for i in range(k)):
                    tot += 1
        out[k] = tot
    return out


print("=" * 78)
print("(A)+(B)  hp odd (Rédei, control) and the PARITY OF hc,  n = 3..%d exhaustive" % NEX)
print("=" * 78)
for n in range(3, NEX + 1):
    m = n * (n - 1) // 2
    hp_odd = True
    par = {0: 0, 1: 0}
    par_strong = {0: 0, 1: 0}
    par_weak = {0: 0, 1: 0}
    mod4 = {}
    for bits in range(1 << m):
        adj = from_bits(bits, n)
        hpv = hp_count(adj, n)
        hcv = hc_count(adj, n)
        if hpv % 2 == 0:
            hp_odd = False
        par[hcv % 2] += 1
        mod4[hcv % 4] = mod4.get(hcv % 4, 0) + 1
        if is_strong(adj, n):
            par_strong[hcv % 2] += 1
        else:
            par_weak[hcv % 2] += 1
    tot = 1 << m
    print("  n = %d  (%d tournaments)" % (n, tot))
    print("     hp always ODD (Rédei control) : %s" % hp_odd)
    print("     hc parity   : even %d (%.1f%%),  odd %d (%.1f%%)"
          % (par[0], 100.0 * par[0] / tot, par[1], 100.0 * par[1] / tot))
    print("     hc mod 4    : %s" % dict(sorted(mod4.items())))
    print("     among STRONG    : even %d, odd %d" % (par_strong[0], par_strong[1]))
    print("     among NON-strong: even %d, odd %d" % (par_weak[0], par_weak[1]))

print()
print("=" * 78)
print("(D)  SAMPLING at n = 7, 8 -- is any parity law an n <= 6 artifact?")
print("=" * 78)
for n in (7, 8):
    m = n * (n - 1) // 2
    par = {0: 0, 1: 0}
    hp_odd = True
    for _ in range(NSAMP):
        adj = from_bits(random.getrandbits(m), n)
        if hp_count(adj, n) % 2 == 0:
            hp_odd = False
        par[hc_count(adj, n) % 2] += 1
    print("  n = %d  (%d samples) : hp always odd %s ; hc even %d, odd %d"
          % (n, NSAMP, hp_odd, par[0], par[1]))

print()
print("=" * 78)
print("(E)  WHERE THE ODD/EVEN SPLIT SITS: the directed-odd-cycle census")
print("=" * 78)
print("  Omega(T)'s vertices are the directed ODD cycles.  hp = I(Omega, 2), so")
print("  hp = 1 + 2*(number of odd cycles) + 4*(...) -- the 1 is the empty independent")
print("  set and is the whole reason Rédei's count is odd.")
for n in (4, 5, 6):
    m = n * (n - 1) // 2
    samples = [from_bits(random.getrandbits(m), n) for _ in range(3)]
    for adj in samples:
        cen = odd_cycle_census(adj, n)
        hpv = hp_count(adj, n)
        a1 = sum(cen.values())
        print("     n=%d  odd-cycle census %s  sum a_1 = %-4d  hp = %-5d  hp - 1 - 2a_1 = %d"
              % (n, cen, a1, hpv, hpv - 1 - 2 * a1))
print()
print("  The residual hp - 1 - 2*a_1 is 4*alpha_2 + 8*alpha_3 + ... , i.e. the")
print("  higher independent sets of Omega.  Divisible by 4 in every row above is the")
print("  check that the OCF grading is doing what it should.")
