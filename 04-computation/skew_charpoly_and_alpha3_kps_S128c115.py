#!/usr/bin/env python3
"""skew_charpoly_and_alpha3_kps_S128c115.py -- kind-pasteur-2026-07-20-S128c115

TWO THINGS.

(A) THE NAMED NEXT STEP.  Three vertex-disjoint odd cycles need 3+3+3 = 9 vertices,
    so alpha_3 = 0 for every tournament with n <= 8, and therefore

        hp(T) = 1 + 2*alpha_1 + 4*alpha_2   EXACTLY  for all n <= 8,

    extending THM-1445's n <= 5 regime by three.  Proved by the vertex count; checked
    here.

(B) THE OWNER'S POLYNOMIAL, x(x^2+7)(x^4+14x^2+17) = x^7 + 21x^5 + 115x^3 + 119x, AND
    THE ODD/EVEN <-> sin/cos LINK.

    Give a tournament its SKEW adjacency S (S_ij = +1 if i->j, -1 if j->i, 0 on the
    diagonal).  S is a real skew-symmetric matrix, so its eigenvalues are purely
    IMAGINARY and come in +-i*lambda pairs, plus a forced 0 when n is odd.  Hence

        n odd   ->  char(S) is an ODD polynomial   (the sin archetype)
        n even  ->  char(S) is an EVEN polynomial  (the cos archetype)

    and in every case char(S)(iy) is, up to a unit, a REAL-ROOTED polynomial in y --
    exactly the property that makes sin and cos the archetypes of odd and even
    functions with real zeros.  THAT is the sin/cos content of the odd/even theme: it
    is not a metaphor, it is the spectral parity of the skew adjacency.

    The owner's polynomial has degree 7 and is odd, so it is a candidate char(S) at
    n = 7.  Its x^5 coefficient is 21 = C(7,2), which is FORCED (see below), so the
    shape is right.  Compare the Paley tournament T_7, whose skew spectrum is
    {0, +-i*sqrt7 each thrice}, giving char = x(x^2+7)^3 = x^7+21x^5+147x^3+343x.
    The owner's polynomial is the same x(x^2+7)*(x^4+14x^2+c) shape with c = 17 rather
    than 49 -- so if it is realised, it is realised by a NON-Paley tournament.

    This script searches the n = 7 tournaments for it.

WHY THE x^5 COEFFICIENT IS FORCED: for skew S, tr(S) = 0 and tr(S^2) = -sum_{i!=j}
S_ij^2 = -n(n-1) = -2*C(n,2).  Newton then gives e_2 = (tr^2 - tr(S^2))/2 = C(n,2), so
the x^{n-2} coefficient of char(S) is ALWAYS C(n,2) = 21 at n = 7, for every
tournament.  It carries no information; the discriminating coefficients are x^3 and x.
"""
import sys
import random
from itertools import permutations
import numpy as np

random.seed(23)
NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 120000


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


def skew(adj, n):
    S = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i, j] = 1 if (adj[i] >> j) & 1 else -1
    return S


def charpoly(A, n):
    """Faddeev-LeVerrier, exact over the integers. Returns [1, c_{n-1}, ..., c_0]."""
    I = np.eye(n, dtype=np.int64)
    M = I.copy()
    cs = [1]
    for k in range(1, n + 1):
        AM = A.dot(M)
        c = -int(np.trace(AM)) // k
        cs.append(c)
        M = AM + c * I
    return cs


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


def odd_cycle_sets(adj, n):
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


# ------------------------------------------------------------------ (A)
print("=" * 78)
print("(A) alpha_3 = 0 for n <= 8  =>  hp = 1 + 2*alpha_1 + 4*alpha_2 EXACTLY")
print("=" * 78)
print("  Three vertex-disjoint odd cycles need 3+3+3 = 9 vertices.  So the identity is")
print("  a counting fact for n <= 8; the check confirms it lands where predicted.")
for n, trials in ((6, 300), (7, 120)):
    bad = 0
    m = n * (n - 1) // 2
    for _ in range(trials):
        adj = from_bits(random.getrandbits(m), n)
        cyc = odd_cycle_sets(adj, n)
        a1 = len(cyc)
        a2 = sum(1 for i in range(len(cyc)) for j in range(i + 1, len(cyc))
                 if cyc[i] & cyc[j] == 0)
        a3 = 0
        for i in range(len(cyc)):
            for j in range(i + 1, len(cyc)):
                if cyc[i] & cyc[j]:
                    continue
                for k2 in range(j + 1, len(cyc)):
                    if cyc[k2] & (cyc[i] | cyc[j]) == 0:
                        a3 += 1
        if hp_count(adj, n) != 1 + 2 * a1 + 4 * a2 or a3 != 0:
            bad += 1
    print("  n = %d : %d tournaments sampled, violations : %d  -> %s"
          % (n, trials, bad, "HOLDS" if bad == 0 else "FAILS"))

# ------------------------------------------------------------------ (B)
print()
print("=" * 78)
print("(B) SKEW SPECTRAL PARITY: char(S) is ODD for n odd, EVEN for n even")
print("=" * 78)
for n in (4, 5, 6, 7):
    m = n * (n - 1) // 2
    okpar = True
    forced = True
    for _ in range(200):
        adj = from_bits(random.getrandbits(m), n)
        cs = charpoly(skew(adj, n), n)
        # cs[k] is the coefficient of x^{n-k}
        for k in range(1, n + 1):
            if k % 2 == 1 and cs[k] != 0:
                okpar = False
        if cs[2] != n * (n - 1) // 2:
            forced = False
    print("  n = %d : char(S) has only %s powers : %-5s ; x^{n-2} coeff = C(n,2) : %s"
          % (n, "even" if n % 2 == 0 else "odd", okpar, forced))
print("  -> the parity of char(S) is exactly the parity of n.  Since S is real skew-")
print("     symmetric its spectrum is purely imaginary in +-pairs, so char(S)(iy) is a")
print("     REAL-ROOTED odd/even polynomial in y -- the sin/cos archetype, spectrally.")

print()
print("=" * 78)
print("(B2) IS x^7 + 21x^5 + 115x^3 + 119x REALISED AT n = 7?")
print("=" * 78)
target = [1, 0, 21, 0, 115, 0, 119, 0]
paley = [1, 0, 21, 0, 147, 0, 343, 0]
n = 7
m = n * (n - 1) // 2
seen = {}
hit = None
for t in range(NSAMP):
    bits = random.getrandbits(m)
    adj = from_bits(bits, n)
    cs = charpoly(skew(adj, n), n)
    key = tuple(cs)
    if key not in seen:
        seen[key] = bits
    if cs == target and hit is None:
        hit = bits
print("  sampled %d tournaments on 7 vertices; distinct skew char polys : %d"
      % (NSAMP, len(seen)))
print("  target  x^7+21x^5+115x^3+119x realised : %s" % (hit is not None))
if hit is not None:
    adj = from_bits(hit, n)
    arcs = sorted((i, j) for i in range(n) for j in range(n)
                  if i != j and (adj[i] >> j) & 1)
    print("     witness bits = %d" % hit)
    print("     arcs = %s" % (arcs,))
    print("     hp = %d" % hp_count(adj, n))
print("  Paley form x(x^2+7)^3 = x^7+21x^5+147x^3+343x realised : %s"
      % (tuple(paley) in seen))
print()
print("  the x^3 coefficients actually attained (sorted):")
xs = sorted({k[4] for k in seen})
print("     %s" % xs)
print("  the x coefficients actually attained (sorted):")
print("     %s" % sorted({k[6] for k in seen}))
