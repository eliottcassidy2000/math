#!/usr/bin/env python3
"""H_on_the_ladder_kps_S128c129.py -- kind-pasteur-2026-07-20-S128c129

NAMED-NEXT (1) of THM-1750: place H (Hamiltonian-path count / OCF) on the moment-nullcone
ladder.  The trace moments tr(A^m) sit on the RATIONAL floor (Cayley-Hamilton, depth n).
Is H a function of those moments -- i.e. spectral -- or does it live ABOVE the ceiling?

THE TEST: group labeled tournaments by their MOMENT VECTOR (tr(A^1),...,tr(A^n)) -- which by
Newton's identities is exactly the characteristic polynomial / spectrum -- and ask whether H
is CONSTANT on each group.  If two co-spectral tournaments have DIFFERENT H, then H is NOT a
spectral moment: it breaks out of the rational-algebraic-holonomic ladder (H is #P-hard, a
permanent).  THM-133's H = (462 - tr(A^4))/2 holds only for the highly symmetric Z_7
circulants; the general question is whether that spectral reduction survives.
"""
import sys
import numpy as np


def from_bits(bits, n):
    A = np.zeros((n, n), dtype=np.int64)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    return A


def moment_vec(A, n):
    v = []
    P = np.eye(n, dtype=np.int64)
    for k in range(1, n + 1):
        P = P @ A
        v.append(int(np.trace(P)))
    return tuple(v)


def ham_paths(A, n):
    """# directed Hamiltonian paths via bitmask DP.  dp[mask][last]."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if not c:
                continue
            for w in range(n):
                if not (mask >> w & 1) and A[last][w]:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full][last] for last in range(n))


print("=" * 84)
print("IS H SPECTRAL?  H constant on each moment-vector (co-spectral) class?")
print("=" * 84)
for n in (4, 5, 6):
    m = n * (n - 1) // 2
    groups = {}
    for bits in range(1 << m):
        A = from_bits(bits, n)
        mv = moment_vec(A, n)
        Al = A.tolist()
        h = ham_paths(Al, n)
        if mv not in groups:
            groups[mv] = {}
        groups[mv][h] = groups[mv].get(h, 0) + 1
    # find a co-spectral class with >1 distinct H
    split = [(mv, sorted(hs)) for mv, hs in groups.items() if len(hs) > 1]
    print("  n=%d : %d distinct moment-vectors; co-spectral classes with SPLIT H: %d"
          % (n, len(groups), len(split)))
    if split:
        mv, hs = split[0]
        print("     example: moment vector %s  carries H in %s  -> H is NOT spectral at n=%d"
              % (mv, hs, n))
    else:
        print("     H is CONSTANT on every co-spectral class at n=%d (H spectral here)" % n)
    sys.stdout.flush()

print()
print("=" * 84)
print("READING")
print("=" * 84)
print("  If H splits within a co-spectral class, H is NOT determined by the moment vector")
print("  (tr A^1..tr A^n), so it is NOT on the rational/algebraic/holonomic ladder -- it is")
print("  the tournament's #P-hard invariant (a permanent), one rung ABOVE the holonomic")
print("  ceiling.  THM-133's spectral H = (462-tr A^4)/2 is then a SYMMETRY collapse special")
print("  to Z_7 circulants, not the general law.  If H never splits, H IS spectral and sits")
print("  on the rational floor with detection depth n -- the stronger unification.")
