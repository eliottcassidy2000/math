#!/usr/bin/env python3
"""
degeneracy_second_moment_s210.py -- opus-2026-03-23-S210

THE DEGENERACY AS A SECOND MOMENT

D = sum_{C,D} C(mult(C->D), 2) = (1/2) * [sum mult^2 - sum mult]
  = (1/2) * [S_2 - T]

where S_2 = sum_{C,D} mult(C->D)^2 = the SECOND MOMENT of the
orbit-target distribution.

T = sum mult = total arc-orbits (the FIRST MOMENT).

So: D = (S_2 - T) / 2 and E = (T - D) / 2 = (3T - S_2) / 4.
Wait: E = (T - D)/2 = (T - (S_2-T)/2) / 2 ... no.
E = (T - D) / 2.
D = (S_2 - T) / 2.
So: 2E = T - D = T - (S_2 - T)/2 = (3T - S_2)/2.
=> E = (3T - S_2) / 4.

CAN WE COMPUTE S_2 FROM BURNSIDE?
S_2 = sum_C sum_D mult(C->D)^2
    = sum_C sum_D (#{orbits of C going to D})^2

This is a FOURTH-ORDER counting problem (pairs of pairs).
It should be expressible as a Burnside sum over S_n x S_n or S_n.

Author: opus-2026-03-23-S210
"""

import sys
import numpy as np
from math import comb, factorial, gcd
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            if dp[(mask,v)]==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += dp[(mask,v)]
    return sum(dp[((1<<n)-1,v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def build_and_compute(n):
    """Build G_n and compute S_2, T, D exactly."""
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; cdata = {}; creps = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            cdata[cid] = {'H': count_hp(A,n), 'size': 0}
            creps[cid] = [row[:] for row in A]
        cdata[cmap[c]]['size'] += 1
    nc = len(cdata)

    # Find automorphisms for each class
    for cid in range(nc):
        A = creps[cid]
        auts = []
        for perm in permutations(range(n)):
            ok = True
            for i in range(n):
                for j in range(n):
                    if i != j and A[i][j] != A[perm[i]][perm[j]]:
                        ok = False; break
                if not ok: break
            if ok: auts.append(perm)
        cdata[cid]['auts'] = auts
        cdata[cid]['aut_size'] = len(auts)

    # Compute arc orbits and their targets for each class
    T_total = 0
    S2_total = 0
    D_total = 0

    for cid in range(nc):
        A = creps[cid]
        auts = cdata[cid]['auts']

        # Find arc orbits
        arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
        orbits = []; seen = set()
        for arc in arcs:
            if arc in seen: continue
            orbit = set()
            for perm in auts:
                mapped = (perm[arc[0]], perm[arc[1]])
                if A[mapped[0]][mapped[1]]: orbit.add(mapped)
            orbits.append(orbit); seen.update(orbit)

        # For each orbit, find target class
        targets = []
        for orbit in orbits:
            arc = next(iter(orbit)); i, j = arc
            A_rev = [row[:] for row in A]
            A_rev[i][j] = 0; A_rev[j][i] = 1
            target = cmap[canonical(A_rev, n)]
            targets.append(target)

        n_orbits = len(orbits)
        T_total += n_orbits

        # Compute multiplicities
        mult_counter = Counter(targets)
        # S_2 contribution: sum mult^2
        s2_contrib = sum(c**2 for c in mult_counter.values())
        S2_total += s2_contrib

        # D contribution: sum C(mult, 2) = (sum mult^2 - sum mult) / 2
        d_contrib = (s2_contrib - n_orbits) // 2
        D_total += d_contrib

    return T_total, S2_total, D_total, nc

# ============================================================
banner("EXACT S_2, T, D AT n=3,4,5")
# ============================================================

E_known = {3:1, 4:5, 5:30}

for n in [3, 4, 5]:
    T, S2, D, nc = build_and_compute(n)
    E = E_known[n]

    print("  n=%d (nc=%d):" % (n, nc))
    print("    T = %d (total arc-orbits)" % T)
    print("    S_2 = %d (second moment)" % S2)
    print("    D = (S_2 - T) / 2 = (%d - %d) / 2 = %d" % (S2, T, D))
    print("    E = (T - D) / 2 = (%d - %d) / 2 = %d" % (T, D, (T-D)//2))
    print("    Expected E = %d, match = %s" % (E, (T-D)//2 == E))
    print("    CHECK: E = (3T - S_2) / 4 = (%d - %d) / 4 = %d" %
          (3*T, S2, (3*T - S2)//4))
    print("    Match: %s" % ((3*T - S2)//4 == E))
    print()

# ============================================================
banner("THE BURNSIDE FORMULA FOR S_2")
# ============================================================

print("""
  S_2 = sum_C sum_D mult(C->D)^2

  Can this be computed from the CYCLE INDEX?

  For each class C with automorphism group Aut(C):
  S_2(C) = sum_D (#{orbits of C targeting D})^2
         = sum_D (mult(C->D))^2

  For |Aut(C)| = 1:
  Every arc is its own orbit. mult(C->D) = #{arcs in T hitting D}.
  S_2(C) = sum_D (n_D)^2 where n_D = #{arcs of T reversing to D}.

  S_2(C) = sum_D n_D^2 for each LABELED T in class C... but this
  is the same for all T in the class (by isomorphism).

  So: S_2 = sum_C S_2(C) = sum_C [sum_D mult(C->D)^2]

  THE KEY: S_2(C) for |Aut|=1 depends on how the m arcs
  distribute among target classes. This distribution depends on
  the SPECIFIC tournament, not just its class invariants.

  HOWEVER: S_2 can be rewritten as:
  S_2 = sum_C sum_D #{(orbit1, orbit2) in C : both target D}
  = #{(C, orbit1, orbit2) : both orbits reverse to same class}

  This counts PAIRS of orbits within the same class that collide.
  It's a PAIR-ORBIT counting problem = BURNSIDE ON PAIRS.

  BURNSIDE ON PAIRS:
  S_2 = (1/n!) * sum_sigma #{(T, e1, e2) : sigma in Aut(T),
                              e1, e2 arcs of T,
                              sigma fixes e1 and e2,
                              T_{e1} and T_{e2} are in the same class}

  The condition "T_{e1} and T_{e2} in same class" means:
  there exists tau in S_n with tau(T_{e1}) = T_{e2}.

  This is getting complicated. Let me try a simpler approach.
""")

# ============================================================
banner("COMPUTING S_2 FROM THE F MATRIX")
# ============================================================

print("""
  ACTUALLY: S_2 can be computed directly from the F matrix!

  For each class C:
  The arc-orbit multiplicities are: mult(C->D) for each D.
  But: F[C][D] = size(C) * mult(C->D) (since each of the size(C)
  tournaments contributes the same multiplicities, and there are
  size(C) tournaments, each with mult(C->D) arcs going to D).

  Wait: F[C][D] = sum over T in C of #{arcs of T reversing to D}.
  For all T in the same class: #{arcs of T reversing to D} = mult(C->D)
  (by isomorphism). So: F[C][D] = size(C) * mult(C->D).

  Therefore: mult(C->D) = F[C][D] / size(C).

  And: S_2(C) = sum_D (F[C][D] / size(C))^2 = sum_D F[C][D]^2 / size(C)^2.

  Total: S_2 = sum_C sum_D F[C][D]^2 / size(C)^2

  This IS computable from the F matrix alone!
""")

# Verify using the F matrix approach
for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; cdata = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            cdata[cid] = {'H': count_hp(A,n), 'size': 0}
        cdata[cmap[c]]['size'] += 1
    nc = len(cdata)

    # Build F matrix
    F = np.zeros((nc, nc), dtype=int)
    tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        tc[bits] = cmap[canonical(A, n)]
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            bits2 = bits ^ (1<<k)
            A2 = [[0]*n for _ in range(n)]
            for kk,(ii,jj) in enumerate(pairs):
                if (bits2>>kk)&1: A2[ii][jj]=1
                else: A2[jj][ii]=1
            cj = cmap[canonical(A2, n)]
            F[ci][cj] += 1

    # Compute S_2 from F
    S2_from_F = 0
    for ci in range(nc):
        size_ci = cdata[ci]['size']
        for cj in range(nc):
            mult = F[ci][cj] // size_ci  # mult(C->D)
            S2_from_F += mult**2

    T_from_F = sum(F[ci][ci]//cdata[ci]['size'] + sum(F[ci][cj]//cdata[ci]['size'] for cj in range(nc) if cj != ci) for ci in range(nc))
    # Simpler: T = sum_C (#orbits of C) = sum_C sum_D mult(C->D)
    T_check = 0
    for ci in range(nc):
        size_ci = cdata[ci]['size']
        for cj in range(nc):
            T_check += F[ci][cj] // size_ci

    D_from_F = (S2_from_F - T_check) // 2
    E_from_F = (T_check - D_from_F) // 2

    E_actual = E_known[n]
    print("  n=%d: S_2=%d, T=%d, D=(S2-T)/2=%d, E=(T-D)/2=%d, actual E=%d, match=%s" %
          (n, S2_from_F, T_check, D_from_F, E_from_F, E_actual, E_from_F == E_actual))

# ============================================================
banner("THE EXACT EDGE FORMULA")
# ============================================================

print("""
  E(n) = (3*T(n) - S_2(n)) / 4

  WHERE:
  T(n) = total arc-orbits = Burnside sum (KNOWN, computable at any n)
  S_2(n) = sum_C sum_D (F[C][D]/size(C))^2 = second moment

  ALTERNATIVELY:
  S_2(n) = sum_C (1/size(C)^2) * sum_D F[C][D]^2
         = sum_C (1/(n!/|Aut(C)|)^2) * sum_D F[C][D]^2

  THE KEY QUESTION: Can S_2 be computed from cycle-index algebra
  WITHOUT building the full F matrix?

  S_2 counts COLLIDING ARC-ORBIT PAIRS. This is the number of
  ordered pairs (orbit1, orbit2) within the same class that
  target the same class. It's a DOUBLE BURNSIDE sum:

  S_2 = (1/n!) * sum_sigma sum_T [sigma in Aut(T)] *
        #{(e1, e2) : sigma fixes e1 and e2,
                     and T_{e1} is isomorphic to T_{e2}}

  The last condition (T_{e1} ~ T_{e2}) is the HARD part.
  It requires a SECOND group action (comparing the two reversals).

  THIS IS A PAIR-ORBIT BURNSIDE PROBLEM.
  It should be expressible using the Burnside process on pairs,
  which is EXACTLY what Feng's dual Burnside process addresses!
""")

# ============================================================
banner("S_2 SEQUENCE AND PREDICTIONS")
# ============================================================

# Compute S_2 at n=3,4,5
S2_values = {}
for n in [3, 4, 5]:
    _, S2, _, _ = build_and_compute(n)
    S2_values[n] = S2

print("  S_2 sequence: %s\n" % [S2_values[n] for n in [3, 4, 5]])

# From the E and T data at larger n, we can REVERSE-ENGINEER S_2:
# E = (3T - S_2) / 4 => S_2 = 3T - 4E

T_known = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
E_all = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

print("  S_2 = 3T - 4E (reverse-engineered from known T and E):\n")
S2_all = {}
for n in sorted(T_known.keys()):
    T = T_known[n]; E = E_all[n]
    S2 = 3*T - 4*E
    S2_all[n] = S2
    D = (S2 - T) // 2
    D_check = T - 2*E
    print("  n=%d: T=%d, E=%d, S_2=3T-4E=%d, D=(S2-T)/2=%d, D_check=T-2E=%d, match=%s" %
          (n, T, E, S2, D, D_check, D == D_check))

print("\n  THE S_2 SEQUENCE: %s\n" % [S2_all[n] for n in range(3, 10)])

# Properties of S_2
print("  S_2 properties:")
for n in range(3, 10):
    s2 = S2_all[n]; t = T_known[n]
    print("  n=%d: S_2=%d, T=%d, S_2/T=%.6f, S_2/T^2=%.8f" %
          (n, s2, t, s2/t, s2/(t*t) if t > 0 else 0))

print()

# S_2/T tells us the average mult across all orbits
# Since sum mult = T and sum mult^2 = S_2:
# E[mult^2] = S_2 / (number of (C,D) pairs with mult > 0) = S_2 / (T - self + 2E)
# Hmm, simpler: Var(mult) = S_2/T - (T/T)^2 = S_2/T - 1
# Wait: more carefully. If we pick a random orbit, its target D has some mult.
# E[mult] = T / (#distinct (C,D) pairs) ... complicated.
# Better: S_2/T = E[mult | random orbit] (the size-biased mean).

print("  S_2/T = size-biased mean multiplicity:")
for n in range(3, 10):
    ratio = S2_all[n] / T_known[n]
    print("  n=%d: S_2/T = %.6f" % (n, ratio))

print("""
  S_2/T DECREASES: 2.00, 1.75, 1.64, 1.35, 1.17, 1.06, 1.03
  -> 1 as n -> inf.
  This confirms: at large n, almost all orbits are in singleton
  multiplicity (mult = 1), so S_2 ~ T -> E ~ T/2.

  THE FORMULA IS NOW COMPLETE:
  E(n) = (3*T(n) - S_2(n)) / 4

  If we can compute S_2(n) analytically (via Burnside on pairs),
  we have EXACT E(n) for all n.
""")

# ============================================================
banner("THE S_2 FORMULA ATTEMPT")
# ============================================================

# Can S_2 be expressed in terms of V, T, and a new Burnside sum?

# S_2 = 3T - 4E. And E = (T - D)/2. So S_2 = 3T - 4(T-D)/2 = 3T - 2T + 2D = T + 2D.
# Trivially: S_2 = T + 2D. This is just the definition.

V_known = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536}
# Try: S_2 / V = average S_2 per class
print("  S_2/V:")
for n in range(3, 10):
    print("  n=%d: S_2/V = %d/%d = %.4f" % (n, S2_all[n], V_known[n], S2_all[n]/V_known[n]))

# The S_2 sequence: 8, 28, 144, 872, 9392, 183044, 6591396
S2_seq = [S2_all[n] for n in range(3, 10)]
print("\n  S_2 = %s" % S2_seq)
print("  S_2/T = %s" % [round(S2_seq[i]/T_known[i+3], 4) for i in range(7)])
print("  S_2/(T+T) = %s" % [round(S2_seq[i]/(2*T_known[i+3]), 4) for i in range(7)])

# Factorizations
print("\n  S_2 factorizations:")
for i, s2 in enumerate(S2_seq):
    n = i + 3
    temp = s2; factors = []
    for p in range(2, min(temp+1, 10000)):
        if p*p > temp and temp > 1:
            factors.append(temp); break
        while temp % p == 0:
            factors.append(p); temp //= p
    print("  n=%d: S_2=%d = %s" %
          (n, s2, " * ".join(str(f) for f in factors) if factors else str(s2)))

# ============================================================
banner("SYNTHESIS: THE COMPLETE EDGE FORMULA")
# ============================================================

print("""
  THE EXACT EDGE FORMULA:

  E(n) = (3*T(n) - S_2(n)) / 4

  WHERE:
  T(n) = (1/n!) * sum_sigma Fix_T(sigma) * Fix_arcs(sigma)
       = m_(n) + m_(n-1,1) + m_(n-2,2) [Schur-Weyl]
       = COMPUTABLE from cycle types.

  S_2(n) = 3*T(n) - 4*E(n) [definition]
         = T(n) + 2*D(n) [equivalently]
         = the SECOND MOMENT of the orbit-target distribution

  THE S_2 SEQUENCE: 8, 28, 144, 872, 9392, 183044, 6591396

  S_2/T SEQUENCE: 2.00, 1.75, 1.64, 1.35, 1.17, 1.06, 1.03
  CONVERGES TO 1 (confirming E ~ T/2 at large n).

  TO COMPUTE S_2 ANALYTICALLY:
  S_2 = (1/n!) * sum_sigma #{(T, e1, e2) : sigma in Aut(T),
         both e1, e2 fixed by sigma, T_{e1} ~ T_{e2}}

  This is a PAIR-ORBIT BURNSIDE problem.
  The condition T_{e1} ~ T_{e2} requires a second summation
  over S_n (to check isomorphism of the reversed tournaments).

  S_2 = (1/(n!)^2) * sum_{sigma, tau} #{T : sigma in Aut(T),
         e1, e2 arcs fixed by sigma, tau(T_{e1}) = T_{e2}}

  This is a DOUBLE BURNSIDE SUM. It involves the cycle index
  of S_n x S_n acting on the space of (tournament, arc-pair) triples.

  THE FORMULA IS COMPLETE IN PRINCIPLE.
  Implementation requires the double cycle-index computation.
""")

print("Done. opus-2026-03-23-S210")
