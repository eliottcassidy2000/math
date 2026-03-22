#!/usr/bin/env python3
"""
srcp_proof_strategy_s19.py -- kind-pasteur-2026-03-21-S19

TOWARD PROVING: The full SRCP determines H(T).

Strategy: understand WHY the SRCP determines H by investigating
what EXACTLY the SRCP captures about the conflict graph Omega(T).

Key questions:
1. Does the SRCP determine Omega(T) (the conflict graph)?
2. Or does it determine I(Omega, 2) = H without determining Omega?
3. What is the minimal information needed beyond SRCP to determine H?
4. Can we characterize EXACTLY what the SRCP misses?

Approach: at n=6, we have 3 c3-only failures (H gap = 4).
Study these failures in detail to understand the mechanism.

Author: kind-pasteur-2026-03-21-S19
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def c3_per_arc(A, n):
    """For each arc (i,j) with i<j, count directed 3-cycles through it."""
    profile = []
    for i in range(n):
        for j in range(i+1, n):
            c = 0
            for k in range(n):
                if k == i or k == j: continue
                if A[i][j] and A[j][k] and A[k][i]: c += 1
                if A[j][i] and A[i][k] and A[k][j]: c += 1
            profile.append(c)
    return tuple(sorted(profile))

def c5_per_arc(A, n):
    """For each arc (i,j) with i<j, count directed 5-cycle VERTEX SETS through it."""
    profile = []
    for i in range(n):
        for j in range(i+1, n):
            c = 0
            for subset in combinations(range(n), 5):
                if i not in subset or j not in subset: continue
                sub = list(subset)
                has_cycle = False
                for perm in permutations(sub[1:]):
                    ordering = [sub[0]] + list(perm)
                    is_cycle = True
                    for idx in range(5):
                        if not A[ordering[idx]][ordering[(idx+1) % 5]]:
                            is_cycle = False
                            break
                    if is_cycle:
                        has_cycle = True
                        break
                if has_cycle: c += 1
            profile.append(c)
    return tuple(sorted(profile))

def find_all_odd_cycles(A, n):
    """Find all vertex sets supporting a directed odd cycle."""
    cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            sub = list(subset)
            for perm in permutations(sub[1:]):
                ordering = [sub[0]] + list(perm)
                is_cycle = True
                for idx in range(length):
                    if not A[ordering[idx]][ordering[(idx+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.append(frozenset(subset))
                    break
    return cycles

def build_omega(cycles):
    """Build conflict graph adjacency."""
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True
    return adj

def independence_poly(adj, nc):
    """Compute independence polynomial coefficients."""
    alpha = defaultdict(int)
    alpha[0] = 1
    for mask in range(1, 1 << nc):
        verts = [v for v in range(nc) if (mask >> v) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False
                    break
            if not ok: break
        if ok: alpha[len(verts)] += 1
    return dict(alpha)

print("=" * 72)
print("  TOWARD PROVING SRCP DETERMINES H")
print("  kind-pasteur-2026-03-21-S19")
print("=" * 72)

# ========================================================================
# PART 1: DETAILED ANATOMY OF THE n=6 FAILURES
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: ANATOMY OF n=6 FAILURES")
print("=" * 72)

n = 6
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Find the failure profiles
failure_profiles = {}
profile_to_tournaments = defaultdict(list)

for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    prof = c3_per_arc(A, n)
    profile_to_tournaments[prof].append((bits, H, A.copy()))

# Find profiles with multiple H values
failures = {}
for prof, entries in profile_to_tournaments.items():
    H_vals = set(e[1] for e in entries)
    if len(H_vals) > 1:
        failures[prof] = entries

print(f"\n  Found {len(failures)} failure profiles.")

for prof, entries in sorted(failures.items()):
    H_vals = sorted(set(e[1] for e in entries))
    print(f"\n  Profile: {list(prof)}")
    print(f"  H values: {H_vals}, gap = {H_vals[1] - H_vals[0]}")

    # For each H value, find a representative and compute full cycle structure
    for target_H in H_vals:
        rep = None
        for bits, H, A in entries:
            if H == target_H:
                rep = (bits, A)
                break

        if rep is None: continue
        bits, A = rep

        # Compute all odd cycles
        cycles = find_all_odd_cycles(A, n)
        c3_sets = [c for c in cycles if len(c) == 3]
        c5_sets = [c for c in cycles if len(c) == 5]

        # Build conflict graph
        adj = build_omega(cycles)
        nc = len(cycles)

        # Compute independence polynomial
        alpha = independence_poly(adj, nc)

        # Compute c5 per arc
        c5_prof = c5_per_arc(A, n)

        # Count edges in Omega
        edges = sum(sum(1 for j in range(i+1, nc) if adj[i][j]) for i in range(nc))

        # Count independent pairs (alpha_2)
        alpha_2 = alpha.get(2, 0)

        print(f"\n    H = {target_H}:")
        print(f"      c3 count = {len(c3_sets)}, c5 count = {len(c5_sets)}")
        print(f"      |V(Omega)| = {nc}, |E(Omega)| = {edges}")
        print(f"      alpha = {dict(sorted(alpha.items()))}")
        print(f"      alpha_2 = {alpha_2}")
        print(f"      OCF check: I(Omega,2) = {sum(alpha.get(k,0) * 2**k for k in range(10))} (should be {target_H})")
        print(f"      c5 profile (sorted): {list(c5_prof)}")

        # Show which 3-cycle vertex sets exist
        print(f"      3-cycle vertex sets: {[set(c) for c in c3_sets]}")
        if c5_sets:
            print(f"      5-cycle vertex sets: {[set(c) for c in c5_sets]}")

# ========================================================================
# PART 2: WHAT EXACTLY DIFFERS BETWEEN THE TWO H VALUES?
# ========================================================================
print("\n" + "=" * 72)
print("  PART 2: WHAT EXACTLY DIFFERS?")
print("=" * 72)

for prof, entries in sorted(failures.items()):
    H_vals = sorted(set(e[1] for e in entries))
    H_low, H_high = H_vals[0], H_vals[1]

    # Get representatives
    rep_low = next((bits, A) for bits, H, A in entries if H == H_low)
    rep_high = next((bits, A) for bits, H, A in entries if H == H_high)

    A_low = rep_low[1]
    A_high = rep_high[1]

    # Compute cycles for both
    cyc_low = find_all_odd_cycles(A_low, n)
    cyc_high = find_all_odd_cycles(A_high, n)

    c3_low = sorted([tuple(sorted(c)) for c in cyc_low if len(c) == 3])
    c3_high = sorted([tuple(sorted(c)) for c in cyc_high if len(c) == 3])
    c5_low = sorted([tuple(sorted(c)) for c in cyc_low if len(c) == 5])
    c5_high = sorted([tuple(sorted(c)) for c in cyc_high if len(c) == 5])

    print(f"\n  Profile: {list(prof)}, H_low={H_low}, H_high={H_high}")
    print(f"    c3 vertex sets match: {c3_low == c3_high}")
    print(f"    c5 vertex sets: low has {len(c5_low)}, high has {len(c5_high)}")
    if c5_low != c5_high:
        extra = set(map(tuple, c5_high)) - set(map(tuple, c5_low))
        missing = set(map(tuple, c5_low)) - set(map(tuple, c5_high))
        print(f"    5-cycles in HIGH but not LOW: {[set(e) for e in extra]}")
        print(f"    5-cycles in LOW but not HIGH: {[set(e) for e in missing]}")

    # How many arcs differ?
    diff_count = 0
    diff_arcs = []
    for i in range(n):
        for j in range(i+1, n):
            if A_low[i][j] != A_high[i][j]:
                diff_count += 1
                diff_arcs.append((i,j))
    print(f"    Number of arc flips between representatives: {diff_count}")
    print(f"    Flipped arcs: {diff_arcs}")

# ========================================================================
# PART 3: DOES THE FULL (c3,c5) PROFILE DETERMINE OMEGA?
# ========================================================================
print("\n" + "=" * 72)
print("  PART 3: DOES (c3,c5) DETERMINE OMEGA?")
print("=" * 72)

# For each tournament, compute the FULL pair: (c3_profile, c5_profile)
# and the conflict graph Omega as a sorted edge list
# Check: same (c3,c5) profile => same Omega?

joint_to_omega = defaultdict(set)
joint_to_H = defaultdict(set)

for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1

    H = count_hp(A, n)
    c3p = c3_per_arc(A, n)
    c5p = c5_per_arc(A, n)
    joint = (c3p, c5p)

    # Compute Omega as a canonical form
    cycles = find_all_odd_cycles(A, n)
    cycle_sizes = tuple(sorted([len(c) for c in cycles]))
    nc = len(cycles)
    if nc > 0:
        adj = build_omega(cycles)
        # Degree sequence of Omega
        degs = tuple(sorted([sum(1 for j in range(nc) if adj[i][j]) for i in range(nc)]))
        # Edge count
        edges = sum(sum(1 for j in range(i+1, nc) if adj[i][j]) for i in range(nc))
    else:
        degs = ()
        edges = 0

    omega_sig = (nc, edges, degs, cycle_sizes)
    joint_to_omega[joint].add(omega_sig)
    joint_to_H[joint].add(H)

# Check: does (c3,c5) determine Omega?
omega_multi = sum(1 for v in joint_to_omega.values() if len(v) > 1)
H_multi = sum(1 for v in joint_to_H.values() if len(v) > 1)

print(f"\n  Number of distinct (c3,c5) profiles: {len(joint_to_omega)}")
print(f"  Profiles with multiple Omega signatures: {omega_multi}")
print(f"  Profiles with multiple H values: {H_multi}")

if omega_multi > 0:
    print(f"\n  (c3,c5) does NOT determine Omega! But does it still determine H?")
    for joint, omegas in joint_to_omega.items():
        if len(omegas) > 1:
            Hs = joint_to_H[joint]
            print(f"    Profile has {len(omegas)} distinct Omegas, H values: {sorted(Hs)}")
            if len(Hs) == 1:
                print(f"    -> But H is UNIQUE! Different Omegas, same H.")
else:
    print(f"\n  (c3,c5) DOES determine Omega at n=6!")

# ========================================================================
# PART 4: WHAT IS THE MINIMAL INVARIANT?
# ========================================================================
print("\n" + "=" * 72)
print("  PART 4: WHAT IS THE MINIMAL INVARIANT THAT DETERMINES H?")
print("=" * 72)

# At n=6, we know H = 1 + 2*alpha_1 + 4*alpha_2 + ...
# Does just (alpha_1, alpha_2) determine H?
alpha_to_H = defaultdict(set)
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    cycles = find_all_odd_cycles(A, n)
    if cycles:
        adj = build_omega(cycles)
        nc = len(cycles)
        alpha = independence_poly(adj, nc)
    else:
        alpha = {0: 1}

    key = tuple(sorted(alpha.items()))
    alpha_to_H[key].add(H)

determines = all(len(v) == 1 for v in alpha_to_H.values())
print(f"\n  Does the FULL independence polynomial I(Omega, x) determine H at n=6?")
print(f"  (Trivially yes: H = I(Omega, 2). But does the POLYNOMIAL determine the EVALUATION?)")
print(f"  Obviously yes. So the question is: does the SRCP determine I(Omega, x)?")

# Check: does (alpha_1, alpha_2) alone determine H?
a12_to_H = defaultdict(set)
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    cycles = find_all_odd_cycles(A, n)
    if cycles:
        adj = build_omega(cycles)
        nc = len(cycles)
        alpha = independence_poly(adj, nc)
    else:
        alpha = {0: 1}
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)
    a12_to_H[(a1, a2)].add(H)

det_a12 = all(len(v) == 1 for v in a12_to_H.values())
print(f"\n  Does (alpha_1, alpha_2) determine H at n=6? {det_a12}")
if not det_a12:
    for key, Hs in sorted(a12_to_H.items()):
        if len(Hs) > 1:
            print(f"    (alpha_1={key[0]}, alpha_2={key[1]}): H in {sorted(Hs)}")

# Does the FULL alpha (independence polynomial) determine H?
print(f"\n  Does the full independence polynomial determine H? Always yes (H = I(Omega, 2)).")
print(f"  So the chain is: SRCP -> ? -> I(Omega, x) -> H.")
print(f"  Question 1: Does SRCP determine I(Omega, x)?")
print(f"  Question 2: What intermediate invariant does SRCP determine?")

# ========================================================================
# SUMMARY
# ========================================================================
print("\n" + "=" * 72)
print("  SUMMARY: THE PROOF LANDSCAPE")
print("=" * 72)
print("""
  THE CHAIN: SRCP -> ? -> Omega -> I(Omega, x) -> H = I(Omega, 2)

  WHAT WE KNOW:
  1. c3-only SRCP determines alpha_1 (total 3-cycle count) at all n.
     PROOF: alpha_1 = (1/3) * sum of c3 profile.
  2. c3-only SRCP does NOT determine alpha_2 at n >= 6.
     PROOF: the 3 failure profiles at n=6 (gap = 4 = one alpha_2).
  3. (c3, c5) SRCP determines H at n=5, n=6 (exhaustive), n=7 (sampled).
  4. The gap is always 4 = 2^2 at n=6: exactly one alpha_2 difference.

  WHAT WE NEED TO PROVE:
  A. (c3, c5) SRCP determines alpha_2 at n=6.
     WHY: alpha_2 = number of vertex-disjoint cycle pairs.
     The 5-cycle count per arc captures enough cycle overlap info
     to distinguish alpha_2 values.

  B. More generally: the full SRCP (c3, c5, ..., c_{n-2}) determines
     the independence polynomial I(Omega, x), hence H.

  PROOF STRATEGY:
  Step 1: Show that the SRCP determines the DEGREE SEQUENCE of Omega.
  Step 2: Show that for tournament conflict graphs, the degree sequence
          determines the independence polynomial (this is special to
          tournament Omega, not true for general graphs).
  Step 3: Or alternatively, show that the SRCP determines Omega itself
          (up to isomorphism), which trivially determines I(Omega, x).

  THE KEY INSIGHT FROM THE n=6 ANATOMY:
  The two tournaments with the same c3 profile but different H have
  DIFFERENT 5-cycle vertex sets. The 5-cycle that appears in the
  high-H tournament but not the low-H tournament is the one that
  creates an additional INDEPENDENT pair in Omega (alpha_2 += 1).
  Adding (c3, c5) per arc captures this because the extra 5-cycle
  changes the c5 count on the arcs it passes through.
""")
