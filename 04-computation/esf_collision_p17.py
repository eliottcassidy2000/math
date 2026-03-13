#!/usr/bin/env python3
"""
esf_collision_p17.py — opus-2026-03-13-S70

Check whether any two non-isomorphic Z_p*-orbits at p=17 share the same
elementary symmetric functions e_j(Q). This tests whether THM-145 has content
beyond orbit classification.

At p=17, m=8: there are 2^8 = 256 total orientations, divided into Z_p*-orbits.
Each orbit has 16 elements (if stabilizer is trivial) or fewer.

We compute Q and e_j for one representative per orbit, then check for collisions.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict

def compute_Q(p, S):
    """Compute Q_k = |S_hat(k)|^2 for k=1,...,m."""
    m = (p-1)//2
    omega = np.exp(2j * np.pi / p)
    Q = []
    for k in range(1, m+1):
        s_hat = sum(omega**(s*k) for s in S)
        Q.append(abs(s_hat)**2)
    return Q

def esf(Q, j):
    """Elementary symmetric function e_j of Q values."""
    from itertools import combinations as C
    return sum(np.prod([Q[i] for i in combo]) for combo in C(range(len(Q)), j))

def all_esf(Q):
    """All elementary symmetric functions."""
    m = len(Q)
    return tuple(round(esf(Q, j), 8) for j in range(1, m+1))

def get_orbits(p):
    """Compute Z_p*-orbits of tournament orientations on Z_p."""
    m = (p-1)//2
    pairs = []
    seen = set()
    for s in range(1, p):
        if s not in seen:
            pairs.append((s, p-s))
            seen.add(s)
            seen.add(p-s)

    # Each orientation: choose one from each pair
    orientations = []
    for bits in range(2**m):
        S = frozenset(pairs[i][0] if (bits >> i) & 1 else pairs[i][1] for i in range(m))
        orientations.append(S)

    # Group into Z_p*-orbits
    orbit_map = {}  # orientation -> orbit representative
    orbits = {}     # representative -> list of orientations

    for S in orientations:
        if S in orbit_map:
            continue
        # Generate orbit: c*S mod p for all c in Z_p*
        orbit = set()
        for c in range(1, p):
            cS = frozenset((c*s) % p for s in S)
            # Also need to handle that cS might use both elements of a pair
            # cS is valid iff for each pair (a, p-a), exactly one is in cS
            valid = True
            for a, b in pairs:
                if a in cS and b in cS:
                    valid = False
                    break
                if a not in cS and b not in cS:
                    valid = False
                    break
            if valid:
                orbit.add(cS)
                orbit_map[cS] = S

        if orbit:
            orbits[S] = list(orbit)

    return orbits

print("="*70)
print("ESF COLLISION CHECK AT VARIOUS PRIMES")
print("="*70)

for p in [5, 7, 11, 13, 17]:
    m = (p-1)//2
    print(f"\n{'─'*60}")
    print(f"p = {p}, m = {m}, total orientations = {2**m}")
    print(f"{'─'*60}")

    orbits = get_orbits(p)
    print(f"  Number of Z_p*-orbits: {len(orbits)}")

    # Compute ESF for each orbit representative
    esf_to_orbits = defaultdict(list)
    orbit_data = []

    for rep, members in orbits.items():
        Q = compute_Q(p, rep)
        esfs = all_esf(Q)
        esf_to_orbits[esfs].append((rep, Q))
        orbit_data.append((rep, len(members), Q, esfs))

    orbit_data.sort(key=lambda x: x[3])  # sort by ESF

    # Check for collisions
    collisions = {k: v for k, v in esf_to_orbits.items() if len(v) > 1}

    if collisions:
        print(f"  *** ESF COLLISIONS FOUND: {len(collisions)} ***")
        for esfs, reps in collisions.items():
            print(f"    e_j = {esfs}")
            for rep, Q in reps:
                print(f"      S = {sorted(rep)}, Q = [{', '.join(f'{q:.4f}' for q in Q)}]")
    else:
        print(f"  No ESF collisions — all {len(orbits)} orbits have distinct e_j(Q)")

    # Print orbit summary
    if p <= 13:
        print(f"\n  Orbit details:")
        for rep, size, Q, esfs in orbit_data:
            print(f"    S={str(sorted(rep)):30s} |orbit|={size:2d}  e_1={esfs[0]:8.4f}")

# For p=17, also check orbit sizes
if True:
    p = 17
    m = 8
    print(f"\n{'='*70}")
    print(f"DETAILED p=17 ANALYSIS")
    print(f"{'='*70}")

    orbits = get_orbits(p)
    print(f"Number of orbits: {len(orbits)}")

    # Size distribution
    sizes = [len(v) for v in orbits.values()]
    from collections import Counter
    size_counts = Counter(sizes)
    print(f"Orbit size distribution: {dict(sorted(size_counts.items()))}")
    print(f"Sum of orbit sizes: {sum(sizes)} (should be {2**m})")

    # ESF data
    esf_to_orbits = defaultdict(list)
    orbit_data = []

    for rep, members in orbits.items():
        Q = compute_Q(p, rep)
        esfs = all_esf(Q)
        esf_to_orbits[esfs].append((sorted(rep), Q, len(members)))
        orbit_data.append((sorted(rep), len(members), esfs[0], esfs))

    collisions = {k: v for k, v in esf_to_orbits.items() if len(v) > 1}
    print(f"\nESF collisions: {len(collisions)}")

    if collisions:
        for esfs, reps in collisions.items():
            print(f"\n  Collision at e_1={esfs[0]:.4f}, e_2={esfs[1]:.4f}:")
            for rep, Q, size in reps:
                print(f"    S = {rep}, |orbit| = {size}")
                print(f"    Q = [{', '.join(f'{q:.4f}' for q in Q)}]")

    # Sort by e_1 and display
    orbit_data.sort(key=lambda x: x[2])
    print(f"\n  Orbits sorted by e_1 (total number of orbits: {len(orbit_data)}):")
    for rep, size, e1, esfs in orbit_data:
        print(f"    |orbit|={size:2d}  e_1={e1:8.4f}  S={rep}")

print("\nDONE.")
