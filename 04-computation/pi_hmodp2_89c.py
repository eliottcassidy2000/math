#!/usr/bin/env python3
"""
pi_hmodp2_89c.py — Deep analysis of H(P_p) mod p² and Wilson connection
opus-2026-03-14-S89c

Key question: What is H(P_p)/p mod p?
Data so far: 1, 6, 10, 17 for p=3,7,11,19
That's p-2, p-1, p-1, p-2. Mixed pattern!

Let's compute more and look for the real pattern.
"""

from itertools import permutations
from math import gcd, factorial
from fractions import Fraction

def paley_tournament(p):
    """Build adjacency for Paley tournament on Z/pZ."""
    # QR = quadratic residues mod p
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    adj = {}
    for i in range(p):
        adj[i] = set()
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i].add(j)
    return adj

def count_hamiltonian_paths(adj, n):
    """Count Hamiltonian paths by DP on bitmask."""
    # dp[mask][v] = number of Hamiltonian paths ending at v using vertices in mask
    dp = [{} for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for u in adj[v]:
                if mask & (1 << u) == 0:
                    new_mask = mask | (1 << u)
                    dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]

    full = (1 << n) - 1
    return sum(dp[full].values())

print("=" * 70)
print("H(P_p) mod p² — Extended Analysis")
print("=" * 70)

# Compute for p = 3, 5, 7, 11, 13
# p=5 is NOT ≡ 3 mod 4, so no Paley. Skip.
# p=13 is ≡ 1 mod 4, skip.
# Paley primes ≡ 3 mod 4: 3, 7, 11, 19, 23, 31, 43, ...

paley_primes = [3, 7, 11, 19, 23]

for p in paley_primes:
    if p > 19:
        print(f"\n  p={p}: SKIPPING (bitmask DP needs 2^{p} states)")
        continue

    adj = paley_tournament(p)
    H = count_hamiltonian_paths(adj, p)

    print(f"\n  p={p}:")
    print(f"    H(P_p) = {H}")
    print(f"    H mod p = {H % p}")
    print(f"    H/p = {H // p}")
    print(f"    H/p mod p = {(H // p) % p}")
    print(f"    p-1 = {p-1}, p-2 = {p-2}")
    print(f"    H mod p² = {H % (p*p)}")

    # Check: is H ≡ -2 mod p² ?
    print(f"    H mod p² vs p(p-2) = {p*(p-2)}: {'✓' if H % (p*p) == p*(p-2) else '✗'}")
    print(f"    H mod p² vs p(p-1) = {p*(p-1)}: {'✓' if H % (p*p) == p*(p-1) else '✗'}")

    # Wilson connection: (p-1)! ≡ -1 mod p
    wils = factorial(p-1)
    print(f"    (p-1)! = {wils}")
    print(f"    H / ((p-1)!/2^{(p-1)//2}) = {Fraction(H, wils // (2**((p-1)//2)))}")

    # Look at H mod various things
    print(f"    H mod 2p = {H % (2*p)}")
    print(f"    H mod 4p = {H % (4*p)}")

print()
print("=" * 70)
print("Deeper: H(P_p) ≡ ? (mod p²)")
print("=" * 70)

# Let's check H mod p² exactly
for p in [3, 7, 11, 19]:
    adj = paley_tournament(p)
    H = count_hamiltonian_paths(adj, p)
    r = H % (p*p)
    print(f"  p={p}: H = {H}, H mod p² = {r}")
    # What is r / p mod p?
    print(f"    r/p = {r // p} (should be H/p mod p = {(H//p) % p})")
    # Is r = p * (p - 2)?
    print(f"    p*(p-2) = {p*(p-2)}, match: {r == p*(p-2)}")
    # Is r = p * ((p-1)! mod p²)/p ?
    # Wilson: (p-1)! ≡ -1 mod p. What about mod p²?
    w = factorial(p-1) % (p*p)
    print(f"    (p-1)! mod p² = {w}")
    # Wolstenholme: For p≥5, (p-1)! ≡ -1 mod p² iff p is Wolstenholme prime
    # Actually: (p-1)! ≡ -1 + a_p * p (mod p²) where a_p = ...
    # Standard result: (p-1)! ≡ -1 (mod p) always
    # For p ≥ 5: (p-1)! ≡ -1 (mod p²) iff p is Wolstenholme prime (only known: 16843, 2124679)
    print(f"    ((p-1)! + 1) / p mod p = {((w + 1) // p) % p if (w + 1) % p == 0 else 'N/A'}")

print()
print("=" * 70)
print("The REAL pattern: H(P_p) as fraction of (p-1)!")
print("=" * 70)

for p in [3, 7, 11, 19]:
    adj = paley_tournament(p)
    H = count_hamiltonian_paths(adj, p)
    f = Fraction(H, factorial(p-1))
    print(f"  p={p}: H/(p-1)! = {f} = {float(f):.10f}")
    # H / (mean H) = H / (p!/2^{p-1}) = H * 2^{p-1} / p!
    ratio = Fraction(H * 2**(p-1), factorial(p))
    print(f"    H/mean = {ratio} = {float(ratio):.6f}")

print()
print("=" * 70)
print("H/p mod p: looking for the REAL pattern")
print("=" * 70)

vals = []
for p in [3, 7, 11, 19]:
    adj = paley_tournament(p)
    H = count_hamiltonian_paths(adj, p)
    r = (H // p) % p
    vals.append((p, r))
    print(f"  p={p}: H/p mod p = {r}")

print(f"\n  Values: {[v[1] for v in vals]}")
print(f"  p-1:    {[v[0]-1 for v in vals]}")
print(f"  p-2:    {[v[0]-2 for v in vals]}")

# Check: H ≡ p(p-2) mod p² for p=3,19 but H ≡ p(p-1) for p=7,11
# That means H/p ≡ -2 mod p for p=3,19 and H/p ≡ -1 mod p for p=7,11
# Hmm. -2 = p-2 and -1 = p-1.
# p=3: H/p ≡ -2 mod 3
# p=7: H/p ≡ -1 mod 7
# p=11: H/p ≡ -1 mod 11
# p=19: H/p ≡ -2 mod 19
# Pattern by p mod 8?
# p=3 ≡ 3 mod 8 → -2
# p=7 ≡ 7 mod 8 → -1
# p=11 ≡ 3 mod 8 → -1  (breaks p mod 8 pattern)
# p=19 ≡ 3 mod 8 → -2
# Hmm, not clean.

# Try p mod 12:
print(f"\n  p mod 12: {[(p, p%12) for p,_ in vals]}")
# p=3 ≡ 3, p=7 ≡ 7, p=11 ≡ 11, p=19 ≡ 7
# r: -2, -1, -1, -2
# Nope.

# Let's look at (-1)^{(p-3)/4}
print(f"\n  (-1)^((p-3)/4):")
for p, r in vals:
    exp = (p - 3) // 4
    sign = (-1)**exp
    print(f"    p={p}: exp={exp}, sign={sign}, r mod p = {r}, -1-sign = {-1-sign}")
    # p=3: exp=0, sign=+1, r=1=-2, -1-sign=-2 ✓
    # p=7: exp=1, sign=-1, r=6=-1, -1-sign=0 ✗

# Try class number h(-p)
print(f"\n  Class numbers h(-p):")
# h(-3)=1, h(-7)=1, h(-11)=1, h(-19)=1 (all class number 1)
# So that's not it.

# Let's just see if there's a Legendre symbol pattern
# (2/p) = (-1)^{(p²-1)/8}
print(f"\n  Legendre symbol (2/p):")
for p, r in vals:
    leg2 = pow(2, (p-1)//2, p)
    if leg2 == p-1:
        leg2 = -1
    print(f"    p={p}: (2/p) = {leg2}, H/p mod p = {r}")
    # p=3: (2/p) = -1, r = -2
    # p=7: (2/p) = +1, r = -1
    # p=11: (2/p) = -1, r = -1
    # p=19: (2/p) = -1, r = -2

# Maybe: H/p ≡ -1 - (2/p) mod p?
print(f"\n  Test: H/p ≡ -1 - (2/p) mod p?")
for p, r in vals:
    leg2 = pow(2, (p-1)//2, p)
    if leg2 > p//2:
        leg2 = leg2 - p  # make it -1
    pred = (-1 - leg2) % p
    match = "✓" if pred == r else "✗"
    print(f"    p={p}: -1-(2/p) = -1-({leg2}) = {-1-leg2} ≡ {pred} mod p, actual = {r} {match}")

# Maybe H/p ≡ (2/p) - 2 mod p?
print(f"\n  Test: H/p ≡ (2/p) - 2 mod p?")
for p, r in vals:
    leg2 = pow(2, (p-1)//2, p)
    if leg2 > p//2:
        leg2_signed = leg2 - p
    else:
        leg2_signed = leg2
    pred = (leg2_signed - 2) % p
    match = "✓" if pred == r else "✗"
    print(f"    p={p}: (2/p)-2 = {leg2_signed}-2 = {leg2_signed-2} ≡ {pred} mod p, actual = {r} {match}")

print()
print("=" * 70)
print("Cycle counts mod p for Paley")
print("=" * 70)

# Count directed odd cycles in P_p by length
def count_directed_cycles(adj, n, max_len=None):
    """Count directed cycles by length using backtracking."""
    if max_len is None:
        max_len = n
    counts = {}

    for length in range(3, max_len + 1, 2):  # odd cycles only
        count = 0
        # Use DFS from each starting vertex, only count cycles where start = min vertex
        for start in range(n):
            # Find all paths of given length starting at start,
            # where all internal vertices > start, ending at vertex with edge back to start
            stack = [(start, 1, 1 << start)]
            while stack:
                v, depth, mask = stack.pop()
                if depth == length:
                    # Check if there's an edge back to start
                    if start in adj[v]:
                        count += 1
                    continue
                for u in adj[v]:
                    if u == start and depth < length:
                        continue  # don't return early
                    if u < start:
                        continue  # canonical: start is minimum vertex
                    if mask & (1 << u):
                        continue  # already visited
                    stack.append((u, depth + 1, mask | (1 << u)))

        counts[length] = count

    return counts

for p in [3, 7, 11]:
    adj = paley_tournament(p)
    print(f"\n  P_{p}:")
    cycles = count_directed_cycles(adj, p)
    for k in sorted(cycles.keys()):
        c = cycles[k]
        print(f"    t_{k} = {c}, mod p = {c % p}")

print()
print("=" * 70)
print("ALL Hamiltonian path orbits under Z/pZ")
print("=" * 70)

# For small p, let's look at the orbits more carefully
for p in [3, 7]:
    adj = paley_tournament(p)
    # Generate all Hamiltonian paths
    paths = []
    for perm in permutations(range(p)):
        valid = True
        for i in range(p-1):
            if perm[i+1] not in adj[perm[i]]:
                valid = False
                break
        if valid:
            paths.append(perm)

    # Group into orbits under cyclic shift
    seen = set()
    orbits = []
    for path in paths:
        if path in seen:
            continue
        orbit = set()
        shifted = path
        for _ in range(p):
            shifted = tuple((v + 1) % p for v in shifted)
            orbit.add(shifted)
        orbits.append(orbit)
        seen.update(orbit)

    print(f"\n  P_{p}: H = {len(paths)}, orbits = {len(orbits)}")
    print(f"    Orbit sizes: {sorted([len(o) for o in orbits])}")
    if p <= 7:
        for i, orbit in enumerate(orbits[:5]):
            rep = sorted(orbit)[0]
            print(f"    Orbit {i}: representative = {rep}, size = {len(orbit)}")

print()
print("=" * 70)
print("CONGRUENCE EXPLORATION: H mod small primes")
print("=" * 70)

for p in [3, 7, 11, 19]:
    adj = paley_tournament(p)
    H = count_hamiltonian_paths(adj, p)
    print(f"\n  P_{p}: H = {H}")
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
        if q <= 30:
            print(f"    H mod {q:2d} = {H % q}")

print("\n\nDone!")
