#!/usr/bin/env python3
"""
Proof attempt: M = (H/n)*I for circulant tournaments at odd n.

WHAT WE KNOW:
1. M is circulant (rotation automorphism)
2. M is symmetric (THM-030)
3. At odd n: tr(M) = H
4. At odd n: M[a,b] = 0 when T[a,b] = 0 AND T is position-uniform?
   Actually, this last point is only for NONHAM, not general circulant.

KEY QUESTION: What additional property of the transfer matrix forces
the off-diagonal circulant coefficients to vanish?

APPROACH: Look at the EIGENVALUES of M for circulant tournaments.
If M = (H/n)*I, then all eigenvalues are H/n.
A circulant matrix has eigenvalues λ_k = sum_j m_j * ω^{jk}.
For all λ_k to be equal, we need m_j = 0 for j ≠ 0.

Let me look at this from the consecutive-position formula perspective.
M[a,b] = sum_j (-1)^j * N(a,b,j) where N(a,b,j) counts Ham paths
with {a,b} at positions {j, j+1}.

For a circulant tournament: N(a,b,j) depends only on (b-a) mod n and j.
So M[a,b] = f((b-a) mod n) where f(d) = sum_j (-1)^j * N(d, j).

For M to be scalar: f(d) = 0 for d ≠ 0.
This means: sum_j (-1)^j * N(d, j) = 0 for all d ≠ 0.

Since n is odd: n-1 is even, so j ranges from 0 to n-2.
The alternating sum cancels iff N(d, j) is palindromic: N(d, j) = N(d, n-2-j).

We verified palindromicity ⟺ scalar M at n=5.
For circulant tournaments, the palindromicity follows from the rotation symmetry
combined with path reversal properties.

Let me verify this chain of reasoning.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def circulant_tournament(n, gen_set):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in gen_set:
                A[i][j] = 1
    return A

def ham_paths(A):
    n = len(A)
    return [p for p in permutations(range(n))
            if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]

def compute_N_circulant(A, n):
    """N(d, j) = #{Ham paths with vertices at distance d at positions j, j+1}."""
    paths = ham_paths(A)
    N = {}
    for d in range(n):
        for j in range(n-1):
            N[(d, j)] = 0

    for p in paths:
        for j in range(n-1):
            d = (p[j+1] - p[j]) % n
            N[(d, j)] += 1

    return N

# =====================================================================
print("=" * 70)
print("N(d,j) PALINDROMICITY FOR CIRCULANT TOURNAMENTS")
print("=" * 70)

for n in [5, 7]:
    gen_sets = []
    half = list(range(1, (n+1)//2))
    for mask in range(1 << len(half)):
        gs = set()
        for k, d_val in enumerate(half):
            if mask & (1 << k):
                gs.add(d_val)
            else:
                gs.add(n - d_val)
        gen_sets.append(frozenset(gs))

    print(f"\nn={n}:")
    for gs in gen_sets[:4]:  # limit for speed at n=7
        A = circulant_tournament(n, gs)
        H = len(ham_paths(A))
        N = compute_N_circulant(A, n)

        all_palindromic = True
        for d in range(1, n):
            vals = [N[(d, j)] for j in range(n-1)]
            palindromic = all(vals[j] == vals[n-2-j] for j in range(n-1))
            if not palindromic:
                all_palindromic = False
                print(f"  gen={sorted(gs)}: d={d}, N = {vals} NOT palindromic")
            alt_sum = sum((-1)**j * vals[j] for j in range(n-1))
            # For scalar M: alt_sum should be 0 for d > 0,
            # or more precisely alt_sum should be H/n for d=0 and 0 for d>0.

        if all_palindromic:
            # Show one example
            d_sample = 1
            vals = [N[(d_sample, j)] for j in range(n-1)]
            print(f"  gen={sorted(gs)}: H={H}, N(d=1,j) = {vals}, "
                  f"palindromic? Yes, alt_sum_d1={sum((-1)**j * vals[j] for j in range(n-1))}")

# =====================================================================
print()
print("=" * 70)
print("WHY PALINDROMIC? — ROTATION + REVERSAL ARGUMENT")
print("=" * 70)

# For a circulant tournament T on Z/nZ:
# Rotation σ: maps path P = (p_0,...,p_{n-1}) to (p_0+1,...,p_{n-1}+1) mod n.
# This preserves edge structure (circulant automorphism).
#
# So N(d, j) = #{paths with (p_j, p_{j+1}) having p_{j+1} - p_j ≡ d (mod n)}
# is invariant under rotation: every path counted at position j with distance d
# has a rotated image also counted at position j with distance d.
#
# Now consider PATH REVERSAL: reversing P = (p_0,...,p_{n-1}) gives
# P^rev = (p_{n-1},...,p_0). For a general tournament, reversal maps T to T^op.
# But the reversed path uses edges (p_{k+1}, p_k) instead of (p_k, p_{k+1}).
# So P^rev is a Ham path of T^op (the reversed tournament).
#
# For circulant T with gen set S, T^op has gen set {n-d : d in S} = S^c.
# T = T^op iff S is self-complementary: S = {n-d : d in S}.
# This is NOT true in general.
#
# However: maybe there's a COMBINED symmetry (rotation + reversal)
# that forces palindromicity.
#
# Reflection symmetry: r: i ↦ -i (mod n).
# This maps edge i→j to (-i)→(-j). The edge direction is (j-i) mod n → ((-j)-(-i)) mod n = (-(j-i)) mod n = (n-(j-i)) mod n.
# So reflection maps gen_set S to {n-d : d in S} = S^c.
#
# So r is an ISOMORPHISM from T(S) to T(S^c) = T^op.
# This means: T^op ≅ T via the reflection r.
#
# Now: for path P = (p_0,...,p_{n-1}) in T:
# The reflected-reversed path r(P^rev) = (-p_{n-1}, ..., -p_0) is a path in
# T^{op, reflected} = T (since r maps T^op to T).
#
# In r(P^rev): position j has vertex -p_{n-1-j}.
# The distance at position j is: (-p_{n-2-j}) - (-p_{n-1-j}) = p_{n-1-j} - p_{n-2-j}.
# Original distance at position n-2-j was: p_{n-1-j} - p_{n-2-j} (same!).
#
# Wait — original distance at position k is p_{k+1} - p_k.
# Position n-2-j: distance is p_{n-1-j} - p_{n-2-j}.
# In r(P^rev), position j has distance: -p_{n-2-j} - (-p_{n-1-j}) = p_{n-1-j} - p_{n-2-j}.
# This equals the original distance at position n-2-j!
#
# So the map P ↦ r(P^rev) sends a path with distance d at position j
# to a path with distance d at position n-2-j.
# And this map is a BIJECTION (since r and reversal are both involutions).
#
# Therefore: N(d, j) = N(d, n-2-j) for all d, j.
# This is PALINDROMICITY!

print("""
PROOF OF PALINDROMICITY:

For circulant tournament T on Z/nZ with generating set S:

1. The reflection r: i ↦ -i (mod n) is an isomorphism from T to T^op.
   (Because it maps edge distance d to n-d, which maps S to its complement.)

2. Define the map φ on Hamiltonian paths:
   φ(P) = r(P^rev), where P^rev is the path reversal.

   If P = (p_0, p_1, ..., p_{n-1}), then:
   φ(P) = (-p_{n-1}, -p_{n-2}, ..., -p_0)  (mod n)

3. φ maps Ham paths of T to Ham paths of T:
   - P^rev is a Ham path of T^op
   - r maps T^op to T isomorphically
   - So r(P^rev) is a Ham path of T

4. φ is a bijection (it's an involution: φ² = identity).

5. KEY: If P has distance d = p_{j+1} - p_j at position j,
   then φ(P) has distance d at position (n-2-j).

   This is because the distance at position k in φ(P) is:
   (-p_{n-2-k}) - (-p_{n-1-k}) = p_{n-1-k} - p_{n-2-k}
   which is the distance of P at position (n-2-k).

6. Therefore N(d, j) = N(d, n-2-j) for all d and j. QED.

COROLLARY: M = (H/n)*I for circulant tournaments at odd n.

  Since N(a,b,j) = N((b-a) mod n, j) is palindromic in j,
  and n-1 is even at odd n, the alternating sum:
  sum_j (-1)^j * N(d, j) = 0 for all d > 0.
  (Palindromic sequence of even length has zero alternating sum.)

  So M[a,b] = 0 for a ≠ b, and M[a,a] = H/n (from trace).
""")

# =====================================================================
# Verify the bijection φ explicitly
# =====================================================================
print("=" * 70)
print("EXPLICIT VERIFICATION OF φ BIJECTION")
print("=" * 70)

n = 5
for gs in [frozenset({1, 4}), frozenset({1, 2})]:
    A = circulant_tournament(n, gs)
    paths = ham_paths(A)
    H = len(paths)

    print(f"\n  gen={sorted(gs)}: H={H}")

    # Check that φ(P) is always a valid path
    phi_count = 0
    for p in paths:
        # φ(P) = (-p_{n-1}, ..., -p_0) mod n
        phi_p = tuple((-p[n-1-k]) % n for k in range(n))
        # Check if phi_p is a valid Ham path
        valid = all(A[phi_p[i]][phi_p[i+1]] == 1 for i in range(n-1))
        if valid:
            phi_count += 1
        else:
            print(f"    FAIL: P={p}, φ(P)={phi_p} not valid!")

    print(f"    φ maps all {H} paths to valid paths: {phi_count == H}")

    # Check φ is involution
    is_involution = True
    for p in paths:
        phi_p = tuple((-p[n-1-k]) % n for k in range(n))
        phi_phi_p = tuple((-phi_p[n-1-k]) % n for k in range(n))
        if phi_phi_p != p:
            is_involution = False
            break
    print(f"    φ² = id: {is_involution}")

    # Verify distance preservation
    distance_ok = True
    for p in paths:
        phi_p = tuple((-p[n-1-k]) % n for k in range(n))
        for j in range(n-1):
            d_orig = (p[j+1] - p[j]) % n  # distance at position j in P
            d_phi = (phi_p[n-2-j+1] - phi_p[n-2-j]) % n  # distance at position n-2-j in φ(P)
            # Actually: position n-2-j in φ(P) has distance:
            d_phi_at_mirror = (phi_p[(n-2-j)+1] - phi_p[n-2-j]) % n
            if d_orig != d_phi_at_mirror:
                distance_ok = False
                print(f"    DIST FAIL at j={j}: d_orig={d_orig}, d_phi@{n-2-j}={d_phi_at_mirror}")
                break
        if not distance_ok:
            break
    print(f"    Distance preservation (N(d,j) = N(d,n-2-j)): {distance_ok}")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
