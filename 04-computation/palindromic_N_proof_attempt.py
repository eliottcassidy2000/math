#!/usr/bin/env python3
"""
WHY does position-uniform imply palindromic N(a,b,j)?

Key insight: PATH REVERSAL.

If P = (v_0, ..., v_{n-1}) is a Ham path in T, then P_rev = (v_{n-1}, ..., v_0)
is a Ham path in T^op.

N_T(a,b,j) counts paths where {a,b} is at positions {j,j+1}.
N_{T^op}(a,b,j) counts paths where {a,b} is at positions {j,j+1} in T^op.

By reversal: N_T(a,b,j) = N_{T^op}(a,b,n-2-j). (THM-051 proof)

So palindromic N_T(a,b,j) = N_T(a,b,n-2-j) is equivalent to:
    N_T(a,b,j) = N_{T^op}(a,b,j) for all j.

QUESTION: When does N_T = N_{T^op}?

Hypothesis: If T is position-uniform AND T^op is position-uniform
(which it always is, since T^op has the same H and reversed position matrix),
then N_T = N_{T^op}.

But this can't be right in general — two tournaments with same H can have
different N values.

Better hypothesis: If T is SELF-COMPLEMENTARY (T isomorphic to T^op),
then N_T(a,b,j) = N_{T^op}(sigma(a),sigma(b),j) = N_T(a,b,n-2-j)
via path reversal.

Wait — the isomorphism maps vertices. Let sigma be the isomorphism T -> T^op.
Then N_{T^op}(a,b,j) = N_T(sigma^{-1}(a), sigma^{-1}(b), j).
Combined with reversal: N_T(a,b,j) = N_T(sigma^{-1}(a), sigma^{-1}(b), n-2-j).

This gives palindromic N for the PAIR {a,b} only if sigma fixes {a,b}.
In general it gives palindromic N if we SUM over the orbit of {a,b} under sigma.

ALTERNATIVE: For VERTEX-TRANSITIVE T, Aut(T) acts transitively on vertices.
If T^op is isomorphic to T via some sigma in Aut(T), AND sigma acts transitively
on pairs, then palindromic N follows.

Actually, let's think more carefully. For circulant T on Z/nZ:
- T^op is also circulant with generator set {n-d : d in S} where S is T's generator.
- If T is self-complementary: there exists r with r*S = {n-d : d in S} (mod n).

For Paley T_p: S = QR(p), and -1 is a QR iff p = 1 mod 4.
Actually for Paley at p=1 mod 4: S = QR, and S^op = {p-d : d in QR} = {-d : d in QR}.
Since -1 is a QR (p=1 mod 4), -QR = QR. So S^op = S, meaning T = T^op (not just isomorphic, EQUAL).

For p = 3 mod 4: -1 is NOT a QR. Then S^op != S. But the map i -> -i sends S to S^op.
So sigma: i -> -i mod p gives the isomorphism T -> T^op.

So for ALL Paley tournaments: T is self-complementary (either T=T^op or via i->-i).

For general circulant T: T may or may not be self-complementary.

Let me test: among the n=5 position-uniform tournaments, are ALL self-complementary?

kind-pasteur-2026-03-06-S25d
"""

from itertools import permutations
import numpy as np

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T

def opposite_tournament(T, n):
    Tp = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            Tp[(i,j)] = 1 - T[(i,j)]
    return Tp

def are_isomorphic(T1, T2, n):
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if T1[(perm[i], perm[j])] != T2[(i,j)]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True, perm
    return False, None

def ham_paths_list(T, n):
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                ok = False; break
        if ok:
            paths.append(perm)
    return paths

def position_matrix(paths, n):
    P = np.zeros((n, n), dtype=int)
    for perm in paths:
        for k in range(n):
            P[perm[k], k] += 1
    return P

# ============================================================
# n=5: Are ALL position-uniform T self-complementary?
# ============================================================
print("=" * 70)
print("n=5: Position-uniform => Self-complementary?")
print("=" * 70)

n = 5
pairs_list = [(i,j) for i in range(n) for j in range(i+1, n)]
uniform_selfcomp = 0
uniform_not_selfcomp = 0
uniform_count = 0

for bits in range(1 << len(pairs_list)):
    T = tournament_from_bits(n, bits)
    paths = ham_paths_list(T, n)
    H = len(paths)
    if H % n != 0:
        continue
    P = position_matrix(paths, n)
    if not all(P[v,k] == H // n for v in range(n) for k in range(n)):
        continue

    uniform_count += 1
    Top = opposite_tournament(T, n)
    is_iso, sigma = are_isomorphic(T, Top, n)

    if is_iso:
        uniform_selfcomp += 1
    else:
        uniform_not_selfcomp += 1
        if uniform_not_selfcomp <= 3:
            scores = sorted([sum(T[(i,j)] for j in range(n) if j != i) for i in range(n)], reverse=True)
            print(f"  NOT self-comp: bits={bits}, H={H}, scores={scores}")

print(f"\n  Position-uniform: {uniform_count}")
print(f"  Self-complementary: {uniform_selfcomp}")
print(f"  NOT self-comp: {uniform_not_selfcomp}")

if uniform_not_selfcomp == 0:
    print("\n  ==> ALL position-uniform n=5 are self-complementary!")
else:
    print("\n  ==> Some position-uniform are NOT self-complementary")

# ============================================================
# Now check: for self-comp T at odd n, does isomorphism sigma
# force palindromic N via path reversal?
# ============================================================
print("\n" + "=" * 70)
print("Self-complementary + path reversal => palindromic N?")
print("=" * 70)

# Take a position-uniform tournament
for bits in range(1 << len(pairs_list)):
    T = tournament_from_bits(n, bits)
    paths = ham_paths_list(T, n)
    H = len(paths)
    if H != 15:  # Only Paley-like
        continue
    P = position_matrix(paths, n)
    if not all(P[v,k] == H // n for v in range(n) for k in range(n)):
        continue

    Top = opposite_tournament(T, n)
    is_iso, sigma = are_isomorphic(T, Top, n)
    if not is_iso:
        continue

    print(f"\n  bits={bits}, H={H}, sigma={sigma}")

    # N_T(a,b,j) via paths in T
    N_T = [[[0]*(n-1) for _ in range(n)] for _ in range(n)]
    for perm in paths:
        for j in range(n-1):
            a, b = perm[j], perm[j+1]
            N_T[a][b][j] += 1
            N_T[b][a][j] += 1

    # Check N_T(a,b,j) = N_T(sigma^{-1}(a), sigma^{-1}(b), n-2-j)?
    # sigma: T -> T^op, so sigma^{-1}: T^op -> T
    # N_{T^op}(a,b,j) = N_T(sigma^{-1}(a), sigma^{-1}(b), j) by iso
    # N_T(a,b,j) = N_{T^op}(a,b,n-2-j) by reversal
    # => N_T(a,b,j) = N_T(sigma^{-1}(a), sigma^{-1}(b), n-2-j)

    sigma_inv = [0]*n
    for i in range(n):
        sigma_inv[sigma[i]] = i

    holds = True
    for a in range(n):
        for b in range(a+1, n):
            for j in range(n-1):
                sa = sigma_inv[a]
                sb = sigma_inv[b]
                # N_T(a,b,j) should equal N_T({sa,sb}, n-2-j)
                # Since N is symmetrized: N[a][b] = N[b][a]
                expected = N_T[sa][sb][n-2-j]
                actual = N_T[a][b][j]
                if actual != expected:
                    holds = False

    print(f"  N_T(a,b,j) = N_T(sigma^-1(a), sigma^-1(b), n-2-j)? {holds}")

    # If sigma acts transitively on pairs AND the identity holds,
    # then summing over the orbit of {a,b} gives palindromic N.
    # But we need palindromic for EACH pair, not just orbit sums.

    # Actually, for palindromic N for pair {a,b}, we need:
    # N_T(a,b,j) = N_T(a,b,n-2-j)
    # From the identity: N_T(a,b,j) = N_T(sigma^{-1}(a), sigma^{-1}(b), n-2-j)
    # This gives palindromic iff sigma fixes every pair {a,b} or
    # N_T is constant on orbits of pairs under sigma.

    # Check which pairs sigma fixes:
    fixed_pairs = []
    for a in range(n):
        for b in range(a+1, n):
            sa, sb = sigma_inv[a], sigma_inv[b]
            pair = (min(sa,sb), max(sa,sb))
            if pair == (a,b):
                fixed_pairs.append((a,b))

    print(f"  sigma^-1 fixes {len(fixed_pairs)}/{n*(n-1)//2} pairs: {fixed_pairs}")

    # For non-fixed pairs: do orbits give the same N values?
    # Pair {a,b} maps to {sigma^-1(a), sigma^-1(b)}.
    # N_T(a,b,j) = N_T(sigma^-1(a), sigma^-1(b), n-2-j)
    # For palindromic: N_T(a,b,j) = N_T(a,b, n-2-j)
    # This requires sigma^-1({a,b}) = {a,b} for each pair.
    # OR: the orbit of {a,b} under sigma^-1 is size 2, and the
    # two pairs in the orbit have the SAME N values.

    # Check this:
    for a in range(n):
        for b in range(a+1, n):
            sa, sb = sigma_inv[a], sigma_inv[b]
            pair2 = (min(sa,sb), max(sa,sb))
            nab = [N_T[a][b][j] for j in range(n-1)]
            nab2 = [N_T[pair2[0]][pair2[1]][j] for j in range(n-1)]
            if nab != nab2:
                print(f"  ({a},{b}) maps to {pair2}: N differs!")
                print(f"    N({a},{b}) = {nab}")
                print(f"    N{pair2} = {nab2}")

    break  # One example

# ============================================================
# KEY INSIGHT: For vertex-transitive T at odd n,
# sigma acts on pairs. If the orbit structure guarantees
# palindromic N, we're done.
# ============================================================
print("\n" + "=" * 70)
print("Orbit analysis of sigma on vertex pairs")
print("=" * 70)

for bits in range(1 << len(pairs_list)):
    T = tournament_from_bits(n, bits)
    paths = ham_paths_list(T, n)
    H = len(paths)
    if H != 15:
        continue
    P = position_matrix(paths, n)
    if not all(P[v,k] == H // n for v in range(n) for k in range(n)):
        continue

    Top = opposite_tournament(T, n)
    is_iso, sigma = are_isomorphic(T, Top, n)
    if not is_iso:
        continue

    sigma_inv = [0]*n
    for i in range(n):
        sigma_inv[sigma[i]] = i

    print(f"\n  bits={bits}, sigma={sigma}, sigma_inv={sigma_inv}")

    # Find orbits of pairs under sigma_inv
    visited = set()
    orbits = []
    for a in range(n):
        for b in range(a+1, n):
            if (a,b) in visited:
                continue
            orbit = []
            ca, cb = a, b
            while True:
                pair = (min(ca,cb), max(ca,cb))
                if pair in visited:
                    break
                visited.add(pair)
                orbit.append(pair)
                ca, cb = sigma_inv[ca], sigma_inv[cb]
            orbits.append(orbit)

    N_T = [[[0]*(n-1) for _ in range(n)] for _ in range(n)]
    for perm in paths:
        for j in range(n-1):
            a, b = perm[j], perm[j+1]
            N_T[a][b][j] += 1
            N_T[b][a][j] += 1

    for orbit in orbits:
        nvals = []
        for (a,b) in orbit:
            nvals.append([N_T[a][b][j] for j in range(n-1)])
        print(f"  Orbit {orbit}: N values = {nvals}")
        # Check: are these the REVERSE of each other (via sigma chain)?
        if len(orbit) == 1:
            print(f"    Fixed pair: palindromic = {nvals[0] == nvals[0][::-1]}")
        elif len(orbit) == 2:
            print(f"    N[0] reversed = {nvals[0][::-1]}")
            print(f"    N[1] = {nvals[1]}")
            print(f"    Match? {nvals[0][::-1] == nvals[1]}")

    break

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
