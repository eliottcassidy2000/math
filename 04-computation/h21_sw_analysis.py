#!/usr/bin/env python3
"""
H=21 gap analysis through size-weighted independence polynomial.

H(T) = sw(0) + sw(1) + sw(2) + ... = 1 + sw(1) + sw(2) + ...

For H=21: sw(1) + sw(2) + ... = 20

sw(j) = sum over independent sets S with sum (|C_i|-1)/2 = j of 2^|S|

At each j, sw(j) is a sum of terms 2^k where k = |S| (number of cycles).

Key constraint: sw(j) is always even (since 2^k >= 2 for k >= 1).
But sw(0) = 1 (always).

So H = 1 + (even number) = odd. ✓ (This is Redei's theorem!)

For H = 21: sum of even numbers = 20.

What are the possible sw(j) values at each j?

At j=1: only 3-cycles contribute. sw(1) = 2 * (# 3-cycles) = 2*t_3.
  So sw(1) is always even, and sw(1)/2 = t_3 (# directed 3-cycles).

At j=2: 5-cycles (each weight 2) and disjoint 3-cycle pairs (each weight 4).
  sw(2) = 2*t_5 + 4*p_{33}  where p_{33} = # vertex-disjoint 3-cycle pairs.

At j=3: 7-cycles (weight 2), (5,3)-pairs (weight 4), and (3,3,3)-triples (weight 8).
  sw(3) = 2*t_7 + 4*p_{53} + 8*p_{333}

For H = 21 = 1 + sw(1) + sw(2) + sw(3) + ...
Need: sw(1) + sw(2) + sw(3) + ... = 20.

Since sw(1) = 2*t_3, we need t_3 = 0, 1, 2, ..., 10 (at most).

Case t_3=0: No 3-cycles → tournament is transitive or has no cyclic triples.
  At n>=5, t_3=0 only for transitive tournament (H=1). Dead end.

Wait, actually t_3 can be 0 for non-transitive tournaments?
No — if there's any non-transitive triple, t_3 >= 1. And for n>=3,
the only tournament with t_3=0 is the transitive one, which has H=1.

Case t_3=1: sw(1) = 2. Need sw(2)+sw(3)+... = 18.
  t_3=1 means exactly one cyclic triple. This is very constrained.
  At n=4: transitive + one flip → t_3=1, H=3 or 5. Too small for H=21.
  At n=5: t_3=1 gives small H values.

Let me compute what sw(1) values are achievable at each n and what
range of H values each allows.

opus-2026-03-07-S39
"""
from itertools import combinations, permutations
from collections import defaultdict


def held_karp(n, adj):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            if dp[S][v] == 0: continue
            for u in range(n):
                if S & (1<<u): continue
                if adj[v][u]:
                    dp[S|(1<<u)][u] += dp[S][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))


def make_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A


def count_3cycles(n, A):
    """Count directed 3-cycles."""
    t3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if A[a][b] and A[b][c] and A[c][a]:
                    t3 += 1
                elif A[a][c] and A[c][b] and A[b][a]:
                    t3 += 1
    return t3


# === n=6 exhaustive: H by t_3 ===
print("=== n=6: H values grouped by t_3 (# directed 3-cycles) ===")
n = 6
t3_to_H = defaultdict(set)
H_to_t3 = defaultdict(set)
total = 1 << (n*(n-1)//2)

for bits in range(total):
    A = make_tournament(n, bits)
    H = held_karp(n, A)
    t3 = count_3cycles(n, A)
    t3_to_H[t3].add(H)
    H_to_t3[H].add(t3)

print(f"Total tournaments: {total}")
print(f"\nt_3 -> achievable H values:")
for t3 in sorted(t3_to_H.keys()):
    H_vals = sorted(t3_to_H[t3])
    print(f"  t_3={t3}: H in {H_vals}")

# Focus on H near 21
print(f"\nH values near 21:")
for h in range(17, 28):
    if h in H_to_t3:
        print(f"  H={h}: t_3 in {sorted(H_to_t3[h])}")
    else:
        print(f"  H={h}: NOT ACHIEVABLE")

# === Specifically analyze H=19 and H=23 (flanking the gap) ===
print(f"\n=== H=19 analysis ===")
for bits in range(total):
    A = make_tournament(n, bits)
    H = held_karp(n, A)
    if H == 19:
        t3 = count_3cycles(n, A)
        # Compute sw decomposition
        edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]}
        # Count 5-cycles
        t5 = 0
        for verts in combinations(range(n), 5):
            for p in permutations(verts):
                if all((p[i], p[(i+1)%5]) in edge_set for i in range(5)):
                    t5 += 1
        t5 //= 5  # Each 5-cycle counted 5 times by rotations

        # Count disjoint 3-cycle pairs
        three_cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if A[a][b] and A[b][c] and A[c][a]:
                        three_cycles.append(frozenset({a,b,c}))
                    elif A[a][c] and A[c][b] and A[b][a]:
                        three_cycles.append(frozenset({a,b,c}))
        disjoint_pairs = sum(1 for i in range(len(three_cycles))
                           for j in range(i+1, len(three_cycles))
                           if not (three_cycles[i] & three_cycles[j]))

        print(f"  bits={bits}: t_3={t3}, t_5={t5}, disjoint_3_pairs={disjoint_pairs}")
        print(f"    sw(1)={2*t3}, sw(2)={2*t5+4*disjoint_pairs}")
        print(f"    H = 1 + {2*t3} + {2*t5+4*disjoint_pairs} = {1+2*t3+2*t5+4*disjoint_pairs}")
        break  # Just one example

print(f"\n=== H=23 analysis ===")
for bits in range(total):
    A = make_tournament(n, bits)
    H = held_karp(n, A)
    if H == 23:
        t3 = count_3cycles(n, A)
        edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]}
        t5 = 0
        for verts in combinations(range(n), 5):
            for p in permutations(verts):
                if all((p[i], p[(i+1)%5]) in edge_set for i in range(5)):
                    t5 += 1
        t5 //= 5

        three_cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if A[a][b] and A[b][c] and A[c][a]:
                        three_cycles.append(frozenset({a,b,c}))
                    elif A[a][c] and A[c][b] and A[b][a]:
                        three_cycles.append(frozenset({a,b,c}))
        disjoint_pairs = sum(1 for i in range(len(three_cycles))
                           for j in range(i+1, len(three_cycles))
                           if not (three_cycles[i] & three_cycles[j]))

        print(f"  bits={bits}: t_3={t3}, t_5={t5}, disjoint_3_pairs={disjoint_pairs}")
        print(f"    sw(1)={2*t3}, sw(2)={2*t5+4*disjoint_pairs}")
        print(f"    H = 1 + {2*t3} + {2*t5+4*disjoint_pairs} = {1+2*t3+2*t5+4*disjoint_pairs}")
        break

# === Key insight: what (sw(1), sw(2)) pairs are achievable? ===
print(f"\n=== Achievable (sw(1), sw(2)) at n=6 ===")
sw_pairs = defaultdict(int)
H_by_sw = defaultdict(set)

for bits in range(total):
    A = make_tournament(n, bits)
    H = held_karp(n, A)
    t3 = count_3cycles(n, A)
    edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]}
    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all((p[i], p[(i+1)%5]) in edge_set for i in range(5)):
                t5 += 1
    t5 //= 5
    three_cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if A[a][b] and A[b][c] and A[c][a]:
                    three_cycles.append(frozenset({a,b,c}))
                elif A[a][c] and A[c][b] and A[b][a]:
                    three_cycles.append(frozenset({a,b,c}))
    dp = sum(1 for i in range(len(three_cycles))
            for j in range(i+1, len(three_cycles))
            if not (three_cycles[i] & three_cycles[j]))

    sw1 = 2*t3
    sw2 = 2*t5 + 4*dp
    sw_pairs[(sw1, sw2)] += 1
    H_by_sw[(sw1, sw2)].add(H)

# Which (sw1, sw2) can give H=21?
# H = 1 + sw1 + sw2 = 21 => sw1 + sw2 = 20
print(f"Pairs (sw1, sw2) with sw1 + sw2 = 20:")
found = False
for (s1, s2), cnt in sorted(sw_pairs.items()):
    if s1 + s2 == 20:
        print(f"  ({s1}, {s2}): {cnt} tournaments, H = {sorted(H_by_sw[(s1, s2)])}")
        found = True
if not found:
    print(f"  NONE! No tournament achieves sw1 + sw2 = 20.")

print(f"\nPairs with sw1 + sw2 near 20:")
for (s1, s2), cnt in sorted(sw_pairs.items()):
    if 16 <= s1 + s2 <= 24:
        Hs = sorted(H_by_sw[(s1, s2)])
        print(f"  ({s1}, {s2}): sum={s1+s2}, H={Hs}")
