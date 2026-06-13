#!/usr/bin/env python3
"""
Investigate the parity constraint: t_3=5 forces t_5 even at n=6.

This is the mechanism behind the H=21 gap. Let's understand WHY.

For a tournament on n=6 vertices:
  t_3 = C(6,3) - sum_v C(d_v^+, 2) = 20 - sum_v C(d_v, 2)

where d_v^+ is the out-degree of vertex v. Since sum d_v = C(6,2) = 15,
and each d_v is in {0,1,2,3,4,5}.

t_3 = 20 - sum C(d_v, 2)

For t_3=5: sum C(d_v, 2) = 15.
The score sequence d = (d_0,...,d_5) with sum=15 and sum C(d_v,2)=15.

Let's enumerate all score sequences with t_3=5 and compute t_5 for each.

Also: is there a parity law relating t_3 and t_5?

Known: t_3 + t_5 + t_7 + ... = alpha_1 (total odd cycles in Omega)
       I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...

But t_3, t_5, t_7 are individual cycle counts by length.

opus-2026-03-07-S39
"""
from itertools import combinations, permutations
from collections import defaultdict
from math import comb


def score_sequences_n6():
    """Generate all possible score sequences for n=6 tournaments."""
    # Score sequence: (d_0,...,d_5) with sum=15, 0<=d_i<=5
    seqs = []
    for d0 in range(6):
        for d1 in range(6):
            for d2 in range(6):
                for d3 in range(6):
                    for d4 in range(6):
                        d5 = 15 - d0 - d1 - d2 - d3 - d4
                        if 0 <= d5 <= 5:
                            seqs.append(tuple(sorted([d0,d1,d2,d3,d4,d5])))
    return list(set(seqs))


def t3_from_scores(scores):
    """Compute t_3 from score sequence."""
    return 20 - sum(comb(d, 2) for d in scores)


# === Score sequences with t_3=5 ===
print("=== Score sequences with t_3=5 at n=6 ===")
all_scores = score_sequences_n6()
scores_t3_5 = [s for s in all_scores if t3_from_scores(s) == 5]

print(f"Total distinct score sequences: {len(all_scores)}")
print(f"Score sequences with t_3=5: {len(scores_t3_5)}")
for s in sorted(scores_t3_5):
    print(f"  {s}: sum C(d,2) = {sum(comb(d,2) for d in s)}")

# === Now exhaustive: for each tournament with t_3=5, what is t_5? ===
print("\n=== Exhaustive t_5 for t_3=5 tournaments at n=6 ===")
n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

t3_5_t5_counts = defaultdict(int)

for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    adj_bits = [0]*n
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[j][i] = 1
            adj_bits[j] |= (1 << i)
        else:
            adj[i][j] = 1
            adj_bits[i] |= (1 << j)

    # Count t3
    t3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if adj[a][b] and adj[b][c] and adj[c][a]:
                    t3 += 1
                elif adj[a][c] and adj[c][b] and adj[b][a]:
                    t3 += 1

    if t3 != 5:
        continue

    # Count t5
    edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and adj[i][j]}
    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all((p[i], p[(i+1)%5]) in edge_set for i in range(5)):
                t5 += 1
    t5 //= 5  # divide by rotations

    t3_5_t5_counts[t5] += 1

    # Also compute score sequence
    scores = tuple(sorted(bin(adj_bits[v]).count('1') for v in range(n)))

print(f"t_5 distribution for t_3=5 tournaments:")
for t5 in sorted(t3_5_t5_counts.keys()):
    parity = "EVEN" if t5 % 2 == 0 else "ODD"
    print(f"  t_5={t5} ({parity}): {t3_5_t5_counts[t5]} tournaments")

all_even = all(t5 % 2 == 0 for t5 in t3_5_t5_counts.keys())
print(f"\nAll t_5 values are even: {all_even}")

# === Check parity of t_5 for ALL t_3 values ===
print("\n=== Parity of t_5 by t_3 at n=6 ===")
t3_to_t5_parity = defaultdict(lambda: {'even': 0, 'odd': 0})

for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[j][i] = 1
        else:
            adj[i][j] = 1

    t3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if adj[a][b] and adj[b][c] and adj[c][a]:
                    t3 += 1
                elif adj[a][c] and adj[c][b] and adj[b][a]:
                    t3 += 1

    edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and adj[i][j]}
    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all((p[i], p[(i+1)%5]) in edge_set for i in range(5)):
                t5 += 1
    t5 //= 5

    if t5 % 2 == 0:
        t3_to_t5_parity[t3]['even'] += 1
    else:
        t3_to_t5_parity[t3]['odd'] += 1

print(f"{'t3':>4s} | {'t5 even':>10s} | {'t5 odd':>10s} | {'t3+t5 parity':>12s}")
print("-" * 50)
for t3 in sorted(t3_to_t5_parity.keys()):
    ev = t3_to_t5_parity[t3]['even']
    od = t3_to_t5_parity[t3]['odd']
    # When t5 is always even, t3+t5 has same parity as t3
    t3_parity = "even" if t3 % 2 == 0 else "odd"
    forced = "ALWAYS" if ev == 0 or od == 0 else "MIXED"
    t5_parity = "even" if od == 0 else ("odd" if ev == 0 else "mixed")
    print(f"{t3:4d} | {ev:10d} | {od:10d} | t3={t3_parity}, t5={t5_parity} -> {forced}")

# Key question: is t3+t5 always even/odd for a given n=6 tournament?
print("\n=== Is t3+t5 always the same parity at n=6? ===")
t3_plus_t5_parities = set()
for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[j][i] = 1
        else:
            adj[i][j] = 1

    t3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if adj[a][b] and adj[b][c] and adj[c][a]:
                    t3 += 1
                elif adj[a][c] and adj[c][b] and adj[b][a]:
                    t3 += 1

    edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and adj[i][j]}
    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all((p[i], p[(i+1)%5]) in edge_set for i in range(5)):
                t5 += 1
    t5 //= 5

    t3_plus_t5_parities.add((t3 + t5) % 2)

print(f"Achievable (t3+t5) mod 2 values: {t3_plus_t5_parities}")
if len(t3_plus_t5_parities) == 1:
    p = list(t3_plus_t5_parities)[0]
    print(f"t3 + t5 is ALWAYS {'even' if p==0 else 'odd'} at n=6!")
    print(f"This means sw(1) + sw(2) = 2*(t3+t5) + 4*p33 is always ≡ {2*p} mod 4")
    if p == 0:
        print(f"So H = 1 + sw(1) + sw(2) ≡ 1 mod 4 always at n=6!")
    else:
        print(f"So H = 1 + sw(1) + sw(2) ≡ 3 mod 4 always at n=6!")
