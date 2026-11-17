#!/usr/bin/env python3
"""heavy_ties_ic.py: (a) for every first creator of orders 4-9, are ALL its deletion classes of maximum count first
creators (not only the tie-broken one)?  (b) Moon's interval chains IC(p_1..p_k) (THM-4744): vertices x_0..x_N, all arcs
forward except the back arc of each block from its last vertex to its first; check H(IC(p)) = prod (2^(p_i - 1) + 1) for
every composition of N <= 8 into parts >= 2."""
import itertools
from collections import defaultdict
exec(open('nearjoin_check.py').read().split('# (1)')[0])
T = '/private/tmp/claude-501/-Users-e-Documents-GitHub-ephrepos-math/d63fb198-a11f-4364-a0ae-e072089f7b73/scratchpad/work/tooth/'
H = defaultdict(list)
for line in open(T + 'counts_z.txt'):
    n, h = map(int, line.split()); H[n].append(h)
seen = set(); first = {}
for n in range(1, 10):
    new = set(H[n]) - seen; first[n] = [h in new for h in H[n]]; seen |= set(H[n])
for n in range(4, 10):
    parents = [list(map(int, l.split())) for l in open(T + f'parents_{n}.txt')]
    fc = [i for i in range(len(H[n])) if first[n][i]]; bad = 0; ties = 0; badc = set()
    for i in fc:
        mx = max(H[n-1][p] for p in parents[i]); tops = [p for p in parents[i] if H[n-1][p] == mx]
        if len(tops) > 1: ties += 1
        if not all(first[n-1][p] for p in tops): bad += 1; badc.add(H[n][i])
    print(f"order {n}: first creators {len(fc)}, with a tie for the heaviest deletion {ties}, with some heaviest deletion not a first creator {bad} {sorted(badc) if badc else ''}")
def IC(parts):
    N = sum(parts); n = N + 1; A = [[1 if i < j else 0 for j in range(n)] for i in range(n)]; s = 0
    for p in parts:
        a, b = s, s + p; A[a][b] = 0; A[b][a] = 1; s += p
    return A
bad = tot = 0
for N in range(2, 9):
    for k in range(1, N // 2 + 1):
        for parts in itertools.product(range(2, N + 1), repeat=k):
            if sum(parts) != N: continue
            tot += 1; h = ham_all(IC(parts))[1][-1]; pred = 1
            for p in parts: pred *= 2 ** (p - 1) + 1
            if h != pred: bad += 1; print("  mismatch", parts, h, pred)
print(f"Moon's interval chains: {tot} compositions of N <= 8 into parts >= 2; H(IC(p)) = prod (2^(p_i-1) + 1) fails {bad} times")
