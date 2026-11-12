#!/usr/bin/env python3
"""Verifier, claim M0: brute force (permutations only) test of the brief's formula
   H(T+v) = H(T) + sum_{a->b, a in A = N+(v), b in B = N-(v)} mu(a->b)      [claimed false by the lane]
and of the corrected single-in-neighbour law H(T + v_{u}) = 2H(T) - s(u) + D(u),
D(u) = #{ordered covers (P1 u)(y P2), both nonempty, with y -> u}."""
import itertools
from vlib import pairs_of, adjm_from_mask

def paths(n, A):
    return [p for p in itertools.permutations(range(n)) if all(A[p[i]][p[i+1]] for i in range(n-1))]

for c in range(2, 6):
    na = c*(c-1)//2
    tot = 0; bad_brief = 0; bad_single = 0; nsingle = 0; maxdiff = 0
    for m in range(1 << na):
        A = adjm_from_mask(c, m)
        P = paths(c, A)
        H = len(P)
        # arc usage
        mu = [[0]*c for _ in range(c)]
        for p in P:
            for i in range(c-1):
                mu[p[i]][p[i+1]] += 1
        s = [sum(1 for p in P if p[0] == x) for x in range(c)]
        # D(u): ordered covers (P1 ending u)(P2 starting y), y -> u
        D = [0]*c
        for perm in itertools.permutations(range(c)):
            for cut in range(1, c):
                if all(A[perm[i]][perm[i+1]] for i in range(cut-1)) and all(A[perm[i]][perm[i+1]] for i in range(cut, c-1)):
                    u = perm[cut-1]; y = perm[cut]
                    if A[y][u]:
                        D[u] += 1
        for B in range(1 << c):
            # T + v_B : x -> v for x in B, v -> y otherwise
            A2 = [row[:] + [0] for row in A] + [[0]*(c+1)]
            for x in range(c):
                if (B >> x) & 1: A2[x][c] = 1
                else: A2[c][x] = 1
            Hv = len(paths(c+1, A2))
            brief = H + sum(mu[a][b] for a in range(c) for b in range(c) if a != b and A[a][b] and not (B>>a)&1 and (B>>b)&1)
            tot += 1
            if brief != Hv:
                bad_brief += 1; maxdiff = max(maxdiff, Hv - brief)
                assert Hv > brief
            if bin(B).count('1') == 1:
                u = B.bit_length() - 1
                nsingle += 1
                if Hv != 2*H - s[u] + D[u]:
                    bad_single += 1
    print(f"core order {c}: {tot} (core, cut) pairs; brief's formula wrong in {bad_brief} (always an undercount; max deficit {maxdiff}); "
          f"single-in-neighbour law 2H - s(u) + D(u): {nsingle} cases, failures {bad_single}")
