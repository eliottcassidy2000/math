#!/usr/bin/env python3
"""Audit A: Monte Carlo on actual integers for the switch-coverage model of THM-4640 (8).
Random odd n with B bits; a family (head h, cmin, D') applies if n's word is 1^r h c ... with c >= cmin and r >= D'.
Every applicable switch is checked on the integer: U^j(n) = U^j(m), m = (n+1)/2^D' - 1, j = r + |h| + 1.
Families: canonical (D' = D) and twin (D' = D - 1) from canonical primitive pairs with exponents <= E, and the root
(reset) switch (h = (), c >= 3, D' = 1).  Usage: python3 a7b_montecarlo.py [E=16] [NS=200000] [B=400]"""
import sys, random, itertools, math
from collections import defaultdict

E = int(sys.argv[1]) if len(sys.argv) > 1 else 16
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
B = int(sys.argv[3]) if len(sys.argv) > 3 else 400
random.seed(4640)


def K_of(s):
    K = -1
    for e in s:
        K = 3 * K + (1 << e)
    return K


def head_word(rep):
    if not rep:
        return ()
    return (rep[0] + 1,) + tuple(rep[j] - rep[j - 1] for j in range(1, len(rep)))


reps = defaultdict(list)
for nn in range(E + 1):
    for s in itertools.combinations(range(2, E + 1), nn):
        reps[K_of(s)].append(s)
fam = defaultdict(list)   # head -> list of (cmin, D', kind)
for k, rs in reps.items():
    if k == -1:
        continue
    for s, t in itertools.combinations(rs, 2):
        if s[-1] == t[-1]:
            continue
        if len(s) > len(t):
            s, t = t, s
        D = len(t) - len(s); cmin = max(1, 1 + t[-1] - s[-1])
        fam[head_word(s)].append((cmin, D, "canonical"))
        fam[head_word((1,) + s)].append((cmin, D - 1, "twin"))
maxlen = max(len(h) for h in fam) + 1


def U(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n


def Uj(n, j):
    for _ in range(j):
        n = U(n)
    return n


cnt = defaultdict(int); verified = 0
for _ in range(NS):
    n = random.getrandbits(B) | 1 | (1 << (B - 1))
    r = ((n + 1) & -(n + 1)).bit_length() - 2          # v2(n+1) - 1
    x = n; w = []
    for _ in range(r + maxlen):
        y = 3 * x + 1; a = (y & -y).bit_length() - 1; x = y >> a; w.append(a)
    u = w[r:]
    kinds = set(); first = None
    if r >= 1 and u[0] >= 3:
        kinds.add("root"); first = first or (1, 0)
    for L in range(1, len(u)):
        for (cmin, Dp, kind) in fam.get(tuple(u[:L]), ()):
            if u[L] >= cmin and r >= Dp:
                kinds.add(kind)
                if first is None or kind != "root":
                    first = (Dp, L)
    for kd in kinds:
        cnt[kd] += 1
    if kinds - {"root"}:
        cnt["canonical or twin"] += 1
    if kinds:
        cnt["any"] += 1
        Dp, L = first
        m = (n + 1) // 2 ** Dp - 1
        j = r + L + 1
        assert Uj(n, j) == Uj(m, j)
        verified += 1
print(f"E <= {E}, {NS} random odd {B}-bit n (switch merges verified on all {verified} covered integers):")
for kd in ("canonical", "twin", "canonical or twin", "root", "any"):
    p = cnt[kd] / NS
    print(f"   {kd:18s} {p:.5f} +- {math.sqrt(p * (1 - p) / NS):.5f}")
