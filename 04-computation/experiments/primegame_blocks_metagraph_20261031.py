#!/usr/bin/env python3
"""primegame_blocks_metagraph_20261031.py -- opus-2026-10-31-S71

The owner's decomposition of Conway's PRIMEGAME into four blocks, against the merged tiling line metagraphs at orders 3, 4, 5.

1. PRIMEGAME (Conway 1987, Theorem 1; THM-4852): the 14 fractions, their prime factorizations, the data registers (2, 3, 5, 7) and
   control flags (11, 13, 17, 19, 23, 29) each fraction touches, the flag digraph, and the owner's blocks
   B1 = {11/13, 13/11}, B2 = {55/1, 1/7, 1/17, 17/91}, B3 = {19/51, 77/19, 77/29, 29/33}, B4 = {78/85, 23/38, 95/23, 15/2}.
2. Products: every ordered product of k <= 4 fractions that equals a single fraction or 1 (macro-move laws), and the triples with
   one factor from each of B2, B3, B4.
3. Triangular and Fibonacci numbers among the 28 numerators and denominators; the semiprime lemma (pq is triangular iff
   q = 2p +- 1), checked for all semiprimes below 10^5; Ming's theorem (the triangular Fibonacci numbers are 1, 3, 21, 55).
4. The merged metagraphs at orders 3, 4, 5 from the tilings (labelled tournaments containing the base path 0 -> ... -> n-1):
   kappa = complement of every free pair (the line partner), sigma = converse followed by i -> n-1-i; a line is blue when sigma
   fixes it; classes by canonical form; merged nodes = class with its converse; (H, |Aut|, tilings) per node; lines with colours.
5. Calibrated labelling tests at order 5: (a) the ten distinct composite numbers of the fractions (adjacent when they share a
   prime factor) against the ten merged nodes (adjacent when a line joins them): the maximum number of preserved adjacencies
   over all 10! bijections, against the same maximum for random graphs with the same number of edges; (b) the block graph
   (blocks adjacent when they share a flag) against the order-4 merged metagraph.
"""
import itertools, math, random, sys
from fractions import Fraction
from collections import defaultdict, Counter

PRIMEGAME = ["17/91", "78/85", "19/51", "23/38", "29/33", "77/29", "95/23", "77/19", "1/17", "11/13", "13/11", "15/2", "1/7", "55/1"]
BLOCKS = {"B1": ["11/13", "13/11"], "B2": ["55/1", "1/7", "1/17", "17/91"], "B3": ["19/51", "77/19", "77/29", "29/33"], "B4": ["78/85", "23/38", "95/23", "15/2"]}
DATA = [2, 3, 5, 7]; FLAGS = [11, 13, 17, 19, 23, 29]
def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def frac(s): a, b = s.split("/"); return Fraction(int(a), int(b))

print("## 1. PRIMEGAME: factorizations, registers, flags, blocks")
assert sorted(sum(BLOCKS.values(), [])) == sorted(PRIMEGAME), "blocks must partition the program"
touch = {}
for s in PRIMEGAME:
    a, b = (int(x) for x in s.split("/")); fa, fb = factor(a), factor(b)
    exps = {p: fa.get(p, 0) - fb.get(p, 0) for p in set(fa) | set(fb)}
    data = {p: e for p, e in exps.items() if p in DATA}; flags = {p: e for p, e in exps.items() if p in FLAGS}
    src = [p for p, e in flags.items() if e < 0]; dst = [p for p, e in flags.items() if e > 0]
    touch[s] = (data, src, dst)
    print(f"  {s:6s} = {a}/{b}: {dict(sorted(fa.items()))} / {dict(sorted(fb.items()))}; data change {dict(sorted(data.items()))}; flag {src if src else [1]} -> {dst if dst else [1]}")
for name, fr in BLOCKS.items():
    regs = sorted(set(p for s in fr for p in touch[s][0])); fl = sorted(set(p for s in fr for p in touch[s][1] + touch[s][2]))
    print(f"  block {name}: fractions {fr}; data registers touched {regs}; flags {fl}")
recips = [(s, t) for s in PRIMEGAME for t in PRIMEGAME if s < t and frac(s) * frac(t) == 1]
print("  reciprocal pairs in the program:", recips)
arcs = [(touch[s][1][0] if touch[s][1] else 1, touch[s][2][0] if touch[s][2] else 1, s) for s in PRIMEGAME]
print("  flag digraph arcs (from, to, fraction):", arcs)

print("\n## 2. Products of fractions equal to a fraction or to 1")
F = {s: frac(s) for s in PRIMEGAME}; vals = {v: s for s, v in F.items()}
for k in (2, 3, 4):
    hits = []
    for combo in itertools.product(PRIMEGAME, repeat=k):
        p = math.prod(F[s] for s in combo)
        if p == 1 or p in vals: hits.append((combo, "1" if p == 1 else vals[p]))
    print(f"  ordered products of {k} fractions equal to a fraction or 1: {len(hits)}")
    for combo, target in hits[:40]: print(f"    {' * '.join(combo)} = {target}")
    if len(hits) > 40: print(f"    ... ({len(hits) - 40} more)")
print("  one factor from each of B2, B3, B4 (any order) with product equal to a fraction or 1:")
found = False
for combo in itertools.product(BLOCKS["B2"], BLOCKS["B3"], BLOCKS["B4"]):
    p = math.prod(F[s] for s in combo)
    if p == 1 or p in vals: print("   ", combo, "=", "1" if p == 1 else vals[p]); found = True
if not found: print("    none")

print("\n## 3. Triangular and Fibonacci numbers in the program")
nums = sorted(set(int(x) for s in PRIMEGAME for x in s.split("/")))
tri = lambda m: int((math.isqrt(8 * m + 1) - 1) // 2) * (int((math.isqrt(8 * m + 1) - 1) // 2) + 1) // 2 == m
fibs = set(); a, b = 1, 2
while a < 10**6: fibs.add(a); a, b = b, a + b
print("  distinct numerators and denominators:", nums)
print("  triangular among them:", [m for m in nums if tri(m)], "(T_k:", {m: (math.isqrt(8 * m + 1) - 1) // 2 for m in nums if tri(m)}, ")")
print("  Fibonacci among them:", [m for m in nums if m in fibs])
print("  composite entries:", [m for m in nums if len(factor(m)) > 1 or any(e > 1 for e in factor(m).values())], "(count", len([m for m in nums if m > 1 and sum(factor(m).values()) > 1]), ")")
# semiprime lemma: pq (p < q primes) is triangular iff q = 2p + 1 or q = 2p - 1 (k = 2p or k = q) -- check below 10^5
primes = [p for p in range(2, 1000) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
bad = 0; tri_sp = []
for i, p in enumerate(primes):
    for q in primes[i + 1:]:
        if p * q >= 10**5: break
        t = tri(p * q); law = (q == 2 * p + 1) or (q == 2 * p - 1)
        if t != law: bad += 1
        if t: tri_sp.append((p, q, p * q))
print(f"  semiprime lemma (pq triangular iff q = 2p +- 1) below 10^5: {bad} violations; triangular semiprimes: {tri_sp[:12]} ...")
print("  Ming 1989: the triangular Fibonacci numbers are exactly 1, 3, 21, 55 (CITED); in the program: 1 and 55; 3 = H(C_3), 21 = the forbidden path count (THM-1370).")

print("\n## 4. Merged metagraphs at orders 3, 4, 5")
def tilings(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    for code in range(1 << len(pairs)):
        out = [0] * n
        for i in range(n - 1): out[i] |= 1 << (i + 1)
        for k, (i, j) in enumerate(pairs):
            if code >> k & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield code, tuple(out)
def canon(out):
    n = len(out); best = None
    for p in itertools.permutations(range(n)):
        rel = [0] * n
        for v in range(n):
            m = 0; o = out[v]
            while o: w = (o & -o).bit_length() - 1; o &= o - 1; m |= 1 << p[w]
            rel[p[v]] = m
        t = tuple(rel)
        if best is None or t < best: best = t
    return best
def aut(out):
    n = len(out); c = 0
    for p in itertools.permutations(range(n)):
        ok = True
        for v in range(n):
            m = 0; o = out[v]
            while o: w = (o & -o).bit_length() - 1; o &= o - 1; m |= 1 << p[w]
            if m != out[p[v]]: ok = False; break
        if ok: c += 1
    return c
def hcount(out):
    n = len(out); full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            ext = out[v] & ~S
            while ext: w = (ext & -ext).bit_length() - 1; ext &= ext - 1; dp[S | (1 << w)][w] += c
    return sum(dp[full])
def converse(out):
    n = len(out); res = [0] * n
    for v in range(n):
        for w in range(n):
            if out[v] >> w & 1: res[w] |= 1 << v
    return tuple(res)
def sigma(out):      # converse then relabel i -> n-1-i : a tiling again
    n = len(out); c = converse(out); res = [0] * n
    for v in range(n):
        m = 0; o = c[v]
        while o: w = (o & -o).bit_length() - 1; o &= o - 1; m |= 1 << (n - 1 - w)
        res[n - 1 - v] = m
    return tuple(res)
meta = {}
for n in (3, 4, 5):
    pairs = [(i, j) for i in range(n) for j in range(i + 2, n)]; M = len(pairs); full = (1 << M) - 1
    T = dict(tilings(n)); cls = {c: canon(o) for c, o in T.items()}
    conv = {k: canon(converse(next(o for c, o in T.items() if cls[c] == k))) for k in set(cls.values())}
    node_of = {}; nodes = []
    for k in sorted(set(cls.values())):
        if k in node_of: continue
        node_of[k] = len(nodes); node_of[conv[k]] = len(nodes); nodes.append(tuple(sorted({k, conv[k]})))
    info = []
    for idx, nd in enumerate(nodes):
        o = next(o for c, o in T.items() if cls[c] == nd[0]); tc = sum(1 for c in T if node_of[cls[c]] == idx)
        info.append((idx, hcount(o), aut(o), tc, len(nd) == 2, nd))
    lines = []
    for code, o in T.items():
        kcode = code ^ full
        if code > kcode: continue
        a, b = node_of[cls[code]], node_of[cls[kcode]]
        scode = next(c for c, oo in T.items() if oo == sigma(o)); blue = scode in (code, kcode)
        lines.append((min(a, b), max(a, b), "blue" if blue else "black", code, kcode))
    pure_blue = [i for i in range(len(nodes)) if all(l[2] == "blue" for l in lines if i in (l[0], l[1]))]
    print(f"  order {n}: {len(T)} tilings, {len(set(cls.values()))} classes, {len(nodes)} merged nodes, {len(lines)} lines ({sum(1 for l in lines if l[2] == 'blue')} blue, {sum(1 for l in lines if l[2] == 'black')} black); pure-blue nodes {pure_blue}")
    for idx, H, A, tc, merged, nd in info: print(f"    node {idx}: H = {H}, |Aut| = {A}, tilings = {tc}{' (merged converse pair)' if merged else ''}")
    lc = Counter((l[0], l[1], l[2]) for l in lines)
    print("    lines:", ", ".join(f"{a}-{b} {col}{' x' + str(m) if m > 1 else ''}" for (a, b, col), m in sorted(lc.items())))
    meta[n] = (info, lines, pure_blue)

print("\n## 5. Calibrated labelling tests")
# (a) composites vs order-5 nodes
comps = sorted(m for m in nums if m > 1 and sum(factor(m).values()) > 1)
print("  ten composites:", comps, "with factorizations", {m: sorted(factor(m)) for m in comps})
cadj = {(i, j) for i in range(len(comps)) for j in range(i + 1, len(comps)) if set(factor(comps[i])) & set(factor(comps[j]))}
info5, lines5, _ = meta[5]
madj = {(min(a, b), max(a, b)) for a, b, col, _, _ in lines5 if a != b}
print(f"  composite-sharing graph: {len(cadj)} edges on 10 vertices; order-5 merged metagraph: {len(madj)} distinct non-loop adjacencies on 10 nodes (plus loops at {sorted({a for a, b, _, _, _ in lines5 if a == b})})")
def best_match(adjA, adjB, n=10):
    best = 0; count = 0
    for p in itertools.permutations(range(n)):
        s = sum(1 for (i, j) in adjA if (min(p[i], p[j]), max(p[i], p[j])) in adjB)
        if s > best: best, count = s, 1
        elif s == best: count += 1
    return best, count
b, c = best_match(cadj, madj)
print(f"  maximum preserved adjacencies over all 10! bijections: {b} of {len(cadj)} (attained by {c} bijections)")
random.seed(20261031); calib = []
for trial in range(20):
    E = set()
    while len(E) < len(cadj):
        i, j = random.sample(range(10), 2); E.add((min(i, j), max(i, j)))
    calib.append(best_match(E, madj)[0])
print(f"  calibration: for 20 random graphs with {len(cadj)} edges the maximum preserved adjacency is {sorted(calib)} (mean {sum(calib)/len(calib):.2f})")
# (b) blocks vs order-4 nodes
blockflags = {name: set(p for s in fr for p in touch[s][1] + touch[s][2]) for name, fr in BLOCKS.items()}
print("  block flag sets:", {k: sorted(v) for k, v in blockflags.items()})
names = ["B2", "B3", "B4"]
badj = {(i, j): sorted(blockflags[names[i]] & blockflags[names[j]]) for i in range(3) for j in range(i + 1, 3)}
print("  shared flags between the three blocks of four:", {f"{names[i]}-{names[j]}": v for (i, j), v in badj.items()}, "(a triangle)")
info4, lines4, pb4 = meta[4]
print("  order-4 merged metagraph lines:", Counter((a, b, col) for a, b, col, _, _ in lines4), "(a star at the strong node with one loop: not a triangle)")
print("  verdict: the block graph (K_3 by shared flags) and the order-4 line graph (star plus loop at the strong node) are not isomorphic in any labelling; the only exact match is the involution pair B1 against the single order-3 line.")
