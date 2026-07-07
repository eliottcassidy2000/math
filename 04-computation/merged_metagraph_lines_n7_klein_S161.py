#!/usr/bin/env python3
"""
klein-2026-07-07-S161 -- BLUE/BLACK LINE STRUCTURE of the merged metagraph, completed to n=7,
with the parity laws now THEOREMS (proofs in the reflection; verified here to n=7).

THEOREMS (one-liners from Redei's odd-H + the grid-reflection involution rho, where
rho(t)'s tournament = (relabel v -> n+1-v)(T)^op, so rho maps fiber([T]) -> fiber([T^op])):
 T1. Every unmerged iso-class fiber is ODD (tilings x |Aut| = H odd; |Aut| odd).
 T2. SC class <=> its fiber contains a gridsym tiling (rho-stable fiber, odd size => odd
     # fixed points >= 1); NS class <=> zero gridsym (a fixed tiling would force [T]=[T^op]).
     Hence: SC nodes are pure-blue or mixed, NEVER pure black; NS-merged nodes are ALWAYS
     pure black. (The n<=7 'verified' claim of CLAUDE.md, now proved for all n.)
 T3. At an SC node: #gridsym = bx + 2bs === fiber === 1 (mod 2) => blue-cross bx ODD;
     #non-gridsym = kx + 2ks === 0 => black-cross kx EVEN. At an NS-merged node: blue
     absent, fiber = 2|C| EVEN = kx + 2ks. ('Blues contribute odd amounts, blacks even.')
 T4. Tripartite: a blue line's endpoints are gridsym => both endpoint classes contain
     gridsym => never pure-black. A black line's endpoints are non-gridsym => both endpoint
     classes contain non-gridsym => never pure-blue.
 T5. #blue lines = 2^{(m+f)/2 - 1}, f = floor((n-1)/2) (= half the rho-fixed tilings;
     flip preserves gridsym and is fixed-point-free).

THIS SCRIPT: refinement-based canonicalizer (validated vs brute force at n<=5); full census
n=4..7: node types, line allocation (blue/black x self/cross) x endpoint-type pairs, fiber
decompositions, gridsym-per-fiber; tests: pure-blue-never-self-loops; self-loop scaling;
node-type counts for formula hunting.
"""
import itertools
from collections import defaultdict, Counter

def build(n):
    tiles = [(x, y) for x in range(1, n+1) for y in range(1, x) if x - y >= 2]
    tidx = {t: k for k, t in enumerate(tiles)}
    sigma = [tidx[(n+1-y, n+1-x)] for (x, y) in tiles]
    return tiles, tidx, sigma

def adj_from_tiling(tv, n, tiles):
    A = [[0]*n for _ in range(n)]
    for a in range(2, n+1): A[a-1][a-2] = 1
    for b, (x, y) in enumerate(tiles):
        if (tv >> b) & 1: A[x-1][y-1] = 1
        else: A[y-1][x-1] = 1
    return A

def canon_tournament(A, n):
    """canonical form via iterated refinement + brute force within cells; returns a bytes key."""
    # initial colors by out-degree
    col = [sum(A[i]) for i in range(n)]
    for _ in range(n):
        prof = []
        for i in range(n):
            outs = sorted(col[j] for j in range(n) if A[i][j])
            ins  = sorted(col[j] for j in range(n) if A[j][i])
            prof.append((col[i], tuple(outs), tuple(ins)))
        ranks = {p: r for r, p in enumerate(sorted(set(prof)))}
        newcol = [ranks[p] for p in prof]
        if newcol == col: break
        col = newcol
    # cells by color
    cells = defaultdict(list)
    for i in range(n): cells[col[i]].append(i)
    ordered_cells = [cells[c] for c in sorted(cells)]
    # candidate orderings: product of per-cell permutations (small in practice)
    best = None
    def rec(prefix, rest):
        nonlocal best
        if not rest:
            perm = prefix
            key = bytearray()
            for i in range(n):
                v = 0
                for j in range(n):
                    v = (v << 1) | A[perm[i]][perm[j]]
                key += v.to_bytes((n+7)//8, 'big')
            kb = bytes(key)
            if best is None or kb < best: best = kb
            return
        cell, tail = rest[0], rest[1:]
        for p in itertools.permutations(cell):
            rec(prefix + list(p), tail)
    total = 1
    for c in ordered_cells:
        f = 1
        for i in range(2, len(c)+1): f *= i
        total *= f
    if total > 40320:   # bail to coarse safety: full brute force (rare; n<=7 regular-ish)
        for perm in itertools.permutations(range(n)):
            key = bytearray()
            for i in range(n):
                v = 0
                for j in range(n):
                    v = (v << 1) | A[perm[i]][perm[j]]
                key += v.to_bytes((n+7)//8, 'big')
            kb = bytes(key)
            if best is None or kb < best: best = kb
        return best
    rec([], ordered_cells)
    return best

def analyze(n, validate=False):
    tiles, tidx, sigma = build(n)
    m = len(tiles); full = (1 << m) - 1
    canon_key = {}; sc_flag = {}
    node = [None]*(1 << m)
    gridsym = [False]*(1 << m)
    for tv in range(1 << m):
        A = adj_from_tiling(tv, n, tiles)
        c = canon_tournament(A, n)
        Aop = [[A[j][i] for j in range(n)] for i in range(n)]
        co = canon_tournament(Aop, n)
        key = min(c, co)
        node[tv] = key
        sc_flag[key] = (c == co)
        gridsym[tv] = all(((tv >> b) & 1) == ((tv >> sigma[b]) & 1) for b in range(m))
    fiber = Counter(node)
    gcount = Counter()
    for tv in range(1 << m):
        if gridsym[tv]: gcount[node[tv]] += 1
    def category(k):
        if not sc_flag[k]: return "pure-black"
        g = gcount[k]; f = fiber[k]
        return "pure-blue" if g == f else ("mixed" if g > 0 else "SC-allblack!")
    inc = defaultdict(Counter); pairstats = Counter()
    nblue = nblack = 0
    for tv in range(1 << m):
        ftv = tv ^ full
        if tv > ftv: continue
        blue = gridsym[tv]
        if blue: nblue += 1
        else: nblack += 1
        a, b = node[tv], node[ftv]
        ca, cb = category(a), category(b)
        if a == b:
            inc[a]['bs' if blue else 'ks'] += 1
            pairstats[('blue' if blue else 'black', 'SELF', ca)] += 1
        else:
            inc[a]['bx' if blue else 'kx'] += 1
            inc[b]['bx' if blue else 'kx'] += 1
            pairstats[('blue' if blue else 'black', 'CROSS', tuple(sorted((ca, cb))))] += 1
    cats = Counter(category(k) for k in fiber)
    return dict(n=n, m=m, nblue=nblue, nblack=nblack, fiber=fiber, gcount=gcount,
                category=category, inc=inc, cats=cats, pairstats=pairstats, sc=sc_flag)

if __name__ == "__main__":
    for n in (4, 5, 6, 7):
        R = analyze(n)
        m = R['m']; fx = (n-1)//2
        print(f"\n===== n={n}: m={m}; lines={1<<(m-1)} = blue {R['nblue']} + black {R['nblack']}"
              f"  [T5 predicts blue = {1 << ((m+fx)//2 - 1)}: {'OK' if R['nblue'] == 1 << ((m+fx)//2-1) else 'FAIL'}] =====")
        print(f"  merged nodes = {len(R['fiber'])}; categories: {dict(R['cats'])}")
        # theorem checks T2/T3
        viol = 0
        sl = Counter();
        for k in R['fiber']:
            c = R['category'](k); f = R['fiber'][k]; g = R['gcount'][k]
            I = R['inc'][k]; bx, kx, bs, ks = I['bx'], I['kx'], I['bs'], I['ks']
            sl[c] += bs + ks
            if c == "pure-black" and not (g == 0 and bx == 0 and bs == 0 and f % 2 == 0 and kx % 2 == 0): viol += 1
            if c in ("pure-blue", "mixed") and not (g % 2 == 1 and bx % 2 == 1 and kx % 2 == 0 and f % 2 == 1): viol += 1
            if not (f == bx + kx + 2*(bs + ks)): viol += 1
        print(f"  T2/T3 violations: {viol}   self-loop lines by cat: PB={sl['pure-black']}, MX={sl['mixed']}, PBlue={sl['pure-blue']}")
        print(f"  line-type x endpoint-type allocation:")
        for key in sorted(R['pairstats']):
            print(f"    {key}: {R['pairstats'][key]}")
        # fiber stats per category
        for c in ("pure-blue", "mixed", "pure-black"):
            fs = sorted(R['fiber'][k] for k in R['fiber'] if R['category'](k) == c)
            if fs:
                print(f"  {c:>10}: count={len(fs)}, fiber sizes min..max = {fs[0]}..{fs[-1]}, total tilings = {sum(fs)}")
        gtot = sum(R['gcount'].values())
        print(f"  gridsym tilings total = {gtot} (predict 2^((m+f)/2) = {1 << ((m+fx)//2)})")
