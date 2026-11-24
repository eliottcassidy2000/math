#!/usr/bin/env python3
"""Exact exhaustion for two-sided order-3, 5, and 7 sharp-cover symmetries."""
from collections import Counter
from itertools import combinations
from pathlib import Path
from runpy import run_path

model = run_path(str(Path(__file__).with_name("gewirtz-model.py")))
G = model["group"](model["G"])
I = tuple(range(56))
n = len(I)
check = model["check"]
check(len(G) == 80640, "group order")
derange = model["derange"]
D = [g for g in G if derange(g)]

def mul(a, b):
    return tuple(a[i] for i in b)

def power(g, p):
    out = I
    for b in bin(p)[2:]:
        out = mul(out, out)
        if b == "1":
            out = mul(g, out)
    return out

def bits(x):
    while x:
        b = x & -x
        yield b.bit_length()-1
        x ^= b

def solve(num_cells, footprints, initial):
    """Algorithm X: branch on the smallest list of candidates for an uncovered cell."""
    masks = [sum(1 << j for j in row) for row in footprints]
    incidence = [0] * num_cells
    for i, row in enumerate(footprints):
        for j in row:
            incidence[j] |= 1 << i

    def blocked(row):
        out = 0
        for j in row:
            out |= incidence[j]
        return out

    kills = [blocked(row) for row in footprints]
    ALL = (1 << num_cells) - 1
    usable = (1 << len(footprints)) - 1
    counts = Counter()
    def search(C, U, depth):
        if not U:
            return True
        counts[depth] += 1
        m = len(footprints)+1
        for b in bits(U):
            options = C & incidence[b]
            size = options.bit_count()
            if not size:
                return False
            if size < m:
                m, choices = size, options
                if size == 1:
                    break
        return any(search(C & ~kills[i], U & ~masks[i], depth+1) for i in bits(choices))

    starts = 0
    for used in initial:
        starts += 1
        if search(usable & ~blocked(used), ALL & ~sum(1 << j for j in used), 0):
            raise RuntimeError("A symmetry-invariant cover was found.")
    return starts, tuple(counts[i] for i in range(max(counts, default=-1)+1))

for p, expected in [(3,(2240,36,2,1048)),(5,(8064,10,1,628)),(7,(5760,14,0,448))]:
    elements = [g for g in G if g != I and power(g,p)==I]
    h = elements[0]
    ih = model["inverse"](h)
    H = [power(h,i) for i in range(p)]
    cent = [g for g in G if mul(g,h)==mul(h,g)]
    fixed_labels = sum(i==j for i,j in enumerate(h))
    cell_names = [min(n*k[i]+k[j] for k in H) for i in I for j in I]
    ix = {a:i for i,a in enumerate(sorted(set(cell_names)))}
    cells = [[ix[cell_names[n*i+j]] for j in I] for i in I]
    N = len(ix)
    check((len(elements),len(cent),fixed_labels,N)==expected, f"order-{p} inventory")
    check(len(elements)*len(cent)==len(G), f"single order-{p} conjugacy class")

    def orbit_footprints(pool, fix_id):
        footprints = set()
        seen = set()
        for g in pool:
            if g in seen:
                continue
            hgh = tuple(h[g[i]] for i in ih)
            if hgh == g:
                continue
            seen.add(g)
            while hgh != g:
                seen.add(hgh)
                hgh=tuple(h[hgh[i]] for i in ih)
            row = tuple(sorted(cells[i][j] for i,j in enumerate(g)))
            # A full p-orbit of permutations must use p distinct cells in each
            # quotient slot it hits (fixed slots therefore cannot be used).
            if len(set(row)) == n and all(not (fix_id and h[i]==i and h[j]==j)
                                           for i,j in enumerate(g)):
                footprints.add(row)
        return sorted(footprints)

    def fixed_cov(q):
        return {cells[i][j] for i,j in enumerate(q)}

    ID = fixed_cov(I)
    cd = [g for g in cent if derange(g)]
    def initial():
        for r in range(1, len(cd)+2):
            if (n-r)%p or (fixed_labels and r!=fixed_labels):
                continue
            for chosen in combinations(cd,r-1):
                W = set(ID)
                for q in chosen:
                    f = fixed_cov(q)
                    if not W.isdisjoint(f):
                        break
                    W.update(f)
                else:
                    yield W

    feet = orbit_footprints(D, True)
    s, trace = solve(N, feet, initial())
    expected_rows = {3:(5250,15),5:(1106,1),7:(507,64)}
    check((len(feet),s)==expected_rows[p], f"order-{p} normalized row inventory")
    print(f"p={p}: {N} cells, {len(feet)} valid nonfixed orbits, {s} normalized fixed cases; nodes {trace}")
    if p == 7:
        feet = orbit_footprints(G, False)
        check(len(feet)==1102, "order-seven no-fixed list")
        _, trace = solve(N, feet, [set()])
        print(f"p={p}: {len(feet)} nonfixed orbits without normalization; nodes {trace}")
