#!/usr/bin/env python3
"""THM-484 / HYP-2415 (claudebox-2026-06-11-S4): the involution modulus 24, the eQR extremal
ladder and its first failure at 72, the digit-square puzzle (answer 20) and its base-dependence.
Pure python, exact. Depends on tournament_barba_flag_cbx1.py + qr_ladder helpers (inlined)."""
import sys, itertools
from math import gcd
sys.path.insert(0, '/home/claude/math-research/04-computation')
from tournament_barba_flag_cbx1 import paley_S, skewH_from_drt

def binrows(H): return [int(''.join('1' if x < 0 else '0' for x in row), 2) for row in H]
def basis_of(vecs):
    b = []
    for v in vecs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v); b.sort(reverse=True)
    return b
def aset(b):
    w = [0]
    for v in b: w += [x ^ v for x in w]
    return set(w)
def mindist(b):
    md = 1 << 30; w = 0; gp = 0
    for i in range(1, 1 << len(b)):
        g = i ^ (i >> 1); bit = (g ^ gp).bit_length() - 1; gp = g
        w ^= b[bit]; c = w.bit_count()
        if 0 < c < md: md = c
    return md
def rm(r, m):
    pts = list(itertools.product([0, 1], repeat=m)); gens = []
    for d in range(r + 1):
        for mo in itertools.combinations(range(m), d):
            v = 0
            for k, p in enumerate(pts):
                val = 1
                for i in mo: val &= p[i]
                if val: v |= (1 << k)
            gens.append(v)
    return basis_of(gens)
def equiv_perm(C1, C2, n):
    t = aset(C2)
    for perm in itertools.permutations(range(n)):
        if all((sum((1 << perm[i]) for i in range(n) if (w >> i) & 1)) in t for w in aset(C1)):
            return True
    return False

# digit-square map
def Sb(n, b):
    s = 0
    while n: s += (n % b) ** 2; n //= b
    return s
def base_cycles(b):
    cap = 3 * (b - 1) ** 2 + 5
    seen = {}; cyc = set()
    for start in range(1, cap):
        path = []; x = start
        while x not in path and x not in seen:
            path.append(x); x = Sb(x, b)
        if x in path: cyc.add(tuple(sorted(path[path.index(x):])))
        for p in path: seen[p] = True
    return sorted(len(c) for c in cyc)

if __name__ == '__main__':
    print("== Part 1: maximal involution modulus ==")
    units = lambda n: [a for a in range(1, n) if gcd(a, n) == 1]
    exp2 = lambda n: all((a * a) % n == 1 for a in units(n))
    mods = [n for n in range(1, 1000) if exp2(n)]
    print(f"  exponent-2 moduli (Z/n)*: {mods}  -> max = {max(mods)}")
    assert mods == [1, 2, 3, 4, 6, 8, 12, 24]
    print(f"  (Z/24)* = {units(24)}, phi=8, all squares -> {{ {','.join(str(a*a%24) for a in units(24))} }}")
    assert all((a * a) % 24 == 1 for a in units(24))

    print("\n== Part 2: e8 = RM(1,3) = (Z/24)*-indexed (the first Gleason generator) ==")
    e8 = basis_of(binrows(skewH_from_drt(paley_S(7))))
    print(f"  C(I+S(Paley7)) ~ RM(1,3): {equiv_perm(e8, rm(1, 3), 8)}  (8 coords = F_2^3 = (Z/24)*)")
    assert equiv_perm(e8, rm(1, 3), 8)

    print("\n== Part 3/HYP-2415: eQR extremal ladder, first failure at 72 ==")
    ext = lambda n: 4 * (n // 24) + 4
    for q in (7, 23, 31, 47):
        d = mindist(basis_of(binrows(skewH_from_drt(paley_S(q)))))
        print(f"  q={q}: eQR({q+1}) d={d}, extremal {ext(q+1)}, EXTREMAL={d==ext(q+1)}")
        assert d == ext(q + 1)
    print(f"  q=71: eQR(72) d=12 (known, QR(71)); extremal {ext(72)} => NOT extremal = the open [72,36,16]")

    print("\n== Part 4: the digit-square puzzle (answer) + base dependence ==")
    cyc = [4]
    while True:
        nx = Sb(cyc[-1], 10)
        if nx == cyc[0]: break
        cyc.append(nx)
    print(f"  base-10 unhappy cycle: {cyc}  (len {len(cyc)})")
    given = [4, 16, 37, 58, 89, 145, 42]
    ans = [v for v in cyc if v not in given]
    print(f"  PUZZLE: given {given} -> unknown start node = {ans}")
    assert ans == [20] and len(cyc) == 8
    print("  cycle lengths by base 2..12:", {b: base_cycles(b) for b in range(2, 13)})
    print("  => unique 8-cycle at bases 6 and 10 only; '8 = phi(24)' is coincidence (honest)")
    print("\nall checks passed")
