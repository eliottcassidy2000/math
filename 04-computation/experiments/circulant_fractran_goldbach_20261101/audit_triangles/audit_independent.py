#!/usr/bin/env python3
"""Independent audit of 05-knowledge/results/quadratic_2916_odd_pair_triangles_n4_20261102.md.

Own code, Python 3 stdlib only; imports nothing from the audited scripts. The only file read is the audited
quadratic_triangle_n4.out, to compare its printed triangles with an independent recomputation.

A  x -> x^2 - 29/16 on x = a/4, a odd (cycle, oddness, escape lemma, image chain, fibres, depths, mod-16 descent,
   automorphisms of f on B).
B  n = 4 tiling metagraph (owner gauge: path 4->3->2->1, tiles (4,1),(3,1),(4,2)): classes, fibres, converse merge,
   flip lines, class line multigraph vs f on f(B), natural Klein four-groups on both sides.
C  P(m) = least prime = +-1 (mod m), even m in [4, 20000]; residue pairs R(m) vs Goldbach pairs G(m); witness types.
D  clique-carrier formula of chart_clique_norm_20261103 at r = 3, t = (3,3,5).
E  odd-pair triangle and PRIMEGAME trial-division rows, compared with the audited .out.

usage: python3 audit_independent.py
"""
import itertools
import math
import os
import sys
from bisect import bisect_left, bisect_right
from collections import Counter
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
BAD = []


def check(name, ok, detail=''):
    print(f'  [{"OK" if ok else "XX"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        BAD.append(name)


# ---------------------------------------------------------------- A
def fq(x):
    return x * x - Fr(29, 16)


def g(a):
    return (a * a - 29) // 4


def part_A():
    print('A  x -> x^2 - 29/16 on odd/4')
    cyc = [Fr(-7, 4), Fr(5, 4), Fr(-1, 4)]
    check('A1 3-cycle -7/4 -> 5/4 -> -1/4 -> -7/4', all(fq(cyc[i]) == cyc[(i + 1) % 3] for i in range(3)))
    odd_ok = all((a * a - 29) % 4 == 0 and ((a * a - 29) // 4) % 2 == 1 for a in range(-20001, 20002, 2))
    check('A2 g(a) = (a^2-29)/4 is an odd integer for odd a (|a| <= 20001; proof: a^2 = 1 mod 8)', odd_ok)
    # escape lemma: g(a) - |a| = ((|a|-2)^2 - 33)/4 >= 4 for |a| >= 9
    lem = all(4 * (g(a) - abs(a)) == (abs(a) - 2) ** 2 - 33 for a in range(-1001, 1002, 2))
    lem_bound = all(g(a) - abs(a) >= 4 for a in range(-200001, 200002, 2) if abs(a) >= 9)
    check('A3 escape lemma: 4(g(a)-|a|) = (|a|-2)^2 - 33, so g(a) > |a| (>= 9) for odd |a| >= 9', lem and lem_bound,
          'then the orbit is strictly increasing: unbounded')
    Bn = [-7, -5, -3, -1, 1, 3, 5, 7]
    check('A4 closure: g maps {+-1,+-3,+-5,+-7} into itself (so those 8 orbits are bounded)',
          all(g(a) in Bn for a in Bn))
    # brute force: iterate until |a| >= 9 (escapes by A3) or a repeat
    brute = []
    for a in range(-200001, 200002, 2):
        seen, x = set(), a
        while abs(x) <= 7 and x not in seen:
            seen.add(x)
            x = g(x)
        if abs(x) <= 7:
            brute.append(a)
    check('A5 bounded iff |a| <= 7 (odd |a| <= 200001, via A3)', brute == Bn)
    B = [Fr(a, 4) for a in Bn]
    fB = sorted(set(map(fq, B)))
    f2B = sorted(set(map(fq, fB)))
    f3B = sorted(set(map(fq, f2B)))
    sh1 = sorted(Counter(fq(x) for x in B).values(), reverse=True)
    sh2 = sorted(Counter(fq(x) for x in fB).values(), reverse=True)
    sh12 = sorted(Counter(fq(fq(x)) for x in B).values(), reverse=True)
    pairs_ok = all(fq(x) == fq(-x) for x in B) and len(fB) == 4
    check('A6 |B|,|f B|,|f^2 B|,|f^3 B| = 8,4,3,3; fibres (2,2,2,2) = {x,-x}, then (2,1,1)',
          (len(B), len(fB), len(f2B), len(f3B)) == (8, 4, 3, 3) and sh1 == [2, 2, 2, 2] and sh2 == [2, 1, 1]
          and pairs_ok and fB == sorted([Fr(-7, 4), Fr(-5, 4), Fr(-1, 4), Fr(5, 4)]) and set(f2B) == set(cyc),
          f'composite f^2 fibres on B: {sh12}')
    depth = {}
    for x in B:
        d, y = 0, x
        while y not in cyc:
            y, d = fq(y), d + 1
        depth[x] = d
    want = {Fr(-7, 4): 0, Fr(5, 4): 0, Fr(-1, 4): 0, Fr(-5, 4): 1, Fr(1, 4): 1, Fr(7, 4): 1, Fr(3, 4): 2, Fr(-3, 4): 2}
    check('A7 depths 0 {-7/4,5/4,-1/4}, 1 {-5/4,1/4,7/4}, 2 {+-3/4}', depth == want)
    # mod 16: B is a full set of odd residues (true of ANY 8 consecutive odd integers); g does not descend mod 16
    full16 = sorted(a % 16 for a in Bn) == [1, 3, 5, 7, 9, 11, 13, 15]
    generic16 = all(sorted((s + 2 * i) % 16 for i in range(8)) == [1, 3, 5, 7, 9, 11, 13, 15]
                    for s in range(-101, 101, 2))
    nodesc = all((g(a + 16) - g(a)) % 16 == 8 for a in range(-999, 1000, 2))
    check('A8 B = units mod 16 (GENERIC: any 8 consecutive odd integers); g(a+16) = g(a) + 8 (mod 16), so f does '
          'not act on residues mod 16', full16 and generic16 and nodesc)
    # canonical 2-adic cube: b = (a-1)/2 mod 8; negation a -> -a is bitwise complement b -> 7 - b
    code = {a: ((a - 1) // 2) % 8 for a in Bn}
    neg_is_flip = all(code[-a] == 7 - code[a] for a in Bn)
    tri = all(g(2 * b + 1) == 2 * (b * (b + 1) // 2 - 4) + 1 for b in range(-500, 500))
    f2fib = {}
    for a in Bn:
        f2fib.setdefault(g(g(a)), set()).add(code[a])
    bit = lambda c, i: (c >> i) & 1
    sw02 = lambda c: (c & 2) | (bit(c, 0) << 2) | bit(c, 2)
    v4orb = []
    for c in range(8):
        o = frozenset({c, 7 - c, sw02(c), 7 - sw02(c)})
        if o not in v4orb:
            v4orb.append(o)
    check('A9 2-adic cube on B (b = (a-1)/2 mod 8): -a <-> complement; g is b -> b(b+1)/2 - 4 (triangular); '
          'f^2-fibres = orbits of <complement, swap(b0,b2)> (an ad hoc swap)',
          neg_is_flip and tri and sorted(map(sorted, f2fib.values())) == sorted(map(sorted, v4orb)))
    # automorphisms of the functional graph of f on B
    idx = {x: i for i, x in enumerate(B)}
    fmap = [idx[fq(x)] for x in B]
    auts = [p for p in itertools.permutations(range(8)) if all(p[fmap[i]] == fmap[p[i]] for i in range(8))]
    tau = [p for p in auts if list(p) != list(range(8))]
    moved = [B[i] for i in range(8) if tau and tau[0][i] != i]
    check('A10 Aut(f on B) = Z/2, the swap 3/4 <-> -3/4 (negation is a deck symmetry f(-x) = f(x), not an automorphism)',
          len(auts) == 2 and sorted(moved) == [Fr(-3, 4), Fr(3, 4)])
    # credit check: c = -29/16 is the tau = 1 member of the classical (Walde-Russo 1994) rational 3-cycle family
    def wr(t):
        c = -(t ** 6 + 2 * t ** 5 + 4 * t ** 4 + 8 * t ** 3 + 9 * t ** 2 + 4 * t + 1) / (4 * t ** 2 * (t + 1) ** 2)
        x = (t ** 3 + 2 * t ** 2 + t + 1) / (2 * t * (t + 1))
        return c, x
    fam_ok = True
    for t in (Fr(1), Fr(2), Fr(1, 2), Fr(3), Fr(-2, 3), Fr(5, 7), Fr(-3), Fr(7, 2)):
        c, x = wr(t)
        y = x
        for _ in range(3):
            y = y * y + c
        fam_ok &= (y == x) and (x * x + c != x)
    check('A11 3-cycle family c(t) = -(t^6+2t^5+4t^4+8t^3+9t^2+4t+1)/(4t^2(t+1)^2): period 3 at 8 sample t; '
          't = 1 gives c = -29/16 with cycle point 5/4', fam_ok and wr(Fr(1)) == (Fr(-29, 16), Fr(5, 4)))
    # in the 2-adic cube, multiplication by the 2-torsion units {1,7,9,15} mod 16 is translation by {0,3,4,7}
    shifts = {}
    for u in (1, 7, 9, 15):
        s = {code[a] ^ (((u * a - 1) // 2) % 8) for a in Bn}
        shifts[u] = s.pop() if len(s) == 1 else None
    check('A12 2-adic cube: x7, x9, x15 (mod 16) act on B as translations by 011, 100, 111 (free; f-incompatible by A8)',
          shifts == {1: 0, 7: 3, 9: 4, 15: 7}, f'{shifts}')
    neg = [idx[-x] for x in B]
    t0 = list(tau[0])
    comp = [neg[t0[i]] for i in range(8)]
    fix = lambda p: sum(1 for i in range(8) if p[i] == i)
    mult = {}
    for u in (7, 9, 15):
        perm = []
        for a in Bn:
            r = (u * a) % 16
            r = r - 16 if r > 7 else r
            perm.append(Bn.index(r))
        mult[u] = fix(perm)
    return dict(fB=fB, f2B=f2B, sh1=sh1, sh2=sh2, sh12=sh12, fixB=sorted([fix(neg), fix(t0), fix(comp)]),
                fix16=sorted(mult.values()))


# ---------------------------------------------------------------- B
def part_B(A):
    print('B  n = 4 tiling metagraph (path 4->3->2->1; tiles (4,1),(3,1),(4,2); bit 1 = arc y -> x)')
    n = 4
    tiles = [(4, 1), (3, 1), (4, 2)]
    path = [(4, 3), (3, 2), (2, 1)]
    perms = list(itertools.permutations(range(1, n + 1)))

    def arcs(mask):
        s = set(path)
        for i, (x, y) in enumerate(tiles):
            s.add((y, x) if (mask >> i) & 1 else (x, y))
        return s

    def canon(s):
        best = None
        for p in perms:
            rl = {v: p[v - 1] for v in range(1, n + 1)}
            key = tuple(sorted((rl[u], rl[v]) for (u, v) in s))
            if best is None or key < best:
                best = key
        return best

    def scores(s):
        return tuple(sorted(sum(1 for (u, v) in s if u == w) for w in range(1, n + 1)))

    T = {m: arcs(m) for m in range(8)}
    cls = {m: canon(T[m]) for m in range(8)}
    keys = sorted(set(cls.values()))
    ci = {m: keys.index(cls[m]) for m in range(8)}
    fib = Counter(ci.values())
    sc = {ci[m]: scores(T[m]) for m in range(8)}
    conv = {ci[m]: keys.index(canon({(v, u) for (u, v) in T[m]})) for m in range(8)}
    merged = []
    for c in range(len(keys)):
        node = frozenset({c, conv[c]})
        if node not in merged:
            merged.append(node)
    mfib = sorted((sum(fib[c] for c in node) for node in merged), reverse=True)
    check('B1 8 tilings -> 4 classes, fibres {5,1,1,1}; converse merge -> 3 nodes, fibres {5,2,1}',
          len(keys) == 4 and sorted(fib.values(), reverse=True) == [5, 1, 1, 1] and mfib == [5, 2, 1],
          f'class: scores/fibre {[(c, sc[c], fib[c]) for c in sorted(fib)]}; converse {conv}')
    strong = [c for c in sc if sc[c] == (1, 1, 2, 2)][0]
    lines = [(m, m ^ 7) for m in range(8) if m < (m ^ 7)]
    cl_edges = Counter(tuple(sorted((ci[a], ci[b]))) for a, b in lines)
    star = all(strong in e for e in cl_edges) and cl_edges[(strong, strong)] == 1 and len(lines) == 4
    check('B2 4 flip lines {t, t+111}; class line multigraph = star at the strong class with one loop there', star,
          f'{dict(cl_edges)}')
    fB = A['fB']
    fe = Counter(tuple(sorted((fB.index(x), fB.index(fq(x))))) for x in fB)
    iso = []
    for p in itertools.permutations(range(4)):
        img = Counter(tuple(sorted((p[u], p[v]))) for (u, v), k in cl_edges.items() for _ in range(k))
        if img == fe:
            iso.append(p)
    check('B3 no bijection classes -> f(B) carries the line multigraph onto undirected f on f(B) (loops kept)',
          not iso, f'f on f(B): {dict(fe)} (triangle + pendant, no loop) on {[str(x) for x in fB]}')
    mnode = {c: [i for i, node in enumerate(merged) if c in node][0] for c in range(4)}
    me = Counter(tuple(sorted((mnode[ci[a]], mnode[ci[b]]))) for a, b in lines)
    f2B = A['f2B']
    f2e = Counter(tuple(sorted((f2B.index(x), f2B.index(fq(x))))) for x in f2B)
    iso3 = any(Counter(tuple(sorted((p[u], p[v]))) for (u, v), k in me.items() for _ in range(k)) == f2e
               for p in itertools.permutations(range(3)))
    check('B4 merged level: line multigraph on 3 nodes vs f on the 3-cycle: not isomorphic', not iso3,
          f'merged lines {dict(me)}; cycle {dict(f2e)}')
    # grid transpose (x,y) -> (n-y+1, n-x+1) and relabel v -> n+1-v: the converse, on tilings
    tr = [tiles.index((n - y + 1, n - x + 1)) for (x, y) in tiles]
    trans = lambda m: sum(((m >> i) & 1) << tr[i] for i in range(3))
    tr_ok = all(ci[trans(m)] == conv[ci[m]] for m in range(8))
    check('B5 grid transpose = converse on tilings (swaps tiles (3,1),(4,2), fixes (4,1))',
          tr_ok and tr == [0, 2, 1])
    # Klein four-group <flip, transpose> on tilings
    orbits = []
    for m in range(8):
        o = frozenset({m, m ^ 7, trans(m), trans(m) ^ 7})
        if o not in orbits:
            orbits.append(o)
    osh = sorted(map(len, orbits), reverse=True)
    fixT = sorted([sum(1 for m in range(8) if m ^ 7 == m), sum(1 for m in range(8) if trans(m) == m),
                   sum(1 for m in range(8) if trans(m) ^ 7 == m)])
    # line chain: lines -> converse orbits of lines
    lorb = []
    for a, b in lines:
        o = frozenset({frozenset({a, b}), frozenset({trans(a), trans(b)})})
        if o not in lorb:
            lorb.append(o)
    lsh = sorted(map(len, lorb), reverse=True)
    check('B6 line chain 8 -> 4 lines -> 3 converse-orbits: shapes (2,2,2,2), (2,1,1); composite (4,2,2) = f^2 '
          'fibres (4,2,2) [GENERIC: equal-shape partition chains are always isomorphic]',
          lsh == [2, 1, 1] and osh == [4, 2, 2] and A['sh12'] == [4, 2, 2])
    cfib = sorted(mfib, reverse=True)
    check('B7 owner chain composite 8 -> 3 has fibres {5,2,1} vs f^2 on B {4,2,2}: differ', cfib != A['sh12'],
          f'{cfib} vs {A["sh12"]}')
    check('B8 natural Klein four-groups are not conjugate in S_8: fixed points <flip,transpose> on tilings vs '
          '<negation, Aut(f|B)> on B vs mod-16 units {1,7,9,15} on B', fixT != A['fixB'] and fixT != A['fix16'],
          f'{fixT} vs {A["fixB"]} vs {A["fix16"]}')


# ---------------------------------------------------------------- C
def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return s


def part_C():
    print('C  primes modulo even m until +-1')
    LIM = 600000
    S = sieve(LIM)
    primes = [p for p in range(3, LIM + 1) if S[p]]  # odd primes
    N = 20000
    P = {}
    for m in range(4, N + 1, 2):
        k = 1
        while True:
            c = [x for x in (k * m - 1, k * m + 1) if x <= LIM and S[x]]
            if c:
                P[m] = min(c)
                break
            assert (k + 1) * m + 1 <= LIM
            k += 1
    below = all(P[m] < m * m for m in P)
    near4 = sum(1 for m in P if P[m] in (m - 1, m + 1))
    near6 = sum(1 for m in P if m >= 6 and P[m] in (m - 1, m + 1))
    mr = max(P, key=lambda m: Fr(P[m], m))
    kmax = max(round(P[m] / m) for m in P)
    kargs = [m for m in P if round(P[m] / m) == kmax]
    check('C1 P(m) < m^2 for all even m in [4, 20000]', below, f'max P = {max(P.values())}')
    check('C2 P(m) in {m-1, m+1}: 4179 of 9999 (m >= 4), 4178 of 9998 (m >= 6); 41.8%',
          near4 == 4179 and len(P) == 9999 and near6 == 4178, f'{near4}/{len(P)} = {near4 / len(P):.4f}')
    check('C3 largest P(m)/m at m = 16012 (P = 17m - 1, ratio 16.99994, not exactly 17)',
          mr == 16012 and kmax == 17 and kargs == [16012], f'P(16012) = {P[16012]} = 17*16012 {P[16012] - 17 * 16012:+d}')
    sub_ok = eq_ok = True
    strict, equal_nonadj = [], []
    nonadj = 0
    wit_i = 0
    wit_other = []
    for m in range(4, N + 1, 2):
        idx = bisect_left(primes, P[m])
        res = bytearray(m)
        for p in primes[:idx]:
            res[p % m] = 1
        R = {(r, m - r) for r in range(3, m // 2 + 1, 2) if res[r] and res[m - r]}
        G = {(p, m - p) for p in range(3, m // 2 + 1, 2) if S[p] and S[m - p]}
        if not G <= R:
            sub_ok = False
        if S[m - 1] or S[m + 1]:
            if R != G:
                eq_ok = False
            continue
        nonadj += 1
        if R == G:
            equal_nonadj.append(m)
            continue
        strict.append(m)
        # every extra pair has a non-prime member, which must be the residue of a prime > m (wrap-around)
        extra = R - G
        assert all(not (S[r] and S[s]) for r, s in extra)
        # witness type (i): Goldbach partition 2m = p + q with q < m < p and p - m composite
        lo, hi = bisect_right(primes, m), bisect_right(primes, 2 * m - 3)
        if any(S[2 * m - p] and not S[p - m] for p in primes[lo:hi]):
            wit_i += 1
        else:
            wit_other.append(m)
    check('C4 G(m) subset of R(m) for all even m in [4, 20000]', sub_ok)
    check('C5 R(m) = G(m) whenever m-1 or m+1 is prime', eq_ok)
    check('C6 no prime neighbour => R(m) strictly contains G(m): 5820 of 5820 (m <= 20000)',
          len(strict) == 5820 and nonadj == 5820 and not equal_nonadj, f'first {strict[:8]}')
    check('C7 witness family (i) {p-m, 2m-p}: 2m = p + q, q < m < p, p - m composite, covers every strict case',
          wit_i == len(strict), f'covered {wit_i}/{len(strict)}; others {wit_other[:10]}')
    # cross-check with the prior prime-wait packet (jacobian_eth_prime_wait_20261102/PRIME-MRB.md, p_*(m) = P(m))
    rows = {120: 239, 144: 431, 208: 1039, 262: 263, 278: 277, 290: 1451, 292: 293, 210: 211, 1688: 15193, 9970: 149551}
    pi = lambda x: bisect_right(primes, x) + (1 if x >= 2 else 0)
    delta_ok = all((pi(P[m]) - pi(m - 2) == 1) == bool(S[m - 1] or S[m + 1]) for m in range(6, N + 1, 2))
    check('C8 P(m) = PRIME-MRB p_*(m) on its rows (120..9970); "prime neighbour" <=> its delta(m) = pi(p_*) - pi(m-2) = 1',
          all(P[m] == v for m, v in rows.items()) and delta_ok)
    return P


# ---------------------------------------------------------------- D
def part_D():
    print('D  clique-carrier dimensions (chart_clique_norm_20261103 formula)')

    def dim(n, r=3, t=(3, 3, 5)):
        return (1 + (r - 1) * math.comb(n, 2)
                + sum(math.comb(n, k) for ti in t for k in range(3, min(n, ti - 1) + 1))
                + sum(1 for ti in t if n >= ti))

    def note_display(n):  # the note's specialisation 1 + 2C(n,2) + sum_{k=3}^{4} C(n,k) + #{t_i <= n}
        return 1 + 2 * math.comb(n, 2) + sum(math.comb(n, k) for k in (3, 4)) + sum(1 for ti in (3, 3, 5) if ti <= n)

    v = [dim(n) for n in (3, 4, 5, 6)]
    check('D1 r = 3, t = (3,3,5): n = 3,4,5,6 -> 10, 20, 39, 69; the note\'s display agrees for n = 3..60',
          v == [10, 20, 39, 69] and all(dim(n) == note_display(n) for n in range(3, 61)), f'{v}')


# ---------------------------------------------------------------- E
def part_E():
    print('E  printed triangles vs recomputation')
    S = sieve(100)
    out = open(os.path.join(HERE, '..', 'quadratic_triangle_n4.out')).read().splitlines()
    tri = {}
    for line in out:
        s = line.strip()
        if s.startswith('m=') and ':' in s and 'P=' not in s:
            mm, row = s[2:].split(':')
            tri[int(mm)] = row.strip()
    mine = {m: ''.join('*' if S[a] and S[m - a] else '.' for a in range(1, m // 2 + 1, 2)) for m in range(6, 41, 2)}
    check('E1 odd-pair triangle rows m = 6..40 match', tri == mine, f'{len(tri)} rows')
    pg = {}
    for line in out:
        s = line.strip()
        if s.startswith('n=') and 'trials' in s:
            nn = int(s[2:s.index(':')])
            d_hi, d_lo = s.split('d = ')[1].split()[0].split('..')
            pg[nn] = (int(d_hi), int(d_lo), int(s.split('(')[1].split()[0]), 'prime' in s)
    ok = True
    for n in range(2, 22):
        D = max(d for d in range(1, n) if n % d == 0)
        ok &= pg.get(n) == (n - 1, D, n - D, D == 1)
        ok &= (D == 1) == (all(n % q for q in range(2, n)))
    check('E2 PRIMEGAME rows n = 2..21: d from n-1 down to D(n), n - D(n) trials, D(n) = 1 iff n prime', ok)
    # shift identity behind the CA reading: chi(m+2-a) = chi(m-(a-2)) (trivial)
    check('E3 shift identity chi(m+2-a) = chi(m-(a-2)) (trivially true; GENERIC)', True)


def main():
    A = part_A()
    part_B(A)
    part_C()
    part_D()
    part_E()
    print()
    print('ALL INDEPENDENT CHECKS AGREE' if not BAD else f'DISAGREEMENTS: {BAD}')
    return 0 if not BAD else 1


if __name__ == '__main__':
    sys.exit(main())
