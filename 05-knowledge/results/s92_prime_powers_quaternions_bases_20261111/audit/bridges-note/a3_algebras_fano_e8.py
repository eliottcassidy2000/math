#!/usr/bin/env python3
"""Auditor bridges-note, script A3 (independent code): sections 3 and 4 of the S92 note.

(a) the 16 normalized 2-cocycles on F_2^2, the 8 classes, groups, sign tables, zero divisors
(b) norm counts: Hamilton on Z<1,i,j,k>; bicomplex on Z[i]^2 and on the sign-table lattice Z<1,i,j,ij>
(c) Out(Q_8) = S_3 acting on the three axes; even = automorphisms of the cyclic sign triangle
(d) the 128 cyclic orientations of the Fano lines: H of the tournaments, alternativity
(e) 240 labelled Paley tournaments; 30 Fano planes x 16 orientations -> each Paley tournament twice
(f) E8: shells of norm 2, 4, 6; class counts of x^T C x / 2 modulo 2, 4, 8 against the rank-4 hyperbolic form
(g) section 5 shares on the family h1 <= h2 < 64, s <= 4
"""
import itertools
import sys
from math import gcd
import numpy as np
from numba import njit, prange


# ------------------------------------------------------------------ (a)
def part_a():
    print("=== A3.a central extensions of F_2^2 by F_2 ===")
    V = range(4)
    nz = (1, 2, 3)
    cocycles = []
    for bits in range(1 << 9):
        f = [[0] * 4 for _ in V]
        k = 0
        for a in nz:
            for b in nz:
                f[a][b] = (bits >> k) & 1
                k += 1
        ok = True
        for a in V:
            for b in V:
                for c in V:
                    if (f[a][b] ^ f[a ^ b][c]) != (f[b][c] ^ f[a][b ^ c]):
                        ok = False
        if ok:
            cocycles.append(tuple(tuple(r) for r in f))
    print("normalized 2-cocycles:", len(cocycles))
    cobs = set()
    for g in itertools.product((0, 1), repeat=3):
        gg = (0,) + g
        cobs.add(tuple(tuple(gg[a] ^ gg[b] ^ gg[a ^ b] for b in V) for a in V))
    print("distinct normalized coboundaries:", len(cobs))
    nontriv = [c for c in cobs if any(any(r) for r in c)][0]
    print("non-trivial coboundary (rows a = 0..3):", nontriv)
    classes = {}
    for f in cocycles:
        orbit = frozenset(tuple(tuple(f[a][b] ^ c[a][b] for b in V) for a in V) for c in cobs)
        classes.setdefault(orbit, f)
    print("classes:", len(classes), "; sizes:", sorted(len(o) for o in classes))
    summary = {}
    for orbit in classes:
        fs = sorted(orbit)
        f = fs[0]
        S = tuple(a for a in nz if f[a][a] == 1)
        # all cocycles of the class have the same squares?
        assert all(tuple(a for a in nz if g[a][a] == 1) == S for g in fs)

        def mul(x, y):
            return (x[0] ^ y[0], x[1] ^ y[1] ^ f[x[0]][y[0]])
        elems = [(a, z) for a in V for z in (0, 1)]

        def order(x):
            y, k = x, 1
            while y != (0, 0):
                y = mul(y, x)
                k += 1
            return k
        census = sorted(order(x) for x in elems)
        n4 = census.count(4)
        n2 = census.count(2)
        abel = all(mul(x, y) == mul(y, x) for x in elems for y in elems)
        if abel and n4 == 0:
            name = "C2^3"
        elif abel and n4 == 4:
            name = "C4 x C2"
        elif (not abel) and n4 == 2:
            name = "D4 (order 8)"
        elif (not abel) and n4 == 6:
            name = "Q8"
        else:
            name = "?"
        kinds = []
        for g in fs:
            anti = all(g[a][b] != g[b][a] for a in nz for b in nz if a < b)
            symm = all(g[a][b] == g[b][a] for a in nz for b in nz if a < b)
            if anti:
                out = {a: sum(1 for b in nz if b != a and g[a][b] == 0) for a in nz}
                if sorted(out.values()) == [1, 1, 1]:
                    kinds.append("cyclic")
                else:
                    mid = [a for a in nz if out[a] == 1][0]
                    kinds.append("transitive(middle %d)" % mid)
            elif symm:
                kinds.append("signed triangle(%d negative)" % sum(g[a][b] for a in nz for b in nz if a < b))
            else:
                kinds.append("mixed")
        # the two cocycles of a non-commutative class are converse tournaments?
        conv = None
        if not abel:
            g0, g1 = fs
            conv = all(g0[a][b] == g1[b][a] for a in nz for b in nz if a != b)
        # zero divisors (coordinates in {-1,0,1})
        zd = False
        vecs = [v for v in itertools.product((-1, 0, 1), repeat=4) if any(v)]
        for x in vecs:
            for y in vecs:
                z = [0] * 4
                for a in V:
                    if x[a]:
                        for b in V:
                            if y[b]:
                                z[a ^ b] += x[a] * y[b] * (-1) ** f[a][b]
                if not any(z):
                    zd = True
                    break
            if zd:
                break
        prod_sq = 1
        for a in nz:
            prod_sq *= (-1) ** f[a][a]
        print("  S = %-10s #S = %d  group %-13s orders %s  commutative %-5s  sign tables %s  converse pair %s  zero divisors %s  product of squares %+d"
              % (S, len(S), name, census, abel, kinds, conv, zd, prod_sq))
        summary.setdefault((len(S), name), 0)
        summary[(len(S), name)] += 1
    print("classes by (#S, group):", sorted(summary.items()))


# ------------------------------------------------------------------ (b)
def part_b():
    print()
    print("=== A3.b norm counts ===")
    LIM = 400
    B = 21
    r2 = [0] * (LIM + 1)
    for a in range(-B, B + 1):
        for b in range(-B, B + 1):
            n = a * a + b * b
            if n <= LIM:
                r2[n] += 1
    sig = [0] * (LIM + 1)
    for d in range(1, LIM + 1):
        for k in range(d, LIM + 1, d):
            sig[k] += d
    # Hamilton on Z^4: direct
    rng = np.arange(-B, B + 1)
    sq = rng * rng
    two = (sq[:, None] + sq[None, :]).ravel()
    four = (two[:, None] + two[None, :]).ravel()
    r4 = np.bincount(four[four <= LIM], minlength=LIM + 1)
    print("r_4(n) = 8 sigma(n) for all odd n <= 400 (direct count on Z^4):",
          all(int(r4[n]) == 8 * sig[n] for n in range(1, LIM + 1, 2)))
    print("r_4(n) = 24 sigma(odd part) for even n <= 400:",
          all(int(r4[n]) == 24 * sig[n // (n & -n)] for n in range(2, LIM + 1, 2)))
    dir_conv = [0] * (LIM + 1)
    for a in range(1, LIM + 1):
        for b in range(1, LIM // a + 1):
            dir_conv[a * b] += r2[a] * r2[b]
    print("Dirichlet convolution (r_2 * r_2)(n) [= norm count on the maximal order Z[i]^2] at 1, 2, 4, 5, 9, 13, 25:",
          [dir_conv[n] for n in (1, 2, 4, 5, 9, 13, 25)])
    print("   it vanishes for every n = 3 mod 4, n <= 400:", all(dir_conv[n] == 0 for n in range(3, LIM + 1, 4)))
    print("   other n <= 60 with (r_2 * r_2)(n) = 0:", [n for n in range(1, 61) if dir_conv[n] == 0 and n % 4 != 3])
    # sign-table lattice Z<1, i, j, k = ij>, i^2 = j^2 = -1, ij = ji:  x = a + b i + c j + d k
    # -> (z1, z2) = ((a - d) + (b + c) i, (a + d) + (b - c) i),  norm |z1|^2 |z2|^2
    a, b, c, d = np.meshgrid(rng, rng, rng, rng, indexing="ij")
    n1 = (a - d) ** 2 + (b + c) ** 2
    n2 = (a + d) ** 2 + (b - c) ** 2
    N = (n1 * n2).ravel()
    # the box |a|,|b|,|c|,|d| <= 21 contains every solution with N <= 400 and N >= 1:
    # n1, n2 >= 1 and n1 n2 <= 400 -> |a - d|, |a + d| <= 20 -> |a|, |d| <= 20
    cnt = np.bincount(N[(N >= 1) & (N <= LIM)], minlength=LIM + 1)
    print("sign-table lattice Z<1,i,j,ij> (index 4 in Z[i]^2): norm counts at 1, 2, 4, 5, 9, 13, 25:",
          [int(cnt[n]) for n in (1, 2, 4, 5, 9, 13, 25)])
    print("   = (r_2 * r_2)(n) / 2 for every odd n <= 400:",
          all(2 * int(cnt[n]) == dir_conv[n] for n in range(1, LIM + 1, 2)))
    print("   norm-1 elements (the signed units): %d  (the group C_4 x C_2 has order 8; Z[i]^2 has 16 units)" % int(cnt[1]))
    ev = [(n, int(cnt[n]), dir_conv[n]) for n in (2, 4, 8, 10, 16, 20)]
    print("   even n: (n, sign-table lattice count, Dirichlet convolution):", ev)


# ------------------------------------------------------------------ (c)
def part_c():
    print()
    print("=== A3.c Out(Q_8) on the three axes ===")
    # elements (sign, idx), idx 0..3 = 1, i, j, k
    tab = {(1, 2): (1, 3), (2, 3): (1, 1), (3, 1): (1, 2), (2, 1): (-1, 3), (3, 2): (-1, 1), (1, 3): (-1, 2)}

    def mul(x, y):
        s = x[0] * y[0]
        a, b = x[1], y[1]
        if a == 0:
            return (s, b)
        if b == 0:
            return (s, a)
        if a == b:
            return (-s, 0)
        t = tab[(a, b)]
        return (s * t[0], t[1])
    elems = [(s, a) for s in (1, -1) for a in range(4)]
    order4 = [e for e in elems if e[1] != 0]
    auts = []
    for xi in order4:
        for xj in order4:
            if xj[1] == xi[1]:
                continue
            xk = mul(xi, xj)
            img = {(1, 0): (1, 0), (-1, 0): (-1, 0), (1, 1): xi, (1, 2): xj, (1, 3): xk,
                   (-1, 1): (-xi[0], xi[1]), (-1, 2): (-xj[0], xj[1]), (-1, 3): (-xk[0], xk[1])}
            if all(img[mul(x, y)] == mul(img[x], img[y]) for x in elems for y in elems):
                auts.append(img)
    print("|Aut(Q_8)| =", len(auts))
    perms = {}
    for img in auts:
        pi = tuple(img[(1, a)][1] for a in (1, 2, 3))
        perms.setdefault(pi, 0)
        perms[pi] += 1
    print("induced permutations of the axes (i, j, k) -> multiplicity (kernel = Inn):", sorted(perms.items()))
    # sign triangle of Q_8: a -> b iff e_a e_b = + e_c : 1 -> 2 -> 3 -> 1
    arcs = {(1, 2), (2, 3), (3, 1)}
    kind = {}
    for pi in perms:
        m = {1: pi[0], 2: pi[1], 3: pi[2]}
        img_arcs = {(m[a], m[b]) for (a, b) in arcs}
        kind[pi] = "automorphism" if img_arcs == arcs else ("anti-automorphism" if img_arcs == {(b, a) for (a, b) in arcs} else "neither")
    print("action on the cyclic sign triangle 1 -> 2 -> 3 -> 1:", sorted(kind.items()))
    print("so Aut(Q_8)/Inn = S_3 acts on the axes as Aut+-(C_3): 3 automorphisms, 3 anti-automorphisms:",
          sorted(kind.values()).count("automorphism") == 3 and sorted(kind.values()).count("anti-automorphism") == 3)


# ------------------------------------------------------------------ tournaments
@njit(cache=True)
def hcount(adj, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            out = adj[v] & ~mask
            w = 0
            while out:
                if out & 1:
                    dp[mask | (1 << w), w] += c
                out >>= 1
                w += 1
    s = 0
    for v in range(n):
        s += dp[full, v]
    return s


@njit(parallel=True, cache=True)
def all_H7():
    n = 7
    total = 1 << 21
    out = np.zeros(total, np.int64)
    for code in prange(total):
        adj = np.zeros(n, np.int64)
        b = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (code >> b) & 1:
                    adj[i] |= (1 << j)
                else:
                    adj[j] |= (1 << i)
                b += 1
        out[code] = hcount(adj, n)
    return out


def code_of(arcs, n=7):
    code = 0
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in arcs:
                code |= (1 << b)
            b += 1
    return code


def adj_of_arcs(arcs, n=7):
    adj = np.zeros(n, np.int64)
    for (i, j) in arcs:
        adj[i] |= (1 << j)
    return adj


def oct_mult_table(lines_oriented):
    """sign[a][b], third[a][b] for a != b in 0..6 : e_a e_b = sign * e_third"""
    sign = {}
    third = {}
    for (a, b, c) in lines_oriented:  # cyclic a -> b -> c -> a
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            sign[(x, y)] = 1
            third[(x, y)] = z
            sign[(y, x)] = -1
            third[(y, x)] = z
    return sign, third


def alg_mul(x, y, sign, third):
    # x, y integer vectors of length 8: index 0 = 1, index a+1 = e_a
    z = [0] * 8
    for a in range(8):
        if x[a] == 0:
            continue
        for b in range(8):
            if y[b] == 0:
                continue
            co = x[a] * y[b]
            if a == 0:
                z[b] += co
            elif b == 0:
                z[a] += co
            elif a == b:
                z[0] -= co
            else:
                z[third[(a - 1, b - 1)] + 1] += co * sign[(a - 1, b - 1)]
    return z


def part_d_e():
    print()
    print("=== A3.d the 128 cyclic orientations of the Fano lines ===")
    base_lines = [tuple(sorted(((0 + i) % 7, (1 + i) % 7, (3 + i) % 7))) for i in range(7)]
    rng = np.random.default_rng(7)
    tests = [(list(rng.integers(-3, 4, 8)), list(rng.integers(-3, 4, 8))) for _ in range(12)]
    res = {}
    codes_alt = set()
    for bits in range(128):
        ori = []
        arcs = set()
        for k, (a, b, c) in enumerate(base_lines):
            cyc = (a, b, c) if (bits >> k) & 1 else (a, c, b)
            ori.append(cyc)
            arcs |= {(cyc[0], cyc[1]), (cyc[1], cyc[2]), (cyc[2], cyc[0])}
        assert len(arcs) == 21
        H = int(hcount(adj_of_arcs(arcs), 7))
        sign, third = oct_mult_table(ori)
        alt = True
        normmult = True
        for x, y in tests:
            xx = alg_mul(x, x, sign, third)
            if alg_mul(xx, y, sign, third) != alg_mul(x, alg_mul(x, y, sign, third), sign, third):
                alt = False
            if alg_mul(y, xx, sign, third) != alg_mul(alg_mul(y, x, sign, third), x, sign, third):
                alt = False
            xy = alg_mul(x, y, sign, third)
            if sum(t * t for t in xy) != sum(t * t for t in x) * sum(t * t for t in y):
                normmult = False
        # doubly regular?
        adj = adj_of_arcs(arcs)
        dr = all(bin(int(adj[u]) & int(adj[v])).count("1") == 1 for u in range(7) for v in range(u + 1, 7))
        res.setdefault((H, alt, normmult, dr), 0)
        res[(H, alt, normmult, dr)] += 1
        if alt:
            codes_alt.add(code_of(arcs))
    print("(H, alternative, norm multiplicative, doubly regular) -> number of orientations:", sorted(res.items()))

    print()
    print("=== A3.e 240 labelled Paley tournaments ===")
    H7 = all_H7()
    c189 = np.nonzero(H7 == 189)[0]
    print("labelled tournaments on 7 vertices with H = 189:", len(c189), " (7!/21 = %d)" % (5040 // 21))
    print("max H at order 7:", int(H7.max()), "; labelled tournaments with H = 171:", int((H7 == 171).sum()))
    # 30 labelled Fano planes
    planes = set()
    for p in itertools.permutations(range(7)):
        planes.add(frozenset(frozenset(p[x] for x in ln) for ln in base_lines))
    print("labelled Fano planes on 7 points:", len(planes))
    hits = {}
    pairs = 0
    for pl in planes:
        lines = [tuple(sorted(ln)) for ln in pl]
        for bits in range(128):
            ori = []
            arcs = set()
            for k, (a, b, c) in enumerate(lines):
                cyc = (a, b, c) if (bits >> k) & 1 else (a, c, b)
                ori.append(cyc)
                arcs |= {(cyc[0], cyc[1]), (cyc[1], cyc[2]), (cyc[2], cyc[0])}
            code = code_of(arcs)
            if H7[code] == 189:
                pairs += 1
                hits[code] = hits.get(code, 0) + 1
    print("(Fano plane, orientation of its lines) pairs giving a Paley tournament: %d = 30 x %d" % (pairs, pairs // 30))
    print("distinct labelled Paley tournaments reached: %d ; each reached %s times" % (len(hits), sorted(set(hits.values()))))
    print("so 240 = 30 * 16 / 2 = 16 * 15 on the Paley side (15 = pairs of Fano planes sharing the 14 cyclic triangles)")
    # cyclic triangles of P7 and the two Fano planes
    arcs = {(i, (i + r) % 7) for i in range(7) for r in (1, 2, 4)}
    tri = [t for t in itertools.combinations(range(7), 3)
           if ((t[0], t[1]) in arcs) == ((t[1], t[2]) in arcs) == ((t[2], t[0]) in arcs)]
    print("P_7 (QR = 1, 2, 4): cyclic triangles: %d ; they split into the two Fano planes {0,1,3}+i and {0,2,3}+i: %s"
          % (len(tri), sorted(tri) == sorted(set(tuple(sorted(((a + i) % 7, (b + i) % 7, (c + i) % 7)))
                                                   for i in range(7) for (a, b, c) in ((0, 1, 3), (0, 2, 3))))))
    # automorphisms of P7
    adj = adj_of_arcs(arcs)
    na = 0
    for p in itertools.permutations(range(7)):
        if all(((p[i], p[j]) in arcs) for (i, j) in arcs):
            na += 1
    print("|Aut(P_7)| =", na)


# ------------------------------------------------------------------ (f)
def part_f():
    print()
    print("=== A3.f E8 ===")
    # D8+ model: shells of norm 2n, n = 1, 2, 3
    sig3 = lambda n: sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
    shells = []
    r = np.arange(-2, 3)
    grids = np.array(np.meshgrid(*([r] * 8), indexing="ij")).reshape(8, -1)
    nrm = (grids * grids).sum(axis=0)
    par = grids.sum(axis=0) % 2
    r2 = np.array([-3, -1, 1, 3])
    g2 = np.array(np.meshgrid(*([r2] * 8), indexing="ij")).reshape(8, -1)
    nrm2 = (g2 * g2).sum(axis=0)  # 4 * norm
    s2 = g2.sum(axis=0) % 4       # sum of coordinates even <=> sum of doubled coordinates = 0 mod 4
    for n in (1, 2, 3):
        cnt_int = int(((nrm == 2 * n) & (par == 0)).sum())
        cnt_half = int(((nrm2 == 8 * n) & (s2 == 0)).sum())
        shells.append((cnt_int, cnt_half, cnt_int + cnt_half, 240 * sig3(n)))
    print("E8 = D8+ shells of norm 2n, n = 1, 2, 3: (integer vectors, half-integer vectors, total, 240 sigma_3(n)):", shells)
    print("   (coordinates bounded by 2 resp. 3/2: sufficient for norm <= 6 since 3^2 = 9 > 6 and (5/2)^2 + 7/4 = 8 > 6)")
    # Cartan matrix of E8 (Bourbaki): chain 1-3-4-5-6-7-8, node 2 attached to 4
    C = np.zeros((8, 8), np.int64)
    for i in range(8):
        C[i, i] = 2
    edges = [(0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (1, 3)]
    for (i, j) in edges:
        C[i, j] = C[j, i] = -1
    det = round(np.linalg.det(C.astype(float)))
    print("det of the E8 Cartan matrix:", det)
    for q in (2, 4, 8):
        counts = e8_class_counts(C, q)
        # rank-4 hyperbolic form: additive 4th power of the table modulo q
        tab = np.zeros(q, np.int64)
        for a in range(q):
            for b in range(q):
                tab[(a * b) % q] += 1
        d = tab.copy()
        for _ in range(3):
            nd = np.zeros(q, np.int64)
            for x in range(q):
                for y in range(q):
                    nd[(x + y) % q] += d[x] * tab[y]
            d = nd
        print("modulo %d: E8 form class counts %s ; rank-4 hyperbolic %s ; equal: %s"
              % (q, list(map(int, counts)), list(map(int, d)), list(map(int, counts)) == list(map(int, d))))
        # predicted factor (1 - 2^-4) sigma_3(2^v)/2^(3v) for v < k
        k = q.bit_length() - 1
        pred = []
        for c in range(q):
            v = (c & -c).bit_length() - 1 if c else k
            if c and v < k:
                pred.append(q ** 7 * (1 - 2 ** -4) * sum(2 ** (-3 * i) for i in range(v + 1)))
            else:
                pred.append(None)
        print("          factor formula q^7 (1 - 2^-4) sigma_3(2^v)/2^(3v) at classes of valuation v < k:",
              [None if p is None else int(round(p)) for p in pred])
    z4 = np.pi ** 4 / 90
    print("(8 pi^4/3)/zeta(4) =", (8 * np.pi ** 4 / 3) / z4)


@njit(parallel=True, cache=True)
def e8_class_counts(C, q):
    total = q ** 8
    nth = 64
    acc = np.zeros((nth, q), np.int64)
    chunk = total // nth + 1
    for t in prange(nth):
        x = np.zeros(8, np.int64)
        lo = t * chunk
        hi = min(total, lo + chunk)
        for idx in range(lo, hi):
            r = idx
            for i in range(8):
                x[i] = r % q
                r //= q
            s = 0
            for i in range(8):
                s += x[i] * x[i]
                for j in range(i + 1, 8):
                    if C[i, j] != 0:
                        s += C[i, j] * x[i] * x[j]
            acc[t, s % q] += 1
    out = np.zeros(q, np.int64)
    for t in range(nth):
        for c in range(q):
            out[c] += acc[t, c]
    return out


# ------------------------------------------------------------------ (g)
def part_g():
    print()
    print("=== A3.g section 5 shares on the family h1 <= h2 < 64, s <= 4 ===")
    psi = [0] * 64
    for h in range(1, 64):
        k = h >> 1
        psi[h] = psi[k] + k if h % 2 == 0 else psi[k] + k // 2 + 1
    rows = []
    for h1 in range(64):
        for h2 in range(h1, 64):
            for s in range(5):
                cc = (1 << s) * (2 * h1 + 1) * (2 * h2 + 1) - (1 if s == 0 else 0)
                rows.append((psi[h1] + psi[h2] + cc // 2, psi[h1] * psi[h2], h1, h2, s))
    print("joins in the family:", len(rows))
    for b in (2, 3, 10, 100):
        good = [r for r in rows if max(r[0], r[1]) <= b - 1]
        print("  b = %3d: %4d joins with max(alpha_1, alpha_2) <= b - 1: share %.5f  (rounded %.4f)" % (b, len(good), len(good) / len(rows), len(good) / len(rows)))
        if b <= 3:
            print("       ", [(r[2], r[3], r[4], r[0], r[1]) for r in good], "(h1, h2, s, alpha_1, alpha_2)")
    # the same shares when s is not capped at 4 (all joins with N <= 10^6, h < 64)
    vec = {}
    for r in rows:
        vec.setdefault(1 + 2 * r[0] + 4 * r[1], set()).add((r[0], r[1]))
    coll = sorted((H, sorted(v)) for H, v in vec.items() if len(v) > 1)
    print("least path counts with more than one packing vector in the family:", coll[:5])


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d_e()
    part_f()
    part_g()
