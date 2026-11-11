"""cd_pg42.py - PG(4,2): the Cayley-Dickson sign tournament CD_5 is not extremal (S23, after the audit request).

Pullback construction (any n): W a subspace of codimension 3, p: Z_2^n -> Z_2^n / W = Z_2^3 the projection.
Arcs between points with distinct nonzero projections are copied from an octonion tournament on Z_2^3; then every
plane complementary to W (there are 8^(n-3)) is a copy of that tournament, hence odd with cyclic lines; all other
lines have their three arcs free and are oriented cyclically.
For n = 5 (W a line) the seven planes containing W lie in the free part and can be made odd as well: 64 + 7 = 71.

modes:
  python cd_pg42.py verify <certificate>      direct recount from the arc word (no linear algebra, no imports)
  python cd_pg42.py build <out file>          F_2 system: the 64 planes disjoint from the line {1,2,3} and the 7
                                              planes containing it odd, all 155 lines cyclic; writes a certificate
  python cd_pg42.py compare <random orders>   CD_5 (50 odd planes; maximal under inclusion?) and greedy lower bounds
"""
import sys, random, itertools


def verify(fn):
    word = None
    for l in open(fn):
        if l.startswith("x ="):
            word = int(l.split("=")[1])
    N = 32
    k = 0
    out = [[False] * N for _ in range(N)]
    for u in range(1, N):
        for v in range(u + 1, N):
            if word >> k & 1:
                out[v][u] = True          # bit 1: v -> u
            else:
                out[u][v] = True
            k += 1
    V = range(1, N)

    def cyc(a, b, c):
        return (out[a][b] and out[b][c] and out[c][a]) or (out[b][a] and out[c][b] and out[a][c])
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
    ncyc = sum(1 for L in lines if cyc(*L))
    quads = sorted({tuple(sorted((a, b, c, a ^ b ^ c))) for a in V for b in V for c in V
                    if len({a, b, c, a ^ b ^ c}) == 4 and (a ^ b ^ c) != 0})
    even = set()
    for q in quads:
        c3 = sum(1 for t in itertools.combinations(q, 3) if cyc(*t))
        if c3 % 2 == 0:
            even.add(q)
    pls = set()
    for a in V:
        for b in V:
            for c in V:
                if a < b < c and c != a ^ b:
                    span = {0}
                    for v in (a, b, c):
                        span |= {s ^ v for s in span}
                    pls.add(tuple(sorted(span - {0})))
    defects = {}
    kinds = {}
    for pl in pls:
        s = set(pl)
        d = sum(1 for q in even if s >= set(q))
        defects[d] = defects.get(d, 0) + 1
        if d == 0:
            m = len(s & {1, 2, 3})
            kinds[m] = kinds.get(m, 0) + 1
    scores = sorted(set(sum(1 for w in V if w != u and out[u][w]) for u in V))
    print("lines %d, cyclic %d; zero-sum 4-sets %d, even %d; planes %d, by number of even 4-sets %s; scores %s"
          % (len(lines), ncyc, len(quads), len(even), len(pls), sorted(defects.items()), scores))
    print("odd planes by number of points on the line {1, 2, 3}:", sorted(kinds.items()))


if sys.argv[1] == "verify":
    verify(sys.argv[2])
    sys.exit(0)

from cd_more import F2System, planes, c3_of, table

n = 5
S = F2System(n)
P = planes(n)
rows = []
for pl in P:
    s = sorted(pl)
    rr = []
    for i in range(7):
        for j in range(i + 1, 7):
            for k in range(j + 1, 7):
                d = s[i] ^ s[j] ^ s[k]
                if d in pl and d > s[k]:
                    rr += S.quad_rows((s[i], s[j], s[k], d))
    assert len(rr) == 14
    rows.append(rr)
V = range(1, 1 << n)
lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
line_rows = []
for (a, b, c) in lines:
    for (p, q) in (((a, b), (b, c)), ((b, c), (c, a))):
        m1, c1 = S.term(*p)
        m2, c2 = S.term(*q)
        line_rows.append((m1 ^ m2, c1 ^ c2))


def add(piv, mask, const):
    while mask:
        hb = mask.bit_length() - 1
        if hb in piv:
            pm, pc = piv[hb]
            mask ^= pm
            const ^= pc
        else:
            piv[hb] = (mask, const)
            return True
    return const == 0


def solution(piv):
    x = 0
    for hb in sorted(piv):
        mask, const = piv[hb]
        val = const ^ (bin(x & (mask ^ (1 << hb))).count("1") & 1)
        if val:
            x |= 1 << hb
    return x


if sys.argv[1] == "build":
    L = {1, 2, 3}
    chosen = [i for i, pl in enumerate(P) if len(set(pl) & L) in (0, 3)]
    piv = {}
    ok = all(add(piv, m, c) for m, c in line_rows)
    for i in chosen:
        for m, c in rows[i]:
            ok &= add(piv, m, c)
    print("planes disjoint from the line {1,2,3}: %d; containing it: %d; all odd with all lines cyclic: %s (rank %d of %d)"
          % (sum(1 for i in chosen if not (set(P[i]) & L)), sum(1 for i in chosen if len(set(P[i]) & L) == 3), ok, len(piv), len(S.var)))
    more = 0
    for j in range(len(P)):
        if j in chosen:
            continue
        p2 = dict(piv)
        if all(add(p2, m, c) for m, c in rows[j]):
            more += 1
    print("further planes that can be made odd as well:", more)
    with open(sys.argv[2], "w") as f:
        f.write("# tournament on the nonzero vectors 1..31 of Z_2^5; pairs u < v in lexicographic order; bit k of x is 1 iff v -> u\n")
        f.write("# all 155 lines cyclic; the 64 planes disjoint from the line {1, 2, 3} and the 7 planes containing it are odd\n")
        f.write("x = %d\n" % solution(piv))
    verify(sys.argv[2])
elif sys.argv[1] == "compare":
    g = table([-1] * n)
    arc = lambda i, j: g[i][j] == 1
    odd = []
    for idx, pl in enumerate(P):
        s = sorted(pl)
        bad = 0
        for i in range(7):
            for j in range(i + 1, 7):
                for k in range(j + 1, 7):
                    d = s[i] ^ s[j] ^ s[k]
                    if d in pl and d > s[k]:
                        bad += c3_of(arc, (s[i], s[j], s[k], d)) % 2 == 0
        if bad == 0:
            odd.append(idx)
    piv = {}
    assert all(add(piv, m, c) for i in odd for m, c in rows[i])
    ext = 0
    for j in range(len(P)):
        if j not in odd:
            p2 = dict(piv)
            ext += all(add(p2, m, c) for m, c in rows[j])
    print("CD_5: %d odd planes of %d; further planes that could be made odd as well: %d" % (len(odd), len(P), ext))
    rng = random.Random(5)
    hist = {}
    base = {}
    assert all(add(base, m, c) for m, c in line_rows)
    for rep in range(int(sys.argv[2])):
        order = list(range(len(P)))
        rng.shuffle(order)
        piv = dict(base)
        size = 0
        for i in order:
            p2 = dict(piv)
            if all(add(p2, m, c) for m, c in rows[i]):
                piv = p2
                size += 1
        hist[size] = hist.get(size, 0) + 1
    print("greedy sets of odd planes (all lines cyclic) over %s random orders, by size: %s" % (sys.argv[2], sorted(hist.items())))
