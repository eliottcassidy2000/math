"""aud_small_checks.py - independent audit (opus-2026-11-04-S75 audit): small exact checks.

1. The almost transitive tournament AT_n (transitive, arc between the two extreme vertices reversed; Achour-Boudabbous-
   Boussairi, Reid-Thomassen) versus the bineutral tournament A_n of Burzio-Demaria 1988, Definition 2
   (arcs x_i -> x_j iff j < i - 1 or j = i + 1), which THM-4868 identifies with its run-length-2 class.
2. The four order-9 route-(b) witnesses v => C3 => x => y => C3' => v: 2-arc-strong, hc, forced arcs, |Aut|, converses,
   and identity with the four classes printed by the audit's exhaustive order-9 run (aud_census_exh_10.out).
3. Theorem A's sharpness family v => C3 => (transitive chain) => C3 => v at orders 7..12.
"""
from itertools import permutations


def tour_from_arcs(n, arcs):
    out = [0] * n
    for a, b in arcs:
        out[a] |= 1 << b
    for a in range(n):
        for b in range(n):
            if a != b and not (out[a] >> b) & 1 and not (out[b] >> a) & 1:
                raise ValueError("pair %d %d missing" % (a, b))
    return out


def strong(n, out, V=None):
    if V is None:
        V = (1 << n) - 1
    if V == 0:
        return True
    s = (V & -V).bit_length() - 1
    seen, fr = 1 << s, [s]
    while fr:
        x = fr.pop()
        for y in range(n):
            if (V >> y) & 1 and (out[x] >> y) & 1 and not (seen >> y) & 1:
                seen |= 1 << y
                fr.append(y)
    if seen != V:
        return False
    seen, fr = 1 << s, [s]
    while fr:
        x = fr.pop()
        for y in range(n):
            if (V >> y) & 1 and (out[y] >> x) & 1 and not (seen >> y) & 1:
                seen |= 1 << y
                fr.append(y)
    return seen == V


def hc_and_arc_counts(n, out):
    """all Hamiltonian cycles through vertex 0, enumerated explicitly (n <= 12)"""
    cycles = []
    path = [0]

    def rec(used):
        v = path[-1]
        if len(path) == n:
            if out[v] & 1:
                cycles.append(list(path))
            return
        nb = out[v] & ~used
        while nb:
            w = (nb & -nb).bit_length() - 1
            nb &= nb - 1
            path.append(w)
            rec(used | (1 << w))
            path.pop()
    rec(1)
    cnt = {}
    for cyc in cycles:
        for i in range(n):
            a, b = cyc[i], cyc[(i + 1) % n]
            cnt[(a, b)] = cnt.get((a, b), 0) + 1
    return len(cycles), cnt


def lambda2(n, out):
    for a in range(n):
        for b in range(n):
            if (out[a] >> b) & 1:
                o2 = list(out)
                o2[a] &= ~(1 << b)
                if not strong(n, o2):
                    return False
    return True


def noncritical(n, out):
    full = (1 << n) - 1
    return [v for v in range(n) if strong(n, out, full & ~(1 << v))]


def scores(out):
    return [bin(o).count("1") for o in out]


def isos(n, o1, o2, limit=None):
    """isomorphisms o1 -> o2 by score-respecting backtracking"""
    s1, s2 = scores(o1), scores(o2)
    if sorted(s1) != sorted(s2):
        return []
    res = []
    perm = [-1] * n
    used = [False] * n

    def rec(i):
        if limit and len(res) >= limit:
            return
        if i == n:
            res.append(list(perm))
            return
        for j in range(n):
            if used[j] or s2[j] != s1[i]:
                continue
            ok = True
            for k in range(i):
                if ((o1[i] >> k) & 1) != ((o2[j] >> perm[k]) & 1):
                    ok = False
                    break
            if ok:
                perm[i] = j
                used[j] = True
                rec(i + 1)
                used[j] = False
                perm[i] = -1
    rec(0)
    return res


def converse(out):
    n = len(out)
    c = [0] * n
    for a in range(n):
        for b in range(n):
            if (out[a] >> b) & 1:
                c[b] |= 1 << a
    return c


def cut_arcs_on_cycle(n, out, cyc):
    k = 0
    for i in range(n):
        a, b = cyc[i], cyc[(i + 1) % n]
        o2 = list(out)
        o2[a] &= ~(1 << b)
        if not strong(n, o2):
            k += 1
    return k


print("1. almost transitive AT_n versus the bineutral A_n (Burzio-Demaria Definition 2)")
for n in range(4, 11):
    AT = tour_from_arcs(n, [(i, j) for i in range(n) for j in range(i + 1, n) if not (i == 0 and j == n - 1)] + [(n - 1, 0)])
    BD = tour_from_arcs(n, [(i, j) for i in range(n) for j in range(n) if j < i - 1 or j == i + 1])
    hAT, _ = hc_and_arc_counts(n, AT)
    hBD, _ = hc_and_arc_counts(n, BD)
    ncAT, ncBD = noncritical(n, AT), noncritical(n, BD)
    iso = len(isos(n, AT, BD, limit=1)) > 0
    cAT = cut_arcs_on_cycle(n, AT, list(range(n)))
    cBD = cut_arcs_on_cycle(n, BD, list(range(n)))
    print("   n=%2d  AT_n: hc=%d non-critical=%d cut arcs=%d | A_n: hc=%d non-critical=%d cut arcs=%d | isomorphic: %s"
          % (n, hAT, len(ncAT), cAT, hBD, len(ncBD), cBD, iso))

print("2. the four order-9 witnesses v => C3 => x => y => C3' => v")
# labels: v=0, C3 = 1,2,3, x=4, y=5, C3' = 6,7,8
wit = {}
for vx in (0, 1):
    for vy in (0, 1):
        arcs = [(1, 2), (2, 3), (3, 1), (6, 7), (7, 8), (8, 6)]
        arcs += [(0, c) for c in (1, 2, 3)] + [(c, 0) for c in (6, 7, 8)]
        arcs += [(c, t) for c in (1, 2, 3) for t in (4, 5, 6, 7, 8)]
        arcs += [(4, 5)] + [(t, c) for t in (4, 5) for c in (6, 7, 8)]
        arcs += [(0, 4) if vx else (4, 0), (0, 5) if vy else (5, 0)]
        T = tour_from_arcs(9, arcs)
        h, cnt = hc_and_arc_counts(9, T)
        forced = sorted(a for a, c in cnt.items() if c == h)
        aut = len(isos(9, T, T))
        name = "v->x" if vx else "x->v"
        name += ", v->y" if vy else ", y->v"
        wit[name] = T
        print("   %-12s 2-arc-strong=%s  2-strong=%s  hc=%d  forced arcs=%s  |Aut|=%d  self-converse=%s"
              % (name, lambda2(9, T), all(strong(9, T, 511 & ~(1 << v)) for v in range(9)), h, forced, aut,
                 len(isos(9, T, converse(T), limit=1)) > 0))
names = list(wit)
for i in range(4):
    for j in range(i + 1, 4):
        a, b = wit[names[i]], wit[names[j]]
        print("   %s vs %s: isomorphic=%s, converse-isomorphic=%s" % (names[i], names[j], len(isos(9, a, b, limit=1)) > 0,
                                                                     len(isos(9, a, converse(b), limit=1)) > 0))
# the four classes printed by the audit's exhaustive run (rows of out-adjacency)
printed = [
    "010100000 001100000 100100000 000000111 111100000 111110000 111011010 111011001 111011100",
    "010010000 001010000 100010000 111000000 000100111 111110000 111101010 111101001 111101100",
    "010100000 001100000 100100000 000010111 111001000 111100000 111011010 111011001 111011100",
    "010001000 001001000 100001000 111000000 111100000 000110111 111110010 111110001 111110100",
]
matched = []
for rows in printed:
    out = [int(r[::-1], 2) for r in rows.split()]
    hit = [nm for nm in names if len(isos(9, out, wit[nm], limit=1)) > 0]
    matched.append(hit)
print("   exhaustive-run classes matched to the family:", matched)

print("3. sharpness family v => C3 => TT_j => C3' => v")
for j in range(0, 4):
    n = 7 + j
    v, C1, mid, C2 = 0, [1, 2, 3], list(range(4, 4 + j)), [4 + j, 5 + j, 6 + j]
    arcs = [(C1[0], C1[1]), (C1[1], C1[2]), (C1[2], C1[0]), (C2[0], C2[1]), (C2[1], C2[2]), (C2[2], C2[0])]
    arcs += [(v, c) for c in C1] + [(c, v) for c in C2]
    blocks = [C1] + [[m] for m in mid] + [C2]
    for bi in range(len(blocks)):
        for bj in range(bi + 1, len(blocks)):
            arcs += [(a, b) for a in blocks[bi] for b in blocks[bj]]
    arcs += [(v, m) for m in mid]
    T = tour_from_arcs(n, arcs)
    h, _ = hc_and_arc_counts(n, T)
    print("   n=%d: 2-arc-strong=%s, cut vertex=%s, hc=%d" % (n, lambda2(n, T), not all(strong(n, T, ((1 << n) - 1) & ~(1 << x)) for x in range(n)), h))
