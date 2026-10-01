#!/usr/bin/env python3
"""Which tournaments are 'designated special' by oriented Hamiltonian cycles and paths (opus S15, 2026-10-01).

A spanning oriented graph D designates the n-tournaments that avoid it.  Grunbaum's theorem: H_n (Hamiltonian path +
arc first->last; the two-block cycle with blocks 1 and n-1) designates exactly C3 and C3[1,C3,1].
Checks (all iso classes, n <= 8):
  1. two-block cycles D(n; p, n-p) (a source and a sink joined by directed paths with p and n-p arcs): the designated
     tournaments are exactly C3 (n=3, p=1), TT2[1,C3] and TT2[C3,1] (n=4, p=2), C3[1,C3,1] (n=5, p=1) and
     TT2[C3,C3] = 'C3 => C3' (n=6, p=2); nothing for n = 7, 8.  All are compositions of 3-cycles.
  2. all oriented Hamiltonian paths: only the antidirected ones, avoided only by C3, R5 (the regular 5-tournament)
     and P7 (Paley) -- Grunbaum / Rosenfeld / Havet-Thomasse
  3. all non-directed oriented Hamiltonian cycles: the designated classes; in particular the 4-block 6-cycles
     001011 and 001101 are avoided ONLY by G_par (two 3-cycles glued 'in parallel' by a matching), 00101 only by
     R5, 0010101 only by P7; P7 minus a vertex (the antiparallel gluing) is one of three avoiders of 001001
  4. path + 1 or 2 extra arcs, n <= 7: H_6 + (0->3) is avoided only by G_par; H_7 + (0->5) only by P7
  5. the two matching gluings of two 3-cycles: antiparallel = P7 minus a vertex (all arcs odd), parallel = G_par
Reproduce: python 04-computation/experiments/tournament_designations_20261001.py   (about 15 minutes)
"""
from collections import Counter
from itertools import combinations, product

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ---------------- tournament library (same conventions as shaved_tournaments_20261001.py)
def refine(n, out, colors):
    while True:
        sigs = [(colors[v], tuple(sorted(colors[w] for w in range(n) if out[v] >> w & 1))) for v in range(n)]
        keys = sorted(set(sigs))
        idx = {k: i for i, k in enumerate(keys)}
        new = [idx[s] for s in sigs]
        if len(keys) == len(set(colors)):
            return new
        colors = new


def code_of(n, out, order):
    c, b = 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            if out[order[i]] >> order[j] & 1:
                c |= 1 << b
            b += 1
    return c


def canon(n, out):
    best = [None, 0]

    def rec(colors):
        colors = refine(n, out, colors)
        if len(set(colors)) == n:
            c = code_of(n, out, sorted(range(n), key=lambda v: colors[v]))
            if best[0] is None or c < best[0]:
                best[0], best[1] = c, 1
            elif c == best[0]:
                best[1] += 1
            return
        cnt = Counter(colors)
        target = min(c for c in cnt if cnt[c] > 1)
        for v in range(n):
            if colors[v] == target:
                nc = [2 * c for c in colors]
                nc[v] = 2 * target - 1
                rec(nc)
    rec([bin(out[v]).count('1') for v in range(n)])
    return best[0], best[1]


def from_code(n, c):
    out = [0] * n
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if c >> b & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1
    return out


def all_classes(nmax):
    res = {1: [(0, 1)]}
    for n in range(2, nmax + 1):
        seen = {}
        for code, _ in res[n - 1]:
            base = from_code(n - 1, code)
            for S in range(1 << (n - 1)):
                out = base[:] + [0]
                for v in range(n - 1):
                    if S >> v & 1:
                        out[n - 1] |= 1 << v
                    else:
                        out[v] |= 1 << (n - 1)
                c, a = canon(n, out)
                seen.setdefault(c, a)
        res[n] = sorted(seen.items())
    return res


def make_embedder(n, arcs):
    outS, inS = [0] * n, [0] * n
    for a, b in arcs:
        outS[a] |= 1 << b
        inS[b] |= 1 << a
    order, chosen, rest = [], 0, set(range(n))
    while rest:
        v = max(rest, key=lambda x: (bin((outS[x] | inS[x]) & chosen).count('1'), bin(outS[x] | inS[x]).count('1'), -x))
        order.append(v)
        rest.remove(v)
        chosen |= 1 << v
    pos = {v: i for i, v in enumerate(order)}
    req_out = [[pos[w] for w in range(n) if outS[v] >> w & 1 and pos[w] < i] for i, v in enumerate(order)]
    req_in = [[pos[w] for w in range(n) if inS[v] >> w & 1 and pos[w] < i] for i, v in enumerate(order)]

    def emb(outT):
        inT = [0] * n
        for v in range(n):
            o = outT[v]
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                inT[w] |= 1 << v
        img = [0] * n

        def rec(i, used):
            if i == n:
                return True
            cand = ((1 << n) - 1) & ~used
            for j in req_out[i]:
                cand &= inT[img[j]]
            for j in req_in[i]:
                cand &= outT[img[j]]
            while cand:
                y = (cand & -cand).bit_length() - 1
                cand &= cand - 1
                img[i] = y
                if rec(i + 1, used | 1 << y):
                    return True
            return False
        return rec(0, 0)
    return emb


C = all_classes(8)
T = {n: [from_code(n, c) for c, _ in C[n]] for n in range(3, 9)}
IDX = {n: {c: i for i, (c, _) in enumerate(C[n])} for n in range(3, 9)}


def cls(out):
    n = len(out)
    return IDX[n][canon(n, out)[0]]


def avoiders(n, arcs):
    emb = make_embedder(n, arcs)
    return [i for i, t in enumerate(T[n]) if not emb(t)]


# ---------------- named tournaments
def compose(frame, parts):
    """lexicographic product frame[parts]: frame = out-masks on k vertices, parts = list of out-mask lists"""
    sizes = [len(p) for p in parts]
    offs = [sum(sizes[:i]) for i in range(len(parts))]
    n = sum(sizes)
    out = [0] * n
    for i, P in enumerate(parts):
        for a in range(len(P)):
            v = offs[i] + a
            for b in range(len(P)):
                if P[a] >> b & 1:
                    out[v] |= 1 << (offs[i] + b)
            for j in range(len(parts)):
                if frame[i] >> j & 1:
                    for b in range(sizes[j]):
                        out[v] |= 1 << (offs[j] + b)
    return out


ONE = [0]
C3 = [0b010, 0b100, 0b001]
TT2 = [0b10, 0b00]
circ = lambda m, S: [sum(1 << ((x + s) % m) for s in S) for x in range(m)]
P7 = circ(7, (1, 2, 4))
R5 = circ(5, (1, 2))
QR6 = [sum(1 << (y - 1) for y in range(1, 7) if P7[x] >> y & 1) for x in range(1, 7)]


def glue(sign):
    """two 3-cycles A = {a_i}, B = {b_i}: a_i -> a_(i+1), b_i -> b_(i+sign), matching a_i -> b_i, all other
    cross arcs b_j -> a_i; sign = +1 'parallel', -1 'antiparallel'"""
    out = [0] * 6
    for i in range(3):
        out[i] |= 1 << ((i + 1) % 3)
        out[3 + i] |= 1 << (3 + (i + sign) % 3)
        out[i] |= 1 << (3 + i)
        for j in range(3):
            if j != i:
                out[3 + j] |= 1 << i
    return out


GPAR, GANTI = glue(+1), glue(-1)
named = {
    3: {cls(C3): 'C3'},
    4: {cls(compose(TT2, [ONE, C3])): 'TT2[1,C3]', cls(compose(TT2, [C3, ONE])): 'TT2[C3,1]'},
    5: {cls(compose(C3, [ONE, C3, ONE])): 'C3[1,C3,1]', cls(R5): 'R5'},
    6: {cls(compose(TT2, [C3, C3])): 'TT2[C3,C3]', cls(QR6): 'P7-v', cls(GPAR): 'G_par'},
    7: {cls(P7): 'P7'},
    8: {},
}


def name(n, i):
    return named[n].get(i, f'class{i}')


# ---------------- 1. two-block cycles
print("=== 1. two-block cycles D(n; p, n-p) ===")
blocks = {}
for n in range(3, 9):
    for p in range(1, n // 2 + 1):
        q = n - p
        t = n - 1
        A = [0] + list(range(1, p)) + [t]
        B = [0] + list(range(p, p + q - 1)) + [t]
        arcs = [(A[k], A[k + 1]) for k in range(len(A) - 1)] + [(B[k], B[k + 1]) for k in range(len(B) - 1)]
        av = avoiders(n, arcs)
        blocks[(n, p)] = sorted(name(n, a) for a in av)
        if av:
            print(f"   n={n} (p,q)=({p},{q}): avoided by {blocks[(n, p)]}")
expect = {(3, 1): ['C3'], (4, 2): ['TT2[1,C3]', 'TT2[C3,1]'], (5, 1): ['C3[1,C3,1]'], (6, 2): ['TT2[C3,C3]']}
check(all(blocks[k] == expect.get(k, []) for k in blocks),
      "two-block cycles designate exactly C3, TT2[1,C3], TT2[C3,1], C3[1,C3,1], TT2[C3,C3] (n <= 8): 3-cycle compositions")

# ---------------- 2. oriented Hamiltonian paths
print("=== 2. oriented Hamiltonian paths ===")
path_exc = []
for n in range(3, 9):
    seen = set()
    for w in product((0, 1), repeat=n - 1):
        cw = min([w, w[::-1], tuple(1 - b for b in w), tuple(1 - b for b in w[::-1])])
        if cw in seen:
            continue
        seen.add(cw)
        arcs = [(i, i + 1) if b else (i + 1, i) for i, b in enumerate(cw)]
        av = avoiders(n, arcs)
        if av:
            anti = all(cw[i] != cw[i + 1] for i in range(len(cw) - 1))
            path_exc.append((n, anti, sorted(name(n, a) for a in av)))
print("   ", path_exc)
check(path_exc == [(3, True, ['C3']), (5, True, ['R5']), (7, True, ['P7'])],
      "the only unavoidable-failures among oriented Hamiltonian paths are the antidirected paths in C3, R5, P7")

# ---------------- 3. oriented Hamiltonian cycles
print("=== 3. oriented (non-directed) Hamiltonian cycles ===")


def canon_cycle_word(w):
    n = len(w)
    c = []
    for s in range(n):
        r = w[s:] + w[:s]
        c.append(r)
        c.append(tuple(1 - b for b in reversed(r)))
    return min(c)


cyc = {}
for n in range(3, 9):
    seen = set()
    for w in product((0, 1), repeat=n):
        if sum(w) in (0, n):
            continue
        cw = canon_cycle_word(w)
        if cw in seen:
            continue
        seen.add(cw)
        arcs = [(i, (i + 1) % n) if b else ((i + 1) % n, i) for i, b in enumerate(cw)]
        av = avoiders(n, arcs)
        if av:
            cyc[(n, ''.join(map(str, cw)))] = sorted(name(n, a) for a in av)
for k, v in sorted(cyc.items()):
    print(f"   n={k[0]} word {k[1]}: {len(v)} avoiders" + (f" {v}" if len(v) <= 4 else ''))
check(cyc.get((6, '001011')) == ['G_par'] and cyc.get((6, '001101')) == ['G_par'] and cyc.get((5, '00101')) == ['R5']
      and cyc.get((7, '0010101')) == ['P7'] and cyc.get((6, '000011')) == ['TT2[C3,C3]']
      and cyc.get((5, '00001')) == ['C3[1,C3,1]'],
      "4-block 6-cycles 001011, 001101 avoided only by G_par (parallel gluing of two 3-cycles); 00101 only by R5; "
      "0010101 only by P7")
check('P7-v' in cyc.get((6, '001001'), []) and len(cyc.get((6, '001001'), [])) == 3,
      "P7 minus a vertex (= antiparallel gluing G_anti) is one of 3 avoiders of 001001, never a unique avoider here")
check(sorted(len(v) for k, v in cyc.items() if k[0] == 8) == [19] and all(set(k[1]) and k[1] == '01' * 4 for k in cyc if k[0] == 8),
      "at n = 8 only the antidirected cycle is avoidable (19 classes); nothing for n = 7 except 0010101")

# ---------------- 4. path + 1 or 2 extra arcs
print("=== 4. Hamiltonian path + 1 or 2 extra arcs ===")
uniq = {}
for n in range(4, 8):
    path = [(k, k + 1) for k in range(n - 1)]
    extra = [(i, j) for i in range(n) for j in range(i + 2, n)] + [(j, i) for i in range(n) for j in range(i + 2, n)]
    for k in (1, 2):
        for ex in combinations(extra, k):
            if len({frozenset(e) for e in ex}) < k:
                continue
            av = avoiders(n, path + list(ex))
            if len(av) == 1:
                uniq.setdefault((n, name(n, av[0])), []).append(ex)
for key in sorted(uniq):
    print(f"   n={key[0]} uniquely designated {key[1]}: {len(uniq[key])} graphs, e.g. {uniq[key][0]}")
check(((0, 3), (0, 5)) in uniq.get((6, 'G_par'), []) and ((0, 5), (0, 6)) in uniq.get((7, 'P7'), [])
      and (6, 'P7-v') not in uniq,
      "H_6 + (0->3) is avoided only by G_par; H_7 + (0->5) only by P7; P7 minus a vertex is never uniquely designated")

# ---------------- 5. the two matching gluings of two 3-cycles
print("=== 5. two 3-cycles glued by a matching ===")


def hp_count_and_arcs(out):
    n = len(out)
    cnt = Counter()
    H = 0
    from itertools import permutations as perms
    for p in perms(range(n)):
        if all(out[p[i]] >> p[i + 1] & 1 for i in range(n - 1)):
            H += 1
            for i in range(n - 1):
                cnt[(p[i], p[i + 1])] += 1
    return H, sum(1 for e in cnt if cnt[e] % 2 == 1)


def extends_to(out, target):
    n = len(out)
    tc = cls(target)
    for S in range(1 << n):
        o = out[:] + [0]
        for v in range(n):
            if S >> v & 1:
                o[n] |= 1 << v
            else:
                o[v] |= 1 << n
        if cls(o) == tc:
            return True
    return False


hpar, opar = hp_count_and_arcs(GPAR)
hant, oant = hp_count_and_arcs(GANTI)
check(cls(GANTI) == cls(QR6) and cls(GPAR) != cls(GANTI) and hpar == hant == 45 and oant == 15 and opar == 9
      and extends_to(GANTI, P7) and not extends_to(GPAR, P7),
      "antiparallel gluing = P7 minus a vertex (all 15 arcs on an odd number of HPs, THM-4524); parallel gluing G_par "
      "has the same scores, |Aut| = 3, H = 45, but 9/15 odd arcs and no extension to P7")

print()
print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
