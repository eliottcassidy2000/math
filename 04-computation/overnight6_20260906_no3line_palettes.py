"""Exact small controls for the palette-conditioned no-three-line theorem.

Standalone: no imports from the repository and no external census input.
Geometry paths are literal determinants and canonical pair-line counting.
Every check remains live under -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
import hashlib
import json
from pathlib import Path
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(predicate, label):
    global GATES
    GATES += 1
    if not predicate:
        raise RuntimeError(label)


def falling(n, r):
    return factorial(n) // factorial(n-r) if 0 <= r <= n else 0


def detzero(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) == (c[0]-a[0])*(b[1]-a[1])


def x_det(board):
    return sum(detzero(*t) for t in combinations(sorted(board), 3))


def x_lines(board):
    lines = defaultdict(set)
    for (x, y), (u, v) in combinations(sorted(board), 2):
        a, b = v-y, x-u
        d = gcd(a, b)
        a, b = a//d, b//d
        if a < 0 or (a == 0 and b < 0):
            a, b = -a, -b
        lines[(a, b, a*x+b*y)].update(((x, y), (u, v)))
    return sum(comb(len(points), 3) for points in lines.values())


def grid_events(n):
    return tuple(t for t in combinations(tuple(product(range(n), repeat=2)), 3)
                 if len({p[0] for p in t}) == 3
                 and len({p[1] for p in t}) == 3 and detzero(*t))


def stats(xs):
    moments = [sum((F(x)**j for x in xs), F(0))/len(xs) for j in (1, 2, 3)]
    a, b, c = moments
    return a, b-a*a, c-3*a*b+2*a*a*a


def factorial_k3(xs):
    mean = sum(map(F, xs), F(0))/len(xs)
    m2 = sum((F(x*(x-1)) for x in xs), F(0))/len(xs)
    m3 = sum((F(x*(x-1)*(x-2)) for x in xs), F(0))/len(xs)
    return m3-3*mean*m2+2*mean**3


def cycle(k):
    return frozenset((i, (i+j) % k) for i in range(k) for j in (0, 1))


def skeleton(parts):
    ans, start = set(), 0
    for k in parts:
        ans.update((r+start, c+start) for r, c in cycle(k))
        start += k
    return frozenset(ans)


def normalize(edges):
    rows = {x: i for i, x in enumerate(sorted({r for r, _ in edges}))}
    cols = {x: i for i, x in enumerate(sorted({c for _, c in edges}))}
    return tuple(sorted((rows[r], cols[c]) for r, c in edges))


@lru_cache(None)
def local_probability(k, edges):
    if not edges:
        return F(1)
    nr = 1+max(r for r, _ in edges)
    nc = 1+max(c for _, c in edges)
    if max(nr, nc) > k:
        return F(0)
    target = cycle(k)
    count = sum(all((rr[r], cc[c]) in target for r, c in edges)
                for rr in permutations(range(k), nr)
                for cc in permutations(range(k), nc))
    return F(count, falling(k, nr)*falling(k, nc))


def palette_probability(edges, rowsets, colsets):
    r_owner = {r: i for i, rs in enumerate(rowsets) for r in rs}
    c_owner = {c: i for i, cs in enumerate(colsets) for c in cs}
    groups = [set() for _ in rowsets]
    for r, c in edges:
        if r_owner[r] != c_owner[c]:
            return F(0)
        groups[r_owner[r]].add((r, c))
    value = F(1)
    for rs, edges_i in zip(rowsets, groups):
        value *= local_probability(len(rs), normalize(edges_i))
    return value


def matching_probability(k, s):
    if s > k:
        return F(0)
    if s == 0:
        return F(1)
    if s == 1:
        return F(2, k)
    if s == 2:
        return F(2*(2*k-3), k*(k-1)**2)
    if s == 3:
        return F(4*(2*k-5), k*(k-1)**2*(k-2))
    raise ValueError(s)


def palette_mean(events, rowsets, colsets):
    ro = {r: i for i, rs in enumerate(rowsets) for r in rs}
    co = {c: i for i, cs in enumerate(colsets) for c in cs}
    total = F(0)
    for event in events:
        if any(ro[r] != co[c] for r, c in event):
            continue
        sizes = Counter(ro[r] for r, _ in event)
        value = F(1)
        for i, s in sizes.items():
            value *= matching_probability(len(rowsets[i]), s)
        total += value
    return total


def family_boards(n, rr, cc):
    rs = tuple(r for r in range(n) if r not in rr)
    cs = tuple(c for c in range(n) if c not in cc)
    base = {(r, c) for r in rr for c in cc}
    if n == 4:
        boards = [frozenset(base | {(r, c) for r in rs for c in cs})]
    else:
        boards = [frozenset(base | {(r, c) for r in rs for c in cs
                                   if c != pp[rs.index(r)]})
                  for pp in permutations(cs)]
    return boards, (tuple(rr), rs), (tuple(cc), cs)


def rectangle_event_decomposition(rowsets, colsets):
    blocks=[frozenset(product(rs,cs)) for rs,cs in zip(rowsets,colsets)]
    paired=0
    for i,(rs,cs) in enumerate(zip(rowsets,colsets)):
        diagonals=(((rs[0],cs[0]),(rs[1],cs[1])),
                   ((rs[0],cs[1]),(rs[1],cs[0])))
        for j,block in enumerate(blocks):
            if i != j:
                paired += sum(detzero(a,b,p) for a,b in diagonals for p in block)
    three=sum(detzero(*points) for inds in combinations(range(len(blocks)),3)
              for points in product(*(blocks[i] for i in inds)))
    return paired+three


def small_palette_model(n):
    events = grid_events(n)
    data, all_x, board_counts = {}, [], Counter()
    for rr in combinations(range(n), 2):
        for cc in combinations(range(n), 2):
            boards, rowsets, colsets = family_boards(n, rr, cc)
            values = [x_det(board) for board in boards]
            board_counts.update(boards)
            all_x.extend(values)
            data[rr, cc] = values
            need(palette_mean(events, rowsets, colsets) == stats(values)[0],
                 f"all-palette mean {n,rr,cc}")
            for board, value in zip(boards, values):
                need(x_lines(board) == value, "two independent geometries")
                need(value == sum(set(event) <= board for event in events),
                     "nonaxis restriction equals all determinant triples")
                if n == 4:
                    need(rectangle_event_decomposition(rowsets,colsets) == value,
                         "all-C4 event occupancy decomposition with diagonal multiplicity")
            for event in events:
                exact = F(sum(set(event) <= board for board in boards), len(boards))
                need(palette_probability(event, rowsets, colsets) == exact,
                     "every singleton event local injection factor")

    # Independent sampling path: two full shore permutations of a fixed cycle graph.
    base = skeleton((2, n-2))
    full_by_palette = defaultdict(Counter)
    full_boards = Counter()
    geometry_cache = {}
    for rp in permutations(range(n)):
        for cp in permutations(range(n)):
            board = frozenset((rp[r], cp[c]) for r, c in base)
            if board not in geometry_cache:
                geometry_cache[board] = x_lines(board)
            key = tuple(sorted(rp[:2])), tuple(sorted(cp[:2]))
            full_by_palette[key][geometry_cache[board]] += 1
            full_boards[board] += 1
    need(set(full_boards) == set(board_counts), "independent board universes")
    expected_board_mult = 32 if n == 4 else 24
    need(set(full_boards.values()) == {expected_board_mult}, "full-label automorphisms")
    multiplier = 16 if n == 4 else 24
    for key, values in data.items():
        need(full_by_palette[key] == Counter({x: v*multiplier for x, v in Counter(values).items()}),
             "every full-label conditional histogram")

    conditionals = [stats(xs) for xs in data.values()]
    mus, variances, thirds = zip(*conditionals)
    avg = lambda seq: sum(seq, F(0))/len(seq)
    terms = (avg(thirds), stats(mus)[2],
             3*(avg([m*v for m, v in zip(mus, variances)])-avg(mus)*avg(variances)))
    mu, var, third = stats(all_x)
    need(sum(terms) == third, "third law of total cumulance")
    need(var == avg(variances)+stats(mus)[1], "total variance")
    if n == 4:
        need(Counter(all_x) == Counter({0:18, 2:8, 4:4, 6:4, 8:2}), "n4 exact palette histogram")
        need((mu, var, third) == (2, F(56, 9), 16), "n4 ordinary cumulants")
        need(terms == (0, 16, 0), "n4 palette-only third cumulant")
        need(factorial_k3(all_x) == F(4, 3), "n4 factorial third cumulant")
        for xs in data.values():
            need(stats(xs)[1:] == (0, 0), "C4 palette board determinism")
            need(factorial_k3(xs) == 2*xs[0], "factorial conditional diagonal survives")
    else:
        need((mu, var, third) == (F(13, 3), F(1769, 225), F(26159, 675)),
             "n5 exact ordinary cumulants")
        need(terms == (F(-871, 675), F(6467, 360), F(7949, 360)),
             "n5 three nonzero cumulance terms")
        # Several unions of events, including impossible degree-three and palette crossing cases.
        for rr, cc in [((0,1),(0,1)), ((0,2),(1,3)), ((1,4),(0,3))]:
            boards, rowsets, colsets = family_boards(n, rr, cc)
            sample = events[::max(1, len(events)//12)]
            for size in (1, 2, 3):
                for family in combinations(sample, size):
                    union = frozenset(p for event in family for p in event)
                    literal = F(sum(union <= board for board in boards), len(boards))
                    need(palette_probability(union, rowsets, colsets) == literal,
                         "whole-event union factorization")
    print(f"n={n}: palettes={len(data)}, distinct boards={len(board_counts)}, full labels={factorial(n)**2}")
    print("  histogram="+str(sorted(Counter(all_x).items())))
    print("  mean/variance/kappa3="+str(tuple(map(str, (mu, var, third)))))
    print("  total-cumulance terms="+str(tuple(map(str, terms))))
    return {"n":n, "hist":sorted(Counter(all_x).items()), "terms":list(map(str, terms))}


def conditional_cluster_controls():
    rs = ((0,1,2), (3,4,5))
    boards = []
    for p in permutations(rs[0]):
        for q in permutations(rs[1]):
            boards.append(frozenset((r,c) for i, pp in enumerate((p,q))
                                    for r in rs[i] for c in rs[i]
                                    if c != pp[rs[i].index(r)]))
    aa = frozenset((i,i) for i in (0,1,2))
    bb = frozenset((i,i) for i in (3,4,5))
    cc = frozenset((i,i) for i in (1,2,3))
    need(len(set(boards)) == 36, "complete conditional C6+C6 universe")
    need(all(detzero(*tuple(event)) for event in (aa,bb,cc)), "three real geometric events")
    def moment(family):
        union = frozenset(p for event in family for p in event)
        actual = F(sum(union <= board for board in boards), len(boards))
        need(actual == palette_probability(union, rs, rs), "n6 conditional joint injection probability")
        return actual
    def cumulant3(a,b,c):
        return (moment((a,b,c))-moment((a,b))*moment((c,))
                -moment((a,c))*moment((b,))-moment((b,c))*moment((a,))
                +2*moment((a,))*moment((b,))*moment((c,)))
    probs = [moment(fam) for fam in ((aa,), (bb,), (cc,), (aa,bb), (aa,cc), (bb,cc), (aa,bb,cc))]
    need(probs == [F(1,3),F(1,3),F(1,3),F(1,9),F(2,9),F(1,6),F(1,9)], "seven connected-kernel probabilities")
    need(cumulant3(aa, aa, bb) == 0, "disconnected repeated-event cumulant")
    need(cumulant3(aa, bb, cc) == F(1,54), "connected cross-cycle cumulant")
    need(moment((aa,bb)) == moment((aa,))*moment((bb,)), "two nonconstant independent indicators")
    need(palette_probability(((0,0),(0,1),(0,2)), rs, rs) == 0,
         "degree-three local hostile")
    need(palette_probability(((0,3),), rs, rs) == 0, "cross-palette edge hostile")
    print("n=6 conditional C6+C6: disconnected kappa(A,A,B)=0; connected kappa(A,B,C)=1/54")


def partial_maps(n):
    return tuple(tuple(zip(dom, vals))
                 for size in range(n+1) for dom in combinations(range(n), size)
                 for vals in permutations(range(n), size))


def merge(a,b):
    d = dict(a)
    for x,y in b:
        if (x in d and d[x] != y) or any(z != x and w == y for z,w in d.items()):
            return None
        d[x] = y
    return tuple(sorted(d.items()))


def cylinder_checks():
    n=3
    maps = partial_maps(n)
    perms = tuple(permutations(range(n)))
    evaluate = lambda a,p: all(p[x] == y for x,y in a)
    for a in maps:
        need(F(sum(evaluate(a,p) for p in perms), factorial(n)) == F(1,falling(n,len(a))),
             "cylinder falling-factorial state")
        for b in maps:
            c = merge(a,b)
            for p in perms:
                need(evaluate(a,p)*evaluate(b,p) == (0 if c is None else evaluate(c,p)),
                     "collision-sensitive cylinder product")
        if len(a)<n:
            x=next(x for x in range(n) if x not in dict(a))
            extensions = [tuple(sorted(a+((x,y),))) for y in range(n) if y not in dict(a).values()]
            for p in perms:
                need(sum(evaluate(b,p) for b in extensions) == evaluate(a,p), "full-assignment extension relation")
    need(merge(((0,0),), ((1,0),)) is None, "disjoint-domain duplicate-label hostile")
    # Two-shore square state, reconstructed both by the multiplication table and actual assignments.
    cylinders = [(maps[(7*j)%len(maps)], maps[(11*j+3)%len(maps)]) for j in range(24)]
    for seed in range(12):
        coeffs = [((seed+3)*(j+5) % 13)-6 for j in range(len(cylinders))]
        square = F(0)
        for a, (ra,ca) in zip(coeffs,cylinders):
            for b, (rb,cb) in zip(coeffs,cylinders):
                r, c = merge(ra,rb), merge(ca,cb)
                if r is not None and c is not None:
                    square += F(a*b,falling(n,len(r))*falling(n,len(c)))
        direct = sum(F(sum(a*evaluate(r,rp)*evaluate(c,cp)
                           for a,(r,c) in zip(coeffs,cylinders))**2)
                     for rp in perms for cp in perms)/factorial(n)**2
        need(square == direct and square >= 0, "positive genuine two-shore square")
    print("injective cylinder algebra: all 34 one-shore partial maps at n=3; positive two-shore squares pass")


def all_c4_and_matching_controls():
    for k in range(2,8):
        for s in range(4):
            matchings=sum(len({r for r,_ in ee})==s and len({c for _,c in ee})==s
                          for ee in combinations(cycle(k),s))
            prob=F(matchings*factorial(s),falling(k,s)**2) if s<=k else F(0)
            need(prob == matching_probability(k,s), "local cycle matching formula")
            need(prob == local_probability(k,tuple((i,i) for i in range(s))),
                 "matching counts versus independent injection sum")
    for m in range(1,9):
        n=2*m
        pairs = factorial(n)//(2**m*factorial(m))
        boards = pairs*pairs*factorial(m)
        need(boards*2**(2*m)*factorial(m) == factorial(n)**2,
             "all-C4 board multiplicity formula")
    for rowsets,colsets in [(((0,2),(1,4),(3,5)),((0,2),(1,4),(3,5))),
                            (((0,1),(2,4),(3,5)),((1,5),(0,3),(2,4)))]:
        board=frozenset(p for rs,cs in zip(rowsets,colsets) for p in product(rs,cs))
        need(rectangle_event_decomposition(rowsets,colsets) == x_det(board),
             "three-C4 two-block plus three-block exact decomposition")
    print("all-C4 board count and s<=3 local matching kernels checked")


def main():
    all_c4_and_matching_controls()
    cylinder_checks()
    data = [small_palette_model(n) for n in (4,5)]
    conditional_cluster_controls()
    digest=hashlib.sha256(json.dumps(data,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    print("semantic_sha256="+digest)
    print(f"PASS: {GATES} optimization-live exact gates")


if __name__ == "__main__":
    main()
