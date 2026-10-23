"""Exact cyclic attachment lattice and shared-boundary tests; stdlib only."""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd, prod
from pathlib import Path
from random import Random
import json
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def require(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def forms(m, n):
    return [(h, s, v) for h in range(1, m+1) if m % h == 0
            for v in range(1, n+1) if n % v == 0
            for s in range(h) if s*(n//v) % h == 0]


def member(x, y, f):
    h, s, v = f
    return y % v == 0 and (v*x-s*y) % (h*v) == 0


def subgroup(m, n, f):
    return frozenset((x, y) for x in range(m) for y in range(n)
                     if member(x, y, f))


def translate(H, o, m, n):
    return frozenset(((x+o[0]) % m, (y+o[1]) % n) for x, y in H)


def normal_form(H, m, n):
    h = next(x for x in range(1, m+1) if (x % m, 0) in H)
    v = next(y for y in range(1, n+1)
             if any((x, y % n) in H for x in range(m)))
    s = next(x for x in range(h) if (x, v % n) in H)
    return h, s, v


def reduce_pair(x, y, f):
    h, s, v = f
    q, r = divmod(y, v)
    return (x-s*q) % h, r


def all_cosets(m, n):
    output = []
    for f in forms(m, n):
        H = subgroup(m, n, f)
        remaining = set(product(range(m), range(n)))
        while remaining:
            o = min(remaining)
            R = translate(H, o, m, n)
            output.append((f, o, R))
            remaining.difference_update(R)
    return output


def determinant(A):
    """Independent integer Bareiss elimination for minor certificates."""
    n = len(A)
    if not n:
        return 1
    A = [list(row) for row in A]
    old = 1
    sign = 1
    for k in range(n-1):
        p = next((r for r in range(k, n) if A[r][k]), None)
        if p is None:
            return 0
        if p != k:
            A[k], A[p] = A[p], A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = A[i][j]*pivot-A[i][k]*A[k][j]
                if numerator % old:
                    raise RuntimeError("nonintegral Bareiss step")
                A[i][j] = numerator//old
        for i in range(k+1, n):
            A[i][k] = 0
        old = pivot
    return sign*A[-1][-1]


def minor_gcd(M, k):
    g = 0
    for rr in combinations(range(len(M)), k):
        for cc in combinations(range(len(M[0])), k):
            g = gcd(g, determinant([[M[r][c] for c in cc] for r in rr]))
            if g == 1:
                return 1
    return abs(g)


def smith_factors(M):
    divisors = [1]+[minor_gcd(M, k) for k in range(1, len(M)+1)]
    return [divisors[k]//divisors[k-1] for k in range(1, len(divisors))]


def presentation(vertices, edges, fs):
    size = 2*len(edges)
    M = [[0]*(size+vertices) for _ in range(size)]
    for e, ((i, j), (h, s, v)) in enumerate(zip(edges, fs)):
        M[2*e][2*e] = h
        M[2*e][2*e+1] = s
        M[2*e+1][2*e+1] = v
        M[2*e][size+i] = 1
        M[2*e+1][size+j] = 1
    return M


def network(orders, edges, fs):
    alphabets = [tuple(product(range(h), range(v))) for h, s, v in fs]
    configurations = set(product(*alphabets))
    image = Counter()
    for a in product(*(range(m) for m in orders)):
        image[tuple(reduce_pair(a[i], a[j], f)
                    for (i, j), f in zip(edges, fs))] += 1
    kernel = image[tuple((0, 0) for _ in edges)]
    require(set(image.values()) == {kernel}, "global witnesses uniform")
    remaining = set(configurations)
    orbit_sizes = Counter()
    while remaining:
        seed = min(remaining)
        orbit = {tuple(reduce_pair(x+u, y+v, f)
                       for (x, y), (u, v), f in zip(seed, action, fs))
                 for action in image}
        remaining.difference_update(orbit)
        orbit_sizes[len(orbit)] += 1
    M = presentation(len(orders), edges, fs)
    factors = smith_factors(M)
    index = prod(factors)
    require(index == sum(orbit_sizes.values()), "Smith/orbit class count")
    require(kernel*prod(h*v for h, s, v in fs) == prod(orders)*index,
            "global kernel formula")
    for offset in configurations:
        flat = [x for pair in offset for x in pair]
        augmented = [row+[x] for row, x in zip(M, flat)]
        require((minor_gcd(augmented, len(M)) == index) == (offset in image),
                "integer lattice membership vs actual shared charts")
    return dict(orders=orders, edges=edges, forms=fs,
                configurations=len(configurations), image=len(image),
                witness_kernel=kernel, smith=factors,
                orbit_sizes={str(k): v for k, v in sorted(orbit_sizes.items())})


def cyclic_tournament(n):
    return tuple(tuple(int(i != j and 1 <= (j-i) % n <= (n-1)//2)
                       for j in range(n)) for i in range(n))


def three_blowup(r):
    R = cyclic_tournament(r)
    return tuple(tuple(R[i//3][j//3] if i % 3 == j % 3
                       else int((j-i) % 3 == 1)
                       for j in range(3*r)) for i in range(3*r))


def preserves(A, p):
    return all(A[i][j] == A[p[i]][p[j]]
               for i in range(len(A)) for j in range(len(A)))


def equality_block(m, n):
    return tuple(tuple(int(i % 3 == j % 3) for j in range(n))
                 for i in range(m))


def shifted(B, a, b):
    m, n = len(B), len(B[0])
    return tuple(tuple(B[(i-a) % m][(j-b) % n] for j in range(n))
                 for i in range(m))


def check_forms():
    total_forms = total_cosets = 0
    for m in range(1, 10):
        for n in range(1, 10):
            found = set()
            for f in forms(m, n):
                h, s, v = f
                H = subgroup(m, n, f)
                require(H not in found, "distinct subgroup normal forms")
                found.add(H)
                require(normal_form(H, m, n) == f, "normal form recovery")
                require(len(H)*h*v == m*n, "pair count")
                left = sum(y == 0 for x, y in H)
                right = sum(x == 0 for x, y in H)
                require(left == m//h, "left kernel")
                require(right*h*v == n*gcd(h, s), "right kernel")
                require(all(((x+a) % m, (y+b) % n) in H
                            for x, y in H for a, b in H), "subgroup closure")
                total_forms += 1
            for f, o, R in all_cosets(m, n):
                require(all(((x, y) in R) == member(x-o[0], y-o[1], f)
                            for x in range(m) for y in range(n)),
                        "affine two-congruence test")
                total_cosets += 1
    return dict(orders_through=9, subgroups=total_forms, cosets=total_cosets)


def check_composition():
    trials = empty = multivalued = 0
    for a, b, c in product(range(1, 4), repeat=3):
        for f, o, R in all_cosets(a, b):
            right = b*gcd(f[0], f[1])//(f[0]*f[2])
            for g, p, S in all_cosets(b, c):
                left = b//g[0]
                counts = Counter((x, z) for x, y in R for yy, z in S if y == yy)
                trials += 1
                if not counts:
                    empty += 1
                    continue
                require(set(counts.values()) == {gcd(right, left)},
                        "hidden-middle fibre multiplicity")
                multivalued += gcd(right, left) > 1
                offset = min(counts)
                H = translate(frozenset(counts), (-offset[0], -offset[1]), a, c)
                hnf = normal_form(H, a, c)
                require(subgroup(a, c, hnf) == H, "composition affine subgroup")
    return dict(orders_through=3, trials=trials, empty=empty,
                with_multiple_middle_witnesses=multivalued)


def exact_elimination(orders, relations, domains=None):
    """Set-based hostile checker for the compact CRT elimination theorem."""
    U = {i: set(range(m)) for i, m in enumerate(orders)} if domains is None else {
        i: set(d) for i, d in enumerate(domains)}
    R = {e: set(r) for e, r in relations.items()}
    weight = 1
    factors = []
    peak_edges = len(R)
    for y in range(len(orders)):
        if any(not d for d in U.values()) or any(not r for r in R.values()):
            return 0, factors, peak_edges
        neighbors = sorted({j if i == y else i for i, j in R if y in (i, j)})
        if not neighbors:
            weight *= len(U[y])
            factors.append(len(U[y]))
            del U[y]
            continue
        owned = {}
        kernel_sizes = []
        for x in neighbors:
            e = (min(x, y), max(x, y))
            oriented = {(a, b) if y < x else (b, a) for a, b in R[e]}
            owned[x] = {(a, b) for a, b in oriented if a in U[y]}
            if not owned[x]:
                return 0, factors, peak_edges
            counts = Counter(b for a, b in owned[x])
            require(len(set(counts.values())) == 1, "constant incident cyclic kernel")
            kernel_sizes.append(next(iter(counts.values())))
            U[x].intersection_update(counts)
        factor = 0
        for k in kernel_sizes:
            factor = gcd(factor, k)
        weight *= factor
        factors.append(factor)
        for x, z in combinations(neighbors, 2):
            composed = {(a, c) for b, a in owned[x] for bb, c in owned[z] if b == bb}
            e = (x, z)
            R[e] = R[e] & composed if e in R else composed
        for e in list(R):
            if y in e:
                del R[e]
        del U[y]
        peak_edges = max(peak_edges, len(R))
    return weight, factors, peak_edges


def check_fill_in():
    helly_trials = 0
    for b in range(1, 13):
        cosets = [frozenset(y for x, y in R) for f, o, R in all_cosets(1, b)]
        for C in product(cosets, repeat=3):
            pairwise = all(C[i] & C[j] for i, j in combinations(range(3), 2))
            joint = C[0] & C[1] & C[2]
            require(bool(joint) == pairwise, "cyclic CRT Helly two")
            if joint:
                require(len(joint) == gcd(gcd(len(C[0]), len(C[1])), len(C[2])),
                        "cyclic triple multiplicity")
            helly_trials += 1
    # Odd-order abelian boundary: pairwise intersections do not suffice.
    grid = set(product(range(3), repeat=2))
    bad = [{(x, y) for x, y in grid if x == 0},
           {(x, y) for x, y in grid if y == 0},
           {(x, y) for x, y in grid if (x+y) % 3 == 1}]
    require(all(bad[i] & bad[j] for i, j in combinations(range(3), 2))
            and not set.intersection(*bad), "C3 squared Helly hostile")
    rng = Random(231023)
    trials = 600
    positive = multiple = fill_in = 0
    max_count = 0
    for test in range(trials):
        n = rng.randrange(2, 7)
        orders = [rng.choice((1, 2, 3, 4, 5, 6, 9)) for _ in range(n)]
        relations = {}
        for i, j in combinations(range(n), 2):
            if rng.randrange(3) == 0:
                continue
            m, k = orders[i], orders[j]
            f = rng.choice(forms(m, k))
            o = (rng.randrange(m), rng.randrange(k))
            relations[(i, j)] = translate(subgroup(m, k, f), o, m, k)
        domains = []
        for m in orders:
            f = rng.choice(forms(1, m))
            o = (0, rng.randrange(m))
            domains.append({y for x, y in translate(subgroup(1, m, f), o, 1, m)})
        direct = sum(all((a[i], a[j]) in r for (i, j), r in relations.items())
                     for a in product(*domains))
        count, factors, peak = exact_elimination(orders, relations, domains)
        require(count == direct, "weighted fill-in vs direct shared choices")
        positive += direct > 0
        multiple += direct > 1
        max_count = max(max_count, direct)
        fill_in += peak > len(relations)
    # A deliberately dense fill-in star, without cycle shortcuts.
    orders = [9, 3, 9, 15, 3]
    rels = {(0, j): subgroup(9, orders[j], (3, 1, 1)) for j in range(1, 5)}
    count, factors, peak = exact_elimination(orders, rels)
    require(count == 135 and factors[0] == 3 and peak == 6,
            "four-port cyclic contraction retains middle multiplicity")
    return dict(helly_triples=helly_trials, seed=231023, network_trials=trials,
                positive=positive, multiple=multiple, max_count=max_count,
                trials_with_added_fill_in=fill_in,
                four_port_star=dict(orders=orders, count=count,
                                    factors=factors, peak_edges=peak),
                noncyclic_hostile="C3^2 cosets x=0,y=0,x+y=1")


def actual_controls():
    cells = {3*r: three_blowup(r) for r in (1, 3, 5)}
    for m, A in cells.items():
        require(all(sum(row) == (m-1)//2 for row in A), "regular actual cell")
        require(all(A[i][j]+A[j][i] == (i != j)
                    for i in range(m) for j in range(m)), "actual tournament")
        require(preserves(A, tuple((i+1) % m for i in range(m))),
                "authenticated cyclic generator")
    rows = []
    for m, n in product(cells, repeat=2):
        B = equality_block(m, n)
        require({sum(row) for row in B} == {n//3}, "row regularity")
        require({sum(B[i][j] for i in range(m)) for j in range(n)} == {m//3},
                "column regularity")
        H = frozenset((a, b) for a in range(m) for b in range(n)
                      if shifted(B, a, b) == B)
        require(normal_form(H, m, n) == (3, 1, 1), "actual block HNF")
        for o in ((0, 0), (1, 0), (0, 1)):
            target = shifted(B, *o)
            R = frozenset((a, b) for a in range(m) for b in range(n)
                          if shifted(B, a, b) == target)
            require(R == translate(H, o, m, n), "actual comparison relation")
        rows.append(dict(m=m, n=n, compatible=len(H),
                         target_above_source=sum(a == 0 for a, b in H)))

    R9 = cells[9]
    autos = [p for p in permutations(range(9)) if preserves(R9, p)]
    require(len(autos) == 81, "full R9 group independently enumerated")
    B = equality_block(9, 9)
    pair_count = sum(all(B[i][j] == B[p[i]][q[j]]
                         for i in range(9) for j in range(9))
                     for p in autos for q in autos)
    require(pair_count == 2187, "full R9 compatible pairs")
    one_sided = sum(all(B[i][j] == B[i][q[j]]
                        for i in range(9) for j in range(9)) for q in autos)
    require(one_sided == 27, "full R9 right kernel")
    return dict(cyclic_subcarrier=rows, full_R9_group=len(autos),
                full_R9_pairs=pair_count, full_R9_right_kernel=one_sided)


def nonabelian_control():
    squares = {1, 2, 4}
    A = tuple(tuple(int((j-i) % 7 in squares) for j in range(7))
              for i in range(7))
    charts = [tuple((pow(2, k, 7)*x+b) % 7 for x in range(7))
              for k in range(3) for b in range(7)]
    require(all(preserves(A, p) for p in charts), "Paley affine charts")
    def compose(p, q):
        return tuple(p[q[x]] for x in range(len(p)))
    t = charts[1]
    d = charts[7]
    inverse_d = tuple(d.index(x) for x in range(7))
    require(compose(compose(d, t), inverse_d) == compose(t, t), "nonabelian conjugacy")
    require(compose(d, t) != compose(t, d), "noncommuting actual charts")
    # B'= (d,id) I, so relation is a translated diagonal stabilizer.
    B = tuple(tuple(int(i == j) for j in range(7)) for i in range(7))
    target = tuple(tuple(B[inverse_d[i]][j] for j in range(7)) for i in range(7))
    R = frozenset((i, j) for i, p in enumerate(charts) for j, q in enumerate(charts)
                  if all(B[x][y] == target[p[x]][q[y]]
                         for x in range(7) for y in range(7)))
    require(len(R) == 21, "actual translated matching relation")
    o = min(R)
    H = translate(R, (-o[0], -o[1]), 21, 21)
    hostile = next(((a, b), (c, d), ((a+c) % 21, (b+d) % 21))
                   for a, b in H for c, d in H
                   if ((a+c) % 21, (b+d) % 21) not in H)
    require(hostile[2] not in H, "numbered nonabelian charts not affine subgroup")
    return dict(group_order=21, natural_label="7*k+b for x -> 2^k*x+b",
                relation_offset=o, subgroup_closure_failure=hostile)


def main():
    result = dict(status="FINITE-EXACT; general scoped proofs in README.md",
                  normal_forms=check_forms(), composition=check_composition(),
                  cyclic_fill_in=check_fill_in(),
                  actual_tournaments=actual_controls(),
                  nonabelian=nonabelian_control())
    triangle = [(0, 1), (1, 2), (0, 2)]
    result["networks"] = [
        network([3, 9, 15], triangle, [(3, 1, 1)]*3),
        network([9, 9, 9], triangle, [(9, 1, 1)]*3),
        network([3, 9, 15], [(0, 1), (1, 2)], [(3, 1, 1)]*2),
        network([9, 9, 9], triangle, [(3, 1, 1), (9, 1, 1), (9, 1, 1)]),
    ]
    result["checks"] = CHECKS
    result = json.loads(json.dumps(result))
    target = HERE / "receipt.json"
    if "--write" in sys.argv:
        target.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    else:
        require(json.loads(target.read_text(encoding="utf-8")) == result,
                "frozen receipt exact replay")
    print(json.dumps(dict(status="PASS", checks=CHECKS,
                          normal_forms=result["normal_forms"],
                          composition=result["composition"],
                          network_smith=[n["smith"] for n in result["networks"]])))


if __name__ == "__main__":
    main()
