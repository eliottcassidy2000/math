"""Exact controls for dual pair-gcd LRC closure; standard library only.

The analytic theorem is in the companion report.  Small integer boxes are
exhaustive controls; the explicitly listed physical rows use actual atlas
components and both orientations of every mixed triple.  No producer import.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, prod, isqrt
import sys

sys.stdout.reconfigure(newline='\n')
Q = 91**6
H = 3*Q//28
G = {7: 90, 8: 30, 9: 9, 10: 4, 11: 2, 12: 1}
gates = 0


def check(ok, message):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(message)


def ceildiv(a, b):
    return -((-a)//b)


def box(B, p, q, x):
    """Return (r,s) with rp+sq=x, or None. Parametrize s modulo p."""
    s = pow(q, -1, p)*x % p if p > 1 else 0
    r = (x-q*s)//p
    lo = max(ceildiv(-B-s, p), ceildiv(r-B, q))
    hi = min((B-s)//p, (r+B)//q)
    if lo > hi:
        return None
    return r-lo*q, s+lo*p


def atlas(x, y):
    d = gcd(x, y)
    p, q = sorted((x//d, y//d))
    n = p+q
    if n > 356:
        return False
    k = 2
    while k*k <= n:
        e = 0
        while n % k == 0:
            e += 1
            n //= k
        if e and (k % 3 != 2 or e > 2):
            return False
        k += 1
    return n == 1 or n % 3 == 2


def components(v):
    todo = set(range(len(v)))
    out = []
    while todo:
        seen = {min(todo)}
        todo -= seen
        queue = list(seen)
        while queue:
            i = queue.pop()
            for j in list(todo):
                if atlas(v[i], v[j]):
                    todo.remove(j)
                    seen.add(j)
                    queue.append(j)
        out.append(tuple(sorted(seen)))
    return out


def mixed(X, Y, Z):
    """Minimal outside coefficient and exact box; inherited height is gated."""
    d = gcd(X, Y)
    p, q = sorted((X//d, Y//d))
    check(q <= Q, 'internal pair height')
    delta = gcd(d, Z)
    c, x = d//delta, Z//delta
    if c > Q:
        return None
    z = box(Q, p, q, x)
    if z is None:
        return None
    r, s = z
    check(c*Z == d*(p*r+q*s) and max(c, abs(r), abs(s)) <= Q,
          'literal mixed witness')
    return c, r, s, p, q, x


def circle(x):
    x %= 1
    return min(x, 1-x)


def actual_row(name, V, U, t, g=1, expect_crossings=0):
    V, U = tuple(sorted(V)), tuple(sorted(U))
    a, b = len(V), len(U)
    row = tuple(t*v for v in V)+tuple(g*u for u in U)
    check(a+b == 13 and a <= b, 'split scope')
    check(gcd(*V) == gcd(*U) == gcd(t, g) == gcd(*row) == 1, 'primitive row')
    check(len(set(row)) == 13 and min(row) > 0, 'positive distinct row')
    check(sum(row) <= Q**2, 'physical box')
    check(components(row) == [tuple(range(a)), tuple(range(a, 13))], 'actual components')
    for W in (V, U):
        for x, y in combinations(W, 2):
            check(max(x, y)//gcd(x, y) <= Q, 'all internal heights')
    crossings = []
    triples = 0
    for A, B in ((row[:a], row[a:]), (row[a:], row[:a])):
        for x, y in combinations(A, 2):
            for z in B:
                triples += 1
                hit = mixed(x, y, z)
                if hit:
                    crossings.append((x, y, z, hit))
    check(triples == 11*a*b//2, 'complete two orientations')
    check(len(crossings) == expect_crossings, 'expected crossing count')
    subsetmax = {k: max(gcd(*s) for s in combinations(row, k)) for k in G}
    check(all(subsetmax[k] <= G[k] for k in G), 'genuine hereditary caps')
    D, A, B = min((gcd(x, y), x, y) for x, y in combinations(V, 2))
    selected = 7*(b+1)*D*max(U) <= a*Q and 7*(b+1)*B >= a*G[b]
    print('ACTUAL', name, 't', t, 'g', g, 'split', (a, b), 'U', U,
          'sum', sum(row), 'mixed', len(crossings), 'caps', tuple(subsetmax.items()),
          'pair_gcd', D, 'new_gate', selected)
    if crossings:
        print('FIRST_CROSSING', name, crossings[0])
    else:
        # Separate domination certificate, independent of box membership.
        check(t*D > 2*Q*g*max(U), 'direct no-crossing lower bound')
        check(t*min(V) > 2*Q*g*max(U), 'single-label no-crossing lower bound')
        check(t % 10 == 5, 'literal phase normalization')
        phase = F(1, 5)+F(1, 2*t)
        check(all(circle(v*F(1, 2)) == F(1, 2) for v in V), 'smaller phase')
        check(all(circle(u*F(1, 5)) >= F(1, 5) for u in U), 'larger phase')
        clearance = min(circle(x*phase) for x in row)
        check(clearance > F(1, 14), 'literal strictly safe physical phase')
        print('PHASE', name, phase, 'clearance', clearance)
    return selected


def choose_t(D, L):
    t = (2*Q*L//D//10+1)*10+5
    while t % 3 == 0:
        t += 10
    return t


def main():
    print('Q', Q, 'H6', H)
    # Independent multiplicative construction of the exact inert-sum atlas.
    sums = {1}
    inert = [p for p in range(2, 357) if p % 3 == 2 and
             all(p % d for d in range(2, isqrt(p)+1))]
    for p in inert:
        sums = {x*p**e for x in sums for e in range(3) if x*p**e <= 356}
    pairs_atlas = 0
    for p in range(1, 178):
        for q in range(p+1, 357-p):
            if gcd(p, q) == 1:
                expected = p+q in sums
                check(atlas(p, q) == expected, 'strict multiplicative atlas')
                pairs_atlas += expected
    check(pairs_atlas == 5855, 'exact atlas count')
    check(not atlas(1, 2) and not atlas(1, 6) and not atlas(1, 7),
          'ramified prime, split prime, and inert exponent-three hostiles')
    print('STRICT_ATLAS', pairs_atlas)
    # Complete small signed boxes, including every point beyond the first gap.
    pairs = points = 0
    for B in range(2, 17):
        for p in range(1, B):
            for q in range(p+1, B+1):
                if gcd(p, q) != 1:
                    continue
                pairs += 1
                literal = {p*r+q*s for r in range(-B, B+1) for s in range(-B, B+1)}
                R = B*(p+q)-(p-1)*(q-1)
                check(B*q < R and 2*R > B*(p+q), 'central interval bounds')
                for x in range(-B*(p+q)-1, B*(p+q)+2):
                    points += 1
                    w = box(B, p, q, x)
                    check((w is not None) == (x in literal), 'complete box')
                    if abs(x) <= R:
                        check(w is not None, 'central coverage')
                    if abs(x) == R+1:
                        check(w is None, 'first gap')
    print('BOX_UNIVERSE', pairs, 'pairs', points, 'points')

    # All small physical scalings: literal cleared coefficient is never erased.
    physical = 0
    for p, q in ((1, 2), (2, 3), (3, 5), (7, 8)):
        B = 11
        R = B*(p+q)-(p-1)*(q-1)
        for t, D, g, u in product(range(1, 8), repeat=4):
            if gcd(t, g) != 1:
                continue
            delta = gcd(t*D, g*u)
            c, x = t*D//delta, g*u//delta
            if c <= B and x <= R:
                r, s = box(B, p, q, x)
                check(c*(g*u) == r*(t*D*p)+s*(t*D*q), 'physical normalization')
                check(max(c, abs(r), abs(s)) <= B, 'physical coefficient height')
                physical += 1
    print('PHYSICAL_TOY_CROSSINGS', physical)

    # Coprime grid identity and nearest-point construction, with wraparound.
    grid = 0
    for t in range(1, 22):
        for g in range(1, 22):
            if gcd(t, g) != 1:
                continue
            for eta, zeta in product((F(1, 3), F(1, 2), F(4, 5)), repeat=2):
                points = [(g*(eta+j)/t) % 1 for j in range(t)]
                check(len(set(points)) == t, 'complete grid')
                check(min(circle(x-zeta) for x in points) <= F(1, 2*t), 'nearest-grid bound')
                grid += 1
    print('EXACT_GRID_CONTROLS', grid)

    table = []
    for a in range(2, 7):
        b = 13-a
        bound = Q*a//(7*(b+1)*355**(a-2))
        threshold = ceildiv(a*G[b], 7*(b+1))
        table.append((a, b, bound, threshold))
    check(table == [(2, 11, 13520696477, 1), (3, 10, 62323312, 1),
                    (4, 9, 257485, 1), (5, 8, 1007, 3), (6, 7, 3, 10)], 'automatic table')
    print('AUTOMATIC_NO_UNIT', table)

    # Star gcd identity tested on all coprime numerator choices for one
    # five-denominator bank. This is a finite algebra control, not an atlas census.
    den = (2, 3, 5, 7, 11)
    P = prod(den)
    stars = 0
    for num in product(range(1, 6), repeat=5):
        if any(gcd(p, q) != 1 for p, q in zip(num, den)):
            continue
        leaves = [p*(P//q) for p, q in zip(num, den)]
        check(gcd(P, *leaves) == 1, 'primitive star')
        for i, j in combinations(range(5), 2):
            check(gcd(leaves[i], leaves[j]) == P//(den[i]*den[j])*gcd(num[i], num[j]), 'star pair gcd')
        stars += 1
    print('STAR_IDENTITY_CONTROLS', stars, 'uniform_L', 3*Q//(28*355**3))

    den = (179, 181, 183, 185, 187)
    P = prod(den)
    V = tuple(sorted([P]+[(356-q)*(P//q) for q in den]))
    D, A, B = min((gcd(x, y), x, y) for x, y in combinations(V, 2))
    check(gcd(*V) == 1 and len(set(V)) == 6, 'actual star shape')
    check(all(atlas(P, (356-q)*(P//q)) for q in den), 'actual star atlas')
    check((D, A, B) == (5929017, 185370716505, 189592176609), 'selected star pair')
    check(min(V) > H and H//D == 10261, 'high-minimum improvement')
    print('SIX_STAR', V, 'pair', (D, A, B), 'reduced', (A//D, B//D), 'Lcut', H//D)
    for name, U in (
        ('nonunit18', (2, 3, 4, 6, 9, 12, 18)),
        ('nonunit9234', (2, 3, 4, 6, 18, 486, 9234)),
    ):
        check(1 not in U, 'no larger unit')
        t = choose_t(D, max(U))
        check(actual_row(name, V, U, t), 'positive theorem gate')
        check(t > 31950, 'whole-six scale-cap hostile')

    # Genuine actual-graph row satisfying all subset caps but lacking entry.
    check(actual_row('graph_not_entry', V, (2, 3, 4, 6, 9, 12, 18), 1,
                     expect_crossings=231), 'same selected gate before entry')

    # Connected primitive cofactor clique: pair gcds stay large.  A safe actual
    # equality row remains outside this sufficient closure condition.
    primes = (127, 139, 151, 157, 163, 193)
    P = prod(primes)
    V = tuple(sorted(P//p for p in primes))
    D = min(gcd(x, y) for x, y in combinations(V, 2))
    check(D == 418499671 and min(V) > H, 'residual shape')
    U = (2, 3, 4, 6, 18, 486, 9234)
    check(not actual_row('safe_gate_failure', V, U, choose_t(D, max(U))), 'gate failure survives')
    print('RESIDUAL_SIX_SHAPE', V, 'pair_min', D, 'Lcut', H//D)

    # Arithmetic-only boundaries, not alleged actual decoder entries.
    check(box(2, 1, 6, 3) is None and box(2, 1, 6, 6) is not None, 'height hostile')
    L = H
    check(28*L <= 3*Q and 90*L > 3*Q, 'small B loses automatic target containment')
    points = {(2*F(j, 4)) % 1 for j in range(4)}
    check(points == {F(0), F(1, 2)} and not any(F(1, 10) <= x <= F(2, 5) for x in points), 'noncoprime grid hostile')
    print('BOUNDARIES', 'no whole-six M6 cap; no target containment without B; no grid without coprimality')
    print('PASS', gates, 'always-active exact gates')


if __name__ == '__main__':
    main()
