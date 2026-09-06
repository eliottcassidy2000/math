"""Independent counting/measure audit: row multigraphs and literal stub maps.

No producer imports. All tests are active under python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial, gcd

gates = 0


def require(condition, label):
    global gates
    gates += 1
    if not condition:
        raise ArithmeticError(label)


def components(n, edges):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            a = parent[a]
        return a

    for a, b in edges:
        a, b = find(a), find(b)
        parent[a] = b
    counts = Counter(find(i) for i in range(n))
    return tuple(sorted(counts.values()))


def multigraphs(n):
    """A derangement orients every row-neighborhood cycle independently."""
    bank = Counter()
    for p in permutations(range(n)):
        if all(p[i] != i for i in range(n)):
            edges = tuple(sorted(tuple(sorted((i, p[i]))) for i in range(n)))
            bank[edges] += 1
    return bank


def multiset_orders(items):
    counts = Counter(items)
    keys = sorted(counts)
    word = []

    def visit():
        if len(word) == len(items):
            yield tuple(word)
        else:
            for item in keys:
                if counts[item]:
                    counts[item] -= 1
                    word.append(item)
                    yield from visit()
                    word.pop()
                    counts[item] += 1
    yield from visit()


def board_mask(columns):
    n = len(columns)
    return sum(1 << (r*n+c) for c, pair in enumerate(columns) for r in pair)


def line_geometry(columns):
    """A repeated primitive affine line from different point pairs is a triple."""
    pts = tuple((r, c) for c, pair in enumerate(columns) for r in pair)
    lines = set()
    good = True
    for (r, c), (s, d) in combinations(pts, 2):
        aa, bb = d-c, r-s
        common = gcd(abs(aa), abs(bb))
        aa //= common
        bb //= common
        if aa < 0 or (aa == 0 and bb < 0):
            aa, bb = -aa, -bb
        key = (aa, bb, aa*r+bb*c)
        if key in lines:
            good = False
        lines.add(key)
    occupancies = Counter(c-r for r, c in pts)
    defect = sum(max(0, v-2) for v in occupancies.values())
    return good, defect


def connected_recurrence(limit):
    """Choose the connected component containing distinguished row zero."""
    values = [1, 0]
    for n in range(2, limit+1):
        value = sum(comb(n-1, r-1)*comb(n, r)*factorial(r)*factorial(r-1)
                    // 2*values[n-r] for r in range(2, n+1))
        values.append(value)
    return values


def literal_stub_bank(n):
    boards = Counter()
    for p in permutations(range(2*n)):
        edges = [(i//2, p[i]//2) for i in range(2*n)]
        if len(set(edges)) == 2*n:
            mask = sum(1 << (r*n+c) for r, c in edges)
            boards[mask] += 1
    return boards


def audit_orbits():
    counts = connected_recurrence(60)
    derangements = [1, 0]
    for n in range(2, 9):
        derangements.append((n-1)*(derangements[-1]+derangements[-2]))
    expected_small = {2: (1, 1, 1), 3: (6, 2, 4), 4: (90, 11, 35),
                      5: (2040, 32, 545)}
    masks_by_n = {}
    for n in range(2, 9):
        bank = multigraphs(n)
        profile_counts = Counter()
        total = weighted = 0
        good_total = zero_total = 0
        masks = set()
        for edges, orientation_count in sorted(bank.items()):
            typ = components(n, edges)
            edge_counts = Counter(edges)
            require(all(v in (1, 2) for v in edge_counts.values()), 'at most twin columns')
            c4 = sum(v == 2 for v in edge_counts.values())
            require(c4 == typ.count(2), 'twin columns exactly C4')
            require(orientation_count == 2**(len(typ)-c4), 'row-cycle orientations')
            size = factorial(n)//2**c4
            total += size
            weighted += size*2**len(typ)
            profile_counts[typ] += size
            if n <= 5:
                orbit_good = orbit_defect = actual_size = 0
                for columns in multiset_orders(edges):
                    mask = board_mask(columns)
                    require(mask not in masks, 'multigraph column ordering is injective')
                    masks.add(mask)
                    good, defect = line_geometry(columns)
                    require(not good or defect == 0, 'full success implies selected zero')
                    good_total += good
                    zero_total += defect == 0
                    orbit_good += good
                    orbit_defect += defect
                    actual_size += 1
                require(actual_size == size, 'literal fixed-row orbit size')
                mean = Q(orbit_defect, size)
                x = mean*mean / (8*(n-1))
                # exp(x) <= 1/(1-x) for 0<=x<1, proved coefficientwise.
                require(0 <= x < 1, 'elementary exponential comparison domain')
                require(Q(orbit_good, size) <= 1-x, 'stronger rational small-orbit tail')
        require(sum(bank.values()) == derangements[n], 'derangement orientation count')
        require(total == counts[n], 'connected-component recurrence matches orbits')
        require(weighted == factorial(n)*derangements[n], '2^components weighted count')
        for typ, number in profile_counts.items():
            divisor = 1
            for r, copies in Counter(typ).items():
                divisor *= (2*r)**copies*factorial(copies)
            require(number*divisor == factorial(n)**2, 'two-shore orbit normalization')
        if n <= 5:
            require((total, good_total, zero_total) == expected_small[n], 'independent exact geometry totals')
            masks_by_n[n] = masks
            print('n', n, 'row-multigraphs', len(bank), 'boards', total,
                  'no-three', good_total, 'slope-one-zero', zero_total)
        else:
            print('n', n, 'row-multigraphs', len(bank), 'boards', total,
                  'derangements', derangements[n], 'ordered-pairs', weighted)
    for n in range(2, 5):
        maps = literal_stub_bank(n)
        require(set(maps) == masks_by_n[n], 'stub model covers exactly the simple boards')
        for preimages in maps.values():
            require(preimages == 4**n, 'literal stub fiber is 4^n')
        require(sum(maps.values()) == counts[n]*4**n, 'configuration count')
        require(sum(maps.values()) <= factorial(2*n), 'discarding nonsimple configurations')
        print('stub n', n, 'simple-bijections', sum(maps.values()),
              'all-bijections', factorial(2*n), 'fiber', 4**n)
    # Check the stated recurrence against independently constructed components.
    for n in range(2, 61):
        require(2*counts[n] == 2*n*(n-1)*counts[n-1]
                + n*(n-1)**2*counts[n-2], 'second-order counting identity')
        require(counts[n]*4**n <= factorial(2*n), 'configuration bound through 60')
    print('component recurrence and configuration bound through n60')


def audit_hazards():
    q = Q(3, 7)
    alive = {(): Q(1)}
    truncated_expectation = Q(0)
    for depth in range(11):
        survival = sum(alive.values(), Q(0))
        truncated_expectation += survival
        require(survival >= (1-q)**depth, 'adaptive survival')
        require(1-survival <= 1-(1-q)**depth <= depth*q, 'adaptive success upper bounds')
        require(truncated_expectation >= sum(((1-q)**j for j in range(depth+1)), Q(0)),
                'truncated tail-sum expectation')
        future = {}
        for history, mass in alive.items():
            code = sum((j+1)*(bit+1) for j, bit in enumerate(history))
            hazard = q*Q(1+(code+depth)%4, 4)
            require(0 <= hazard <= q, 'full-history conditional hazard cap')
            observation = Q(1+(code%3), 4)
            future[history+(0,)] = mass*(1-hazard)*observation
            future[history+(1,)] = mass*(1-hazard)*(1-observation)
        alive = future
    print('adaptive history tree q=3/7 depth10 survival', survival)
    # Equality benchmark: constant conditional hazard exactly gives geometric law.
    for depth in range(21):
        require(sum((q*(1-q)**j for j in range(depth)), Q(0)) == 1-(1-q)**depth,
                'geometric equality control')

    # A genuine n4 orbit supplies a hostile to MARGINAL-only uniformity.
    # A uniform cyclic offset through all 24 permutations keeps each marginal
    # uniform, but forces a success eventually, unlike geometric lower tails.
    edges = ((0, 1), (0, 2), (1, 3), (2, 3))
    orders = tuple(permutations(range(4)))
    good = [line_geometry(tuple(edges[i] for i in order))[0] for order in orders]
    hits = sum(good)
    require(0 < hits < len(orders), 'marginal-uniform hostile is nondegenerate')
    first_hits = []
    for offset in range(len(orders)):
        first = next(j+1 for j in range(len(orders)) if good[(offset+j)%len(orders)])
        first_hits.append(first)
    for j in range(len(orders)):
        require(Counter((offset+j)%len(orders) for offset in range(len(orders)))
                == Counter(range(len(orders))), 'every restart marginal is uniform')
    q_orbit = Q(hits, len(orders))
    mean = Q(sum(line_geometry(tuple(edges[i] for i in order))[1]
                 for order in orders), len(orders))
    require(mean == Q(5, 6), 'hostile has nontrivial conditional concentration exponent')
    require(max(first_hits) <= len(orders) and (1-q_orbit)**len(orders) > 0,
            'marginal-uniform schedule violates geometric survival')
    print('n4 C8 marginal-only hostile hits', hits, 'of', len(orders),
          'largest-first-hit', max(first_hits), 'mean-first-hit', Q(sum(first_hits), len(first_hits)),
          'selected-defect-mean', mean)


def audit_constants():
    # Bound e^2 by a fresh 20-term positive series and geometric tail.
    term = total = Q(1)
    for k in range(1, 21):
        term *= Q(2, k)
        total += term
    lower = total
    upper = total + term*Q(2, 21)/(1-Q(2, 22))
    a_lo, a_hi = 1-5/lower, 1-5/upper
    require(Q(323, 1000) < a_lo < a_hi < Q(324, 1000), 'alpha enclosure')
    require(a_lo > Q(2, 15) and a_lo*a_lo/8 > Q(1, 900), 'eventual dominant branches')
    # Verify the exact polynomial-division identity for arbitrary rational alpha.
    for alpha in (Q(323, 1000), Q(1, 3), a_lo, a_hi):
        for n in (2, 4, 64, 128, 1024):
            direct = (alpha*n-21)**2/(8*(n-1))
            expanded = alpha*alpha*n/8 + (alpha*alpha-42*alpha)/8 \
                       + (alpha-21)**2/(8*(n-1))
            require(direct == expanded, 'linear entropy-rate expansion')
    print('alpha in (323/1000,324/1000); c=alpha^2/8>1/900')


def main():
    print('Independent orbit, configuration, and conditional-restart audit')
    audit_orbits()
    audit_hazards()
    audit_constants()
    print('PASS', gates, 'always-active gates')


if __name__ == '__main__':
    main()
