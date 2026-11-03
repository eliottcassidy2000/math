#!/usr/bin/env python3
"""Finite-exact F5/F125 Frobenius-cube bridge; Python standard library only."""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import isqrt
from pathlib import Path
import json

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def polyadd(a, b):
    c = [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
         for i in range(max(len(a), len(b)))]
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return tuple(c)


def polymul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return polyadd(tuple(c), (0,))


def mmul(a, b):
    return [[polyadd(polymul(a[i][0], b[0][j]), polymul(a[i][1], b[1][j]))
             for j in range(2)] for i in range(2)]


F = list(product(range(5), repeat=3))
ZERO = (0, 0, 0)
ONE = (1, 0, 0)


def add(a, b):
    return tuple((a[i] + b[i]) % 5 for i in range(3))


def mul(a, b):
    c = [0] * 5
    for i in range(3):
        for j in range(3):
            c[i+j] += a[i] * b[j]
    # z^3 = -z-1; descending reduction also handles the z^4 term.
    for k in (4, 3):
        c[k-3] -= c[k]
        c[k-2] -= c[k]
    return tuple(x % 5 for x in c[:3])


def power(a, n):
    b = ONE
    while n:
        if n % 2:
            b = mul(b, a)
        a = mul(a, a)
        n //= 2
    return b


def canonical_curve(a, b):
    return min((pow(u, 4, 5)*a % 5, pow(u, 6, 5)*b % 5) for u in range(1, 5))


def count_f25(a, b):
    # F5[w]/(w^2-2): 2 is not a square in F5.
    f25 = list(product(range(5), repeat=2))
    def m(x,y):
        return ((x[0]*y[0]+2*x[1]*y[1]) % 5, (x[0]*y[1]+x[1]*y[0]) % 5)
    square_counts = Counter(m(y,y) for y in f25)
    count = 1
    for x in f25:
        x3 = m(m(x,x),x)
        rhs = ((x3[0]+a*x[0]+b) % 5, (x3[1]+a*x[1]) % 5)
        count += square_counts[rhs]
    return count


def lines(grid):
    return grid + [[grid[i][j] for i in range(3)] for j in range(3)] + [
        [grid[i][i] for i in range(3)], [grid[i][2-i] for i in range(3)]]


def magic_arrangements(values):
    """Complete search via forced center and opposite pairs; inputs distinct."""
    require(len(values) == len(set(values)) == 9, 'magic search requires nine distinct values')
    total = sum(values)
    if total % 9:
        return []
    center = total // 9
    if center not in values:
        return []
    centered = set(v-center for v in values)
    out = []
    for u, v in permutations(centered - {0}, 2):
        g = [[u, -u-v, v], [v-u, 0, u-v], [-v, u+v, -u]]
        if sorted(x for row in g for x in row) == sorted(centered):
            out.append([[x+center for x in row] for row in g])
    return out


def tournaments():
    """Independent complete census on five vertices, only 1024 labeled masks."""
    pairs = list(combinations(range(5), 2))
    perms = list(permutations(range(5)))

    def matrix(mask):
        a = [[0]*5 for _ in range(5)]
        for k, (i, j) in enumerate(pairs):
            a[i][j] = (mask >> k) & 1
            a[j][i] = 1-a[i][j]
        return a

    def canon(a):
        return min(sum(a[p[i]][p[j]] << k for k, (i, j) in enumerate(pairs)) for p in perms)

    reps = sorted({canon(matrix(mask)) for mask in range(1024)})
    spectra = Counter()
    records = []
    for mask in reps:
        a = matrix(mask)
        # Leibniz determinant of xI-A, coefficients in increasing order.
        cp = (0,)
        for p in perms:
            sign = (-1)**sum(p[i] > p[j] for i in range(5) for j in range(i+1, 5))
            term = (sign,)
            for i in range(5):
                term = polymul(term, (0, 1) if p[i] == i else (-a[i][p[i]],))
            cp = polyadd(cp, term)
        spectra[cp] += 1
        dual = canon([list(row) for row in zip(*a)])
        records.append({'mask': mask, 'characteristic_desc': list(reversed(cp)),
                        'dual_mask': dual, 'self_dual': dual == mask})
    return records, sorted(spectra.values(), reverse=True)


def main():
    roots = [(x**3+x+1) % 5 for x in range(5)]
    require(roots == [1, 3, 1, 1, 4], 'irreducible modulus root evaluations')
    require(0 not in roots, 'degree-three modulus is irreducible')
    require(len(F) == 125, 'field universe')
    require(mul((0, 0, 1), (0, 1, 0)) == (4, 4, 0), 'z^3 reduction')
    for a in F:
        require(mul(a, ONE) == a and add(a, ZERO) == a, 'field identities')
        require(power(a, 125) == a, 'F125 Frobenius')
        if a != ZERO:
            require(mul(a, power(a, 123)) == ONE, 'nonzero field inverse')

    # Two square tests: a literal histogram and Euler's criterion.
    squares = Counter(mul(y, y) for y in F)
    for a in F:
        expected = 1 if a == ZERO else (2 if power(a, 62) == ONE else 0)
        require(squares[a] == expected, 'square histogram / Euler cross-check')

    models = []
    singular = []
    classes = defaultdict(list)
    expected_n5 = {(0,1):6, (0,2):6, (0,3):6, (0,4):6,
                   (1,0):4, (1,1):9, (1,2):4, (1,3):4, (1,4):9,
                   (2,0):2, (2,1):7, (2,4):7, (3,0):10,
                   (3,2):5, (3,3):5, (4,0):8, (4,1):8,
                   (4,2):3, (4,3):3, (4,4):8}
    for a, b in product(range(5), repeat=2):
        if (4*a**3 + 27*b*b) % 5 == 0:
            singular.append([a, b])
            continue
        affine5 = [(x,y) for x,y in product(range(5),repeat=2)
                   if (y*y-x**3-a*x-b) % 5 == 0]
        n5 = 1+len(affine5)
        require(n5 == expected_n5[(a,b)], 'direct F5 census')
        t = 6-n5
        rhs = [add(add(mul(mul(x,x),x), mul((a,0,0),x)), (b,0,0)) for x in F]
        n125 = 1+sum(squares[v] for v in rhs)
        # A second direct count exhausts all 125^2 coordinate pairs.
        n125_pairs = 1+sum(mul(y,y) == v for v in rhs for y in F)
        require(n125 == n125_pairs, 'direct extension point counts agree')
        require(n125 == 126-t**3+15*t, 'Frobenius-cube formula')
        n25 = count_f25(a,b)
        require(n25 == 36-t*t, 'direct F25 quadratic-extension control')
        rep = canonical_curve(a,b)
        classes[rep].append((a,b))
        models.append({'a':a, 'b':b, 'n5':n5, 't':t, 'n25':n25, 'n125':n125,
                       'class_representative':list(rep), 'affine_points_f5':affine5})
    require(len(models) == 20 and len(singular) == 5, 'model universe split')
    require(set(m['n5'] for m in models) == set(range(2,11)), 'cycle-label support')
    require(len(classes) == 12, 'isomorphism classes')
    by_pair = {(m['a'],m['b']):m for m in models}
    class_records = []
    for (a,b), members in sorted(classes.items()):
        t = by_pair[(a,b)]['t']
        require({by_pair[p]['t'] for p in members} == {t}, 'scaling preserves counts')
        twisted = canonical_curve(4*a % 5, 3*b % 5)
        require(twisted != (a,b), 'twist acts freely on classes')
        require(canonical_curve(4*twisted[0] % 5, 3*twisted[1] % 5) == (a,b), 'twist involution')
        require(by_pair[twisted]['t'] == -t, 'twist negates trace')
        class_records.append({'representative':[a,b], 'models':members,
                              't':t, 'n5':6-t, 'n125':126-t**3+15*t,
                              'twist_representative':list(twisted)})
    fibers = Counter(row['t'] for row in class_records)
    require(fibers == Counter({-4:1,-3:1,-2:2,-1:1,0:2,1:1,2:2,3:1,4:1}), 'trace fibers')

    # Symbolic companion-matrix proof, not merely nine numerical tests.
    companion = [[(0,),(-5,)],[(1,),(0,1)]]
    m2 = mmul(companion, companion)
    m3 = mmul(m2, companion)
    for i,j in product(range(2), repeat=2):
        require(polyadd(m2[i][j], polyadd(polymul((0,-1), companion[i][j]),
                                        (5,) if i == j else (0,))) == (0,), 'Cayley-Hamilton')
    trace3 = polyadd(m3[0][0],m3[1][1])
    require(trace3 == (0,-15,0,1), 'symbolic Frobenius cube trace')
    hecke0, hecke1 = (1,), (0,1)
    hecke2 = polyadd(polymul((0,1), hecke1), polymul((-5,), hecke0))
    hecke3 = polyadd(polymul((0,1), hecke2), polymul((-5,), hecke1))
    require(hecke3 == (0,-10,0,1), 'Hecke / symmetric-cube coefficient')
    require(polyadd(hecke3, tuple(-x for x in trace3)) == (0,5), 'Hecke differs by 5t')

    counts = [126-t**3+15*t for t in range(-4,5)]
    witness = [[108,148,122],[140,126,112],[130,104,144]]
    require(sorted(x for row in witness for x in row) == sorted(counts), 'magic witness entries')
    require([sum(line) for line in lines(witness)] == [378]*8, 'all eight magic lines')
    arrangements = magic_arrangements(counts)
    require(len(arrangements) == 8 and witness in arrangements, 'all magic arrangements')
    hecke_values = [t**3-10*t for t in range(-4,5)]
    require(magic_arrangements(hecke_values) == [], 'Hecke cubic is not magic-arrangeable')
    require(magic_arrangements([t**3 for t in range(-4,5)]) == [], 'literal cubes not magic-arrangeable')
    require(len(magic_arrangements(list(range(2,11)))) == 8, 'affine positive control')
    square_entries = [n for n in counts if isqrt(n)**2 == n]
    require(square_entries == [144], 'not a square-valued magic square')
    counts25 = [36-t*t for t in range(-4,5)]
    require(sum(counts25) == 264 and sum(counts25) % 9 != 0,
            'F25 indexed multiset cannot be integer magic: nonintegral center')
    cycle = [4,8,3,7,2,6,10,5,9,4]
    transformed_cycle = [126-(6-n)**3+15*(6-n) for n in cycle]
    require(len(set(transformed_cycle[:-1])) == 9 and transformed_cycle[0] == transformed_cycle[-1], 'cycle support transport')

    trecords, tfibers = tournaments()
    inherited = json.loads((Path(__file__).parent/'sources/tournament_classification.json').read_text())
    wanted = [{k:r[k] for k in ('mask','characteristic_desc','dual_mask','self_dual')}
              for r in inherited['classes']]
    require(trecords == wanted, 'independent tournament census agrees with source')
    require(tfibers == [3,2,1,1,1,1,1,1,1], 'tournament spectrum fibers')
    selfdual = sum(r['self_dual'] for r in trecords)
    require(selfdual == 8 and (len(trecords)+selfdual)//2 == 10, 'converse involution')
    require(tfibers != sorted(fibers.values(),reverse=True), 'no quotient-preserving bijection')
    require(selfdual != 0, 'no twist/converse-equivariant bijection')

    result = {'status':'FINITE-EXACT; proposed precise lambda-cubed route, not confirmed user referent',
              'checks_passed':CHECKS, 'universe':{'short_models':25,'nonsingular':20,'singular':singular,
              'field_modulus':'z^3+z+1 over F5','modulus_root_values':roots},
              'models':models, 'isomorphism_classes':class_records,
              'class_trace_fibers':dict(sorted(fibers.items())), 'twist_orbits':6,
              'symbolic_trace3_coefficients_ascending':trace3,
              'symbolic_hecke3_coefficients_ascending':hecke3,
              'trace_to_extension_count':[[t,126-t**3+15*t] for t in range(-4,5)],
              'f25_control':{'indexed_values':counts25,'sum':sum(counts25),
                             'nonintegral_required_center':'88/3'},
              'magic_square':witness,'line_sums':[378]*8,'magic_arrangement_count':len(arrangements),
              'square_entries':square_entries,'cycle':cycle,'transformed_cycle':transformed_cycle,
              'tournament_classes':trecords,'tournament_spectrum_fibers':tfibers,
              'tournament_selfdual_classes':selfdual,'tournament_converse_orbits':10,
              'hostile_controls':{'hecke_cubic_not_magic':True,'literal_cube_not_magic':True,
              'incompatible_12_to_9_fibers':True,'incompatible_involutions':True}}
    print(json.dumps(result,sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
