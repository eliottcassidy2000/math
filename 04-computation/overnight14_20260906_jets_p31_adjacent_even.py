"""Two adjacent even ideals and the all-depth p31 factor pair; exact arithmetic."""
from fractions import Fraction as Q
from itertools import combinations
from math import comb, factorial

gates = 0


def check(ok, message):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(message)


def vp(n, p=31):
    if n == 0:
        return 10**9
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def factorial_vp(n, p=31):
    answer = 0
    while n:
        n //= p
        answer += n
    return answer


def shift_scalar(m, s):
    h = m-s
    value = Q(1)
    for c in range(m, m+2*s):
        value *= Q(factorial(c), factorial(c-h))
    for j in range(h, m):
        value /= Q(factorial(j), factorial(j-h))**2
    check(value.denominator == 1, 'integral even shift scalar')
    return value.numerator


def det_bareiss(matrix):
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(len(a)-1):
        hit = next((r for r in range(k, len(a)) if a[r][k]), None)
        if hit is None:
            return 0
        if hit != k:
            a[k], a[hit] = a[hit], a[k]
            sign = -sign
        pivot = a[k][k]
        for r in range(k+1, len(a)):
            for c in range(k+1, len(a)):
                value = pivot*a[r][c]-a[r][k]*a[k][c]
                check(value % previous == 0, 'fraction-free exact division')
                a[r][c] = value//previous
        previous = pivot
    return sign*a[-1][-1]


def even_minor(m, s, a):
    h = m-s
    return [[comb(c, j)*x**(c-j) for c in range(m, m+2*s)]
            for x in (1, a) for j in range(h, m)]


def scalar_controls():
    count = 0
    for s in range(1, 7):
        for m in (s, s+1, s+3):
            scalar = shift_scalar(m, s)
            for a in (2, 3):
                actual = det_bareiss(even_minor(m, s, a))
                expected = scalar*(a*(a-1))**(s*s)
                check(actual == expected, 'universal signed even minimal determinant')
                count += 1
    c28, c30 = shift_scalar(16, 14), shift_scalar(16, 15)
    closed28 = Q(factorial(43)*factorial(42), factorial(15)**3*factorial(14)**3)
    closed30 = Q(factorial(45), factorial(15)**3)
    check(c28 == closed28 and c30 == closed30, 'two factorial scalar expressions')
    check(vp(c28) == 2 and vp(c30) == 1, 'literal scalar contents')
    check(factorial_vp(43)+factorial_vp(42)-3*factorial_vp(15)-3*factorial_vp(14) == 2,
          'Legendre scalar28 content')
    check(factorial_vp(45)-3*factorial_vp(15) == 1, 'Legendre scalar30 content')
    for s, scalar in ((14, c28), (15, c30)):
        for a in (2, 3, 4, -29):
            actual = det_bareiss(even_minor(16, s, a))
            check(actual == scalar*(a*(a-1))**(s*s), 'literal large normalized even determinant')
            check(vp(actual) == vp(scalar), 'unit monomial causes no valuation change')
            count += 1
    print('signed literal even determinants', count)
    print('C28', c28, 'v31', vp(c28))
    print('C30', c30, 'v31', vp(c30))


def complete_bands(r):
    allrows = [(x, j) for x in (1, 2) for j in range(16)]
    allcols = list(range(16, 48))
    s, omitted_count = r//2, 32-r
    rowmax = 2*sum(range(16-s, 16))
    colmin = sum(range(16, 16+r))
    rowsets, colsets = {0: [], 1: []}, {0: [], 1: []}
    index_count = 0
    for omitted in combinations(range(32), omitted_count):
        index_count += 1
        rowloss = rowmax-240+sum(allrows[i][1] for i in omitted)
        colloss = 1008-sum(allcols[i] for i in omitted)-colmin
        check(rowloss >= 0 and colloss >= 0, 'all row/column losses nonnegative')
        if rowloss <= 1:
            rowsets[rowloss].append([x for i, x in enumerate(allrows) if i not in omitted])
        if colloss <= 1:
            colsets[colloss].append([c for i, c in enumerate(allcols) if i not in omitted])
    check(index_count == comb(32, r), 'entire index-set universe')
    check(len(rowsets[0]) == len(colsets[0]) == 1, 'unique minimum even minor')
    check(len(rowsets[1]) == 4 and len(colsets[1]) == 1, 'complete next row and column bands')
    check(colmin-rowmax == 3*s*s, 'even minimum weight')
    next_count = 0
    for rowloss in (0, 1):
        for rows in rowsets[rowloss]:
            for cols in colsets[1-rowloss]:
                next_count += 1
                check(sum(cols)-sum(j for _, j in rows) == 3*s*s+1, 'next weight')
                if r == 28:
                    check(all(j >= 1 for _, j in rows), 'rank28 next band retains no zero jets')
                    check(31 in cols and all(comb(31, j) % 31 == 0 for _, j in rows),
                          'identically zero column31 in every rank28 next minor')
    check(next_count == 5, 'five next-band even minors')
    if r == 30:
        # The next coefficient need not be divisible by31; the weight gap suffices.
        rows = next(rows for rows in rowsets[1] if any(j == 0 for _, j in rows))
        check(any(comb(31, j) % 31 for _, j in rows), 'do not import rank28 zero-column argument')
    print('rank', r, 'index-pair universe', index_count**2,
          'minimum weight', colmin-rowmax, 'minimal minors1 next minors5')


def smith_global_pivots(a, e):
    """Literal original matrix, fixed precision, least-valuation DVR pivots."""
    precision = 48*e+3
    modulus = 31**precision
    nodes = (0, 31**e, 31**e*a)
    matrix = [[comb(c, j)*x**(c-j) % modulus if c >= j else 0
               for c in range(48)] for x in nodes for j in range(16)]
    exponents = []
    while matrix:
        n = len(matrix)
        d, r, c = min((vp(matrix[r][c]), r, c) for r in range(n) for c in range(n))
        check(d < precision, 'visible Smith pivot')
        matrix[0], matrix[r] = matrix[r], matrix[0]
        for row in matrix:
            row[0], row[c] = row[c], row[0]
        scale = 31**d
        inverse = pow(matrix[0][0]//scale, -1, modulus)
        pivot = [(v*inverse) % modulus for v in matrix[0]]
        check(pivot[0] == scale, 'unit-only normalization')
        reduced = []
        for row in matrix[1:]:
            check(row[0] % scale == 0, 'integral elimination multiplier')
            multiplier = row[0]//scale
            reduced.append([(row[j]-multiplier*pivot[j]) % modulus for j in range(1, n)])
        exponents.append(d)
        matrix = reduced
    check(exponents == sorted(exponents) and sum(exponents) == 768*e, 'Smith ordering and determinant')
    return tuple(exponents)


def full_controls_and_kernel_window():
    roots = {3, 11, 15, 17, 21, 29}
    bank = [(3, e) for e in (0, 1, 2, 5)] + [(4, e) for e in (0, 1, 2, 5)]
    bank += [(879, 3), (-20, 4)]
    for a, e in bank:
        parts = smith_global_pivots(a, e)
        kappa = int(a % 31 in roots)
        if e == 0:
            check(parts == (0,)*48, 'depth-zero boundary')
        else:
            check(sum(parts[:44]) == 588*e+2, 'literal full D44 / residual E28')
            check(sum(parts[:46]) == 675*e+1, 'literal full D46 / residual E30')
            check(sum(parts[:45]) == 631*e+1+kappa, 'inherited odd ideal control')
            check(parts[44:46] == (43*e-1+kappa, 44*e-kappa), 'actual all-depth factor-pair formula')
        print('a', a, 'e', e, 'E28,E29,E30',
              tuple(sum(parts[:16+r]) for r in (28, 29, 30)), 'pair', parts[44:46])
    for e in range(1, 21):
        check(589*e+1 >= 588*e+2 and 590*e >= 588*e+2,
              'rank28 next/higher lower bounds')
        check(676*e >= 675*e+1, 'rank30 higher lower bound')
        ordinary, exceptional = (43*e-1, 44*e), (43*e, 44*e-1)
        check(sum(ordinary) == sum(exceptional) == 87*e-1, 'determinant-blind pair sum')
        check((exceptional[0] == exceptional[1]) == (e == 1), 'exact equality boundary')
        for b in range(0, 45*e+2):
            difference = sum(min(b, v) for v in exceptional)-sum(min(b, v) for v in ordinary)
            check(difference == int(43*e <= b <= 44*e-1), 'exact two-factor kernel window')
    print('full original matrices', len(bank), 'pair kernel windows through e20')


def main():
    import sys
    sys.stdout.reconfigure(newline='\n')
    print('p31 adjacent even ideals E28=588e+2, E30=675e+1 and residual factors29/30')
    scalar_controls()
    complete_bands(28)
    complete_bands(30)
    full_controls_and_kernel_window()
    print('PASS', gates, 'always-active gates')


if __name__ == '__main__':
    main()
