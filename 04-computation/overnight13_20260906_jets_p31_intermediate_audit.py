"""Independent universal cofactor/shift and p31 intermediate-ideal audit."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import comb, factorial, gcd

gates = 0


def require(ok, label):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(label)


def polynomial(s):
    """Expand the two-binomial Hermite coefficient by ordinary convolution."""
    out = [0]*(s+1)
    for j in range(s+1):
        scalar = comb(s+j, j)*comb(2*s-j, s-j)
        for k in range(j+1):
            out[s-j+k] += scalar*comb(j, k)*(-1)**(j-k)
    content = reduce(gcd, out)
    return content, tuple(x//content for x in out)


def evaluate(coefficients, x):
    return sum(c*x**j for j, c in enumerate(coefficients))


def valuation(x, p=31):
    if x == 0:
        return 10**9
    if isinstance(x, F):
        return valuation(x.numerator, p)-valuation(x.denominator, p)
    result = 0
    while x % p == 0:
        result += 1
        x //= p
    return result


def determinant_q(matrix):
    """Rational Gaussian determinant, unlike the producer's Bareiss path."""
    a = [[F(x) for x in row] for row in matrix]
    ans = F(1)
    for k in range(len(a)):
        choices = [r for r in range(k, len(a)) if a[r][k]]
        if not choices:
            return F(0)
        r = min(choices, key=lambda r: abs(a[r][k].numerator))
        if r != k:
            a[k], a[r] = a[r], a[k]
            ans = -ans
        pivot = a[k][k]
        ans *= pivot
        for r in range(k+1, len(a)):
            ratio = a[r][k]/pivot
            for c in range(k+1, len(a)):
                a[r][c] -= ratio*a[k][c]
    return ans


def rows_for(m, s, which):
    d = m-s-1
    return [(x, j) for x in (1, 2) for j in range(d, m)
            if j > d or x == which]


def literal_minor(m, s, a, which):
    return [[comb(c, j)*(1 if x == 1 else a)**(c-j)
             for c in range(m, m+2*s+1)] for x, j in rows_for(m, s, which)]


def derivative_ratio(m, s):
    d = m-s-1
    ratio = F(1)
    for c in range(m, m+2*s+1):
        ratio *= F(factorial(c), factorial(c-d))
    for _, j in rows_for(m, s, 1):
        ratio /= F(factorial(j), factorial(j-d))
    return ratio


def inverse_coefficient(s, x, a):
    """Multiply the two reciprocal Taylor series directly over Q."""
    others = [b for b in (0, 1, a) if b != x]
    series = [F(1)]+[F(0)]*s
    for b in others:
        factor = [F((-1)**j*comb(s+j, j), (x-b)**(s+1+j))
                  for j in range(s+1)]
        series = [sum((series[k]*factor[j-k] for k in range(j+1)), F(0))
                  for j in range(s+1)]
    return series[s]


def universal_controls():
    determinants = 0
    for s in range(7):
        content, q = polynomial(s)
        for a in (2, 3, 5):
            t = content*evaluate(q, a)
            u = content*evaluate(tuple(reversed(q)), a)
            require(inverse_coefficient(s, a, a)
                    == F((-1)**s*t, (a*(a-1))**(2*s+1)), 'Hermite reciprocal at a')
            require(inverse_coefficient(s, 1, a)
                    == F((-1)**(s+1)*u, (a-1)**(2*s+1)), 'Hermite reciprocal at one')
        for m in (s+1, s+2, s+4):
            k = content*derivative_ratio(m, s)
            require(k.denominator == 1, 'integral content after rational shift')
            for a in (2, 3):
                for which in (1, 2):
                    expected = k*(a-1)**(s*s)
                    if which == 1:
                        expected *= a**(s*s)*evaluate(q, a)
                    else:
                        expected *= (-1)**s*a**((s+1)**2)*evaluate(tuple(reversed(q)), a)
                    actual = determinant_q(literal_minor(m, s, a, which))
                    require(actual == expected, 'signed universal minimal determinant')
                    determinants += 1
    require(derivative_ratio(4, 2) == F(560, 3), 'rational-shift hostile')
    print('universal signed rational determinants', determinants, 's0..6; m=s+1,s+2,s+4')
    return determinants


def packet_controls():
    content, q = polynomial(14)
    q2 = tuple(reversed(q))
    require(content == 19380 and gcd(content, 31) == 1, 'T14 primitive content')
    k = content*derivative_ratio(16, 14)
    require(k == 23038504627568008043520 and valuation(k) == 1, 'shifted packet content')
    companion = []
    for i in range(17):
        numerator = (q2[i] if i < 15 else 0)-(q[i-2] if 2 <= i < 17 else 0)
        require(numerator % 31 == 0, 'global divided companion identity')
        companion.append(numerator//31)
    roots = []
    for a in range(2, 31):
        pair = (evaluate(q, a) % 31, evaluate(q2, a) % 31)
        require(pair[1] == a*a*pair[0] % 31, 'packet proportionality')
        if pair == (0, 0):
            roots.append(a)
            require(evaluate(companion, a) % 31 != 0, 'common valuation cannot reach two')
    require(roots == [3, 11, 15, 17, 21, 29], 'complete field root bank')
    units = [evaluate(companion, a) % 31 for a in roots]
    require(units == [11, 3, 17, 20, 14, 28], 'divided-unit values')
    # Rational units, negative values, and large higher digits test the p-adic scope.
    bank = [F(a+31*d, 1+31*k) for a in range(2, 31)
            for d, k in ((0, 0), (1, 0), (-7, 2), (31**3+5, 31**2+2))]
    for a in bank:
        residue = a.numerator*pow(a.denominator, -1, 31) % 31
        require(min(valuation(evaluate(q, a)), valuation(evaluate(q2, a)))
                == int(residue in roots), 'rational and higher-digit packet loss')
    # The six exceptional classes form a single anharmonic orbit of 3.
    orbit = {3, (1-3) % 31, pow(3, -1, 31), (1-pow(3, -1, 31)) % 31,
             3*pow(2, -1, 31) % 31, pow(-2, -1, 31)}
    require(orbit == set(roots), 'node-relabeling orbit')
    for a in roots:
        require(sum(comb(15, j)**2*a**j for j in range(16)) % 31 != 0
                and a not in (2, 16, 30), 'all six classes are ordinary and non-AP')
    require((valuation(evaluate(q, 879)), valuation(evaluate(q2, 879))) == (3, 1),
            'deep individual cancellation does not change the common ideal')
    for a in (2, 3, 4):
        for which in (1, 2):
            expected = k*a**(196 if which == 1 else 225)*(a-1)**196 \
                       * evaluate(q if which == 1 else q2, a)
            require(determinant_q(literal_minor(16, 14, a, which)) == expected,
                    'literal p31 29-minor by rational elimination')
    print('T14 content', content, 'K', k)
    print('exceptional residues', roots, 'divided units', units)
    print('rational/higher-digit packet controls', len(bank), 'literal large minors', 6)
    return set(roots)


def weight_controls():
    row_indices = [(x, j) for x in (1, 2) for j in range(16)]
    column_indices = list(range(16, 48))
    row_bands, column_bands = {}, {}
    for omit in combinations(range(32), 3):
        rows = [r for i, r in enumerate(row_indices) if i not in omit]
        columns = [c for i, c in enumerate(column_indices) if i not in omit]
        row_loss = 239-sum(j for _, j in rows)
        column_loss = sum(columns)-870
        require(row_loss >= 0 and column_loss >= 0, 'global least weight')
        row_bands.setdefault(row_loss, []).append(rows)
        column_bands.setdefault(column_loss, []).append(columns)
    require(len(row_bands[0]) == 2 and len(row_bands[1]) == 4, 'two lowest row bands complete')
    require(len(column_bands[0]) == 1 and len(column_bands[1]) == 1, 'two lowest column bands complete')
    next_count = 0
    for rloss in (0, 1):
        for rows in row_bands[rloss]:
            for columns in column_bands[1-rloss]:
                next_count += 1
                require(sum(columns)-sum(j for _, j in rows) == 632, 'next weight exactly')
                if min(j for _, j in rows) >= 1:
                    require(31 in columns and all(comb(31, j) % 31 == 0 for _, j in rows),
                            'single identically zero mod31 column')
                else:
                    require(sum(j == 0 for _, j in rows) == 1
                            and all(j == 0 or j >= 2 for _, j in rows), 'exceptional row type')
                    require(31 in columns and 32 in columns
                            and all(comb(c, j) % 31 == 0 for c in (31, 32)
                                    for _, j in rows if j >= 2), 'two columns with one-row support')
    require(next_count == 6, 'six next-band minors')
    all_minors = sum(map(len, row_bands.values()))*sum(map(len, column_bands.values()))
    require(all_minors == comb(32, 29)**2 == 24601600, 'full index universe covered by bands')
    for e in (1, 2, 3, 100):
        require(632*e+1 >= 631*e+2 and 633*e >= 631*e+2, 'higher-band all-depth lower bounds')
    print('all minor index pairs', all_minors, 'minimal', 2, 'next', next_count)


def smith_global_pivots(a, e):
    """Fixed-precision DVR diagonalization using globally least-valued pivots.

    No common-layer division and no import of the producer's peeling code.
    Materialize the original full integer Hasse matrix before reducing it.
    """
    p, m = 31, 16
    precision = max(1, 48*e+3)
    modulus = p**precision
    nodes = (0, p**e, p**e*a)
    matrix = [[comb(c, j)*x**(c-j) if c >= j else 0
               for c in range(48)] for x in nodes for j in range(m)]
    matrix = [[v % modulus for v in row] for row in matrix]
    exponents = []
    while matrix:
        n = len(matrix)
        d, r, c = min((valuation(matrix[r][c]), r, c)
                      for r in range(n) for c in range(n))
        require(d < precision, 'finite-precision pivot visible')
        matrix[0], matrix[r] = matrix[r], matrix[0]
        for row in matrix:
            row[0], row[c] = row[c], row[0]
        scale = p**d
        unit = matrix[0][0]//scale
        inverse = pow(unit, -1, modulus)
        pivot = [(v*inverse) % modulus for v in matrix[0]]
        require(pivot[0] == scale, 'unit-only pivot normalization')
        reduced = []
        for row in matrix[1:]:
            require(row[0] % scale == 0, 'integral DVR elimination multiplier')
            multiplier = row[0]//scale
            reduced.append([(row[j]-multiplier*pivot[j]) % modulus for j in range(1, n)])
        matrix = reduced
        exponents.append(d)
    require(exponents == sorted(exponents), 'nondecreasing exact valuations')
    require(sum(exponents) == 768*e, 'original confluent determinant')
    return tuple(exponents)


def full_matrix_controls(roots):
    bank = [(a, e) for e in (0, 1, 2, 3) for a in (2, 3, 4, 11, 16, 29, 30)]
    bank += [(3+31*d, e) for d, e in ((1, 1), (2, 2), (-7, 3), (31**2+4, 1))]
    profiles = {}
    for a, e in bank:
        parts = smith_global_pivots(a, e)
        kappa = int(a % 31 in roots)
        require(parts[:16] == (0,)*16, 'zero-node identity block')
        require(sum(parts[:45]) == (631*e+1+kappa if e >= 1 else 0), 'full D45 equals residual E29')
        profiles[a, e] = parts
    for a in (3, 4):
        h = sum(comb(15, j)**2*a**j for j in range(16)) % 31
        require(h != 0 and a not in (2, 16, 30), 'same ordinary non-AP class')
        parts = profiles[a, 1]
        require(parts[-1] == 47, 'same largest exponent')
        print('a', a, 'full Hasse e1 Smith valuations', parts)
        print('a', a, 'D45', sum(parts[:45]), 'kernel exponent mod31^43',
              sum(min(43, v) for v in parts))
    require(profiles[3, 1][44:46] == (43, 43)
            and profiles[4, 1][44:46] == (42, 44), 'new intermediate partition hostile')
    print('fresh full original Hasse controls', len(bank), 'depths0..3 and higher lifts')


def main():
    print('Independent universal Hermite packet and all-depth p31 E29 audit')
    universal_controls()
    roots = packet_controls()
    weight_controls()
    full_matrix_controls(roots)
    print('PASS', gates, 'always-active gates')


if __name__ == '__main__':
    main()
