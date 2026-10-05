"""Exact 5/8/9 connections and the boundary between clocks and paid routes.

Standard library; python -B and -O give identical output. No import-time census.
"""
from fractions import Fraction as F
from itertools import permutations, product


def need(ok, message):
    if not ok:
        raise ValueError(message)


def mul(x, y, modulus=None):
    a, b = x
    c, d = y
    z = a*c+b*d, a*d+b*c+b*d
    return tuple(v % modulus for v in z) if modulus else z


def power(x, k, modulus=None):
    result = (1, 0)
    while k:
        if k & 1:
            result = mul(result, x, modulus)
        x = mul(x, x, modulus)
        k //= 2
    return result


def norm(x):
    a, b = x
    return a*a+a*b-b*b


def sign(x):
    a, b = map(F, x)
    if not b:
        return (a > 0)-(a < 0)
    if b < 0:
        return -sign((-a, -b))
    if a >= 0:
        return 1
    t = -a/b
    if t <= 1:
        return 1
    if t >= 2:
        return -1
    return 1 if t*t-t-1 < 0 else -1


def vp(n, p):
    need(type(n) is int and n != 0, 'nonzero exact integer valuation')
    n, k = abs(n), 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def actual(n, word):
    need(type(n) is int and n > 0 and n % 2, 'positive odd source')
    for a in word:
        need(n != 1 and vp(3*n+1, 2) == a, 'first-hit exact valuation guard')
        n = (3*n+1) >> a
    return n


def carrier(word):
    p, q, b = 1, 1, 0
    for a in word:
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def field_controls():
    orbit = tuple(power((0, 1), k, 3) for k in range(8))
    need(len(set(orbit)) == 8 and set(orbit) == set(product(range(3), repeat=2))-{(0, 0)},
         'golden F9 nonzero orbit')
    need(power((0, 1), 4, 3) == (2, 0), 'golden half-turn')
    need(tuple(norm(x) % 3 for x in orbit) == (1, 2)*4, 'alternating norm')
    need(mul((2, 1), (2, 1), 5) == (0, 0), 'ramified5 nilpotent hostile')
    need(power((0, 1), 20, 5) == (1, 0) and
         all(power((0, 1), d, 5) != (1, 0) for d in (4, 10)), 'mod5 clock20')
    print('Discriminant5, field F9, golden clock8:', orbit)
    print('Norm alternates1,2; mod5 phi-3 is nonzero square-zero, clock20.')
    total = 0
    for a in range(1, 5):
        m = 3**a
        states = tuple(product(range(m), repeat=2))
        need(len({mul((0, 1), v, m) for v in states}) == m*m, 'golden permutation')
        q = (a+1)//2
        factor = pow(9*pow(8, -1, m), q, m)
        need(all((factor*x % m, factor*y % m) == (0, 0) for x, y in states),
             'shifted G is nilpotent on the ternary register')
        period = 8*3**(a-1)
        need(power((0, 1), period, m) == (1, 0) and
             power((0, 1), period//2, m) != (1, 0), 'golden order2 factor')
        if a > 1:
            need(power((0, 1), period//3, m) != (1, 0), 'golden order3 factor')
        total += len(states)
    # Complete linear intertwiners F3^2 -> F3^2 for A*M=zero*A.
    intertwiners = []
    for a, b, c, d in product(range(3), repeat=4):
        if (b % 3, (a+b) % 3, d % 3, (c+d) % 3) == (0, 0, 0, 0):
            intertwiners.append((a, b, c, d))
    need(intertwiners == [(0, 0, 0, 0)], 'only zero linear intertwiner')
    print('All ternary ring states at precisions1..4:', total,
          '; all81 linear maps: only zero intertwines the two dynamics.')


def golden_return_control():
    x = (F(1, 3), F(1, 3))
    start, digits, phases = x, [], []
    for _ in range(8):
        phases.append(tuple(int(3*v) % 3 for v in x))
        y = mul((0, 1), x)
        digit = int(sign((y[0]-1, y[1])) > 0)
        digits.append(digit)
        x = y[0]-digit, y[1]
        need(sign(x) > 0 and sign((x[0]-1, x[1])) < 0, 'interior golden phase')
    need(x == start and ''.join(map(str, digits)) == '10100000' and len(set(phases)) == 8,
         'golden8-cycle and its ordered parity word')
    p, q, b = carrier((1, 5))
    need((p, q, b) == (9, 64, 5), 'rational realization carry')
    value, orbit = F(1, 11), []
    for digit in digits:
        orbit.append(value)
        need(value.denominator % 2 and value.numerator % 2 == digit, 'rational parity gate')
        value = 3*value+1 if digit else value/2
    need(value == F(1, 11) and len(set(orbit)) == 8, 'rational exact period')
    need(F(b, q-p) == F(1, 11) and b % (q-p) != 0, 'integer-realization hostile')
    print('Golden theta=(1+phi)/3: parity10100000; arithmetic return(9n+5)/64, fixed1/11.')
    print('Same carry, different budget: word12 fixes-5; word15 is legal35mod128 and35->53->5.')
    need(carrier((1, 2)) == (9, 8, 5) and actual(35, (1, 5)) == 5, 'cost-switch controls')
    for a in range(1, 9):
        p, q, b = carrier((1, a))
        need(p == 9 and b == 5 and q == 2**(a+1), 'all final-exponent costs')
        r = ((q-b)*pow(p, -1, 2*q)) % (2*q)
        for n in (r, r+2*q, r+8*q):
            need(actual(n, (1, a)) == F(p*n+b, q), 'independent legal source')


def guard_controls():
    count = 0
    for n in range(1, 1024, 2):
        fuel = (vp(n+5, 2)-1)//3
        for q in range(fuel+1):
            y = actual(n, (1, 2)*q)
            need(F(9**q*(n+5), 8**q)-5 == y, 'signed anchor formula')
            need(vp(y+5, 2) == vp(n+5, 2)-3*q and
                 vp(y+5, 3) == vp(n+5, 3)+2*q, 'dual-prime fuel transport')
            need(F(8**q*(y+5), 9**q)-5 == n, 'inverse returns original source')
            if q:
                need(y % 9 == 4 and (8*y-5)//9 < y, 'mod9 smaller predecessor')
            count += 1
    need(F(9*3+5, 8) == 4 and (vp(3+5, 2)-1)//3 == 0, 'integral-even endpoint hostile')
    for t in range(64):
        n = 27+64*t
        need((vp(n+5, 2)-1)//3 == 1, '27cell exactly one12 block')
        y = actual(n, (1, 2))
        need(y == 31+72*t and y % 8 == 7 and y % 2048 != 155, 'H exit exclusion')
        need(actual(n, (1, 2, 1, 1)) == F(81*n+85, 32), 'forced four-letter corridor')
        need((8*y-5)//9 == n, 'mod9 round trip does not pay original')
    print('Odd sources1..1023:', count, 'legal repeat prefixes with exact binary fuel and ternary contraction.')
    print('27mod64: exactly one12; output31+72t is7mod8, blocks immediate H; mod9 predecessor returns source.')
    for t in range(64):
        n, child, smaller = 219+256*t, 139+162*t, 123+144*t
        checkpoint = actual(n, (1, 2, 1, 1))
        need(checkpoint == 4*child+1 and F(81*n+53, 128) == child,
             'quarter-child interface')
        need(F(8*child-5, 9) == smaller == F(9*n-3, 16) and
             actual(smaller, (1, 2)) == child, 'positive inverse-G composition')
        need(vp(child+5, 3) == 2 and vp(smaller+5, 3) == 0,
             'exactly one inverse block in ternary coordinate')
        a = vp(3*checkpoint+1, 2)
        need(a >= 3 and actual(checkpoint, (a,)) == actual(child, (a-2,)),
             'source and child actual common future')
        ratio = F(smaller+5, n+5)*F(9, 8)**4
        need(smaller < child < n and ratio <= F(6561, 7168), 'four-credit interface control')
    print('Positive composition219mod256: K=139+162t, inverse-G child L=123+144t; four-credit ratio<=6561/7168.')


def graph(mask):
    flips = ((0, 2), (1, 3), (0, 3))
    arcs = {(i, j) for i in range(4) for j in range(i+1, 4)}
    for k, edge in enumerate(flips):
        if mask >> k & 1:
            arcs.remove(edge)
            arcs.add(edge[::-1])
    return arcs


def canon(arcs):
    pairs = tuple((i, j) for i in range(4) for j in range(i+1, 4))
    return min(tuple(int((p[i], p[j]) in arcs) for i, j in pairs) for p in permutations(range(4)))


def paths(n, arcs):
    dp = {(1 << v, v): 1 for v in range(n)}
    for size in range(1, n):
        for (mask, last), count in tuple(dp.items()):
            if mask.bit_count() != size:
                continue
            for nxt in range(n):
                if not mask >> nxt & 1 and (last, nxt) in arcs:
                    key = mask | 1 << nxt, nxt
                    dp[key] = dp.get(key, 0)+count
    return sum(dp.get(((1 << n)-1, last), 0) for last in range(n))


def scc_sizes(n, arcs):
    reach = [[i == j or (i, j) in arcs for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] |= reach[i][k] and reach[k][j]
    classes = {frozenset(j for j in range(n) if reach[i][j] and reach[j][i]) for i in range(n)}
    return tuple(len(c) for c in sorted(classes, key=lambda c: -sum(reach[min(c)])))


def tournament_controls():
    classes = {}
    for mask in range(8):
        classes.setdefault(canon(graph(mask)), []).append(mask)
    need(sorted(classes.values()) == [[0], [1], [2], [3, 4, 5, 6, 7]], '8 presentations,5 strong')
    need(tuple(paths(4, graph(k)) for k in range(8)) == (1, 3, 3, 5, 5, 5, 5, 5), 'Hamiltonian5')
    need(canon(graph(3)) == canon(graph(4)) and canon(graph(2)) != canon(graph(5)), 'flip memory loss')
    joins = []
    for left, right in ((1, 2), (2, 1)):
        arcs = graph(left) | {(i+4, j+4) for i, j in graph(right)} | {(i, j) for i in range(4) for j in range(4, 8)}
        need(paths(8, arcs) == 9, 'order-join path count9')
        joins.append(scc_sizes(8, arcs))
    need(joins == [(3, 1, 1, 3), (1, 3, 3, 1)], 'intrinsic order sidecar')
    print('Tournament8 masks ->4 classes; strong fiber5; diamond joins each H9, SCC words', joins)


if __name__ == '__main__':
    field_controls()
    golden_return_control()
    guard_controls()
    tournament_controls()
    print('PASS: exact transfer maps and hostile boundaries; no all-source guard coverage claim.')
