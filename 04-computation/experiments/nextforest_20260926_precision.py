"""Exact precision collisions for odd Collatz, with capped actual-source scans."""
from collections import Counter


def need(ok, context):
    if not ok:
        raise RuntimeError(context)


def valuation(n):
    return (n & -n).bit_length() - 1


def oddpart(n):
    return n >> valuation(n)


def U(n):
    return oddpart(3 * n + 1)


def step(state):
    q, r, u = state
    a = valuation(3 * u + 1)
    v = (3 * u + 1) >> a
    if a < r:
        return (3 * q, r - a, v), None, 'consume'
    if a > r:
        return (v, a - r, 3 * q), None, 'swap'
    return None, oddpart(3 * q + v), 'collision'


def value(state):
    q, r, u = state
    return (q << r) + u


def split(y):
    need(y > 1 and y & 1, ('split domain', y))
    r = y.bit_length() - 1
    return 1, r, y - (1 << r)


def main():
    counts = Counter()
    canonical = 0
    for q in range(1, 64, 2):
        for u in range(1, 128, 2):
            for r in range(1, 11):
                state = q, r, u
                nxt, end, tag = step(state)
                actual = U(value(state))
                need((value(nxt) if nxt else end) == actual, ('exact', state))
                counts[tag] += 1
                if u < 1 << r:
                    a = valuation(3 * u + 1)
                    if a >= r:
                        v = (3 * u + 1) >> a
                        need(v == 1, ('canonical v', state))
                        if a == r:
                            need(r % 2 == 0 and actual == U(q) and actual < value(state), ('canonical collision', state))
                        else:
                            need(r % 2 == 1 and a == r + 1 and actual == 3*q+2, ('canonical swap', state))
                            need((actual < value(state)) == (r >= 3), ('canonical sign', state))
                        canonical += 1
    print('exact triple universe: q odd1..63,u odd1..127,r1..10:', dict(sorted(counts.items())))
    print('canonical boundary cases:', canonical)

    # The controller explicitly charges one extra U before a fresh split.
    visits = Counter()
    stopped = capped = 0
    for source in range(3, 2048, 2):
        n = source
        state = None
        for tick in range(96):
            if n == 1:
                stopped += 1
                break
            if state is None:
                n = U(n)
                visits['canonical_reset_U'] += 1
                state = split(n) if n > 1 else None
            else:
                nxt, endpoint, tag = step(state)
                n = U(n)
                need(n == (value(nxt) if nxt else endpoint), ('source prefix', source, tick))
                visits[tag] += 1
                state = nxt
        else:
            capped += 1
    print('actual sources odd3..2047, cap96 U steps:', dict(sorted(visits.items())))
    print('within cap reached1:', stopped, '; capped:', capped, '(no inferred eventual collision or termination)')

    expected = [(1,5,9),(3,3,7),(9,2,11),(27,1,17),(13,1,81),(61,1,39)]
    state = split(U(27))
    for target in expected:
        need(state == target, ('27 state', state, target))
        nxt, end, _ = step(state)
        state = nxt
    need(end == 121 and U(121) == 91 and split(91) == (1,6,27), '27 reentry')
    n = 91
    path = [n]
    for _ in range(4):
        n = U(n)
        path.append(n)
    need(path == [91,137,103,155,233], ('233 path', path))
    print('27 split path:', expected, 'then collision121, reset91=(1,6,27)')
    print('samecore9 higher precision:', path)

    for k in range(2, 129):
        u = (4**k - 1)//3
        nxt, end, tag = step((1,1,u))
        need(tag == 'swap' and nxt == (1,2*k-1,3), ('recreated precision', k))
        need(U(u+2) == 2**(2*k-1)+3, ('hostile actual', k))
    print('unbounded precision-recreation family: k=2..128 checked')

    rows = []
    for r in range(1, 258, 2):
        n = (2**(r+1)+17)//3
        current = n
        orbit = []
        for j in range(1, 7):
            current = U(current)
            orbit.append(current)
        need(orbit[0] == 2**r+9, ('family source', r))
        if r >= 11:
            explicit = [2**r+9,3*2**(r-2)+7,9*2**(r-3)+11,
                        27*2**(r-4)+17,81*2**(r-6)+13,243*2**(r-9)+5]
            need(orbit == explicit, ('six step identity', r))
            need(all(y > n for y in orbit[:5]) and orbit[5] < n, ('first descent6', r))
        if r == 9:
            need(orbit[5] == 31 and all(y > n for y in orbit[:5]), 'r9 boundary')
        if r <= 15:
            rows.append((r,n,orbit))
    print('family first six odd iterates:', rows)
    print('first descent exactly6 for every odd r>=9 proved; finite identities through r257 checked')
    print('PASS: explicit checks active under -O')


if __name__ == '__main__':
    main()
