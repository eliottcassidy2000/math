#!/usr/bin/env python3
"""Exact controls: inherited R clock, state strata, boundary probes and all widths.

Standard library only. No imported research implementation; fail closed even -O.
Outputs are written only when run as a script.
"""
from collections import Counter
from itertools import combinations, product
from math import comb, gcd, lcm
from pathlib import Path
import json

SEEDS = (153, 8345, 8337, 8369, 177, 185)
PAIRS = tuple(combinations(range(6), 2))


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def step(v):
    n = len(v) // 2
    a, b = v[:n], v[n:]
    return tuple(b) + tuple((a[i]+b[i]+b[(i+1) % n]) % 3 for i in range(n))


def lift(v):
    n = len(v) // 2
    a, b = list(v[:n]), list(v[n:])
    for track in (0, 1, 2, 3):
        if track == 0:
            target, values = a, [a[i]+b[i]+b[(i+1) % n] for i in range(n)]
        elif track == 1:
            target, values = a, [a[i]+b[i] for i in range(n)]
        elif track == 2:
            target, values = b, [a[i]-b[i] for i in range(n)]
        else:
            target, values = a, [a[i]-b[i] for i in range(n)]
        for i, value in enumerate(values):
            old = target[i]
            delta = (0, 1, -1)[(value-old) % 3]
            target[i] = (old+delta) % 6
            need((SEEDS[old] ^ SEEDS[target[i]]).bit_count() == bool(delta),
                 'each nontrivial gate is one literal arc flip')
    return tuple(a+b)


def period(v, update, limit):
    w = v
    for t in range(1, limit+1):
        w = update(w)
        if w == v:
            return t
    raise RuntimeError('period search cap exceeded')


def depth(v):
    """t-adic depth of a,b in F3[z]/((z-1)^n), polynomial coefficients."""
    n = len(v)//2
    val = n
    for row in (v[:n], v[n:]):
        for j in range(n):
            if sum(comb(i, j)*row[i] for i in range(j, n)) % 3:
                val = min(val, j)
                break
    return n-val


def pow3ceil(d):
    p = 1
    while p < d:
        p *= 3
    return p


def census(n, base, update):
    unseen = set(product(range(base), repeat=2*n))
    lengths = Counter()
    while unseen:
        v = min(unseen)
        w, size = v, 0
        while w in unseen:
            unseen.remove(w)
            w = update(w)
            size += 1
        need(w == v, 'bijection: no tails or previously visited cycles')
        lengths[size] += 1
    return dict(sorted(lengths.items()))


def rows(mask):
    out = [0]*6
    for k, (i, j) in enumerate(PAIRS):
        u, v = (i, j) if mask >> k & 1 else (j, i)
        out[u] |= 1 << v
    return out


def subsets(r):
    n = len(r)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s in range(1, 1 << n):
        for i, c in enumerate(dp[s]):
            bits = r[i] & ~s
            while bits:
                b = bits & -bits
                bits -= b
                dp[s | b][b.bit_length()-1] += c
    return dp


def birth(r, cut):
    n = len(r)
    out = r.copy()+[0]
    for i in range(n):
        if cut >> i & 1:
            out[i] |= 1 << n
        else:
            out[n] |= 1 << i
    return out


def hp(r):
    return sum(subsets(r)[-1])


def kroot(r, q):
    dp = subsets(r)
    hs = [sum(row) for row in dp]
    hs[0] = 1
    allbits = (1 << len(r))-1
    return sum(dp[s][q]*hs[allbits ^ s] for s in range(1, allbits+1))


def rank(a):
    a = [list(row) for row in a]
    row = 0
    for col in range(len(a[0])):
        piv = next((i for i in range(row, len(a)) if a[i][col] % 3), None)
        if piv is None:
            continue
        a[row], a[piv] = a[piv], a[row]
        inv = a[row][col] % 3  # both nonzero elements are self-inverse
        a[row] = [(inv*x) % 3 for x in a[row]]
        for i in range(len(a)):
            if i != row:
                c = a[i][col]
                a[i] = [(x-c*y) % 3 for x, y in zip(a[i], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def observer(n):
    columns = []
    for i in range(2*n):
        v = tuple(int(i == j) for j in range(2*n))
        col = []
        for _ in range(2*n):
            col.append(v[0])
            v = step(v)
        columns.append(col)
    return [list(row) for row in zip(*columns)]


# Polynomials over F3, coefficients low first; [] is zero.
def trim(p):
    p = [x % 3 for x in p]
    while p and p[-1] == 0:
        p.pop()
    return p


def add(a, b, scale=1):
    return trim([(a[i] if i < len(a) else 0) +
                 scale*(b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def mul(a, b):
    out = [0]*max(0, len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def divmodp(a, b):
    need(bool(b), 'polynomial denominator')
    a = trim(a)
    out = [0]*max(0, len(a)-len(b)+1)
    while a and len(a) >= len(b):
        shift = len(a)-len(b)
        c = a[-1]*b[-1] % 3
        out[shift] = c
        a = add(a, [0]*shift+[c*x for x in b], -1)
    return trim(out), a


def pgcd(a, b):
    while b:
        a, b = b, divmodp(a, b)[1]
    return trim([x*a[-1] for x in a]) if a else []


def ppow(a, n, mod=None):
    out = [1]
    while n:
        if n & 1:
            out = mul(out, a)
            if mod:
                out = divmodp(out, mod)[1]
        n //= 2
        a = mul(a, a)
        if mod:
            a = divmodp(a, mod)[1]
    return out


def factors(n):
    out = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d = 3 if d == 2 else d+2
    if n > 1:
        out.append(n)
    return out


def characteristic(n):
    return add(ppow([2, 2, 1], n), [0]*n+[1], -1)


def clock_order(n):
    m, threes = n, 1
    while m % 3 == 0:
        m //= 3
        threes *= 3
    p = characteristic(m)
    derivative = trim([i*p[i] for i in range(1, len(p))])
    repeated = pgcd(p, derivative)
    need(repeated == ([1, 0, 1] if m % 8 == 0 else [1]),
         'only the order-eight spatial modes ramify')
    squarefree, remainder = divmodp(p, repeated)
    need(not remainder, 'exact squarefree quotient')
    rest, prime_order, groups = squarefree, 1, []
    for d in range(1, len(squarefree)):
        if rest == [1]:
            break
        group = pgcd(rest, add(ppow([0, 1], 3**d, rest), [0, 1], -1))
        if len(group) == 1:
            continue
        candidate = 3**d-1
        for prime in factors(candidate):
            while candidate % prime == 0 and ppow([0, 1], candidate//prime, group) == [1]:
                candidate //= prime
        need(ppow([0, 1], candidate, group) == [1], 'root order upper bound')
        groups.append({'irreducible_degree': d, 'polynomial': group, 'root_lcm': candidate})
        prime_order = lcm(prime_order, candidate)
        rest, rem = divmodp(rest, group)
        need(not rem, 'distinct-degree division')
    need(rest == [1], 'all factors accounted for')
    order = prime_order*threes*(3 if m % 8 == 0 else 1)
    pn = characteristic(n)
    need(ppow([0, 1], order, pn) == [1], 'full polynomial order')
    for prime in factors(order):
        need(ppow([0, 1], order//prime, pn) != [1], 'minimal polynomial order')
    return {'n': n, 'order': order, 'prime_to_3_part': prime_order,
            'ramification': m % 8 == 0, 'groups': groups}


def main():
    probes = []
    for cut in range(64):
        values = [hp(birth(rows(seed), cut)) for seed in SEEDS]
        if len(set(values)) == 6:
            probes.append({'cut': cut, 'counts': values})
    expected = [109, 115, 101, 77, 73, 81]
    need(probes[0] == {'cut': 2, 'counts': expected}, 'one-vertex six-state reader')
    for seed, value in zip(SEEDS, expected):
        r = rows(seed)
        reverse = [sum(1 << j for j in range(6) if r[j] >> i & 1) for i in range(6)]
        starts = subsets(reverse)[-1][1]
        need(hp(r) == 29 and 29-starts+kroot(r, 1) == value, 'rooted split formula')
    strata = []
    for n in (1, 3):
        bydepth = Counter()
        for v in product(range(3), repeat=2*n):
            d = depth(v)
            p = 1 if d == 0 else 8*pow3ceil(d)
            need(period(v, step, 8*n) == p, 'depth predicts least period')
            bydepth[(d, p)] += 1
        need(all(count == (1 if d == 0 else 8*9**(d-1))
                 for (d, _), count in bydepth.items()), 'depth populations')
        logical = census(n, 3, step)
        physical = census(n, 6, lift)
        need(physical == {p: c*4**n for p, c in logical.items()}, 'all-sheet cycle census')
        strata.append({'n': n, 'depth_states': [[d, p, c] for (d, p), c in sorted(bydepth.items())],
                       'logical_cycles': logical, 'physical_cycles': physical})
    large = []
    for n in (9, 27, 81):
        v = (1,)+(0,)*(2*n-1)
        lp, pp = period(v, step, 8*n), period(v, lift, 8*n)
        need(lp == pp == 8*n, 'delta clock inherits exact period')
        w = v
        for _ in range(4*n):
            w = step(w)
        need(w == tuple(-x % 3 for x in v), 'half-lap antisymmetry')
        large.append({'n': n, 'logical_period': lp, 'physical_period': pp})
    observability = []
    for n in list(range(1, 21))+[27]:
        matrix = observer(n)
        need(rank(matrix) == 2*n and rank(matrix[:-1]) == 2*n-1,
             'sharp linear boundary history length')
        observability.append({'n': n, 'rank': 2*n, 'one_fewer_rank': 2*n-1})
    arbitrary = [clock_order(n) for n in list(range(1, 13))+[16]]
    for n in (2, 4, 5, 7, 8):
        entry = next(e for e in arbitrary if e['n'] == n)
        v = (1,)+(0,)*(2*n-1)
        need(period(v, step, entry['order']) == entry['order'], 'cyclic impulse independent orbit')
        entry['delta_physical_period'] = period(v, lift, 2*entry['order'])
    carry_census = []
    for n in (1, 2, 3, 4):
        states = [(v[0]+3*b,)+v[1:] for v in product(range(3), repeat=2*n) for b in (0, 1)]
        histories = [() for _ in states]
        counts = []
        for length in range(1, 2*n+1):
            histories = [h+(v[0],) for h, v in zip(histories, states)]
            count = len(set(histories))
            expected = 6 if length == 1 else (2*9**n if length == 2*n else 10*3**(length-1))
            need(count == expected, 'exact nonlinear history partition')
            counts.append(count)
            following = [lift(v) for v in states]
            for v, w in zip(states, following):
                logical = tuple(x % 3 for x in v)
                y0, y1, y2 = logical[0], step(logical)[0], step(step(logical))[0]
                signed = lambda x: (0, 1, -1)[x % 3]
                carry = signed(y2-y0)+signed(y1)-signed(y2)
                need((w[0]-v[0]-carry) % 6 == 0, 'three-trit physical carry formula')
            states = following
        carry_census.append({'n': n, 'histories_by_length': counts})
    # An unobserved deck bit is forever invisible at the observed cell.
    v = (1, 0, 0, 0, 0, 0)
    w = (1, 3, 0, 0, 0, 0)
    for _ in range(48):
        need(v[0] == w[0] and v != w, 'remote deck hostile')
        v, w = lift(v), lift(w)
    result = {'scope': 'exact finite controls; uniform claims have proofs in PROOFS.md',
              'selector_bank': probes, 'strata': strata, 'large_basis_clocks': large,
              'boundary_observability': observability, 'arbitrary_widths': arbitrary,
              'physical_boundary_histories': carry_census,
              'remote_deck_hostile': 'same a0 forever; different labelled states'}
    Path(__file__).with_name('results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'selectors': len(probes), 'strata': strata,
                      'orders': [(e['n'], e['order']) for e in arbitrary]}, indent=2))


if __name__ == '__main__':
    main()
