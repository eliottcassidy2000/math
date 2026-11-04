"""Exact controls for T1-T3; independent DP, words, cycles, squarefree products."""
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament(n, code):
    a = [0]*n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> bit & 1:
            a[u] |= 1 << v
        else:
            a[v] |= 1 << u
    return a


def opposite(a):
    return [((1 << len(a))-1) ^ (1 << u) ^ row for u, row in enumerate(a)]


def ends(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1 << v)
                dp[mask][v] = (sum(dp[rest][u] for u in range(n)
                                     if rest >> u & 1 and a[u] >> v & 1)
                              if rest else 1)
    return dp


def h(a):
    return sum(ends(a)[-1]) if a else 1


def literal_h(a):
    return sum(all(a[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(a))))


def extend(a, x, y):
    n = len(a)
    out = a[:] + [1 << (n+1), 0]
    for j, mask in enumerate((x, y)):
        for u in range(n):
            if mask >> u & 1:
                out[u] |= 1 << (n+j)
            else:
                out[n+j] |= 1 << u
    return out


def inflate(q, a):
    r = len(a)
    out = a[:] + [0]*(len(q)-1)
    for u in range(len(q)):
        for v in range(len(q)):
            if q[u] >> v & 1:
                for i in (range(r) if u == 0 else (r+u-1,)):
                    for j in (range(r) if v == 0 else (r+v-1,)):
                        out[i] |= 1 << j
    return out


def sidecar(a):
    n, full = len(a), (1 << len(a))-1
    e, s = ends(a), ends(opposite(a))
    c = [[0]*n for _ in a]
    for mask in range(1, full):
        for u in range(n):
            for v in range(n):
                c[u][v] += e[mask][u]*s[full ^ mask][v]
    return s, e, c


def pxor(*ps):
    out = set()
    for p in ps:
        out.symmetric_difference_update(p)
    return out


def pmul(p, q):
    out = set()
    for u in p:
        for v in q:
            if not u & v:
                key = u | v
                if key in out:
                    out.remove(key)
                else:
                    out.add(key)
    return out


def sector(p, parity):
    return sum(mask.bit_count() % 2 == parity for mask in p) % 2


def path_matrix(a):
    n = len(a)
    g = [[set() for _ in a] for _ in a]
    for start in range(n):
        dp = [0]*(1 << n)
        dp[1 << start] = 1 << start
        for mask in range(1, 1 << n):
            for u in range(n):
                if dp[mask] >> u & 1:
                    g[start][u].add(mask)
                    todo = a[u] & ~mask
                    while todo:
                        bit = todo & -todo
                        todo -= bit
                        dp[mask | bit] ^= bit
    return g


def contrast_formula(a, x, y, g):
    """(F_T-F_Top)/2 mod2 from T2, including the paired two-path term."""
    n = len(a)
    s = [pxor(*row) for row in g]
    e = [pxor(*(g[u][v] for u in range(n))) for v in range(n)]
    p = pxor(*(s[u] for u in range(n) if x >> u & 1))
    q = pxor(*(s[u] for u in range(n) if y >> u & 1))
    X = pxor(*(g[u][v] for u in range(n) if x >> u & 1
               for v in range(n) if y >> v & 1))
    esx = pxor(*(pxor(e[u], s[u]) for u in range(n) if x >> u & 1))
    esy = pxor(*(pxor(e[u], s[u]) for u in range(n) if y >> u & 1))
    middle = pmul(set(range(1, 1 << n)), pmul(pxor(p, q), X))
    return (sector(pxor(esx, esy), 0)
            ^ sector(pxor(esx, pmul(p, q), middle), 1))


def cubic(a, u, v, tables=None):
    s, e = tables if tables is not None else (ends(opposite(a)), ends(a))
    return sum(s[m][u]+e[m][u] for m in range(1 << len(a))
               if not m >> v & 1 and m.bit_count() % 2 == 0) % 2


def direct_cubic(a, u, v):
    value = 0
    ao = opposite(a)
    # Coefficient of a_u a_v c_v. All unselected variables are zero.
    for mask in range(8):
        x = ((mask & 1) << u) | (((mask >> 1) & 1) << v)
        y = ((mask >> 2) & 1) << v
        value ^= ((h(extend(a, x, y))-h(extend(ao, x, y)))//2) & 1
    return value


def odd_cycles(a):
    n = len(a)
    total = 0
    for k in range(3, n+1, 2):
        for vertices in combinations(range(n), k):
            root = vertices[0]
            for tail in permutations(vertices[1:]):
                p = (root,) + tail + (root,)
                total += all(a[u] >> v & 1 for u, v in zip(p, p[1:]))
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = Path(__file__).read_bytes()
    counts = dict(cores=0, path_skew_entries=0, cubic_entries=0,
                  actual_cubic_queries=0, response_queries=0,
                  literal_response_queries=0, ambient_pairs=0, cycle_checks=0,
                  integer_boolean_coefficients=0, nonzero_exact_cubics=0)
    controls = [tournament(n, code) for n in range(1, 6)
                for code in range(1 << (n*(n-1)//2))]
    controls += [tournament(6, 408), [4, 17, 42, 131, 45, 75, 31, 119]]
    witnesses = []
    for a in controls:
        n = len(a)
        s, e, c = sidecar(a)
        g = path_matrix(a)
        counts['cores'] += 1
        for u in range(n):
            for v in range(n):
                if u == v:
                    continue
                odd_path = sector(pxor(g[u][v], g[v][u]), 1)
                need((c[u][v]-c[v][u]) % 2 == 0, 'cover skew even')
                need(((c[u][v]-c[v][u])//2) % 2 == odd_path, 'T3 odd paths')
                counts['path_skew_entries'] += 1
                z = cubic(a, u, v, (s, e))
                if n % 2 == 0:
                    simple = sum(s[m][u] for m in range(1 << n)
                                 if not m >> v & 1 and m.bit_count() % 2 == 0
                                 and m.bit_count() <= n-4) % 2
                    need(z == simple, 'T2 even simplification')
                if n <= 4 or (u, v) in ((0, 1), (1, 0)):
                    need(z == direct_cubic(a, u, v), 'actual cubic coefficient')
                    counts['actual_cubic_queries'] += 8
                counts['cubic_entries'] += 1
        masks = range(1 << n) if n <= 4 else sorted({0, 1, 2, 3, (1 << n)-1})
        difference_table = [0]*(1 << (2*n)) if n <= 4 else None
        for x in masks:
            for y in masks:
                first, second = extend(a, x, y), extend(opposite(a), x, y)
                values = h(first), h(second)
                need((values[0]-values[1]) % 2 == 0, 'odd path totals')
                need(((values[0]-values[1])//2) % 2 == contrast_formula(a, x, y, g),
                     'T2 complete contrast formula')
                counts['response_queries'] += 1
                if difference_table is not None:
                    difference_table[x | (y << n)] = values[0]-values[1]
                if n <= 3 or (n == 4 and x == y == 1):
                    need(values == (literal_h(first), literal_h(second)), 'literal response')
                    counts['literal_response_queries'] += 1
        if difference_table is not None:
            for bit in range(2*n):
                for mask in range(1 << (2*n)):
                    if mask >> bit & 1:
                        difference_table[mask] -= difference_table[mask ^ (1 << bit)]
            for mask, coefficient in enumerate(difference_table):
                x, y = mask & ((1 << n)-1), mask >> n
                if mask.bit_count() > 3 or (mask.bit_count() == 2 and (not x or not y)):
                    need(coefficient == 0, 'T4 exact cubic and pure-quadratic bound')
                counts['integer_boolean_coefficients'] += 1
                counts['nonzero_exact_cubics'] += bool(coefficient and mask.bit_count() == 3)
        if n <= 4 or a in controls[-2:]:
            need(h(a) % 4 == (1+2*odd_cycles(a)) % 4, 'OCF residue control')
            counts['cycle_checks'] += 1
        if n >= 6:
            witnesses.append(dict(rows=a, H=h(a), endpoint_parity_balanced=
                                  all((x-y) % 2 == 0 for x, y in zip(s[-1], e[-1])),
                                  skew_symmetric_mod4=all((c[u][v]-c[v][u]) % 4 == 0
                                                         for u in range(n) for v in range(n)),
                                  nonzero_cubic=[(u, v) for u in range(n) for v in range(n)
                                                 if u != v and cubic(a, u, v, (s, e))]))
    # Ports are inside the substituted module and have its uniform exterior incidences.
    for r in range(1, 4):
        for code in range(1 << (r*(r-1)//2)):
            a = tournament(r, code)
            for b in (a, opposite(a), tournament(r, 0)):
                for qn in range(1, 4):
                    for qc in range(1 << (qn*(qn-1)//2)):
                        q = tournament(qn, qc)
                        for x, y in ((0, 0), (1, 2 if r > 1 else 1), ((1 << r)-1, 1)):
                            aa, bb = extend(a, x, y), extend(b, x, y)
                            difference = h(inflate(q, aa))-h(inflate(q, bb))
                            need((difference-h(q)*(h(aa)-h(bb))) % 4 == 0,
                                 'T1 arbitrary-pair ambient response')
                            counts['ambient_pairs'] += 1
    need(Path(__file__).read_bytes() == source, 'source changed')
    result = dict(status='FINITE-EXACT controls for proved T1-T4',
                  universe='All 1099 labelled cores orders 1..5; named 6/8 hostiles. '
                           'Complete two-port tables through 4; specified mask controls above. '
                           'All labelled modules orders 1..3 and quotients orders 1..3; '
                           'three module comparisons and three two-port masks each.',
                  counts=counts, witnesses=witnesses, source_sha256=sha256(source).hexdigest())
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
