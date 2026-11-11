"""Independent literal audit: no imports from construction/counting lanes.

Vertices have persistent names; all original arcs are checked explicitly.
Endpoint DP is controlled against literal permutations at small orders.
"""
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def source(k):
    names = ('r', 'a', 'b') + tuple('p'+str(i) for i in range(k)) + ('B', 'M', 'A')
    letters = (3,)*k + (2, 3, 1)
    edges = {('a', 'r'), ('b', 'r'), ('a', 'b')}
    chain = names[3:]
    edges.update(('r', x) for x in chain)
    edges.update((x, y) for i, x in enumerate(chain) for y in chain[i+1:])
    for x, letter in zip(chain, letters):
        for j, port in enumerate(('a', 'b')):
            edges.add((x, port) if letter & (1 << j) else (port, x))
    return names, edges


def gadget(names, edges, t):
    new = tuple(str(t)+x for x in ('x', 'y', 'z'))
    x, y, z = new
    found = set(edges)
    for v in names:
        found.add((v, x))
        found.add((v, y))
        found.add((z, v) if v == 'B' else (v, z))
    found.update(((x, y), (y, z), (z, x)))
    return names+new, found


def masks(names, edges):
    rows = tuple(sum(1 << j for j, w in enumerate(names) if (v, w) in edges)
                 for v in names)
    for i in range(len(rows)):
        need(not rows[i] >> i & 1, 'loop')
        for j in range(i):
            need(bool(rows[i] >> j & 1) != bool(rows[j] >> i & 1), 'tournament')
    return rows


def tables(rows):
    n = len(rows)
    end = [[0]*n for _ in range(1 << n)]
    hp = [1] + [0]*((1 << n)-1)
    incoming = [sum(1 << u for u in range(n) if rows[u] >> v & 1) for v in range(n)]
    for mask in range(1, 1 << n):
        bits = mask
        while bits:
            bit = bits & -bits
            bits -= bit
            v = bit.bit_length()-1
            prev = mask ^ bit
            if not prev:
                end[mask][v] = 1
            else:
                allowed = prev & incoming[v]
                count = 0
                while allowed:
                    b = allowed & -allowed
                    allowed -= b
                    count += end[prev][b.bit_length()-1]
                end[mask][v] = count
        hp[mask] = sum(end[mask])
    return hp, end


def strong(rows):
    n = len(rows)
    for initial in range(n):
        seen = 1 << initial
        while True:
            old = seen
            for v in range(n):
                if seen >> v & 1:
                    seen |= rows[v]
            if seen == old:
                break
        if seen != (1 << n)-1:
            return False
    return True


def profile(names, edges, literal=False):
    rows = masks(names, edges)
    hp, ends = tables(rows)
    full = (1 << len(rows))-1
    root = names.index('r')
    @lru_cache(None)
    def starts(mask, v):
        rest = mask ^ (1 << v)
        if not rest:
            return 1
        return sum(starts(rest, w) for w in range(len(rows))
                   if rest >> w & 1 and rows[v] >> w & 1)
    h = hp[full]
    s = starts(full, root)
    f = sum(ends[m][root]*hp[full ^ m] for m in range(1, full+1) if m >> root & 1)
    residual = full
    for x in ('r', 'a', 'b'):
        residual ^= 1 << names.index(x)
    out = dict(order=len(rows), H=h, q=hp[full ^ (1 << root)], s=s,
               defect=f-h, residual_H=hp[residual], strong=strong(rows))
    if literal:
        count = sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(rows)-1))
                    for p in permutations(range(len(rows))))
        need(h == count, 'literal permutation control')
        out['literal_permutation_H'] = count
    return out


def exact(k, t):
    d = 1 << k
    if t == 0:
        return (14*d+3, 5, 2, 1)
    h_num = (57*d+8)*12**t + ((12*t+66)*d+2*t+17)*3**t
    q_num = 24*12**t+(6*t+21)*3**t
    b_num = 16*12**t+4*(t+1)*3**t
    r_num = 5*12**t+(2*t+5)*3**t
    need(all(x % 9 == 0 for x in (h_num, q_num, b_num, r_num)), 'integrality')
    return tuple(x//9 for x in (h_num, q_num, b_num, r_num))


def main():
    records = []
    for k in range(1, 8):
        original, original_edges = source(k)
        src = profile(original, original_edges, literal=k <= 2)
        need((src['H'], src['q'], src['s'], src['defect']) ==
             (34*(1 << k)+5, 9, 9, 4), 'source family')
        names = tuple(x for x in original if x != 'M')
        edges = {(v, w) for v, w in original_edges if v != 'M' and w != 'M'}
        inherited_names, inherited_edges = names, edges.copy()
        max_t = 3 if k == 1 else 2 if k <= 3 else 1
        for t in range(max_t+1):
            if t:
                names, edges = gadget(names, edges, t)
            got = profile(names, edges, literal=(k, t) in {(1, 0), (1, 1), (2, 0)})
            h, q, b, r = exact(k, t)
            need((got['H'], got['q'], got['s'], got['defect'], got['residual_H']) ==
                 (h, q, q, b, r), 'closed response')
            need(got['strong'], 'strong endpoints')
            need({(v, w) for v, w in edges if v in inherited_names and w in inherited_names}
                 == inherited_edges, 'persistent old induced graph')
            if t == 1:
                need(h == 3*src['H']+2, 'uniform ternary branch')
            records.append(dict(k=k, t=t, profile=got))
    # Independent scalar recurrences, not graph evidence for arbitrary t.
    for k in range(1, 31):
        for t in range(1, 51):
            vals = [exact(k, t+j)[0] for j in range(4)]
            need(vals[3] == 18*vals[2]-81*vals[1]+108*vals[0], 'cell recurrence')
            h, q, b, r = exact(k, t)
            need(b-2*r == 2*3**(t-1)*(4**t-1), 'exposure surplus identity')
    receipt = dict(status='PASS', method='Fresh named graph construction, endpoint subset recurrence, literal permutations',
                   graph_cases=len(records), source_scales=7, scalar_controls=1500,
                   no_construction_lane_imports=True, records=records)
    (HERE/'retained_reset.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print('PASS: seven source scales,', len(records), 'composed graphs, persistent arcs and full root responses')


if __name__ == '__main__':
    main()
