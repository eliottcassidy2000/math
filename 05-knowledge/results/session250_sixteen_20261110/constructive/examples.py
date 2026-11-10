#!/usr/bin/env python3
"""Proof-derived exact examples; no external packages or census claims."""
import json
import random
from pathlib import Path


def tournament(n, code):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            a, b = (i, j) if (code >> k) & 1 else (j, i)
            out[a] |= 1 << b
            k += 1
    return out


def join(a, b):
    n, m = len(a), len(b)
    return [x | (((1 << m) - 1) << n) for x in a] + [x << n for x in b]


def signature(out, selected):
    n = len(out)
    size = 1 << n
    # Literal subset endpoint DP, independently of the gluing formulas.
    ends = [[0] * n for _ in range(size)]
    starts = [[0] * n for _ in range(size)]
    for v in range(n):
        ends[1 << v][v] = starts[1 << v][v] = 1
    for mask in range(1, size):
        for v in range(n):
            if not mask & (1 << v):
                continue
            rest = mask ^ (1 << v)
            if rest:
                ends[mask][v] = sum(ends[rest][u] for u in range(n)
                                    if rest & (1 << u) and out[u] & (1 << v))
                starts[mask][v] = sum(starts[rest][u] for u in range(n)
                                      if rest & (1 << u) and out[v] & (1 << u))
    h = [sum(x) for x in ends]
    h[0] = 1
    s = [sum(row[v] for v in range(n) if selected & (1 << v)) for row in starts]
    e = [sum(row[v] for v in range(n) if selected & (1 << v)) for row in ends]
    full = size - 1
    q = sum(h[x] * h[full ^ x] for x in range(size))
    m = sum(h[x] * s[full ^ x] for x in range(size))
    nn = sum(e[x] * h[full ^ x] for x in range(size))
    t = sum(e[x] * s[full ^ x] for x in range(size))
    return (h[full], q, s[full], e[full], m, nn, t)


def glue(a, b):
    d, q, s, e, m, n, t = a
    D, Q, S, E, M, N, T = b
    return (d*D, q*Q, s*D, d*E, m*Q+d*M,
            n*D+q*N, t*D+d*T+e*S+m*N)


def root(out, selected):
    n = len(out)
    return [x | (0 if selected & (1 << i) else 1 << n)
            for i, x in enumerate(out)] + [selected]


def tail(out, selected, t):
    # Root beats precisely selected old vertices; all old deletion vertices
    # beat the new tail, and every tail vertex beats the root.
    n = len(out)
    r = root(out, selected)
    r[:n] = [x | (((1 << t)-1) << (n+1)) for x in r[:n]]
    for i in range(t):
        r.append((1 << n) | sum(1 << (n+1+j) for j in range(i+1, t)))
    return r


def main():
    rng = random.Random(250)
    count = 0
    for _ in range(80):
        n, k = rng.randrange(1, 5), rng.randrange(1, 5)
        a = tournament(n, rng.getrandbits(n*(n-1)//2))
        b = tournament(k, rng.getrandbits(k*(k-1)//2))
        x, y = rng.randrange(1 << n), rng.randrange(1 << k)
        sa, sb = signature(a, x), signature(b, y)
        assert sa[4] >= x.bit_count()*sa[0]
        assert sa[5] >= x.bit_count()*sa[0]
        assert glue(sa, sb) == signature(join(a, b), x | (y << n))
        count += 1
        alpha = sa[0]*sb[0]-sa[6]*sb[0]-sa[0]*sb[6]
        beta = sa[4]*sb[5]
        empty = glue(sa, sb)
        if x and y:
            assert empty[0]-empty[6] <= 0
        assert empty[0]-empty[6] == alpha-beta-sa[3]*sb[2]
        for r in (1, 2):
            p = tournament(r, (1 << (r*(r-1)//2))-1)
            actual = signature(join(join(a, p), b), x | (y << (n+r)))
            predicted = glue(glue(sa, (1, 2**r, 0, 0, 0, 0, 0)), sb)
            assert actual == predicted
            assert actual[0]-actual[6] == alpha-beta*2**r
            assert actual[4] == sa[4]*sb[1]*2**r+sa[0]*sb[4]
    # An inherited small unit seed supplies an actual rooted-ray control.
    cycle = [1 << 1, 1 << 2, 1 << 0]
    seed = signature(cycle, 3)
    assert seed[0]-seed[6] == 1
    prefix = tournament(2, 1)
    suffix = tournament(3, 7)
    padded = signature(join(join(prefix, cycle), suffix), 3 << 2)
    assert padded[0]-padded[6] == 1
    assert padded[4] == 8*seed[4] == 48
    assert signature([0], 1) == (1, 2, 1, 1, 1, 1, 0)
    ray = []
    for t in (1, 2, 3):
        h = signature(tail(cycle, 3, t), 0)[0]
        assert h == 1 + seed[4]*2**t
        ray.append(h)
    # Exhaustive small factor bank, only to exhibit whether this particular
    # seam construction already produces a positive-gap seed in this box.
    bank = []
    for n in range(1, 5):
        for code in range(1 << (n*(n-1)//2)):
            a = tournament(n, code)
            for x in range(1, 1 << n):
                bank.append((n, code, x, signature(a, x)))
    seen = {}
    for entry in bank:
        seen.setdefault(entry[3], entry)
    hit = None
    for aa in seen.values():
        for bb in seen.values():
            a, b = aa[3], bb[3]
            alpha = a[0]*b[0]-a[6]*b[0]-a[0]*b[6]
            beta = a[4]*b[5]
            numerator = alpha-1
            if numerator > 0 and numerator % beta == 0:
                z = numerator // beta
                if z >= 2 and z & (z-1) == 0:
                    hit = {'A': aa[:3], 'C': bb[:3], 'gap': z.bit_length()-1}
                    break
        if hit:
            break
    result = {'status': 'all assertions passed', 'random_join_examples': count,
              'seam_examples': 2*count, 'unit_seed_signature': seed,
              'inherited_unit_seed_ray': ray, 'small_bank_signature_count': len(seen),
              'prefix_2_suffix_3_padded_seed_signature': padded,
              'padding_core_and_padded_odd_slope': [3, 3],
              'positive_gap_seed_in_factors_orders_at_most_4': hit,
              'scope': 'finite examples, not spectrum coverage or a canon audit'}
    print(json.dumps(result, indent=2))
    Path(__file__).with_name('results.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
