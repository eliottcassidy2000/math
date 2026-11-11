#!/usr/bin/env python3
"""Exact dyadic targeting in Moon interval chains. Standard library only."""
import argparse, itertools, json, random

def require(condition, message):
    if not condition:
        raise RuntimeError("Exact check failed: " + message)

def banks(k):
    require(k >= 3, 'k >= 3')
    return [([5], [2, 3]), ([7], [3, 4])] + [([2 * j], [j, j]) for j in range(3, k)]

def weight(parts):
    h = 1
    for p in parts:
        h *= 1 + (1 << p - 1)
    return h

def target(k, r):
    require(r % 2 and k >= 3, 'r % 2 and k >= 3')
    bs = banks(k)
    parts = [a[:] for a, b in bs]
    h = weight(sum(parts, []))
    bits = []
    for j, (a, b) in enumerate(bs, 1):
        bit = h - r >> j & 1
        if bit:
            h = h // weight(a) * weight(b)
            parts[j - 1] = b[:]
        bits.append(bit)
        require((h - r) % (1 << j + 1) == 0, '(h - r) % (1 << j + 1) == 0')
    return {'k': k, 'residue': r % (1 << k), 'bits': bits, 'parts': sum(parts, []), 'order': k * k - k + 7, 'H': h}

def graph(parts):
    n = 1 + sum(parts)
    rows = [sum((1 << j for j in range(i + 1, n))) for i in range(n)]
    start = 0
    for p in parts:
        end = start + p
        rows[start] ^= 1 << end
        rows[end] ^= 1 << start
        start = end
    return rows

def dp(rows):
    n = len(rows)
    d = [{} for _ in range(1 << n)]
    for i in range(n):
        d[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i, v in d[mask].items():
            available = rows[i] & ~mask
            while available:
                bit = available & -available
                available -= bit
                j = bit.bit_length() - 1
                nxt = d[mask | bit]
                nxt[j] = nxt.get(j, 0) + v
    return sum(d[-1].values())

def strong(rows):
    n = len(rows)
    for start in range(n):
        reached = 1 << start
        front = reached
        while front:
            bit = front & -front
            front -= bit
            new = rows[bit.bit_length() - 1] & ~reached
            reached |= new
            front |= new
        if reached != (1 << n) - 1:
            return False
    return True

def run():
    out = {'status': 'PASS', 'scope': 'exact all-target residues through k=15; literal subset DP on all bank blocks through j=7 and 40 small interval chains; fixed-family k=3 graphs', 'full_residue_tests': [], 'literal_tests': 0}
    for k in range(3, 16):
        hs = []
        for r in range(1, 1 << k, 2):
            w = target(k, r)
            h = w['H']
            hs.append(h % (1 << k))
            require(sum(w['parts']) + 1 == w['order'], "sum(w['parts']) + 1 == w['order']")
            require(h % 7, 'h % 7')
        require(len(set(hs)) == 1 << k - 1, 'len(set(hs)) == 1 << k - 1')
        out['full_residue_tests'].append({'k': k, 'targets': len(hs), 'order': k * k - k + 7})
    cases = []
    for j, pair in enumerate(banks(8), 1):
        a, b = pair
        delta = weight(a) - weight(b)
        require(abs(delta) & -abs(delta) == 1 << j, 'abs(delta) & -abs(delta) == 1 << j')
        cases.extend([a, b])
    rng = random.Random(269)
    for _ in range(40):
        parts = [rng.randrange(2, 6) for _ in range(rng.randrange(1, 4))]
        cases.append(parts)
    cases += [target(3, r)['parts'] for r in range(1, 8, 2)]
    for parts in cases:
        rows = graph(parts)
        n = len(rows)
        require(dp(rows) == weight(parts), 'dp(rows) == weight(parts)')
        require(strong(rows), 'strong(rows)')
        scores = sorted((x.bit_count() for x in rows))
        require(scores == [1] + list(range(1, n - 1)) + [n - 2], 'scores == [1] + list(range(1, n - 1)) + [n - 2]')
        triangles = sum((1 for a, b, c in itertools.combinations(range(n), 3) if rows[a] >> b & 1 and rows[b] >> c & 1 and rows[c] >> a & 1 or (rows[b] >> a & 1 and rows[c] >> b & 1 and rows[a] >> c & 1)))
        require(triangles == n - 2, 'triangles == n - 2')
        out['literal_tests'] += 1
    out['target_seven_examples'] = [target(k, 7) for k in [3, 4, 5, 8, 12, 32, 64]]
    out['target_sixty_three_examples'] = [target(k, 63) for k in [6, 8, 12, 32, 64]]
    for k in [32, 64, 128, 256]:
        for _ in range(20):
            r = rng.getrandbits(k) | 1
            w = target(k, r)
            require((w['H'] - r) % (1 << k) == 0, "(w['H'] - r) % (1 << k) == 0")
    out['large_precision_targets'] = 80
    return out
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--target', type=int)
    p.add_argument('--bits', type=int, default=8)
    p.add_argument('--graph', action='store_true')
    a = p.parse_args()
    if a.target is not None:
        result = target(a.bits, a.target)
        if a.graph:
            result['outneighbor_bitrows'] = graph(result['parts'])
    else:
        result = run()
    print(json.dumps(result, indent=2))
