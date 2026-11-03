#!/usr/bin/env python3
"""Exact parity/pairing audit. Standard library only; no producer imports.

Universe: all ordinary tournament classes at orders 0..6, reached via every
fixed-Hamiltonian-spine mask; independently checked against every labelled
tournament mask. No filters except fixed forward spine for the presentation
carrier. Full canonical forms retain converse as a distinct class.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations
from math import factorial, gcd, prod
from pathlib import Path
import argparse
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def odd_partitions(n, least=1):
    if not n:
        yield ()
    for k in range(least, n + 1, 2):
        for tail in odd_partitions(n - k, k):
            yield (k,) + tail


def v2(n):
    require(n > 0, 'v2 requires a positive integer')
    return (n & -n).bit_length() - 1


def burnside(n):
    total = Fraction(0)
    valuations = []
    for part in odd_partitions(n):
        counts = Counter(part)
        exponent = sum((k - 1) // 2 for k in part)
        exponent += sum(gcd(a, b) for i, a in enumerate(part) for b in part[i + 1:])
        denominator = prod(k ** m * factorial(m) for k, m in counts.items())
        total += Fraction(2 ** exponent, denominator)
        valuations.append(exponent - v2(denominator))
    require(total.denominator == 1, f'Burnside integrality n={n}')
    return int(total), min(valuations)


def trail_pairs(fibres, flip):
    """Leave least word in each odd fibre; pair remaining consecutive words."""
    local = {}
    representatives = {}
    for c, words in fibres.items():
        words = sorted(words)
        require(len(words) % 2 == 1, 'fibre must be odd')
        representatives[words[0]] = c
        for a, b in zip(words[1::2], words[2::2]):
            local[a], local[b] = b, a
    pairing = {}
    paths = []
    for start, c in sorted(representatives.items()):
        if c in pairing:
            continue
        word = start
        path = [word]
        seen = {word}
        while True:
            word = flip(word)
            require(word not in seen, 'representative path cannot cycle')
            seen.add(word)
            path.append(word)
            if word in representatives:
                d = representatives[word]
                require(c != d and d not in pairing, 'distinct unpaired endpoints')
                pairing[c], pairing[d] = d, c
                paths.append({'classes': [c, d], 'masks': path})
                break
            word = local[word]
            require(word not in seen, 'local continuation cannot cycle')
            seen.add(word)
            path.append(word)
    require(len(pairing) == len(fibres), 'all classes paired')
    return pairing, paths


def census(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    index = {e: k for k, e in enumerate(edges)}
    chords = [e for e in edges if e[1] != e[0] + 1]
    chord_index = {e: k for k, e in enumerate(chords)}
    q = len(chords)
    full = (1 << len(edges)) - 1
    words = tuple(permutations(range(n)))

    def backward(mask, i, j):
        return ((mask >> index[min(i, j), max(i, j)]) & 1) ^ (i > j)

    def transport(mask, order):
        return sum(backward(mask, order[i], order[j]) << k for k, (i, j) in enumerate(edges))

    def from_spine(mask):
        return sum(((mask >> k) & 1) << index[e] for k, e in enumerate(chords))

    def mirror(mask):
        return sum(((mask >> k) & 1) << chord_index[n - 1 - j, n - 1 - i]
                   for k, (i, j) in enumerate(chords))

    canonical = {}
    aut = {}
    for mask in range(1 << q):
        graph = from_spine(mask)
        if graph in canonical:
            continue
        orbit = {transport(graph, order) for order in words}
        c = min(orbit)
        canonical.update((x, c) for x in orbit)
        aut[c] = factorial(n) // len(orbit)
    require(set(canonical) == set(range(1 << len(edges))), 'all labelled masks reached')

    fibres, blue = defaultdict(list), defaultdict(list)
    for mask in range(1 << q):
        c = canonical[from_spine(mask)]
        fibres[c].append(mask)
        if mirror(mask) == mask:
            blue[c].append(mask)
        require(canonical[from_spine(mirror(mask))] == canonical[full ^ from_spine(mask)],
                'mirror realizes converse')

    h_counts = {}
    for c in fibres:
        h = sum(all(backward(c, p[i], p[i + 1]) == 0 for i in range(n - 1)) for p in words)
        dp = {(1 << i, i): 1 for i in range(n)}
        for subset in range(1, 1 << n):
            for last in range(n):
                count = dp.get((subset, last), 0)
                if count:
                    for nxt in range(n):
                        if not (subset >> nxt) & 1 and not backward(c, last, nxt):
                            key = subset | (1 << nxt), nxt
                            dp[key] = dp.get(key, 0) + count
        h_dp = sum(dp.get(((1 << n) - 1, i), 0) for i in range(n)) if n else 1
        require(h == h_dp, 'permutations agree with subset DP')
        require(h % 2 == aut[c] % 2 == len(fibres[c]) % 2 == 1, 'odd H, aut, fibre')
        require(h == aut[c] * len(fibres[c]), 'free-fibre identity')
        h_counts[c] = h
        self_converse = canonical[full ^ c] == c
        require((c in blue) == self_converse, 'blue iff self-converse')
        if self_converse:
            require(len(blue[c]) % 2 == 1, 'blue fibre odd')

    result = {'n': n, 'classes': len(fibres), 'self_converse': len(blue),
              'presentation_count': 1 << q,
              'F_n': dict(sorted(Counter(h_counts.values()).items())),
              'SC_H': dict(sorted(Counter(h_counts[c] for c in blue).items()))}
    require(result['classes'] == burnside(n)[0], 'independent Burnside count')
    require(result['F_n'][1] == 1, 'unique H=1 class')
    if n >= 3:
        flip = lambda mask: ((1 << q) - 1) ^ mask
        _, all_paths = trail_pairs(fibres, flip)
        blue_pairs, blue_paths = trail_pairs(blue, flip)
        refined_pairs = dict(blue_pairs)
        for c in fibres:
            if c not in blue:
                refined_pairs[c] = canonical[full ^ c]
        require(len(refined_pairs) == len(fibres), 'complete refined pairing')
        require(all(d != c and refined_pairs[d] == c for c, d in refined_pairs.items()),
                'refined pairing is a free involution')
        require(all((c in blue) == (d in blue) for c, d in refined_pairs.items()),
                'refined pairing preserves SC predicate')
        require(all(h_counts[c] == h_counts[d] for c, d in refined_pairs.items() if c not in blue),
                'non-SC pairing preserves H')
        result['full_trail_max_vertices'] = max(len(p['masks']) for p in all_paths)
        result['blue_trails'] = [dict(p, H=[h_counts[c] for c in p['classes']]) for p in blue_paths]
        result['refined_pairs'] = [{'classes': [c, d], 'H': [h_counts[c], h_counts[d]],
                                    'self_converse': c in blue}
                                   for c, d in sorted(refined_pairs.items()) if c < d]
        result['blue_fibres'] = {str(c): blue[c] for c in sorted(blue)}
    return result


def diagram(row, path):
    """Native SVG: exact order-five blue trails, not a schematic class mirror."""
    rows = row['blue_trails']
    width, height = 1080, 160 + 130 * len(rows)
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#f8fafc"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#10223b}.title{font-size:23px;font-weight:bold}.small{font-size:14px}.label{font-size:17px}</style>',
             '<text class="title" x="32" y="38">Order 5: pairing the eight self-converse tournament classes</text>',
             '<text class="small" x="32" y="66">Blue lines: complement every nonspine arc. Orange arcs: pair words inside one class fibre.</text>',
             '<text class="small" x="32" y="90">One least word is left in each odd blue fibre; alternating trails pair the class endpoints.</text>']
    for r, item in enumerate(rows):
        y = 160 + 130 * r
        masks = item['masks']
        step = min(130, 900 / max(1, len(masks) - 1))
        start = (width - step * (len(masks) - 1)) / 2
        for k in range(len(masks) - 1):
            x, x2 = start + step * k, start + step * (k + 1)
            if k % 2 == 0:
                lines.append(f'<line x1="{x}" y1="{y}" x2="{x2}" y2="{y}" stroke="#2979bd" stroke-width="4"/>')
                lines.append(f'<text class="small" x="{(x+x2)/2}" y="{y+28}" text-anchor="middle">τ</text>')
            else:
                lines.append(f'<path d="M{x},{y} Q{(x+x2)/2},{y-62} {x2},{y}" fill="none" stroke="#cb761c" stroke-width="3"/>')
                lines.append(f'<text class="small" x="{(x+x2)/2}" y="{y-39}" text-anchor="middle">same class</text>')
        for k, mask in enumerate(masks):
            x = start + step * k
            endpoint = k in (0, len(masks)-1)
            lines.append(f'<circle cx="{x}" cy="{y}" r="{17 if endpoint else 13}" fill="{"#143d59" if endpoint else "#f8fafc"}" stroke="#143d59" stroke-width="2"/>')
            lines.append(f'<text x="{x}" y="{y+5}" text-anchor="middle" style="font-family:monospace;font-size:14px;fill:{"white" if endpoint else "#143d59"}">{mask}</text>')
        lines.append(f'<text class="label" x="32" y="{y-27}">H = {item["H"][0]} ↔ {item["H"][1]}</text>')
    lines.append(f'<text class="small" x="32" y="{height-28}">The four non-self-converse classes form two converse pairs. Vertex masks are normalized spine words (0–63).</text>')
    lines.append('</svg>')
    path.write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--svg', type=Path)
    args = parser.parse_args()
    rows = [census(n) for n in range(7)]
    divisibility = []
    for n in range(41):
        a, minimum = burnside(n)
        if n >= 3:
            require(minimum >= n // 2, 'termwise stronger divisibility')
            require(a % 2 == 0, 'all n>=3 even')
            if n % 2:
                require(v2(a) == (n - 1) // 2, 'odd-order exact valuation')
            else:
                phi = sum(gcd(k, n) == 1 for k in range(1, n))
                require((a >> (n // 2)) % 2 == (phi // 2) % 2, 'even-order first residue')
        divisibility.append({'n': n, 'a_n': a, 'v2': v2(a), 'minimum_term_v2': minimum})
    result = {'status': 'FINITE-EXACT', 'universe': 'all ordinary classes n=0..6; Burnside partitions n=0..40',
              'checks': 'explicit require, active under python -O', 'census': rows, 'burnside': divisibility}
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end='')
    if args.svg:
        diagram(rows[5], args.svg)


if __name__ == '__main__':
    main()
