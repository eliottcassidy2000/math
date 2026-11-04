"""Independent literal-permutation checks for the incoming integration note.

No producer imports or coefficient banks. Every queried extended tournament
is built from its physical arcs; H counts permutations of all its vertices.
"""
import hashlib
import itertools
import json
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def reverse(T, selected):
    R = list(T)
    for u, v in itertools.combinations(range(len(T)), 2):
        if sum(bool(M & (1 << u) and M & (1 << v)) for M in selected) % 2:
            R[u] ^= 1 << v
            R[v] ^= 1 << u
    return tuple(R)

def module(T, M):
    return all(T[v] & M in (0, M) for v in range(len(T)) if not M >> v & 1)

@lru_cache(None)
def value(T, a, c):
    n = len(T)
    R = list(T) + [1 << (n + 1), 0]
    for v in range(n):
        for w, mask in ((n, a), (n + 1, c)):
            if mask >> v & 1:
                R[v] |= 1 << w
            else:
                R[w] |= 1 << v
    need(all(bool(R[u] >> v & 1) != bool(R[v] >> u & 1)
             for u, v in itertools.combinations(range(n + 2), 2)), 'tournament arcs')
    return sum(all(R[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in itertools.permutations(range(n + 2)))

def run():
    pinned = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest()
              for p in ('INTEGRATION.md', 'check.py', 'receipt.json', 'receipt.optimized.json')}
    expected = json.loads((HERE / 'receipt.json').read_text())
    need((HERE / 'receipt.json').read_bytes() == (HERE / 'receipt.optimized.json').read_bytes(), 'producer modes')
    T = (14, 12, 56, 48, 35, 3)
    Ms = (3, 12, 48)
    need(all(module(T, M) for M in Ms), 'source modules')
    states = [reverse(T, [M for i, M in enumerate(Ms) if bits >> i & 1]) for bits in range(8)]
    decoders = []
    for bits, R in enumerate(states):
        for i, (a, c) in enumerate(((48, 49), (3, 7), (12, 28))):
            F = value(R, a, c)
            need(F == (483 if bits >> i & 1 else 393), 'fixed-coordinate decoder')
            decoders.append([bits, i, a, c, F])
    need(decoders == expected['decoder_queries'], 'all decoder coordinates match')
    mon = 5 | (18 << 6)
    total = 0
    terms = []
    s = mon
    while True:
        a, c = s & 63, s >> 6
        counts = [value(states[k], a, c) for k in range(4)]
        sign = (-1) ** (4 - s.bit_count())
        total += sign * (counts[0] - counts[1] - counts[2] + counts[3])
        terms.append({'a': a, 'c': c, 'four_counts': counts, 'sign': sign})
        if not s:
            break
        s = (s - 1) & mon
    need(total == -2, 'quartic interaction')
    need(terms == expected['mixed_quartic']['sixteen_terms'], 'all finite-difference data match')
    small = []
    for M, N, a, target in ((3, 12, 5, [19, 11, 15, 9]), (3, 7, 9, [15, 9, 5, 9])):
        R = (0, 1, 3, 7)
        need(module(R, M) and module(R, N), 'small original modules')
        counts = [value(reverse(R, ops), a, 0) for ops in ([], [M], [N], [M, N])]
        need(counts == target, 'literal small interaction')
        small.append({'T': list(R), 'modules': [M, N], 'query': [a, 0], 'four_counts': counts,
                      'mixed_difference': counts[0] - counts[1] - counts[2] + counts[3]})
    need(small == expected['small_nonadditivity'], 'all small data match')
    U = reverse(T, [63, *Ms])
    need(U == (50, 48, 11, 3, 44, 12), 'covered source')
    need([value(T, 1, 5), value(U, 1, 5)] == [251, 249], 'covered nonuniform separator')
    uniform = []
    for a_bits, c_bits in itertools.product(range(8), repeat=2):
        a = sum(M for i, M in enumerate(Ms) if a_bits >> i & 1)
        c = sum(M for i, M in enumerate(Ms) if c_bits >> i & 1)
        x, y = value(T, a, c), value(U, a, c)
        need((x - y) % 4 == 0, 'covered uniform mod-four hostile')
        uniform.append([a, c, x, y])
    need(sum(x != y for a, c, x, y in uniform) == 24, 'covered integer mismatch count')
    need(pinned == {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in pinned}, 'frozen sources')
    return {'status': 'PASS', 'source_sha256': pinned,
            'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'method': 'Literal permutations for every physical port-extended tournament. No producer imports.',
            'distinct_literal_counts': value.cache_info().currsize,
            'decoder_queries': decoders, 'quartic_coefficient': total,
            'finite_difference_terms': terms, 'small_nonadditivity': small,
            'covered_uniform_queries': uniform,
            'scope': 'Selected exact incoming controls, not an unrestricted collision theorem or a replay of its entire full-response census.'}

if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
