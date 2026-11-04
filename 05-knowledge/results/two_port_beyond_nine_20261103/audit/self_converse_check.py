"""Literal controls for even self-converse endpoint separation."""
from argparse import ArgumentParser
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).read_bytes()
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def compatible_family(n):
    involution = tuple(u ^ 1 for u in range(n))
    seen, orbits = set(), []
    for u, v in combinations(range(n), 2):
        if (u, v) in seen:
            continue
        other = (involution[v], involution[u])
        orbits.append(((u, v), other))
        seen.add((u, v))
        seen.add(tuple(sorted(other)))
    require(len(orbits) == (n // 2) ** 2, 'free anti-involution edge orbits')
    masks = range(1 << len(orbits)) if n <= 6 else list(range(8)) + [65535 - k for k in range(8)]
    for mask in masks:
        rows = [0] * n
        for bit, orbit in enumerate(orbits):
            for u, v in orbit:
                if mask >> bit & 1:
                    rows[u] |= 1 << v
                else:
                    rows[v] |= 1 << u
        yield rows, involution


def literal_endpoints(rows):
    n = len(rows)
    s, e = [0] * n, [0] * n
    for word in permutations(range(n)):
        if all(rows[u] >> v & 1 for u, v in zip(word, word[1:])):
            s[word[0]] += 1
            e[word[-1]] += 1
    return s, e


def main():
    parser = ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    counts = {}
    for n in (2, 4, 6, 8):
        counts[n] = 0
        for rows, involution in compatible_family(n):
            require(all(((rows[u] >> v) & 1) != ((rows[involution[u]] >> involution[v]) & 1)
                        for u, v in combinations(range(n), 2)), 'actual anti-automorphism')
            s, e = literal_endpoints(rows)
            require(all(s[u] == e[involution[u]] for u in range(n)), 'literal endpoint transfer')
            require(sum(s) % 2 == 1 and any((a-b) % 2 for a, b in zip(s, e)), 'even converse visible modulo two')
            counts[n] += 1
    parity = json.loads((HERE / 'parity.json').read_text())
    for record in parity['even_parity_hostiles']:
        rows = record['rows']
        require(not any(all(((rows[u] >> v) & 1) != ((rows[p[u]] >> p[v]) & 1)
                            for u, v in combinations(range(8), 2))
                        for p in permutations(range(8))), 'parity hostile is not self-converse')
    require(Path(__file__).read_bytes() == SOURCE, 'source stable')
    result = {'status': 'FINITE-EXACT PASS', 'checks': CHECKS,
              'source_sha256': sha256(SOURCE).hexdigest(),
              'literal_controls': counts,
              'order8_masks': list(range(8)) + [65535-k for k in range(8)],
              'scope': 'All fixed-pair anti-involution tournaments through6; sixteen specified order8 cases; both parity hostiles tested against every8! permutation.'}
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(text, encoding='utf8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
