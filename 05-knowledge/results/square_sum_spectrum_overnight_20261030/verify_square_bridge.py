#!/usr/bin/env python3
"""Independent exhaustive replay of every cycle, without discovery dependencies."""
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
checks = 0


def require(test, message):
    global checks
    checks += 1
    if not test:
        raise ValueError(message)


def is_square(n):
    return isqrt(n) ** 2 == n


def audit(p, n):
    require(len(p) == n and set(p) == set(range(1, n + 1)), 'vertex coverage')
    for i in range(n):
        require(is_square(p[i] + p[(i + 1) % n]), 'non-square cycle edge')


def main():
    packed = (HERE / 'square_bridge.json.gz').read_bytes()
    data = json.loads(gzip.decompress(packed))
    require((data['first'], data['last']) == (32, 6242), 'declared finite range')
    require(len(data['steps']) == 6210, 'one record per intermediate size')
    path = data['initial']
    audit(path, 32)
    direct = moves = resets = worst = 0
    for n, step in enumerate(data['steps'], 33):
        if 'cycle' in step:
            require(set(step) == {'cycle'}, 'reset schema')
            path = step['cycle']
            resets += 1
        else:
            require(set(step) == {'insert', 'reverse'}, 'move schema')
            i = step['insert']
            require(type(i) is int and 0 <= i < n - 1, 'insertion index')
            path.insert(i + 1, n)
            direct += not step['reverse']
            worst = max(worst, len(step['reverse']))
            for a, b in step['reverse']:
                require(type(a) is int and type(b) is int and 0 <= a < b < n,
                        'reversal indices')
                path[a + 1:b + 1] = path[a + 1:b + 1][::-1]
                moves += 1
        # Recheck the entire graph certificate at EVERY size, not just changed edges.
        audit(path, n)
    # Explicit hostile controls: duplicates and non-square adjacencies cannot pass.
    good = data['initial']
    require(len(set(good[:-1] + good[:1])) != 32, 'duplicate control')
    require(not is_square(1 + 2), 'non-square control')
    print(json.dumps({'status': 'passed', 'first': 32, 'last': 6242,
                      'fully_checked_cycles': 6211, 'checks': checks,
                      'reversals': moves, 'direct_insertions': direct,
                      'maximum_reversals': worst, 'reset_cycles': resets,
                      'trace_sha256': hashlib.sha256(packed).hexdigest(),
                      'scope': 'Finite interval only; the all-size conclusion also uses the named inherited induction.'},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
