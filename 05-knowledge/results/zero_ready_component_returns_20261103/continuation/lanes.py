"""Exact small-integer lane classification of the selected named tree."""
from pathlib import Path
from math import isqrt, prod
import hashlib
import json

HERE = Path(__file__).resolve().parent


def factors(n):
    result = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            power = 0
            while n % d == 0:
                n //= d
                power += 1
            result.append([d, power])
        d += 1
    if n > 1:
        result.append([n, 1])
    return result


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def main():
    source = HERE / 'named_tree.json'
    tree = json.loads(source.read_text())
    records = []
    for item in [tree['source']] + tree['targets']:
        h = item['H']
        fac = factors(h)
        if prod(p**e for p, e in fac) != h or not all(prime(p) for p, e in fac):
            raise RuntimeError('Factorization or trial-division certificate failed')
        lane = 'prime' if prime(h) else '7p' if h % 7 == 0 and prime(h // 7) else 'auxiliary composite'
        records.append(dict(H=h, branch=item.get('branch', ''), factors=fac, lane=lane))
    print(json.dumps(dict(status='FINITE-EXACT',
        scope='Only the selected source and fourteen targets; no coverage inference',
        primality_method='Trial division by every integer through floor(sqrt(n))',
        records=records, source_sha256={
            str(p.relative_to(HERE.parent.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (source, Path(__file__).resolve())}), indent=2))


if __name__ == '__main__':
    main()
