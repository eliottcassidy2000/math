#!/usr/bin/env python3
"""Physical compiler with an explicitly incomplete profile-certificate oracle.
No missing certificate is inferred from the open all-order lemma.
"""
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise ValueError(message)


def valid(a):
    n = len(a)
    require(n > 0, 'empty output')
    for i, row in enumerate(a):
        require(0 <= row < 1 << n and not (row >> i & 1), 'invalid row')
        for j in range(i):
            require(((row >> j & 1) + (a[j] >> i & 1)) == 1, 'invalid pair')


def count(a):
    """Exact vertex-subset path DP, used only on small certificate graphs."""
    n = len(a)
    require(n <= 20, 'literal DP restricted to at most20 vertices')
    d = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        d[1 << v][v] = 1
    for used in range(1, 1 << n):
        for v, value in enumerate(d[used]):
            if not value:
                continue
            todo = a[v] & ~used
            while todo:
                bit = todo & -todo
                todo -= bit
                d[used | bit][bit.bit_length() - 1] += value
    return sum(d[-1])


def profiles(a, root):
    """Return actual(h,s,L) by subset paths and complementary covers."""
    n = len(a)
    d = [[0] * n for _ in range(1 << n)]
    starts = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for v in range(n):
        d[1 << v][v] = starts[1 << v][v] = 1
    for used in range(1, 1 << n):
        for v in range(n):
            if not (used >> v & 1):
                continue
            rest = used ^ (1 << v)
            if rest:
                for w in range(n):
                    if rest >> w & 1:
                        if a[w] >> v & 1:
                            d[used][v] += d[rest][w]
                        else:
                            starts[used][v] += starts[rest][w]
            h[used] += d[used][v]
    full = (1 << n) - 1
    L = sum(d[used][root] * h[full ^ used] for used in range(1, full + 1))
    return h[-1], starts[-1][root], L


def read_data(directory):
    seeds = {}
    for line in (directory / 'SEEDS.tsv').read_text().splitlines():
        target, n, *a = map(int, line.split())
        require(len(a) == n and target not in seeds, 'invalid seed record')
        seeds[target] = a
    certificates = {}
    for line in (directory / 'PROFILE_CERTIFICATES.tsv').read_text().splitlines():
        B, residue, j, root, n, *a = map(int, line.split())
        require(len(a) == n and (B, residue) not in certificates, 'invalid profile record')
        certificates[B, residue] = (j, root, a)
    return seeds, certificates


def child(a, root):
    """Copy U into core(2,8,3,5), identifying its root with core0."""
    order = [root] + [i for i in range(len(a)) if i != root]
    u = [sum(1 << j for j, w in enumerate(order) if a[v] >> w & 1) for v in order]
    out = [2, 8, 3, 5] + [0] * (len(u) - 1)
    for j in range(1, len(u)):
        x = j + 3
        for i in range(4):
            if i or u[0] >> j & 1:
                out[i] |= 1 << x
            else:
                out[x] |= 1 << i
        for k in range(1, j):
            y = k + 3
            if u[j] >> k & 1:
                out[x] |= 1 << y
            else:
                out[y] |= 1 << x
    valid(out)
    return out


def join_triangle(a):
    n = len(a)
    suffix = ((1 << n) - 1) << 3
    return [2 | suffix, 4 | suffix, 1 | suffix] + [row << 3 for row in a]


def compile_target(target, seeds, certificates):
    require(target >= 23 and target % 2 == 1, 'target must be odd and at least23')
    if target in seeds:
        a = seeds[target]
        valid(a)
        require(count(a) == target, 'seed path count failed')
        return list(a), [{'kind': 'literal_seed', 'target': target}]
    if target % 3 == 0:
        a, trace = compile_target(target // 3, seeds, certificates)
        out = join_triangle(a)
        valid(out)
        return out, trace + [{'kind': 'ordered_triangle_join', 'target': target}]
    residue = 1 if target % 3 == 2 else 2
    B = (target - 2 * residue) // 3
    key = (B, residue)
    if key not in certificates:
        raise ValueError(f'MISSING_OPEN_PROFILE_CERTIFICATE: B={B}, residue={residue}; '
                         'the all-order profile lemma is OPEN')
    j, root, u = certificates[key]
    valid(u)
    require(0 <= j <= 8 and 0 <= root < len(u), 'profile parameter outside scope')
    h, s, L = profiles(u, root)
    require(h + L - s == B - 2 * j and s == residue + 3 * j, 'profile equations failed')
    out = child(u, root)
    require(count(out) == target, 'literal compiled output count failed')
    return out, [{'kind': 'literal_profile_child', 'target': target,
                  'B': B, 'residue': residue, 'j': j,
                  'source_profile': {'h': h, 's': s, 'L': L}}]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('target', type=int)
    p.add_argument('--data-dir', type=Path, default=HERE)
    args = p.parse_args()
    seeds, certificates = read_data(args.data_dir)
    try:
        rows, trace = compile_target(args.target, seeds, certificates)
    except ValueError as exc:
        print(json.dumps({'status': 'NOT_COMPILED', 'target': args.target, 'reason': str(exc)}))
        raise SystemExit(2)
    print(json.dumps({'status': 'CERTIFIED_PHYSICAL_WITNESS', 'target': args.target,
                      'order': len(rows), 'rows': rows, 'trace': trace,
                      'all_order_profile_lemma': 'OPEN'}, indent=2))


if __name__ == '__main__':
    main()
