"""Exact controls for the fixed-valuation negative-Collatz chart."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def valuation(n, prime):
    need(n > 0, 'positive valuation input')
    result = 0
    while n % prime == 0:
        result += 1
        n //= prime
    return result


def phi(a, m):
    return (1 << a) * (2 * m - 1) - 1


def step(h):
    a = valuation(h + 1, 2)
    u = (h + 1) >> a
    return 3 * h + 2 if u % 4 == 1 else (h - 1 - (1 << a)) // 2


def count(rows):
    n = len(rows)
    full = (1 << n) - 1
    for u in range(n):
        need(not (rows[u] >> u & 1) and not rows[u] & ~full, 'row bounds')
        for v in range(u):
            need(((rows[u] >> v) ^ (rows[v] >> u)) & 1, 'tournament pair')
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        dp[1 << u][u] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if dp[mask][last]:
                available = rows[last] & ~mask
                while available:
                    bit = available & -available
                    available -= bit
                    dp[mask | bit][bit.bit_length() - 1] += dp[mask][last]
    return sum(dp[-1])


def ordinal(*blocks):
    total = sum(map(len, blocks))
    result = []
    offset = 0
    for rows in blocks:
        later = ((1 << total) - 1) ^ ((1 << (offset + len(rows))) - 1)
        result.extend((row << offset) | later for row in rows)
        offset += len(rows)
    return tuple(result)


def main():
    controls = 0
    for a in range(1, 13):
        for m in range(1, 4097):
            h = phi(a, m)
            cm = 3 * m - 1 if m % 2 else m // 2
            need(step(h) == phi(a, cm), 'exact conjugacy')
            need(valuation(step(h) + 1, 2) == a, 'preserved stratum')
            need((3 * m - 1) % 2 != m % 2, 'ternary parity flip')
            r = valuation(h + 1, 3)
            base = (h + 1) // 3 ** r - 1
            need(base > 0 and base % 2 == 1 and (base + 1) % 3, 'ternary depth')
            controls += 1
        need(step(phi(a, 1)) == phi(a, 2), 'cycle outward')
        need(step(phi(a, 2)) == phi(a, 1), 'cycle return')
    rows23 = (48, 1, 3, 7, 14, 30)
    triangle = (2, 4, 1)
    four = (2, 4, 9, 3)
    rows45 = ordinal(triangle, triangle, four)
    need(count(rows23) == 23 and count(four) == 5 and count(rows45) == 45,
         'literal positive source counts')
    need(step(23) == 7 and step(45) == 21, 'physical hostile edges')
    k = 17
    ray = [45 * 2 ** k + 7, 135 * 2 ** k + 23, 405 * 2 ** k + 71]
    coords = []
    for j, h in enumerate(ray):
        a = valuation(h + 1, 2)
        coords.append({'H': h, 'a': a, 'r': valuation(h + 1, 3),
                       'm': (((h + 1) >> a) + 1) // 2})
        need(a == 3 and coords[-1]['r'] == j, 'ray valuations')
    need(step(ray[0]) == ray[1] and step(ray[1]) == 135 * 2 ** 16 + 7,
         'ray guarded branches')
    for k in range(3, 60):
        h = 462 * 2 ** k - 5
        need(valuation(h + 1, 2) == 2 and ((h + 1) // 4) % 4 == 3,
             '462 even branch')
    receipt = {
        'status': 'FINITE-EXACT controls; all-input proof in NEGATIVE_COLLATZ_STRATA.md',
        'universe': {'a': [1, 12], 'm': [1, 4096], 'total': controls},
        'literal_sources': [{'H': 23, 'rows': rows23, 'guarded_target': 7},
                            {'H': 45, 'rows': rows45, 'guarded_target': 21}],
        'k17': coords, 'legal_halving_from_middle': step(ray[1]),
        'excluded_target_status': '7 and21 inherited all-order exclusions; not inferred by finite scan',
        'next_ready_target': 3 * ray[-1] + 2,
    }
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if '--write' in sys.argv:
        (HERE / 'negative_strata_receipt.json').write_text(encoded)
    elif (HERE / 'negative_strata_receipt.json').exists():
        need((HERE / 'negative_strata_receipt.json').read_text() == encoded,
             'saved receipt mismatch')
    print(encoded, end='')


if __name__ == '__main__':
    main()
