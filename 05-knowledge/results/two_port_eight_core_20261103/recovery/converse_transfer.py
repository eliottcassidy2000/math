"""Independent transfer of incoming converse identities to full named twelve-core data."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PACKET = HERE.parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def subsets(mask):
    x = mask
    while True:
        yield x
        if not x:
            break
        x = (x - 1) & mask


def path_parities(rows):
    n = len(rows)
    starts = [0] * (1 << n)
    for mask in range(1, 1 << n):
        for u in range(n):
            if mask >> u & 1:
                rest = mask ^ (1 << u)
                if not rest or (rows[u] & starts[rest]).bit_count() % 2:
                    starts[mask] |= 1 << u
    incoming = [sum(1 << u for u in range(n) if rows[u] >> v & 1) for v in range(n)]
    odd = [[0] * n for _ in rows]
    for first in range(n):
        ends = [0] * (1 << n)
        ends[1 << first] = 1 << first
        for mask in range(1, 1 << n):
            if not mask >> first & 1:
                continue
            if mask != 1 << first:
                for last in range(n):
                    if last != first and mask >> last & 1:
                        if (ends[mask ^ (1 << last)] & incoming[last]).bit_count() % 2:
                            ends[mask] |= 1 << last
            if mask.bit_count() % 2:
                for last in range(n):
                    odd[first][last] ^= (ends[mask] >> last) & 1
    return starts, odd


def main():
    source_paths = [HERE / 'converse_transfer.py', PACKET / 'response/normal.json',
                    ROOT / '05-knowledge/results/two_port_twelve_principles_20261103/PROOF.md']
    pins = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source_paths}
    data = json.loads((PACKET / 'response/normal.json').read_text())
    left, right = data['controls']
    n = 12
    masks = [m for m in range(1 << n) if m.bit_count() <= 2]
    delta = {(a, c): left['battery'][i * len(masks) + j] - right['battery'][i * len(masks) + j]
             for i, a in enumerate(masks) for j, c in enumerate(masks)}
    coeff = {}
    for a, c in delta:
        degree = a.bit_count() + c.bit_count()
        value = sum((-1 if (degree - x.bit_count() - y.bit_count()) % 2 else 1) * delta[x, y]
                    for x in subsets(a) for y in subsets(c))
        need(value % 2 == 0, 'integral contrast divided by two')
        coeff[a, c] = value
        if degree == 4 or (degree == 2 and not (a and c)):
            need(value == 0, 'exact canceled quartic or pure quadratic')
    allowed = {(a, c) for a, c in coeff if a.bit_count() + c.bit_count() <= 3
               and not (a.bit_count() + c.bit_count() == 2 and not (a and c))}
    need(len(allowed) == 1753, 'converse-only battery size')
    for a, c in delta:
        need(sum(coeff[x, y] for x in subsets(a) for y in subsets(c) if (x, y) in allowed)
             == delta[a, c], 'reconstruct full battery from reduced contrast coefficients')
    starts, odd = path_parities(left['rows'])
    h = [[0] * n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            for mask in range(1 << n):
                if not mask >> v & 1 and mask.bit_count() in (2, 4, 6, 8):
                    h[u][v] ^= starts[mask] >> u & 1
            a, c = (1 << u) | (1 << v), 1 << v
            need((coeff[a, c] // 2) % 2 == h[u][v], 'shared-label cubic support formula')
            a, c = 1 << v, (1 << u) | (1 << v)
            need((coeff[a, c] // 2) % 2 == h[u][v], 'second shared-label cubic support formula')
            skew = left['C'][u * n + v] - left['C'][v * n + u]
            need(skew % 2 == 0 and (skew // 2) % 2 == (odd[u][v] ^ odd[v][u]),
                 'odd-support directed-path cover-skew formula')
    need(pins == {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in source_paths}, 'source changed')
    print(json.dumps({'status': 'PASS INDEPENDENT TWELVE-CORE CONVERSE TRANSFER',
                      'scope': 'One specified named twelve-core and its converse; inherited all-order proofs remain separately scoped.',
                      'full_battery': len(delta), 'reduced_battery': len(allowed),
                      'exact_coefficient_counts': {str(d): sum(bool(v) for (a, c), v in coeff.items()
                                                             if a.bit_count() + c.bit_count() == d)
                                                   for d in range(5)},
                      'shared_label_cubic_matrix': h,
                      'shared_label_cubic_nonzero': sum(map(sum, h)),
                      'ordered_cubic_checks': 264, 'ordered_path_skew_checks': 132,
                      'source_sha256': pins}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
