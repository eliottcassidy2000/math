"""Exact finite-field character certificates for public ICARM point banks.

No Sage, floating regulators, GRH or BSD. Proves lower bounds, not upper
bounds. The classical 2-descent character lemma is stated in ELLIPTIC.md.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def primes(limit):
    return [p for p in range(3, limit + 1, 2)
            if all(p % d for d in range(3, isqrt(p) + 1, 2))]


def insert(basis, row):
    while row:
        pivot = row.bit_length() - 1
        if pivot in basis:
            row ^= basis[pivot]
        else:
            basis[pivot] = row
            return True
    return False


def binary_rank(rows):
    b = {}
    for row in rows:
        insert(b, row)
    return len(b)


def legendre_bit(v, p):
    t = pow(v % p, (p - 1) // 2, p)
    need(t in (1, p - 1), 'nonzero character argument')
    return int(t == p - 1)


def add(P, Q, p, A, B):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and (y + v) % p == 0:
        return None
    if x == u:
        slope = (3*x*x + 2*A*x + B) * pow(2*y, -1, p) % p
    else:
        slope = (v-y) * pow((u-x) % p, -1, p) % p
    z = (slope*slope-A-x-u) % p
    return z, (-y+slope*(x-z)) % p


def check_finite_group(p, b2, b4, b6, roots):
    """Independent exhaustive addition-table check of every character."""
    A, B, C = b2*pow(4, -1, p)%p, b4*pow(2, -1, p)%p, b6*pow(4, -1, p)%p
    points = [None] + [(x, y) for x in range(p) for y in range(p)
                        if (y*y-x*x*x-A*x*x-B*x-C) % p == 0]
    signatures = {}
    for P in points:
        val = 0
        for i, r in enumerate(roots):
            if P is None:
                bit = 0
            else:
                x, _ = P
                v = (x-r) % p if x != r else (3*r*r+2*A*r+B) % p
                bit = legendre_bit(v, p)
            val |= bit << i
        signatures[P] = val
    for P in points:
        for Q in points:
            need(signatures[add(P, Q, p, A, B)] == signatures[P] ^ signatures[Q],
                 ('character homomorphism', p, P, Q))
    return len(points) ** 2


def certificate(identifier):
    path = HERE / f'icarm_{identifier}.json'
    raw = path.read_bytes()
    data = json.loads(raw)
    a1, a2, a3, a4, a6 = map(int, data['ainvs'])
    points = [tuple(map(Fraction, xy)) for xy in data['points']]
    for x, y in points:
        need(y*y+a1*x*y+a3*y == x*x*x+a2*x*x+a4*x+a6, 'point on curve')
    b2, b4, b6 = a1*a1+4*a2, a1*a3+2*a4, a3*a3+4*a6
    b8 = a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    delta = -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
    need(delta == int(data['discriminant']) and delta != 0, 'discriminant')
    basis, retained, no_two, dependent, group_checks = {}, [], None, 0, []
    for p in primes(2000):
        if delta % p == 0 or any(x.denominator % p == 0 for x, _ in points):
            continue
        roots = [r for r in range(p) if (4*r**3+b2*r*r+2*b4*r+b6) % p == 0]
        if not roots and no_two is None:
            no_two = p
        rows = []
        for r in roots:
            row = 0
            for i, (x, _) in enumerate(points):
                xmod = x.numerator*pow(x.denominator, -1, p) % p
                v = (xmod-r) % p
                if v == 0:
                    v = (12*r*r+2*b2*r+2*b4) % p
                row |= legendre_bit(v, p) << i
            rows.append(row)
            if insert(basis, row):
                retained.append(dict(prime=p, root=r, row=row))
        if len(rows) == 3:
            need(rows[0] ^ rows[1] ^ rows[2] == 0, 'three-root free parity bit')
            dependent += 1
        # Check one split and one single-root complete finite group per curve.
        if roots and len(roots) not in [g['root_count'] for g in group_checks]:
            group_checks.append(dict(prime=p, root_count=len(roots),
                                    pair_checks=check_finite_group(p, b2, b4, b6, roots)))
        if len(basis) == len(points) and no_two is not None:
            break
    need(len(basis) == len(points), 'independent character columns')
    need(no_two is not None, 'no rational two-torsion certificate')
    # A second, row-list elimination and a deliberately duplicated column.
    rows = [r['row'] for r in retained]
    need(binary_rank(rows) == len(points), 'second elimination')
    duplicated = [(r & ~(1 << (len(points)-1))) |
                  ((r & 1) << (len(points)-1)) for r in rows]
    need(binary_rank(duplicated) == len(points)-1, 'duplicate-point hostile')
    return dict(id=identifier, rank_lower_bound=len(points),
                source=f'https://elliptic-rank.icarm.cloud/curve/{identifier}.json',
                source_sha256=hashlib.sha256(raw).hexdigest(),
                no_rational_two_torsion_prime=no_two, retained_rows=retained,
                largest_prime=max(r['prime'] for r in retained),
                three_root_parity_checks=dependent, full_group_checks=group_checks,
                status='FINITE-EXACT rank lower bound using classical character lemma')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    result = [certificate(i) for i in (273, 302, 4239)]
    path = HERE / 'rank_certificate.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        need(json.loads(path.read_text(encoding='utf-8')) == result, 'stored rank receipt')
    print(json.dumps([dict(id=r['id'], rank=r['rank_lower_bound'],
                           maxprime=r['largest_prime'],
                           pair_checks=sum(g['pair_checks'] for g in r['full_group_checks']))
                      for r in result]))


if __name__ == '__main__':
    main()
