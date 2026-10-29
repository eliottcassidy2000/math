#!/usr/bin/env python3
"""Cardinality-optimal good-prime F2 witness banks for ranks 30,31,32.

Extends the August rank-30 observer construction. Uses the audited incoming
point snapshots, but independently recomputes every field root and bit.
Optimality is in this specified local-character observation model only.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

P = Path(__file__).resolve().parent
S = P.parent/'session97_overnight_arithmetic_carriers_20261029'/'ranks'
CHECKS = 0


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(message)


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def add_rows(basis, rows):
    basis = basis.copy()
    for v in rows:
        while v:
            k = v.bit_length()-1
            if k in basis:
                v ^= basis[k]
            else:
                basis[k] = v
                break
    return basis


def independent_rank(rows, width):
    """Column pivot elimination, independent of the high-bit row basis."""
    rows = rows[:]
    rank = 0
    for k in range(width):
        pivot = next((i for i in range(rank, len(rows)) if rows[i] >> k & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(rank+1, len(rows)):
            if rows[i] >> k & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def certify(cid):
    path = S/f'curve{cid}.json'
    c = json.loads(path.read_text())
    a1, a2, a3, a4, a6 = map(int, c['ainvs'])
    pts = [tuple(map(F, xy)) for xy in c['points']]
    width = len(pts)
    for x, y in pts:
        require(y*y+a1*x*y+a3*y == x**3+a2*x*x+a4*x+a6, 'curve incidence')
    b2, b4, b6 = a1*a1+4*a2, 2*a4+a1*a3, a3*a3+4*a6
    disc = int(c['discriminant'])
    b8 = a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    require(disc == -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6,
            'independently computed discriminant')
    basis, banks, no_two = {}, [], None
    for p in range(3, 5000, 2):
        if not prime(p) or disc % p == 0:
            continue
        roots = [r for r in range(p) if (4*r**3+b2*r*r+2*b4*r+b6) % p == 0]
        if not roots and no_two is None:
            no_two = p
        if len(roots) != 3 or any(x.denominator % p == 0 or y.denominator % p == 0
                                 for x, y in pts):
            continue
        residues = [x.numerator*pow(x.denominator, -1, p) % p for x, y in pts]
        if any(x in roots for x in residues):
            continue  # special-order-two values are not needed in this bank
        rows = []
        for r in roots:
            row = sum((pow((x-r) % p, (p-1)//2, p) == p-1) << i
                      for i, x in enumerate(residues))
            rows.append(row)
        require(rows[0] ^ rows[1] ^ rows[2] == 0, 'three-root parity relation')
        trial = add_rows(basis, rows)
        needed = min(2, width-len(basis))
        if len(trial)-len(basis) != needed:
            continue
        banks.append(dict(prime=p, roots=roots, rows=rows,
                          rank_before=len(basis), rank_after=len(trial)))
        basis = trial
        if len(basis) == width and no_two:
            break
    flat = [row for bank in banks for row in bank['rows']]
    require(len(basis) == width, 'full witness rank')
    require(len(banks) == (width+1)//2, 'optimal prime cardinality')
    require(independent_rank(flat, width) == width, 'independent elimination')
    require(no_two is not None, 'odd rational torsion guard')
    for bank in banks:
        without = [row for other in banks if other is not bank for row in other['rows']]
        require(independent_rank(without, width) < width, 'each prime is needed')
    # Duplicate the first point-column into the last; the witness must fail.
    duplicated = [(row & ~(1 << (width-1))) | ((row & 1) << (width-1)) for row in flat]
    require(independent_rank(duplicated, width) < width, 'duplicate-point hostile')
    return dict(curve_id=cid, witness_rank=width, optimal_bank_size=len(banks),
                additional_no_two_torsion_prime=no_two, banks=banks,
                source_sha256=sha256(path.read_bytes()).hexdigest())


def main():
    rows = [certify(cid) for cid in (273, 302, 4239)]
    print(json.dumps(dict(status='PASS; optimal within good-prime binary-character model',
                          checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                          explicit_checks=CHECKS, curves=rows), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
