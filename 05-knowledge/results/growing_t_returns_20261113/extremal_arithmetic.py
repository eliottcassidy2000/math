"""Finite controls for the arithmetic corollaries of the exact ready minimum."""
from pathlib import Path
from math import isqrt
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def mu(M):
    if M < 8:
        return [5, 9, 17, 33][M - 4]
    r = M - 4
    s = (r - 1) // 2
    return 3 * (1 << r) + 3 * (1 << s) + (1 << (r - s)) - 1


def run(output):
    squares = []
    for M in range(4, 1001):
        value = mu(M)
        if M >= 9 and M % 2:
            x = 1 << ((M - 5) // 2)
            need(value == (6 * x - 1) * (x + 1), 'odd-order factorization')
        if M >= 8 and M % 2 == 0:
            q = (M - 4) // 2
            x = 1 << (q - 1)
            need(value == 12 * x * x + 7 * x - 1, 'even-order quadratic')
        if M >= 10 and M % 4 == 2:
            need(value % 3 == 0 and value > 3, 'composite residue class')
        if M >= 11:
            need(value % 8 == 7, 'stabilized residue')
        root = isqrt(value - 1)
        if root * root == value - 1:
            squares.append(dict(order=M, H=value, square_root=root))
    need(squares == [dict(order=4,H=5,square_root=2), dict(order=6,H=17,square_root=4)], 'exact finite near-square controls')
    report = dict(status='FINITE-EXACT controls; infinite consequences proved in EXTREMAL_ARITHMETIC.md',
                  universe='Orders4..1000; formula from the separately proved exact strong unit-ready minimum',
                  prefix=[dict(order=M,H=mu(M)) for M in range(4,79)],
                  square_plus_one_controls=squares,
                  inherited_aliquot_source='05-knowledge/results/uniform_return_metagraph_20261111/arithmetic/PROOF.md')
    (HERE / output).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('PASS:997 orders; square-plus-one controls5 and17 only')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='extremal_arithmetic.json')
    run(parser.parse_args().output)
