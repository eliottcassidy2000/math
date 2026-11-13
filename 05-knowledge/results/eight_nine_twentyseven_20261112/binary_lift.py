"""Explicit binary-residue constructor from the two-switch 8-vertex selector."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'slot_polynomial'))
from check import path_counts

CORE = (56, 117, 129, 22, 100, 76, 141, 59)
BASE_COUNTS = {(0, 0): 513, (1, 0): 515, (0, 1): 517, (1, 1): 495}


def choose(precision, residue):
    if precision < 3 or residue % 2 == 0:
        raise ValueError('precision >=3 and odd residue required')
    modulus = 1 << precision
    residue %= modulus
    switches = next(ab for ab, h in BASE_COUNTS.items() if h % 8 == residue % 8)
    base = BASE_COUNTS[switches]
    target_ratio = residue * pow(base, -1, modulus) % modulus
    t = 0
    for j in range(precision - 3):
        bit_modulus = 1 << (j + 4)
        if pow(9, t, bit_modulus) != target_ratio % bit_modulus:
            t += 1 << j
    if pow(9, t, modulus) != target_ratio:
        raise RuntimeError('digit lift failed')
    return switches, t


def ordinal(left, right):
    n, m = len(left), len(right)
    suffix = ((1 << m) - 1) << n
    return tuple(row | suffix for row in left) + tuple(row << n for row in right)


def transitive(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def realize(precision, residue):
    (a, b), t = choose(precision, residue)
    rows = list(CORE)
    for selected, i, j in ((a, 1, 6), (b, 2, 3)):
        if selected:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    rows = tuple(rows)
    triangle = (2, 4, 1)
    nine = ordinal(triangle, triangle)
    for _ in range(t):
        rows = ordinal(rows, nine)
    maximum = (1 << (precision - 3)) - 1
    rows = ordinal(rows, transitive(6 * (maximum - t)))
    return rows, BASE_COUNTS[(a, b)] * 9 ** t


def controls():
    residue_cases = 0
    for s in range(3, 11):
        seen = set()
        for target in range(1, 1 << s, 2):
            ab, t = choose(s, target)
            address = (ab, t)
            if address in seen or not 0 <= t < (1 << (s - 3)):
                raise RuntimeError('nonunique address')
            seen.add(address)
            if BASE_COUNTS[ab] * pow(9, t, 1 << s) % (1 << s) != target:
                raise RuntimeError('wrong modular output')
            residue_cases += 1
    actual = []
    for target in range(1, 16, 2):
        rows, count = realize(4, target)
        if len(rows) != 14 or path_counts(rows)[-1] != count or count % 16 != target:
            raise RuntimeError('literal common-order construction')
        actual.append({'target': target, 'order': len(rows), 'actual_H': count})
    data = {
        'status': 'FINITE-EXACT controls; all-precision proof in BINARY_LIFT.md',
        'precisions': [3, 10], 'odd_residue_cases': residue_cases,
        'literal_precision4_outputs': actual,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE / 'binary_lift_receipt.json').write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps({'passed': True, 'residue_cases': residue_cases,
                      'literal_outputs': len(actual), 'common_order': 14}))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--precision', type=int)
    ap.add_argument('--residue', type=int)
    args = ap.parse_args()
    if args.precision is None and args.residue is None:
        controls()
    elif args.precision is None or args.residue is None:
        ap.error('supply both --precision and --residue')
    else:
        rows, count = realize(args.precision, args.residue)
        ab, t = choose(args.precision, args.residue)
        print(json.dumps({'precision': args.precision, 'target_residue': args.residue % (1 << args.precision),
                          'switches': ab, 'nine_blocks': t, 'order': len(rows),
                          'H': count, 'outgoing_rows': rows}))


if __name__ == '__main__':
    main()
