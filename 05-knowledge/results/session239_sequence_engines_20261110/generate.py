#!/usr/bin/env python3
"""Convenience CLI for exact sequences from Session239. Standard library only."""
import argparse
import json
from math import comb
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / 'burnside'))
from engines import integer_response
from palette_engines import palette_series


def sequence(name, n):
    if not isinstance(n, int) or n < 0:
        raise ValueError('nonnegative integer order required')
    if name == 'A000568':
        return integer_response(n)[0]
    if name == 'A051337':
        a = integer_response(n)[0]
        strong = [0] * (n + 1)
        for j in range(1, n + 1):
            strong[j] = a[j] - sum(strong[k] * a[j-k] for k in range(1, j))
        strong[0] = 1  # OEIS convention only; empty object is not an SCC.
        return strong
    if name == 'A054946':
        labelled = [1 << comb(j, 2) for j in range(n+1)]
        strong = [0] * (n+1)
        for j in range(1, n+1):
            strong[j] = labelled[j] - sum(comb(j,k)*strong[k]*labelled[j-k] for k in range(1,j))
        return strong[1:]  # Offset 1 in OEIS.
    palettes = {'A000088': (2,2), 'A000273': (4,2), 'A001174': (3,1)}
    if name not in palettes:
        raise ValueError('unsupported sequence')
    result = palette_series(n, *palettes[name])[0]
    return result[1:] if name == 'A001174' else result


def main():
    # These are internally computed exact integers, not untrusted decimal input.
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequence', choices=['A000568','A051337','A054946','A000088','A000273','A001174'], default='A000568')
    parser.add_argument('--through', type=int, default=20)
    parser.add_argument('--binary-precision', type=int, help='A000568 residue mode only')
    args = parser.parse_args()
    if args.binary_precision is not None:
        if args.sequence != 'A000568':
            parser.error('--binary-precision applies only to A000568')
        values = integer_response(args.through, binary_precision=args.binary_precision)[0]
    else:
        values = sequence(args.sequence, args.through)
    print(json.dumps({'sequence': args.sequence, 'offset': 1 if args.sequence in ('A054946','A001174') else 0,
                      'through': args.through, 'binary_precision':args.binary_precision,
                      'terms':[str(x) for x in values]}, indent=2))

if __name__ == '__main__':
    main()
