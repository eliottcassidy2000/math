#!/usr/bin/env python3
"""Is a small quadratic field forced inside every field of definition of EI17's realization scheme? (S40, THM-4665)

If every component field of the zero-dimensional realization scheme contained sqrt(d), then a split-model F_p-point
would need (d/p) = +1, and a twisted-model point (embedding in G_p, complex conjugation acting as Frobenius) would need
(d/p) = +1 for d > 0 and -1 for d < 0. The primes with points are read from the committed S39 outputs.
Usage: python ei17_quadratic_subfield_test_20261017.py
Session opus-2026-10-17-S40.
"""
import os, re
from sympy import jacobi_symbol, factorint

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', '..', '05-knowledge', 'results')


def parse(name):
    hit, miss = [], []
    for line in open(os.path.join(RES, name), encoding='utf-8'):
        m = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
        if m:
            p, f = int(m.group(1)), int(m.group(3))
            (hit if f > 0 else miss).append(p)
    return hit, miss


def leg(d, p):
    return int(jacobi_symbol(d % p, p))


def main():
    S_hit, S_miss = parse('ei17_finite_field_realizations_20261017.out')
    T_hit, T_miss = parse('ei17_twisted_unit_graph_embeddings_20261017.out')
    print('split model: %d primes with points, %d without' % (len(S_hit), len(S_miss)))
    print('twisted model: %d primes with points, %d without' % (len(T_hit), len(T_miss)))
    ds = [d for d in range(-200, 201) if abs(d) > 1 and all(e == 1 for e in factorint(abs(d)).values())]
    split_ok = [d for d in ds if all(leg(d, p) == 1 for p in S_hit if p > 2 and p % abs(d))]
    twist_ok = [d for d in ds if all(leg(d, p) == (1 if d > 0 else -1) for p in T_hit if p > 2 and p % abs(d))]
    print('squarefree d (1 < |d| <= 200) that are squares at every split prime with points:', split_ok)
    print('squarefree d (1 < |d| <= 200) with the forced twisted sign at every twisted prime with points:', twist_ok)
    print('d = -1 (sqrt(-1) in the field): split %s, twisted %s' % (
        all(leg(-1, p) == 1 for p in S_hit), all(leg(-1, p) == -1 for p in T_hit)))
    print('split primes with points:', S_hit)
    print('twisted primes with points:', T_hit)


if __name__ == '__main__':
    main()
