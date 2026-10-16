#!/usr/bin/env python3
"""Exact controls for W2's conductor majorant, not a ROOT regularity test."""
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from math import comb
from pathlib import Path


def need(condition, name):
    if not condition:
        raise RuntimeError(name)


def depth(u, modulus):
    while u and u % 3 == 0:
        u //= 3
        modulus //= 3
    q = 0
    while modulus > 1:
        modulus //= 3
        q += 1
    return q if u else 0


def main():
    branch_checks = 0
    for Q in range(1, 5):
        modulus = 3**Q
        inv2 = pow(2, -1, modulus)
        for initial in range(1, modulus):
            q = depth(initial, modulus)
            # Retain frequency, phase and multiplicity; do not merge phases
            # and pretend the majorant is the actual Fourier norm.
            paths = Counter({(initial, 0): 1})
            for m in range(9):
                survivors = sum(paths.values())
                expected = sum(comb(m, j) for j in range(min(q, m+1)))
                need(survivors == expected, 'exact finite conductor survival')
                branch_checks += 1
                newer = Counter()
                for (u, phase), multiplicity in paths.items():
                    even_u = u*inv2 % modulus
                    newer[(even_u, phase)] += multiplicity
                    odd_u = 3*even_u % modulus
                    if odd_u:
                        need(depth(odd_u, modulus) == depth(u, modulus)-1,
                             'one conductor level lost')
                        newer[(odd_u, (phase+even_u) % modulus)] += multiplicity
                    else:
                        need(depth(u, modulus) == 1, 'rank-one death only at first level')
                paths = newer
    expectations = 0
    for c in (F(1, 3), F(1, 2), F(2, 3)):
        for q in range(1, 13):
            for M in range(25):
                spent = sum(sum(F(comb(m, j))*c**j*(1-c)**(m-j)
                                for j in range(min(q, m+1))) for m in range(M))
                remaining = sum(F(comb(M, j))*c**j*(1-c)**(M-j)*F(q-j, 1)/c
                                for j in range(min(q, M+1)))
                need(spent+remaining == F(q, 1)/c, 'truncated waiting time plus exact tail')
                expectations += 1
    result = dict(status='FINITE-EXACT controls of the proved positive majorant; source regularity OPEN',
                  Fourier_branch_checks=branch_checks, waiting_time_checks=expectations,
                  probabilities=['1/3', '1/2', '2/3'], conductor_levels=[1, 4],
                  branch_depths=[0, 8], waiting_depths=[1, 12], truncations=[0, 24])
    data = json.dumps(result, indent=2, sort_keys=True)+'\n'
    path = Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_critical_conductor_20261015e.json'
    path.write_text(data)
    print(json.dumps(dict(output=str(path), sha256=hashlib.sha256(data.encode()).hexdigest())))


if __name__ == '__main__':
    main()
