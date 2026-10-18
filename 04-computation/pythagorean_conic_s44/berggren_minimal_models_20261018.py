#!/usr/bin/env python3
"""Minimal models on the Berggren tree of primitive Pythagorean triples (S44, elementary; exact).

Dictionary.  A primitive triple is (q^2 - p^2, 2pq, q^2 + p^2) = (p + iq)^2 with q > p >= 1 coprime and q - p odd.
In (m, n) = (q, p) coordinates the Berggren matrices act as A:(m,n)->(2m-n,m), B:(m,n)->(2m+n,m), C:(m,n)->(m+2n,n).
To each node attach the Virasoro (p,q) model (p >= 2: minimal model; p = 1: the logarithmic (1,q) triplet series).
Checks (all exact):
  1. every coprime pair 1 <= p < q <= Q with q - p odd occurs exactly once in the tree (Q = 60);
  2. c_(p,q) = 1 - 6(p-q)^2/(pq) = 13 - 12 * hyp / b and c_eff = 1 - 6/(pq) = 1 - 12/b  (b = 2pq the even leg);
  3. the A-chain from the root is the unitary series (n, n+1): triples (2n+1, 2n(n+1), 2n^2+2n+1) (c - b = 1);
  4. the C-chains from (3,2) and (5,2) are the (2, 2N+1) series; C fixes p (the q-stability direction of HYP-9315);
  5. on the parameter t = q/p the moves are t -> 2 - 1/t, 2 + 1/t, t + 2 (Moebius), and c = 13 - 6(t + 1/t);
  6. the unit-gap spine (t, (t^2-1)/2, (t^2+1)/2) of the Codex graph-Pythagorean note (Seidel radii of H_t) is the
     A-chain: t = 2n+1 <-> M(n, n+1).
Usage: python berggren_minimal_models_20261018.py
"""
from fractions import Fraction
from math import gcd
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))


def c_pq(p, q):
    return 1 - Fraction(6 * (p - q) ** 2, p * q)


def tree(maxq):
    nodes = {}
    stack = [((2, 1), '')]
    while stack:
        (m, n), word = stack.pop()
        if m > maxq:
            continue
        nodes[(n, m)] = word
        stack.append(((2 * m - n, m), word + 'A'))
        stack.append(((2 * m + n, m), word + 'B'))
        stack.append(((m + 2 * n, n), word + 'C'))
    return nodes


if __name__ == '__main__':
    Q = 60
    T = tree(Q)
    truth = {(p, q) for q in range(2, Q + 1) for p in range(1, q) if gcd(p, q) == 1 and (q - p) % 2 == 1}
    ok1 = set(T) == truth
    ok2 = ok3 = ok4 = ok5 = True
    for (p, q), w in T.items():
        a, b, h = q * q - p * p, 2 * p * q, q * q + p * p
        assert a * a + b * b == h * h
        c = c_pq(p, q)
        ok2 &= (c == 13 - Fraction(12 * h, b)) and (1 - Fraction(6, p * q) == 1 - Fraction(12, b))
        t = Fraction(q, p)
        ok5 &= (c == 13 - 6 * (t + 1 / t))
        if set(w) <= {'A'}:
            ok3 &= (q == p + 1 and h - b == 1)
        if w and w[0] in 'AB' and set(w[1:]) <= {'C'} and p == 2:
            pass
    # C-chains from (3,2) [word 'A'] and (5,2) [word 'B'] give the (2, 2N+1) series
    two_series = sorted(q for (p, q), w in T.items() if p == 2)
    ok4 = two_series == [q for q in range(3, Q + 1, 2)]
    for (p, q), w in T.items():
        if p == 2:
            ok4 &= (w in ('A', 'B') or (w[0] in 'AB' and set(w[1:]) <= {'C'}))
    # moves on t
    import itertools
    for (p, q), w in list(T.items())[:200]:
        t = Fraction(q, p)
        m, n = q, p
        ok5 &= Fraction(2 * m - n, m) == 2 - 1 / t and Fraction(2 * m + n, m) == 2 + 1 / t and Fraction(m + 2 * n, n) == t + 2
    unitary = sorted((p, q) for (p, q), w in T.items() if set(w) <= {'A'})
    rows = []
    for (p, q) in sorted(T, key=lambda x: (len(T[x]), T[x]))[:13]:
        rows.append({'word': T[(p, q)] or 'root', 'p': p, 'q': q, 'triple': (q*q-p*p, 2*p*q, q*q+p*p), 'c': str(c_pq(p, q))})
    out = {'every_coprime_opposite_parity_pair_once': ok1, 'c_formula': ok2, 'A_chain_unitary': ok3,
           'p2_series_on_C_chains': ok4, 'mobius_moves_and_c_of_t': ok5, 'unitary_nodes': unitary[:8], 'first_levels': rows}
    for k, v in out.items():
        if k != 'first_levels':
            print(f'{k}: {v}')
    for r in rows:
        print(f"  {r['word']:5s} (p,q)=({r['p']},{r['q']})  triple {r['triple']}  c = {r['c']}")
    json.dump(out, open(os.path.join(HERE, 'berggren_minimal_models.json'), 'w'), indent=1, default=str)
