#!/usr/bin/env python3
"""Audit check (e)(i): Berggren tree <-> Virasoro (p,q) labels, central charges, Moebius moves and branch limits.

Independent re-derivation:
  * tree on Euclid parameters (m, n) = (q, p) from (2, 1) with A:(2m-n, m), B:(2m+n, m), C:(m+2n, n); every coprime
    pair q > p >= 1 with q - p odd appears exactly once (q <= 300);
  * c = 1 - 6(p-q)^2/(pq) equals 13 - 12 h/b (b = 2pq even leg, h = p^2 + q^2), c_eff = 1 - 6/(pq) = 1 - 12/b;
  * branch limits: iterate words A^k, B^k, (CA)^k, (AC)^k, C^k from the root, record t = q/p and c, and compare with
    the claimed limits (A: t -> 1, c -> 1; B: t -> 1 + sqrt2, c -> 13 - 12 sqrt2; C,A alternating: the note claims
    t -> sqrt3, c -> 13 - 8 sqrt3).
"""
import math, json, os
from fractions import Fraction
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))


def c_pq(p, q):
    return 1 - Fraction(6 * (p - q) ** 2, p * q)


def move(w, mn):
    m, n = mn
    return {'A': (2 * m - n, m), 'B': (2 * m + n, m), 'C': (m + 2 * n, n)}[w]


def main():
    Q = 300
    seen = {}
    stack = [((2, 1), '')]
    while stack:
        (m, n), word = stack.pop()
        if m > Q:
            continue
        key = (n, m)
        seen[key] = seen.get(key, 0) + 1
        for w in 'ABC':
            stack.append((move(w, (m, n)), word + w))
    truth = {(p, q) for q in range(2, Q + 1) for p in range(1, q) if gcd(p, q) == 1 and (q - p) % 2}
    bij = set(seen) == truth and all(v == 1 for v in seen.values())
    cform = all(c_pq(p, q) == 13 - Fraction(12 * (p * p + q * q), 2 * p * q) and
                1 - Fraction(6, p * q) == 1 - Fraction(12, 2 * p * q) for (p, q) in truth)
    out = {'bijection_q_le_300': bij, 'c_formula': cform}
    print('tree bijection (q <= 300):', bij, ' c = 13 - 12h/b and c_eff = 1 - 12/b:', cform)
    branches = {}
    for name, word in (('A^k', 'A' * 14), ('B^k', 'B' * 14), ('(CA)^k', 'CA' * 7), ('(AC)^k', 'AC' * 7)):
        mn = (2, 1)
        seq = []
        for k, w in enumerate(word):
            mn = move(w, mn)
            p, q = mn[1], mn[0]
            seq.append({'step': k + 1, 'move': w, 'p': p, 'q': q, 't': float(Fraction(q, p)), 'c': float(c_pq(p, q))})
        branches[name] = seq
        tail = seq[-4:]
        print(f'{name}: last (p,q,t,c): ' + '; '.join(f"({s['p']},{s['q']}, t={s['t']:.6f}, c={s['c']:.6f})" for s in tail))
    r3, r2 = math.sqrt(3), math.sqrt(2)
    print(f'reference: 1+sqrt2 = {1+r2:.6f}, 13-12sqrt2 = {13-12*r2:.6f}; sqrt3 = {r3:.6f}, 2+sqrt3 = {2+r3:.6f}, '
          f'13-8sqrt3 = {13-8*r3:.6f}, c(2+sqrt3) = {13-6*((2+r3)+1/(2+r3)):.6f}')
    ca = branches['(CA)^k']
    after_A = [s for s in ca if s['move'] == 'A']
    after_C = [s for s in ca if s['move'] == 'C']
    out['CA_alternating'] = {
        'after_A_t_limit': after_A[-1]['t'], 'after_A_c_limit': after_A[-1]['c'],
        'after_C_t_limit': after_C[-1]['t'], 'after_C_c_limit': after_C[-1]['c'],
        'claimed_single_limit_t_sqrt3_c_13_minus_8sqrt3': False if abs(after_C[-1]['t'] - r3) > 0.5 else True}
    out['branches'] = branches
    print('C,A alternating: after-A nodes t ->', after_A[-1]['t'], ', c ->', after_A[-1]['c'],
          '; after-C nodes t ->', after_C[-1]['t'], ', c ->', after_C[-1]['c'])
    json.dump(out, open(os.path.join(HERE, 'berggren_audit.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
