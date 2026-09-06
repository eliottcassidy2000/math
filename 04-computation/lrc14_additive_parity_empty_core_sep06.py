#!/usr/bin/env python3
"""Exact finite bases for the all-height additive ternary-unit family.

Fresh normalized ray formulas are compared with the frozen complete literal
interval and modular raw paths. Analytic tails live in the companion note.
"""

from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from math import gcd
from pathlib import Path

R = Q(3, 14)
TARGET = Q(6, 55)
CHECKS = 0
INHERITED_HASH = 'b99af309ff6f8643dedf923f5ee8d86d67b32ff2b0d6510209c565820894f399'


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def inherited_verifier():
    path = Path(__file__).with_name('lrc14_parity_empty_core_sep06.py')
    need(sha256(path.read_bytes()).hexdigest() == INHERITED_HASH,
         'frozen independent literal/raw source hash')
    spec = importlib.util.spec_from_file_location('frozen_parity_verifier', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalized_ray(w):
    a, b, c = w
    t = Q(a, c)
    positive_addresses = [k for k in range(1, (3*c-1)//14+1) if k % 3]
    values = []
    for k in positive_addresses:
        x = Q(k, c)
        values.append((min(2*R, (R*(2-t)-x)/(1-t)),
                       min(2*R, (R*(1+t)-x)/t),
                       min(2*R, (R-x)/(t*(1-t)))))
    projections = tuple(Q(2, c)*sum((row[i] for row in values), Q(0))
                        for i in range(3))
    physical = Q(2, c)*sum((min(row) for row in values), Q(0))
    carriers = {(sign*k, sign*k, -sign*k)
                for k in positive_addresses for sign in (-1, 1)}
    return projections, physical, carriers, tuple(positive_addresses)


def main():
    verifier = inherited_verifier()
    rows = [(a, c-a, c) for c in range(1, 60) if c % 3
            for a in range(1, (c+1)//2)
            if a % 3 and (c-a) % 3 and gcd(a, c) == 1]
    need(len(rows) == 136, 'complete projection finite head')
    need(sum(w[2] <= 33 for w in rows) == 42, 'complete physical finite head')
    need(TARGET-Q(39, 392)-Q(4, 7*60) == Q(1, 12936), 'network tail margin')
    need(TARGET-Q(9, 98)-Q(4, 7*34) == Q(41, 91630), 'physical tail margin')
    equality = {'network': [], 'physical': []}
    digest = sha256()
    controls = {}
    runner_up = (Q(-1), [])
    for w in rows:
        a, b, c = w
        need(a < b < c and a+b == c and sum(v % 2 == 0 for v in w) == 1,
             ('typed additive universe', w))
        E, physical, carriers, addresses = normalized_ray(w)
        literal_E, literal_physical = verifier.native(w)
        raw_E, raw_physical, raw_carriers = verifier.raw(w)
        need((E, physical) == (literal_E, literal_physical), ('native comparison', w))
        need((E, physical, carriers) == (raw_E, raw_physical, raw_carriers),
             ('complete raw comparison', w))
        need(physical <= min(E) and E[0] <= E[1], ('projection and physical order', w))
        t = Q(a, c)
        need(E[0] < (9+3*t)/98+Q(4, 7*c), ('first projection quadrature', w))
        need(E[2] < 12*(1-t+t*t)/98+Q(4, 7*c), ('third projection quadrature', w))
        need(physical < Q(9, 98)+Q(4, 7*c), ('physical quadrature', w))
        for name, value in (('network', min(E)), ('physical', physical)):
            need(value <= TARGET, ('finite target', name, w, value))
            if value == TARGET:
                equality[name].append(w)
        if min(E) < TARGET:
            old, witnesses = runner_up
            if min(E) > old:
                runner_up = min(E), [w]
            elif min(E) == old:
                witnesses.append(w)
        if w in ((1, 4, 5), (2, 5, 7), (1, 10, 11), (28, 31, 59)):
            controls[w] = E, physical, addresses
        digest.update(repr((w, E, physical, tuple(sorted(carriers)))).encode())
    need(equality == {'network': [(1, 10, 11)], 'physical': [(1, 10, 11)]},
         ('exact equality boundary', equality))
    print('universe primitive sorted positive ternary-units, c=a+b, c<=59')
    print('complete_network_head_rows', len(rows))
    print('complete_physical_head_rows_c_le_33', sum(w[2] <= 33 for w in rows))
    print('sharp_bound', TARGET)
    print('unique_equality', equality)
    print('finite_network_runner_up', runner_up)
    print('network_tail c>=60; margin at60', Q(1, 12936))
    print('physical_tail c>=34; margin at34', Q(41, 91630))
    for w, values in sorted(controls.items()):
        print('control', w, values)
    print('complete_semantic_digest', digest.hexdigest())
    print('new_explicit_checks', CHECKS)
    print('inherited_explicit_checks', verifier.CHECKS)
    print('inherited_source_sha256', INHERITED_HASH)
    print('PASS: additive finite bases; all-height tails require companion analytic proof')


if __name__ == '__main__':
    main()
