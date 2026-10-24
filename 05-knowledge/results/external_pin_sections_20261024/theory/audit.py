#!/usr/bin/env python3
"""Reproducible receiver edge cases with an independent state enumerator.

The production compiler supplies routes.  Literal binary tuples independently
define the feasible states, charges, populations, costs, and endpoint bounds.
The production checker's three extra finite controls are also replayed, to
retain the originally reported 82,801-check / 3,924-route receipt.
"""
from collections import defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('receiver_boundary_producer', HERE/'check.py')
producer = module_from_spec(spec)
spec.loader.exec_module(producer)
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def binary_sets(size):
    return [frozenset(j for j, bit in enumerate(bits) if bit)
            for bits in product((0, 1), repeat=size)]


def independent(model):
    R, E = binary_sets(model.h), binary_sets(model.e)
    helper_records = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

    def q(state):
        I, J, y = state
        return (sum(model.R[i] for i in I)+sum(model.E[j] for j in J)
                +(0, model.n-1, 1 % model.n)[y]) % model.n

    def normalized_outside(state, k):
        I, J, _ = state
        J = set(J)
        r, r0 = len(I), max(model.lo, k-model.e)
        if r < model.lo:
            J.difference_update(sorted(J)[:model.lo-r])
            r = model.lo
        elif r > model.hi:
            J.update(sorted(set(range(model.e))-J)[:r-model.hi])
            r = model.hi
        while r > r0:
            J.add(min(set(range(model.e))-J))
            r -= 1
        return J

    fibres = defaultdict(list)
    for I, J, y in product(R, E, range(3)):
        k = len(I)+len(J)
        if model.lo <= k <= model.hi+model.e:
            state = I, J, y
            fibres[k, q(state)].append(state)
    for J, y, r in product(E, range(3), range(model.lo, model.hi+1)):
        charges = {sum(model.R[i] for i in I) % model.n for I in R if len(I) == r}
        need(charges == set(range(model.n)), 'independent complete frame section')

    routes, events = 0, 0
    for (k, total), states in fibres.items():
        for source, target in product(states, repeat=2):
            route = model.route(source, target)
            routes += 1
            need(route[0][0] == source and route[-1][0] == target,
                 'independently enumerated literal endpoints')
            dE = len(normalized_outside(source, k)-normalized_outside(target, k))
            limit = 2*(model.hi-model.lo+1)+dE+abs(source[2]-target[2])+1
            need(len(route)-1 <= limit, 'independent event bound with refresh')
            for state, _ in route:
                I, J, y = state
                need(I in R and J in E and y in range(3), 'declared full record')
                need(len(I)+len(J) == k and q(state) == total,
                     'independent fixed core population and charge')
                need(sum(helper_records[y]) == 1, 'independent helper population')
            for (before, _), (after, kind) in zip(route, route[1:]):
                events += 1
                if kind.startswith('reverse-'):
                    before, after = after, before
                    kind = kind.removeprefix('reverse-')
                I, J, y = before
                I2, J2, y2 = after
                r, s = len(I), len(I2)
                sharp = min(r+s, 2*model.h-r-s)
                costR = sum(model.Rcost[i] for i in I ^ I2)
                costE = sum(model.Ecost[j] for j in J ^ J2)
                costH = sum(c for c, a, b in zip(model.Hcost, helper_records[y],
                                                 helper_records[y2]) if a != b)
                need(len(I ^ I2) <= sharp and
                     costR <= sum(sorted(model.Rcost, reverse=True)[:sharp]),
                     'independent clamped grade cost')
                need(costR+costE+costH == model.event(before, after, kind),
                     'literal occurrence cost agrees with compiler audit')
                need(kind == 'helper' or y == y2, 'explicit helper skeleton only')
                if kind == 'entry-low':
                    ok = r < model.lo and s == model.lo and J2 <= J and len(J-J2) == s-r
                elif kind == 'entry-high':
                    ok = r > model.hi and s == model.hi and J <= J2 and len(J2-J) == r-s
                elif kind == 'adjacent':
                    ok = s == r-1 and J <= J2 and len(J2-J) == 1
                elif kind == 'outside':
                    ok = r == s and len(J ^ J2) == 2
                elif kind == 'helper':
                    ok = r == s and J == J2 and abs(y-y2) == 1
                elif kind == 'refresh':
                    ok = r == s and J == J2 and y == y2
                else:
                    ok = False
                need(ok, 'independent exact event type')
    return routes, events


def main():
    models = [
        producer.ReceiverModel('audit_C2_singleton_band', [0, 1, 0, 1], [0, 1], 2, 2, 2, False),
        producer.ReceiverModel('audit_trivial_empty_receiver', [], [0, 0], 1, 0, 0, False),
        producer.ReceiverModel('audit_trivial_empty_outside_singleton', [0, 0], [], 1, 1, 1, False),
    ]
    receipts = [model.check() for model in models]
    need([row['ordered_endpoint_pairs'] for row in receipts] == [3834, 54, 36],
         'frozen extra-model endpoint counts')
    need(producer.CHECKS == 82801, 'frozen producer boundary check count')
    compiler_checks = producer.CHECKS
    additional = [independent(model) for model in models]
    need([row[0] for row in additional] == [3834, 54, 36],
         'independent literal state enumerator agrees on universe')
    source = frozenset({0}), frozenset(), 0
    target = frozenset({1}), frozenset(), 0
    route = models[-1].route(source, target)
    need(len(route) == 2 and route[-1][1] == 'refresh',
         'same frame distinct receiver needs final plus-one reset')
    print(f'PASS {compiler_checks:,} frozen compiler checks / 3,924 routes; '
          f'{CHECKS:,} independent checks / {sum(x[1] for x in additional):,} events; '
          'literal final refresh confirmed')


if __name__ == '__main__':
    main()
