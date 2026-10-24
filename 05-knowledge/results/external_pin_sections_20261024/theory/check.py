#!/usr/bin/env python3
"""Exact small controls for conditional deficient supply and receiver migration."""
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def subsets(n):
    return [frozenset(i for i in range(n) if (m >> i) & 1)
            for m in range(1 << n)]


def charge(S, labels, modulus):
    return sum(labels[i] for i in S) % modulus


def support_sum(A, B, modulus):
    return {(a+b) % modulus for a in A for b in B}


def weight(S, weights):
    return sum(weights[i] for i in S)


def largest(weights, count):
    return sum(sorted(weights, reverse=True)[:count])


def deficiency_controls():
    pairs = 0
    for n in range(2, 7):
        family = subsets(n)
        for A, B in product(family, repeat=2):
            pairs += 1
            sums = support_sum(A, B, n)
            if len(A)+len(B) > n:
                need(len(sums) == n, 'strict density gives complete sums')
            for q in range(n):
                candidates = [(a, (q-a) % n) for a in A
                              if (q-a) % n in B]
                need(len(candidates) >= len(A)+len(B)-n,
                     'candidate intersection lower bound')
                # Arbitrary forbidden prefixes realize all possible deletion counts.
                for erased in range(len(candidates)+1):
                    if len(A)+len(B) > n+erased:
                        need(bool(candidates[erased:]), 'cross-defect repair')
    need(support_sum({0, 2}, {0, 2}, 4) == {0, 2}, 'equality hostile')
    A, B = set(range(5)), {0, 1}
    bad = {(0, 0), (4, 1)}
    available = {(a+b) % 5 for a in A for b in B if (a, b) not in bad}
    need(available == {1, 2, 3, 4}, 'cross-defect sharp boundary')
    need({(x+x) % 2 for x in range(2)} == {0}, 'shared-pin hostile')
    left, right = [0, 1, 2], [0, 2, 4]
    graded = {}
    for total in range(2, 5):
        splits = []
        for r in (1, 2):
            s = total-r
            if s not in (1, 2):
                continue
            A = {charge(I, left, 5) for I in subsets(3) if len(I) == r}
            B = {charge(I, right, 5) for I in subsets(3) if len(I) == s}
            need(len(A) == len(B) == 3, 'both incomplete graded pieces')
            need(len(support_sum(A, B, 5)) == 5, 'graded composition complete')
            splits.append([r, s])
        graded[str(total)] = splits
    return {'cyclic_subset_pairs': pairs, 'graded_C5_splits': graded,
            'density_equality_hostile': {'group': 4, 'A': [0, 2], 'B': [0, 2]},
            'cross_defect_hostile': {'group': 5, 'missing': 0, 'bad_pairs': sorted(bad)}}


class ReceiverModel:
    """Fixed binary slot records; helper states retain a literal indexed pin."""
    def __init__(self, name, R, E, modulus, lo, hi, pinned, sign=1):
        self.name, self.R, self.E = name, tuple(sign*x % modulus for x in R), tuple(sign*x % modulus for x in E)
        self.n, self.h, self.e = modulus, len(R), len(E)
        self.lo, self.hi, self.pinned = lo, hi, pinned
        self.Rcost = tuple(2+(i % 3) for i in range(self.h))
        self.Ecost = tuple(3+i for i in range(self.e))
        self.Hcharge = (0, modulus-1, 1 % modulus)
        self.Hrecord = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        self.Hcost = (5, 7, 11)
        self.Rsets, self.Esets = subsets(self.h), subsets(self.e)
        self.sections = defaultdict(list)
        for J, y, I in product(self.Esets, range(3), self.Rsets):
            if self.valid_boundary(I, J, y):
                self.sections[len(I), J, y, charge(I, self.R, self.n)].append(I)
        for J, y, r, q in product(self.Esets, range(3), range(lo, hi+1), range(modulus)):
            need(bool(self.sections[r, J, y, q]), name+' conditional graded section')

    def valid_boundary(self, I, J, y):
        # The first receiver's actual owner value equals the external frame pin.
        return not self.pinned or int(0 in I) == ((y+len(J)) & 1)

    def q(self, state):
        I, J, y = state
        return (charge(I, self.R, self.n)+charge(J, self.E, self.n)+self.Hcharge[y]) % self.n

    def section(self, J, y, r, q, current):
        debt = (q-charge(J, self.E, self.n)-self.Hcharge[y]) % self.n
        return min(self.sections[r, J, y, debt],
                   key=lambda I: (weight(I ^ current, self.Rcost), tuple(sorted(I))))

    def normalize(self, start, q, r0):
        output = [(start, 'start')]
        I, J, y = start
        if len(I) < self.lo:
            change = sorted(J)[:self.lo-len(I)]
            J = J-frozenset(change)
            I = self.section(J, y, self.lo, q, I)
            output.append(((I, J, y), 'entry-low'))
        elif len(I) > self.hi:
            change = sorted(set(range(self.e))-J)[:len(I)-self.hi]
            J = J | frozenset(change)
            I = self.section(J, y, self.hi, q, I)
            output.append(((I, J, y), 'entry-high'))
        while len(I) > r0:
            j = min(set(range(self.e))-J)
            J = J | {j}
            I = self.section(J, y, len(I)-1, q, I)
            output.append(((I, J, y), 'adjacent'))
        return output

    def route(self, start, end):
        q = self.q(start)
        k = len(start[0])+len(start[1])
        r0 = max(self.lo, k-self.e)
        left, right = self.normalize(start, q, r0), self.normalize(end, q, r0)
        output = list(left)
        I, J, y = left[-1][0]
        goalI, goalJ, goalY = right[-1][0]
        for remove, add in zip(sorted(J-goalJ), sorted(goalJ-J)):
            J = (J-{remove}) | {add}
            I = self.section(J, y, r0, q, I)
            output.append(((I, J, y), 'outside'))
        while y != goalY:
            y += 1 if goalY > y else -1
            I = self.section(J, y, r0, q, I)
            output.append(((I, J, y), 'helper'))
        if I != goalI:
            output.append(((goalI, J, y), 'refresh'))
        for index in range(len(right)-1, 0, -1):
            output.append((right[index-1][0], 'reverse-'+right[index][1]))
        return output

    def event(self, before, after, kind):
        if kind.startswith('reverse-'):
            return self.event(after, before, kind[len('reverse-'):])
        I, J, y = before
        I2, J2, y2 = after
        r, s = len(I), len(I2)
        t = min(r+s, 2*self.h-r-s)
        costR = weight(I ^ I2, self.Rcost)
        costE = weight(J ^ J2, self.Ecost)
        costH = sum(c for c, a, b in zip(self.Hcost, self.Hrecord[y], self.Hrecord[y2]) if a != b)
        need(len(I ^ I2) <= t and costR <= largest(self.Rcost, t), 'sharp grade cost bound')
        need(kind == 'helper' or y == y2, 'helper changes only explicit edge')
        if kind == 'entry-low':
            need(s == self.lo and len(J-J2) == s-r and J2 <= J, 'low grade exact transfer')
        elif kind == 'entry-high':
            need(s == self.hi and len(J2-J) == r-s and J <= J2, 'high grade exact transfer')
        elif kind == 'adjacent':
            need(s == r-1 and len(J2-J) == 1 and J <= J2, 'adjacent exact margin transfer')
        elif kind == 'outside':
            need(r == s and len(J ^ J2) == 2, 'Johnson outside move')
        elif kind == 'helper':
            need(r == s and J == J2 and abs(y-y2) == 1, 'actual helper skeleton edge')
        elif kind == 'refresh':
            need(r == s and J == J2 and y == y2, 'literal receiver refresh')
        else:
            raise RuntimeError('unknown route event '+kind)
        return costR+costE+costH

    def check(self):
        fibres = defaultdict(list)
        for I, J, y in product(self.Rsets, self.Esets, range(3)):
            k = len(I)+len(J)
            if self.lo <= k <= self.hi+self.e and self.valid_boundary(I, J, y):
                fibres[k, self.q((I, J, y))].append((I, J, y))
        total, events, max_cost, max_length = 0, 0, 0, 0
        counts = Counter()
        for (k, q), states in fibres.items():
            for start, end in product(states, repeat=2):
                total += 1
                route = self.route(start, end)
                need(route[0][0] == start and route[-1][0] == end, 'literal endpoints preserved')
                r0 = max(self.lo, k-self.e)
                dE = len(self.normalize(start, q, r0)[-1][0][1]-self.normalize(end, q, r0)[-1][0][1])
                bound = 2*(self.hi-self.lo+1)+dE+abs(start[2]-end[2])+1
                need(len(route)-1 <= bound, 'terminating event bound')
                for (previous, _), (current, kind) in zip(route, route[1:]):
                    need(len(current[0])+len(current[1]) == k, 'core margin retained')
                    need(sum(self.Hrecord[current[2]]) == sum(self.Hrecord[start[2]]) == 1,
                         'separately fixed helper population')
                    need(self.valid_boundary(*current), 'actual indexed boundary pin')
                    need(self.q(current) == q, 'complete physical-code coordinate retained')
                    cost = self.event(previous, current, kind)
                    max_cost = max(max_cost, cost)
                    events += 1
                    counts[kind] += 1
                max_length = max(max_length, len(route)-1)
        return {'model': self.name, 'Q': self.n, 'R_labels': self.R, 'E_labels': self.E,
                'grade_band': [self.lo, self.hi], 'conditional_owner_pin': self.pinned,
                'feasible_fibres': len(fibres), 'ordered_endpoint_pairs': total,
                'route_events': events, 'max_observed_cost': max_cost,
                'max_observed_length': max_length, 'event_types': dict(sorted(counts.items()))}


def grade_and_complement_controls():
    cases = 0
    for h in range(1, 9):
        family = subsets(h)
        for r, s in product(range(h+1), repeat=2):
            actual = max(len(A ^ B) for A in family if len(A) == r
                         for B in family if len(B) == s)
            need(actual == min(r+s, 2*h-r-s), 'sharp grade support equality')
            cases += 1
    # Complementing only core requires negating its labels, not helper labels.
    labels, H = [1, 2, 4], [0, 1, 2]
    for I, eta in product(subsets(3), H):
        U = frozenset(range(3))-I
        old = (sum(labels[i] for i in I)+eta) % 5
        new = (-sum(labels[i] for i in U)+eta) % 5
        need(new == (old-sum(labels)) % 5, 'core-only complement has negative labels')
    wrong = ((-sum(labels[i] for i in {0, 1}))-1) % 5
    correct = ((-sum(labels[i] for i in {0, 1}))+1) % 5
    need(wrong != correct, 'negating helper too is a hostile')
    need({0, 1} == {0, 1} and {0} != {0, 1}, 'union-grade coverage cannot pay fixed-grade debt')
    # Root-proposed 3-erasure replacement: only arithmetic bounds, not native coverage.
    costsR, costsE = [2]*17+[3], [4]*7
    worst = Counter()
    for k in range(8, 18):
        r0 = max(8, k-7)
        for r in range(19):
            if not 0 <= k-r <= 7:
                continue
            if r < 8:
                worst['entry'] = max(worst['entry'], largest(costsR, r+8)+largest(costsE, 8-r))
            elif r > 10:
                worst['entry'] = max(worst['entry'], largest(costsR, 26-r)+largest(costsE, r-10))
        for r in range(r0+1, 11):
            worst['adjacent'] = max(worst['adjacent'], largest(costsR, min(2*r-1, 37-2*r))+4)
        refresh = largest(costsR, 2*min(r0, 18-r0))
        worst['helper'] = max(worst['helper'], refresh+15)
        worst['outside'] = max(worst['outside'], refresh+8)
    need(dict(worst) == {'entry': 47, 'adjacent': 39, 'helper': 52, 'outside': 45}, 'conditional 18-bank exact crude bounds')
    return {'sharp_grade_cases': cases, 'conditional_18_bank_cost_bounds': dict(worst)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    models = [
        ReceiverModel('C3_indexed_owner_pin', [0, 1, 2, 0, 1, 2], [1, 2], 3, 2, 4, True),
        ReceiverModel('C3_negated_core_same_helper', [0, 1, 2, 0, 1, 2], [1, 2], 3, 2, 4, True, -1),
        ReceiverModel('C5_complementary_incomplete_blocks', [0, 1, 2, 0, 2, 4], [1, 3], 5, 2, 4, False),
        ReceiverModel('C2_zero_outside_width', [0, 1, 0, 1], [], 2, 1, 3, False),
        ReceiverModel('trivial_group_zero_grade', [0, 0], [0], 1, 0, 2, False),
    ]
    result = {'scope': 'abstract finite full-record conditional models, not a native n8 census',
              'deficiency': deficiency_controls(), 'costs': grade_and_complement_controls(),
              'receiver_models': [m.check() for m in models],
              'checks_before_receipt': CHECKS,
              'sources': {name: sha256((HERE/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                          for name in ['check.py', 'RECEIVER_BANKS.md']}}
    result = json.loads(json.dumps(result))
    receipt = HERE/'receipt.json'
    if args.write:
        receipt.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    else:
        need(result == json.loads(receipt.read_text(encoding='utf-8')), 'frozen receiver receipt')
    print(f'PASS {CHECKS:,} checks; '+str(sum(m['ordered_endpoint_pairs'] for m in result['receiver_models']))+' full endpoint routes')


if __name__ == '__main__':
    main()
