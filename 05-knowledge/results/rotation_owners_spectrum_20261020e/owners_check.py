"""Primitive external-port kernels for deficiency-three ownership.

H is a distinguished heavy symbol, L=3 represents a light weight at least3.
The formal kernel budget is j; heavy H is never inflated into triangle atoms.
"""
from collections import Counter
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from fractions import Fraction
from math import comb, factorial, prod
from pathlib import Path

spec = spec_from_file_location("old_owner", Path(__file__).resolve().parents[1] /
                              "deficiency_three_20261020/ownership_check.py")
O = module_from_spec(spec)
spec.loader.exec_module(O)
S, need = O.S, O.need
H = 100


def active_endpoints(defects, weights, gaps):
    t = len(weights)
    adjacency = [set() for _ in range(t)]
    for a, b in defects:
        adjacency[a].add(b)
        adjacency[b].add(a)
    active = {v for v in range(t) if weights[v] != 1 or len(adjacency[v]) != 1}
    for g, value in enumerate(gaps):
        if value != H:
            if g:
                active.add(g - 1)
            if g < t:
                active.add(g)
    todo = list(active)
    while todo:
        v = todo.pop()
        for w in adjacency[v] - active:
            active.add(w)
            todo.append(w)
    return active


def formal_record(rows, weights, order, budget):
    backward = [(a, b) for a, b in combinations(range(len(order)), 2)
                if rows[order[b]] >> order[a] & 1]
    if sum(weights[order[a]] * weights[order[b]] for a, b in backward) > budget:
        return None
    endpoints = sorted({v for e in backward for v in e})
    where = {v: j for j, v in enumerate(endpoints)}
    ew = tuple(weights[order[v]] for v in endpoints)
    edges = tuple((where[a], where[b]) for a, b in backward)
    gaps = []
    for a, b in zip([-1] + endpoints, endpoints + [len(order)]):
        need(b - a <= 2, "formal pair-free core has two unmarked vertices")
        gaps.append(weights[order[a + 1]] if b == a + 2 else 0)
    r, t = len(edges), len(endpoints)
    q, delta, u = budget - r, 2 * r - t, sum(w != H for w in gaps)
    need(2 * q + delta + u == 3, "formal deficiency identity")
    return q, delta, u, ew, tuple(gaps), edges


@lru_cache(None)
def formal_owner(rows, weights, budget):
    heavy = [v for v, w in enumerate(weights) if w == H]
    heavy.sort(key=lambda v: -sum(rows[v] >> z & 1 for z in heavy))
    need(len(heavy) == 2 * budget - 2, "formal heavy count")
    gaps = [[] for _ in range(len(heavy) + 1)]
    for v, w in enumerate(weights):
        if w == H:
            continue
        bits = tuple(rows[z] >> v & 1 for z in heavy)
        i = sum(bits)
        need(bits == (1,) * i + (0,) * (len(heavy) - i), "formal intrinsic gap")
        gaps[i].append(v)
    cross = sum(weights[a] * weights[b]
                for i, j in combinations(range(len(gaps)), 2)
                for a in gaps[i] for b in gaps[j] if rows[b] >> a & 1)
    choices = [O.feasible_gap_orders(S.G.induced(rows, g),
                                    tuple(weights[v] for v in g), budget - cross)
               for g in gaps]
    best = None
    for parts in product(*choices):
        if sum(cost for _, cost in parts) > budget - cross:
            continue
        order = []
        for j, (part, _) in enumerate(parts):
            order.extend(gaps[j][v] for v in part)
            if j < len(heavy):
                order.append(heavy[j])
        record = formal_record(rows, weights, order, budget)
        if record is not None and (best is None or record < best):
            best = record
    need(best is not None, "formal owner absent")
    return best


def kernels(j):
    for q in (0, 1):
        r = j - q
        for delta in range(4 - 2 * q):
            u = 3 - 2 * q - delta
            t = 2 * r - delta
            if t < 2:
                continue
            for defect in combinations(tuple(combinations(range(t), 2)), r):
                if {v for e in defect for v in e} != set(range(t)):
                    continue
                degrees = [sum(v in e for e in defect) for v in range(t)]
                weights = [(1,) * t]
                if q:
                    weights += [tuple(1 + int(v == a) for v in range(t))
                                for a, degree in enumerate(degrees) if degree == 1]
                for ew in weights:
                    for positions in combinations(range(t + 1), u):
                        for values in product(range(4), repeat=u):
                            gaps = [H] * (t + 1)
                            for p, v in zip(positions, values):
                                gaps[p] = v
                            gaps = tuple(gaps)
                            if active_endpoints(defect, ew, gaps) != set(range(t)):
                                continue
                            yield (q, delta, u, ew, gaps, defect)


def is_owner(record, j):
    _, _, _, ew, gaps, defect = record
    rows, weights, edges = O.ordered_template(defect, ew, gaps)
    if len(O.interval_chains(len(rows), edges)) != len(rows):
        return False
    return formal_owner(rows, weights, j) == record


def matchings(vertices):
    if not vertices:
        yield ()
        return
    a = vertices[0]
    for i, b in enumerate(vertices[1:], 1):
        for rest in matchings(vertices[1:i] + vertices[i + 1:]):
            yield ((a, b),) + rest


def compositions(total, boxes):
    if boxes == 1:
        yield (total,)
    elif boxes:
        for first in range(total + 1):
            for rest in compositions(total - first, boxes - 1):
                yield (first,) + rest


def completions(record, j, extra):
    q, delta, u, ew, gaps, defect = record
    slots = [i for i, w in enumerate(gaps) if w == H]
    for allocation in compositions(2 * extra, len(slots)):
        counts = dict(zip(slots, allocation))
        new_ew, new_gaps, labels, bulk = [], [], {}, []
        for i, gap in enumerate(gaps):
            for _ in range(counts.get(i, 0)):
                new_gaps.append(H)
                bulk.append(len(new_ew))
                new_ew.append(1)
            new_gaps.append(gap)
            if i < len(ew):
                labels[i] = len(new_ew)
                new_ew.append(ew[i])
        fixed = tuple((labels[a], labels[b]) for a, b in defect)
        for matching in matchings(tuple(bulk)):
            yield (q, delta, u, tuple(new_ew), tuple(new_gaps),
                   tuple(sorted(fixed + matching)))


def count_kernels(j):
    polynomial, total = Counter(), 0
    examples, sampled = [], Counter()
    for record in kernels(j):
        total += 1
        accepted = is_owner(record, j)
        if accepted:
            polynomial[record[4].count(3)] += 1
        category = (record[:3], record[4].count(3), accepted)
        if sampled[category] < 3:
            examples.append((record, accepted))
            sampled[category] += 1
    print("j=", j, "primitive_records=", total,
          "owner_polynomial_coefficients=", dict(sorted(polynomial.items())), flush=True)
    return polynomial, examples


def bulk_controls(examples):
    tested, large_tested = 0, 0
    for j, records in examples:
        for record, accepted in records:
            for extra in (1, 2):
                for child in completions(record, j, extra):
                    need(is_owner(child, j + extra) == accepted,
                         "bulk insertion changes owner acceptance")
                    tested += 1
            positions = [i for i, w in enumerate(record[4]) if w == 3]
            for values in product((3, 4, 7), repeat=len(positions)):
                gaps = list(record[4])
                for p, value in zip(positions, values):
                    gaps[p] = value
                changed = record[:4] + (tuple(gaps), record[5])
                need(is_owner(changed, j) == accepted, "large light values change owner acceptance")
                large_tested += 1
    print("Exact one/two-edge bulk completion controls:", tested,
          "large-light substitution controls:", large_tested, flush=True)


def main():
    results, examples = {}, []
    for j in (2, 3):
        results[j], sample = count_kernels(j)
        examples.append((j, sample))
    need(results[2][0] == 756, "known k2 D3 count")
    N3 = 3 * sum(results[2].values()) + sum(results[3].values())
    need(N3 == 36473, "known k3 D3 count from finite kernels")
    print("N3(2)=", results[2][0], "N3(3)=", N3, flush=True)
    bulk_controls(examples)
    leading = 0
    for positions in combinations(range(11), 3):
        marked = tuple(sorted({v for p in positions for v in (p, p + 1)}))
        if len(marked) != 6:
            continue
        other = tuple(v for v in range(12) if v not in marked)
        gaps = tuple(3 if g - 1 in positions else H for g in range(13))
        for partners in (other, other[::-1]):
            defects = tuple(sorted(tuple(sorted((a, b))) for a, b in zip(marked, partners)))
            record = (0, 0, 3, (1,) * 12, gaps, defects)
            need(active_endpoints(defects, record[3], gaps) == set(range(12)),
                 "leading kernel activated ports")
            need(is_owner(record, 6), "leading kernel unique-owner claim")
            leading += 1
    need(leading == 168 and comb(9, 3) * factorial(6) == 60480,
         "leading coefficient finite count")
    need(Fraction(16 * 60480, factorial(9)) == Fraction(8, 3),
         "normalized degree-seven leading coefficient")
    for j in range(2, 7):
        for ell in range(4):
            for start in range(5):
                difference = sum((-1) ** (8 - i) * comb(8, i) *
                                 prod(start + i - d for d in range(j - 2)) *
                                 (start + i) ** ell for i in range(9))
                need(difference == 0, "universal eighth-difference basis identity")
    print("Leading-kernel port controls:", leading,
          "; c63=60480, normalized leading coefficient=8/3; 100 recurrence-basis identities.")
    print("PASS: primitive kernel counts and port-preserving bulk ownership controls.")


if __name__ == "__main__":
    main()
