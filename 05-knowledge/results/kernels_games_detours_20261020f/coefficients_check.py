"""Exact D3 kernel coefficients using component grammar and maximal-activation shortcuts."""
from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import comb, factorial, prod
from pathlib import Path
import sys

spec = spec_from_file_location("inherited_port_kernels", Path(__file__).resolve().parents[1] /
                              "rotation_owners_spectrum_20261020e/owners_check.py")
K = module_from_spec(spec)
spec.loader.exec_module(K)
O, S, H, need = K.O, K.S, K.H, K.need


def connected(v, edges):
    seen = {0}
    while True:
        grown = seen | {b for a, b in edges if a in seen} | {a for a, b in edges if b in seen}
        if grown == seen:
            return len(seen) == v
        seen = grown


TYPES = []
for size, edges, charge in ((2, 1, 0), (3, 2, 1), (4, 3, 2), (5, 4, 3), (3, 3, 3)):
    patterns = tuple(e for e in combinations(tuple(combinations(range(size), 2)), edges)
                     if connected(size, e))
    TYPES.append((size, edges, charge, patterns))


def graphs(vertices, edges, charge):
    """Unique component partition: the next component contains the first unused vertex."""
    if not vertices:
        if edges == charge == 0:
            yield ()
        return
    if edges <= 0 or charge < 0:
        return
    for size, count, delta, patterns in TYPES:
        if size > len(vertices) or count > edges or delta > charge:
            continue
        for others in combinations(vertices[1:], size - 1):
            chosen = (vertices[0],) + others
            remaining = tuple(v for v in vertices if v not in chosen)
            for rest in graphs(remaining, edges - count, charge - delta):
                for pattern in patterns:
                    yield tuple(sorted(rest + tuple((chosen[a], chosen[b]) for a, b in pattern)))


def owner_without_global_cache(record, j):
    _, _, _, ew, gaps, defects = record
    rows, weights, edges = O.ordered_template(defects, ew, gaps)
    if len(O.interval_chains(len(rows), edges)) != len(rows):
        return False
    return K.formal_owner.__wrapped__(rows, weights, j) == record


def universe_controls():
    """Independent all-edge-subsets comparison, including absence of duplicates."""
    checked, totals = 0, 0
    for e in range(1, 5):
        for delta in range(4):
            t = 2 * e - delta
            if t < 2:
                continue
            generated = list(graphs(tuple(range(t)), e, delta))
            actual = {edges for edges in combinations(tuple(combinations(range(t), 2)), e)
                      if {v for edge in edges for v in edge} == set(range(t))}
            need(len(generated) == len(set(generated)), "duplicate component-grammar graph")
            need(set(generated) == actual, "component grammar differs from raw edge subsets")
            checked += 1
            totals += len(actual)
    print("independent edge-subset universes=", checked, "graphs=", totals, flush=True)


def exact(j):
    polynomial, records, graph_count = Counter(), 0, 0
    by_signature = Counter()
    for q in (0, 1):
        e = j - q
        for delta in range(4 - 2 * q):
            u, t = 3 - 2 * q - delta, 2 * e - delta
            if t < 2:
                continue
            for defects in graphs(tuple(range(t)), e, delta):
                graph_count += 1
                degree = [sum(v in edge for edge in defects) for v in range(t)]
                endpoint_weights = [(1,) * t]
                if q:
                    endpoint_weights += [tuple(1 + int(v == a) for v in range(t))
                                         for a, d in enumerate(degree) if d == 1]
                for ew in endpoint_weights:
                    for positions in combinations(range(t + 1), u):
                        first = tuple(0 if i in positions else H for i in range(t + 1))
                        if K.active_endpoints(defects, ew, first) != set(range(t)):
                            continue
                        for values in product(range(4), repeat=u):
                            gap = [H] * (t + 1)
                            for p, value in zip(positions, values):
                                gap[p] = value
                            record = (q, delta, u, ew, tuple(gap), defects)
                            records += 1
                            if owner_without_global_cache(record, j):
                                polynomial[values.count(3)] += 1
                                by_signature[q, delta, u, values.count(3)] += 1
    print("exact j=", j, "component_graphs=", graph_count, "primitive_records=", records,
          "polynomial=", tuple(polynomial[i] for i in range(4)), flush=True)
    print("signature_coefficients=", dict(sorted(by_signature.items())), flush=True)
    return tuple(polynomial[i] for i in range(4))


def product_linear(a, b):
    # (3+y)^a (2+y)^b: the b exceptional matched gaps must be occupied.
    p = [1]
    for constant in [3] * a + [2] * b:
        nxt = [0] * (len(p) + 1)
        for i, value in enumerate(p):
            nxt[i] += constant * value
            nxt[i + 1] += value
        p = nxt
    return p + [0] * (4 - len(p))


def shortcut_five():
    counts, polynomial, audits = Counter(), [0] * 4, []
    for delta in range(4):
        e, u, t = 5, 3 - delta, 10 - delta
        for defects in graphs(tuple(range(t)), e, delta):
            degree = [sum(v in edge for edge in defects) for v in range(t)]
            isolated = [edge for edge in defects if all(degree[v] == 1 for v in edge)]
            masks = [sum(1 << v for v in edge) for edge in isolated]
            for positions in combinations(range(t + 1), u):
                marked = 0
                for p in positions:
                    if p:
                        marked |= 1 << (p - 1)
                    if p < t:
                        marked |= 1 << p
                if any(not mask & marked for mask in masks):
                    continue
                bad_positions = [p for p in positions if 0 < p < t and
                                 (p - 1, p) in isolated]
                b = len(bad_positions)
                need(b <= 1, "j5 has at most one internal matched gap")
                counts[delta, u, b] += 1
                contribution = product_linear(u - b, b)
                polynomial = [x + y for x, y in zip(polynomial, contribution)]
                category = (delta, b)
                if sum(item[0] == category for item in audits) < 8:
                    audits.append((category, defects, positions, bad_positions, t, u))
    tested = 0
    for _, defects, positions, bad, t, u in audits:
        for values in product(range(4), repeat=u):
            gap = [H] * (t + 1)
            for p, value in zip(positions, values):
                gap[p] = value
            record = (0, 3 - u, u, (1,) * t, tuple(gap), defects)
            expected = all(gap[p] != 0 for p in bad)
            need(owner_without_global_cache(record, 5) == expected,
                 "j5 shortcut disagrees with full exact owner")
            tested += 1
    print("j5 activation sectors (delta,u,matched_gap)=", dict(sorted(counts.items())), flush=True)
    need(tuple(polynomial) == tuple(210 * x for x in multiply(multiply([7, 2], [19, 6]), [32, 9])),
         "j5 factorization")
    print("shortcut j=5 polynomial=", tuple(polynomial), "full_owner_controls=", tested, flush=True)
    return tuple(polynomial)


def shortcut_six():
    # At maximal activation every marked matching endpoint has an external partner.
    sectors = (comb(9, 3) * factorial(6),
               comb(9, 2) * comb(7, 3) * 3 * factorial(4),
               9 * comb(8, 6) * 90 * factorial(2),
               factorial(9) // (factorial(3) ** 3 * factorial(3)) * 3 ** 3)
    need(sectors == (60480, 90720, 45360, 7560), "j6 placement counts")
    polynomial = [0] * 4
    for u, count in zip((3, 2, 1, 0), sectors):
        polynomial = [x + count * y for x, y in zip(polynomial, product_linear(u, 0))]
    need(polynomial == [7560 * comb(3, i) * 7 ** (3 - i) * 2 ** i for i in range(4)],
         "j6 factorization")
    tested = 0
    for delta in range(4):
        u, t = 3 - delta, 12 - delta
        marked = tuple(range(3 * delta, 3 * delta + 2 * u))
        partners = tuple(range(3 * delta + 2 * u, t))
        positions = tuple(3 * delta + 2 * i + 1 for i in range(u))
        for centers in product(range(3), repeat=delta):
            nonmatching = tuple(sorted((3 * i + min(c, v), 3 * i + max(c, v))
                                       for i, c in enumerate(centers) for v in range(3) if v != c))
            for index, assignment in enumerate(permutations(partners)):
                if delta == 0 and index >= 24:
                    break
                defects = tuple(sorted(nonmatching + tuple(zip(marked, assignment))))
                for values in product(range(4), repeat=u):
                    gap = [H] * (t + 1)
                    for p, value in zip(positions, values):
                        gap[p] = value
                    ew = (1,) * t
                    need(K.active_endpoints(defects, ew, gap) == set(range(t)), "j6 audit activation")
                    record = (0, delta, u, ew, tuple(gap), defects)
                    need(owner_without_global_cache(record, 6), "j6 shortcut disagrees with full owner")
                    tested += 1
    print("shortcut j=6 sector_counts=", sectors, "polynomial=", tuple(polynomial),
          "full_owner_controls=", tested, flush=True)
    return tuple(polynomial)


def multiply(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def double_factorial(n):
    return prod(range(n, 0, -2))


def algebra_controls(results):
    p = [Fraction(0)] * 8
    for j, coefficients in sorted(results.items()):
        falling = [Fraction(1)]
        for a in range(j - 2):
            falling = multiply(falling, [-a, 1])
        factor = Fraction(2 ** (j - 2), factorial(2 * j - 3))
        for i, value in enumerate(multiply(falling, coefficients)):
            p[i] += factor * value
    expected = [Fraction(v, 6) for v in (4536, 15953, 23395, 18378, 8260, 2120, 288, 16)]
    need(p == expected, "expanded P coefficient algebra")
    shifted = [sum(p[i] * comb(i, a) * (-2) ** (i - a) for i in range(a, 8))
               for a in range(8)]
    need(shifted == [Fraction(v, 6) for v in (-110, 61, 135, -22, -140, 8, 64, 16)],
         "shifted N3 polynomial")
    values = []
    for k in range(2, 41):
        x = k - 2
        kernel = sum(sum(c * x ** i for i, c in enumerate(results[j])) *
                     comb(2 * k - 3, 2 * j - 3) * double_factorial(2 * k - 2 * j - 1)
                     for j in range(2, min(6, k) + 1))
        closed = double_factorial(2 * k - 3) * sum(c * k ** i for i, c in enumerate(shifted))
        need(kernel == closed and closed.denominator == 1, "N3 exact kernel/closed algebra")
        if k <= 10:
            values.append(int(closed))
        # Independent all-L matching universe, before primitive extraction.
        top = sum(results[j][3] * comb(2 * k - 3, 2 * j - 3) *
                  double_factorial(2 * k - 2 * j - 1) for j in range(2, min(6, k) + 1))
        need(top == double_factorial(2 * k - 1) * comb(2 * k + 1, 3),
             "all-L leading-stratum decomposition")
    need(values[:2] == [756, 36473], "independent inherited full-class counts")
    need(sum(shifted) == 2 and sum(shifted) != 5, "exceptional k1 boundary retained")
    print("P(x) ascending coefficients=", tuple(str(c) for c in p), flush=True)
    print("6*P(k-2) ascending coefficients=", tuple(int(6 * c) for c in shifted), flush=True)
    print("N3(k), k=1..10=", (5,) + tuple(values), flush=True)
    print("exact kernel/closed and all-L matching controls k=2..40; k1 exception retained", flush=True)


def main():
    if sys.argv[1:] == ["audit-five"]:
        need(exact(5) == shortcut_five(), "complete j5 owner audit")
        print("PASS complete j5 independent full-owner comparison", flush=True)
        return
    selected = list(map(int, sys.argv[1:])) or [2, 3, 4, 5, 6]
    universe_controls()
    results = {}
    for j in selected:
        results[j] = exact(j) if j <= 4 else shortcut_five() if j == 5 else shortcut_six()
    if 2 in results:
        need(results[2] == (756, 792, 270, 30), "inherited C2")
    if 3 in results:
        need(results[3] == (13380, 12974, 4140, 435), "inherited C3")
    if set(results) == set(range(2, 7)):
        algebra_controls(results)
    print("PASS", flush=True)


if __name__ == "__main__":
    main()
