"""Independent top-two-strata template generation for feedback budget k.

The A/B/C generators are specialized to the proved structural classification.
They reuse only the general exact clone/canonicalization routines for checking.
The complete finite compiler at k=1,2 is an independent completeness control.
"""
from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import factorial
from pathlib import Path


def need(value, message):
    if not value:
        raise RuntimeError(message)


spec = spec_from_file_location("next_sparse", Path(__file__).with_name("sparse_count.py"))
S = module_from_spec(spec)
spec.loader.exec_module(S)


def matchings(vertices):
    if not vertices:
        yield ()
        return
    a = vertices[0]
    for b in vertices[1:]:
        rest = tuple(v for v in vertices if v not in (a, b))
        for matching in matchings(rest):
            yield ((a, b),) + matching


def template(defect, gaps, k):
    """Unit-weight endpoints, prescribed gap weights, forward reference order."""
    t = len(gaps) - 1
    slots, weights = {}, []
    for j in range(t + 1):
        if gaps[j]:
            weights.append(gaps[j])
        if j < t:
            slots[j] = len(weights)
            weights.append(1)
    n = len(weights)
    rows = [(1 << n) - (1 << (j + 1)) for j in range(n)]
    for a, b in defect:
        a, b = sorted((slots[a], slots[b]))
        rows[a] ^= 1 << b
        rows[b] ^= 1 << a
    return S.normalize(tuple(rows), weights, k)


def classify(k):
    K = k + 1
    top, A, B, C = set(), set(), set(), set()
    raw = Counter()
    excluded = set()
    for matching in matchings(tuple(range(2 * k))):
        edges = {frozenset(edge) for edge in matching}
        top.add(template(matching, (K,) * (2 * k + 1), k))
        raw["top"] += 1
        for gap in range(2 * k + 1):
            for weight in range(1, k + 1):
                gaps = [K] * (2 * k + 1)
                gaps[gap] = weight
                A.add(template(matching, gaps, k))
                raw["A"] += 1
            gaps = [K] * (2 * k + 1)
            gaps[gap] = 0
            state = template(matching, gaps, k)
            if 0 < gap < 2 * k and frozenset((gap - 1, gap)) in edges:
                need(state[0] == 4 * k - 3, "excluded B canonical core size")
                need(state[2].count(K) == 2 * k - 1, "excluded B heavy count")
                excluded.add(state)
                raw["excluded B"] += 1
            else:
                B.add(state)
                raw["B"] += 1
    vertices = tuple(range(2 * k - 1))
    for triple in combinations(vertices, 3):
        rest = tuple(v for v in vertices if v not in triple)
        for center in triple:
            path = tuple((center, leaf) for leaf in triple if leaf != center)
            for matching in matchings(rest):
                C.add(template(path + matching, (K,) * (2 * k), k))
                raw["C"] += 1
    M = factorial(2 * k) // (2 ** k * factorial(k))
    expected = {"top": M, "A": k * (2 * k + 1) * M,
                "B": 2 * k * M, "C": (k - 1) * M}
    families = {"top": top, "A": A, "B": B, "C": C}
    for name, states in families.items():
        need(len(states) == raw[name] == expected[name], f"{name} injective count")
        q = {"top": 4 * k + 1, "A": 4 * k + 1,
             "B": 4 * k, "C": 4 * k - 1}[name]
        h = 2 * k + int(name == "top")
        need(all(state[0] == q and state[2].count(K) == h for state in states),
             f"{name} canonical q/heavy stratum")
    need(raw["excluded B"] == M, "summed adjacent matching exclusions")
    need(len(top | A | B | C) == sum(map(len, (top, A, B, C))),
         "families overlap")
    need(not excluded & (top | A | B | C), "excluded B remains in top strata")
    L = len(A | B | C)
    need(L == M * (2 * k ** 2 + 4 * k - 1), "next stratum formula")
    shift = 2 * k ** 2 + 5 * k + 1
    need(all(sum(state[2]) == shift for state in top), "top shift")
    leading = Fraction(M, factorial(2 * k))
    next_coefficient = Fraction(L, factorial(2 * k - 1)) + Fraction(
        M * (2 * k + 1 - 2 * shift), 2 * factorial(2 * k - 1))
    need(next_coefficient == Fraction(-3 * k, 2 ** k * factorial(k)),
         "universal next polynomial coefficient")
    if k <= 2:
        complete, _ = S.template_states(k)
        filtered = {state for state in complete if state[2].count(K) >= 2 * k}
        need(filtered == top | A | B | C, "independent complete compiler stratum")
    print("k=", k, "top/A/B/C=", tuple(len(families[name])
                                       for name in ("top", "A", "B", "C")),
          "excluded_raw=", raw["excluded B"], "h=2k total=", L,
          "leading=", str(leading), "next=", str(next_coefficient))


def main():
    for k in range(1, 5):
        classify(k)
    # The symmetric internal A cluster is a positive hostile control: its
    # three rotations survive as automorphisms but still give one cap orbit.
    q, code, weights = template(((0, 1),), (2, 1, 2), 1)
    rows = S.E.decode(q, code)
    best, orders = S.canonical_orders(rows)
    need(best == code, "canonical code in symmetry control")
    stabilizer = [order for order in orders
                  if tuple(weights[v] for v in order) == weights]
    need(len(stabilizer) == 3, "internal A C3 has exactly three rotations")
    print("Internal A directed-triangle symmetry: stabilizer order 3, one template.")
    print("Complete finite compiler strata independently agree at k=1,2.")
    print("PASS (explicit exceptions active under normal and -O Python)")


if __name__ == "__main__":
    main()
