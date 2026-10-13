"""Exact Cayley closing-relation certificates and algebraic positive controls.

No finite-field colouring bound is promoted to a Euclidean colouring bound.
The production finite search classifies a supplied finite abelian Cayley source;
the small coordinate search is explicitly budgeted.  All checks survive -O.
"""
from collections import deque
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, isqrt, lcm, prod
import json


def need(test, message):
    if not test:
        raise ValueError(message)


def natural(n, minimum=0):
    need(type(n) is int and n >= minimum, "exact integer outside domain")


def add(x, y, moduli):
    return tuple((a + b) % m for a, b, m in zip(x, y, moduli))


def scale(k, x, moduli):
    return tuple(k * a % m for a, m in zip(x, moduli))


def vector(x, moduli):
    need(type(x) is tuple and len(x) == len(moduli), "wrong vector type/length")
    need(all(type(a) is int and 0 <= a < m for a, m in zip(x, moduli)),
         "noncanonical exact group coordinates")


@dataclass(frozen=True)
class Cayley:
    moduli: tuple
    directions: tuple  # one representative of each inverse pair

    def __post_init__(self):
        need(type(self.moduli) is tuple and self.moduli, "nonempty moduli tuple")
        for m in self.moduli:
            natural(m, 2)
        need(type(self.directions) is tuple, "directions tuple required")
        for s in self.directions:
            vector(s, self.moduli)
            need(any(s), "zero direction")
            need(s <= scale(-1, s, self.moduli), "canonical inverse-pair representative")
        need(tuple(sorted(set(self.directions))) == self.directions,
             "sorted distinct directions required")

    @property
    def zero(self):
        return (0,) * len(self.moduli)


def make_cayley(moduli, generators):
    need(type(moduli) is tuple and moduli, "moduli tuple required")
    for m in moduli:
        natural(m, 2)
    need(type(generators) is tuple, "generator tuple required")
    for s in generators:
        vector(s, moduli)
        need(any(s), "zero generator")
    return Cayley(moduli, tuple(sorted({min(s, scale(-1, s, moduli))
                                       for s in generators})))


def order(s, moduli):
    return lcm(*(m // gcd(a, m) for a, m in zip(s, moduli)))


def subgroup_words(moduli, generators):
    """BFS words are finite exact membership witnesses; no geometric oracle."""
    zero = (0,) * len(moduli)
    steps = tuple(sorted(set(generators) | {scale(-1, s, moduli) for s in generators}))
    words = {zero: ()}
    todo = deque([zero])
    while todo:
        x = todo.popleft()
        for s in steps:
            y = add(x, s, moduli)
            if y not in words:
                words[y] = words[x] + (s,)
                todo.append(y)
    return words


@dataclass(frozen=True)
class Collapse:
    source: Cayley
    direction: tuple
    order: int
    period: int
    other_word: tuple


def audit_collapse(source, receipt):
    need(type(source) is Cayley and type(receipt) is Collapse, "typed source/receipt")
    need(type(receipt.source) is Cayley and receipt.source == source, "source replacement")
    vector(receipt.direction, source.moduli)
    natural(receipt.order, 2)
    natural(receipt.period, 1)
    need(receipt.direction in source.directions, "missing marked direction")
    need(receipt.order == order(receipt.direction, source.moduli), "wrong cycle order")
    need(receipt.period < receipt.order and receipt.order % receipt.period == 0,
         "proper period dividing cycle order required")
    need(type(receipt.other_word) is tuple, "exact relation word required")
    others = set(source.directions) - {receipt.direction}
    others |= {scale(-1, s, source.moduli) for s in tuple(others)}
    total = source.zero
    for t in receipt.other_word:
        vector(t, source.moduli)
        need(t in others, "relation uses the marked inverse pair")
        # The four source vertices for every commutative square are distinct.
        need(len({source.zero, receipt.direction, t,
                  add(receipt.direction, t, source.moduli)}) == 4,
             "degenerate source square")
        total = add(total, t, source.moduli)
    need(total == scale(receipt.period, receipt.direction, source.moduli),
         "relation does not close the marked period")
    need(total != source.zero, "collapse must identify distinct source vertices")
    return True


@dataclass(frozen=True)
class Product:
    source: Cayley
    factors: tuple
    # Entries: (cyclic factor coordinate, original-source vertex).
    coordinate_map: tuple


def audit_product(source, receipt):
    need(type(source) is Cayley and type(receipt) is Product, "typed source/receipt")
    need(type(receipt.source) is Cayley and receipt.source == source, "source replacement")
    need(type(receipt.factors) is tuple, "factor tuple required")
    for n in receipt.factors:
        natural(n, 2)
    need(receipt.factors == tuple(order(s, source.moduli) for s in source.directions),
         "wrong cyclic factors")
    need(type(receipt.coordinate_map) is tuple, "marked map tuple required")
    expected = []
    for c in product(*(range(n) for n in receipt.factors)):
        x = source.zero
        for a, s in zip(c, source.directions):
            x = add(x, scale(a, s, source.moduli), source.moduli)
        expected.append((c, x))
    # Validate before equality: bool/float tuple aliases are not receipts.
    for item in receipt.coordinate_map:
        need(type(item) is tuple and len(item) == 2, "coordinate entry type")
        vector(item[0], receipt.factors)
        vector(item[1], source.moduli)
    need(receipt.coordinate_map == tuple(expected), "incorrect complete coordinate map")
    need(len({x for _, x in expected}) == len(expected), "cyclic factors overlap")
    need({x for _, x in expected} == set(subgroup_words(source.moduli, source.directions)),
         "generated subgroup not covered")
    return True


def classify(source):
    """Total finite classification; exhaustive subgroup arithmetic only."""
    need(type(source) is Cayley, "typed Cayley source required")
    for s in source.directions:
        words = subgroup_words(source.moduli, tuple(t for t in source.directions if t != s))
        n = order(s, source.moduli)
        d = next(j for j in range(1, n + 1) if scale(j, s, source.moduli) in words)
        if d < n:
            out = Collapse(source, s, n, d, words[scale(d, s, source.moduli)])
            audit_collapse(source, out)
            return out
    factors = tuple(order(s, source.moduli) for s in source.directions)
    pairs = []
    for c in product(*(range(n) for n in factors)):
        x = source.zero
        for a, s in zip(c, source.directions):
            x = add(x, scale(a, s, source.moduli), source.moduli)
        pairs.append((c, x))
    out = Product(source, factors, tuple(pairs))
    audit_product(source, out)
    return out


def independent_basis(source):
    """For elementary abelian sources, retain a maximal generator basis."""
    need(type(source) is Cayley, "typed source required")
    p = source.moduli[0]
    need(all(m == p for m in source.moduli) and prime(p), "elementary abelian group required")
    basis, span = [], {source.zero}
    for s in source.directions:
        if s not in span:
            basis.append(s)
            span = set(subgroup_words(source.moduli, tuple(basis)))
    return make_cayley(source.moduli, tuple(basis))


def prime(p):
    return type(p) is int and p >= 2 and all(p % j for j in range(2, isqrt(p) + 1))


# Exact finite fields.  The small degrees below use independent monic-factor
# irreducibility checking; neither a trusted table nor float arithmetic enters.
def prem(a, b, p):
    a, b = list(a), list(b)
    while len(a) > 1 and not a[-1]:
        a.pop()
    while len(a) >= len(b) and any(a):
        k, c = len(a) - len(b), a[-1] * pow(b[-1], -1, p) % p
        for j, x in enumerate(b):
            a[k + j] = (a[k + j] - c * x) % p
        while len(a) > 1 and not a[-1]:
            a.pop()
    return tuple(a)


def irreducible(p, polynomial):
    n = len(polynomial) - 1
    return all(prem(polynomial, c + (1,), p) != (0,)
               for d in range(1, n // 2 + 1) for c in product(range(p), repeat=d))


@dataclass(frozen=True)
class FiniteField:
    p: int
    polynomial: tuple  # low to high, monic

    def __post_init__(self):
        need(prime(self.p), "prime characteristic required")
        need(type(self.polynomial) is tuple and len(self.polynomial) >= 2,
             "monic polynomial tuple required")
        need(all(type(c) is int and 0 <= c < self.p for c in self.polynomial)
             and self.polynomial[-1] == 1, "canonical monic coefficients required")
        need(irreducible(self.p, self.polynomial), "reducible modulus")

    @property
    def degree(self):
        return len(self.polynomial) - 1

    @property
    def one(self):
        return (1,) + (0,) * (self.degree - 1)

    def times(self, x, y):
        z = [0] * (2 * self.degree - 1)
        for i, a in enumerate(x):
            for j, b in enumerate(y):
                z[i + j] = (z[i + j] + a * b) % self.p
        z = prem(z, self.polynomial, self.p)
        return z + (0,) * (self.degree - len(z))

    def power(self, x, n):
        z = self.one
        while n:
            if n & 1:
                z = self.times(z, x)
            x, n = self.times(x, x), n // 2
        return z


def norm_circle(field):
    need(type(field) is FiniteField and field.degree % 2 == 0, "even degree field required")
    q = field.p ** (field.degree // 2)
    units = tuple(x for x in product(range(field.p), repeat=field.degree)
                  if field.power(x, q + 1) == field.one)
    need(len(units) == q + 1, "wrong norm kernel cardinality")
    return make_cayley((field.p,) * field.degree, units)


# An exact rational cyclotomic quotient for finite positive-coordinate controls.
def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return tuple(a)


def qdiv(a, b):
    r, q = list(map(F, a)), [F(0)] * max(1, len(a) - len(b) + 1)
    while len(r) >= len(b) and any(r):
        k, c = len(r) - len(b), r[-1] / b[-1]
        q[k] = c
        for j, x in enumerate(b):
            r[k + j] -= c * x
        r = list(trim(r))
    return trim(q), trim(r)


def cyclotomic(n):
    f = (-1,) + (0,) * (n - 1) + (1,)
    for d in range(1, n):
        if n % d == 0:
            f, r = qdiv(f, cyclotomic(d))
            need(not any(r), "cyclotomic division")
    return f


class Cyclotomic:
    def __init__(self, n):
        natural(n, 4)
        need(n % 4 == 0, "field must contain i")
        self.n, self.f = n, cyclotomic(n)
        self.degree = len(self.f) - 1
        self.zero, self.one = self.reduce((0,)), self.reduce((1,))
        self.powers = tuple(self.reduce((0,) * j + (1,)) for j in range(n))

    def reduce(self, x):
        r = qdiv(x, self.f)[1]
        return r + (F(0),) * (self.degree - len(r))

    def add(self, x, y):
        return tuple(a + b for a, b in zip(x, y))

    def scale(self, a, x):
        return tuple(a * b for b in x)

    def times(self, x, y):
        z = [F(0)] * (2 * self.degree - 1)
        for i, a in enumerate(x):
            for j, b in enumerate(y):
                z[i + j] += a * b
        return self.reduce(z)

    def conjugate(self, x):
        z = self.zero
        for j, a in enumerate(x):
            z = self.add(z, self.scale(a, self.powers[-j % self.n]))
        return z

    def norm(self, x):
        return self.times(x, self.conjugate(x))

    def difference(self, x, y):
        return self.add(x, self.scale(-1, y))

    def rotation(self, t):
        return self.add(self.scale(F(1 - t*t, 1 + t*t), self.one),
                        self.scale(F(2*t, 1 + t*t), self.powers[self.n // 4]))


def product_edges(labels, factors):
    edges = set()
    for a, b in combinations(range(len(labels)), 2):
        ds = [i for i, (x, y) in enumerate(zip(labels[a], labels[b])) if x != y]
        if len(ds) == 1:
            i = ds[0]
            if (labels[a][i] - labels[b][i]) % factors[i] in (1, factors[i] - 1):
                edges.add((a, b))
    return frozenset(edges)


def unit_edges(field, points):
    need(len(set(points)) == len(points), "colliding geometric source vertices")
    return frozenset((a, b) for a, b in combinations(range(len(points)), 2)
                     if field.norm(field.difference(points[a], points[b])) == field.one)


def realize_product(factors, parameter_cap=24):
    """Budgeted exact positive search; None means cap exhaustion, not impossibility."""
    need(type(factors) is tuple, "factor tuple required")
    for n in factors:
        natural(n, 2)
    natural(parameter_cap)
    field = Cyclotomic(lcm(4, *factors))
    labels, points, phases = ((),), (field.zero,), []
    for i, n in enumerate(factors):
        polygon = [field.zero]
        for j in range(1, n):
            polygon.append(field.add(polygon[-1], field.powers[(j - 1) * field.n // n]))
        next_labels = tuple(c + (j,) for c in labels for j in range(n))
        wanted = product_edges(next_labels, factors[:i + 1])
        for t in range(parameter_cap + 1):
            rotation = field.rotation(t)
            next_points = tuple(field.add(x, field.times(rotation, y))
                                for x in points for y in polygon)
            if len(set(next_points)) == len(next_points) and unit_edges(field, next_points) == wanted:
                labels, points = next_labels, next_points
                phases.append(t)
                break
        else:
            return None
    return field, labels, points, tuple(phases)


def main():
    checks = 0

    def check(ok, why):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(why)

    def rejects(function, *args):
        try:
            function(*args)
        except (ValueError, TypeError):
            check(True, "rejected hostile")
        else:
            check(False, "malformed packet accepted")

    # Exhaust all inverse-pair subsets of selected small abelian groups.
    classified, positives, negatives = 0, 0, 0
    for moduli in ((2,), (3,), (4,), (5,), (6,), (8,), (2, 2), (2, 3), (3, 3), (2, 4)):
        all_vectors = tuple(x for x in product(*(range(m) for m in moduli)) if any(x))
        reps = make_cayley(moduli, all_vectors).directions
        for mask in range(1 << len(reps)):
            source = make_cayley(moduli, tuple(s for i, s in enumerate(reps) if mask >> i & 1))
            out = classify(source)
            full = subgroup_words(moduli, source.directions)
            independent = prod(order(s, moduli) for s in source.directions) == len(full)
            check(isinstance(out, Product) == independent, "independent product-size check")
            if isinstance(out, Collapse):
                check(audit_collapse(source, out), "collapse verifier")
                negatives += 1
            else:
                check(audit_product(source, out), "product verifier")
                positives += 1
            classified += 1

    # Entire norm-circle sources, not preselected direction fragments.
    field_data = ((2, (1, 1, 1)), (3, (1, 0, 1)), (2, (1, 1, 0, 0, 1)),
                  (5, (3, 0, 1)), (7, (1, 0, 1)), (3, (2, 1, 0, 0, 1)),
                  (13, (11, 0, 1)), (17, (14, 0, 1)))
    circles = []
    for p, f in field_data:
        field = FiniteField(p, f)
        source = norm_circle(field)
        q = p ** (field.degree // 2)
        out = classify(source)
        check(isinstance(out, Product) == (q == 3), "norm-circle classification")
        basis = independent_basis(source)
        check(len(basis.directions) == field.degree, "norm units span full field")
        check(isinstance(classify(basis), Product), "retained-basis positive escape")
        circles.append({"q": q, "vertices": q*q, "inverse_pairs": len(source.directions),
                        "retained_basis": len(basis.directions),
                        "deleted_families": len(source.directions) - len(basis.directions),
                        "faithful_unit_realization": q == 3})
        for s in source.directions:
            check(order(s, source.moduli) == p, "actual additive order")

    # A fully 3-colourable 81-vertex hostile is still not a unit point graph.
    frame = make_cayley((3,) * 4, tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
                       + ((2, 2, 2, 2),))
    check(isinstance(classify(frame), Collapse), "ternary frame must fold")
    for x in product(range(3), repeat=4):
        for s in frame.directions:
            check(sum(x) % 3 != sum(add(x, s, frame.moduli)) % 3,
                  "explicit ternary-frame three-colouring")

    positive_coordinates = []
    for factors in ((2,), (3,), (5,), (2, 2, 2), (3, 3), (3, 4), (5, 5)):
        packet = realize_product(factors)
        check(packet is not None, "positive exact coordinate search cap")
        field, labels, points, phases = packet
        expected = product_edges(labels, factors)
        check(unit_edges(field, points) == expected, "all pair distances exact")
        check(len(set(points)) == prod(factors), "all source vertices distinct")
        for a, b in combinations(range(len(points)), 2):
            check((field.norm(field.difference(points[a], points[b])) == field.one)
                  == ((a, b) in expected), "induced edge iff")
        positive_coordinates.append({"factors": factors, "vertices": len(points),
                                     "edges": len(expected), "rational_rotation_parameters": phases,
                                     "cyclotomic_order": field.n})
    check(realize_product((3, 3), 0) is None, "identical rotations lose injectivity")

    sample = make_cayley((5, 5), ((1, 0), (2, 2), (2, 3)))
    collapse = classify(sample)
    check(audit_collapse(sample, collapse), "q5 explicit source")
    rejects(audit_collapse, sample, replace(collapse, period=True))
    rejects(audit_collapse, sample, replace(collapse, order=5.0))
    rejects(audit_collapse, sample, replace(collapse, other_word=(collapse.direction,)))
    rejects(audit_collapse, make_cayley((5, 5), ((1, 0),)), collapse)
    rejects(make_cayley, (True,), ((0,),))
    rejects(make_cayley, (3,), ((1.0,),))
    rejects(make_cayley, (3,), ((0,),))
    rejects(FiniteField, 5, (1, 0, 1))
    rejects(FiniteField, True, (1, 1, 1))
    rejects(realize_product, (3, True))
    rejects(realize_product, (3, 3), 2.0)
    product_source = make_cayley((3, 3), ((1, 0), (0, 1)))
    proof = classify(product_source)
    damaged = (((False, 0), (0, 0)),) + proof.coordinate_map[1:]
    rejects(audit_product, product_source, replace(proof, coordinate_map=damaged))

    print(json.dumps({"checks": checks,
                      "complete_small_source_universe": classified,
                      "product_sources": positives, "collapse_sources": negatives,
                      "norm_circle_controls": circles,
                      "positive_coordinate_controls": positive_coordinates,
                      "q5_collapse": {"direction": collapse.direction, "order": collapse.order,
                                      "period": collapse.period, "other_word": collapse.other_word},
                      "three_colour_nonrealizable_control": "81-vertex ternary five-frame",
                      "scope": "No new six-colour plane construction or Euclidean lower bound."},
                     indent=2))


if __name__ == "__main__":
    main()
