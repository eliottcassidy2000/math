#!/usr/bin/env python3
"""Exact interacting-hub controls beyond positive unit ears. No dependencies."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial
from pathlib import Path
import hashlib
import json
import sys


def need(ok, detail):
    if not ok:
        raise RuntimeError(detail)


@dataclass(frozen=True)
class K:
    # Q(sqrt3,sqrt11), basis 1,sqrt3,sqrt11,sqrt33.
    c: tuple = (F(0), F(0), F(0), F(0))

    def __add__(self, y):
        if not isinstance(y, K):
            y = scalar(y)
        return K(tuple(a + b for a, b in zip(self.c, y.c)))

    __radd__ = __add__

    def __neg__(self):
        return K(tuple(-a for a in self.c))

    def __sub__(self, y):
        return self + (-y if isinstance(y, K) else -F(y))

    def __mul__(self, y):
        if not isinstance(y, K):
            y = scalar(y)
        c = [F(0)] * 4
        for i, a in enumerate(self.c):
            for j, b in enumerate(y.c):
                both = i & j
                c[i ^ j] += a * b * (3 if both & 1 else 1) * (11 if both & 2 else 1)
        return K(tuple(c))

    __rmul__ = __mul__

    def conjugate(self, mask):
        return K(tuple(a * (-1) ** ((j & mask).bit_count()) for j, a in enumerate(self.c)))

    def inv(self):
        n = self.conjugate(1) * self.conjugate(2) * self.conjugate(3)
        d = self * n
        need(d.c[0] != 0 and d.c[1:] == (0, 0, 0), ("field inversion", self, d))
        return K(tuple(a / d.c[0] for a in n.c))

    def __truediv__(self, y):
        if not isinstance(y, K):
            y = scalar(y)
        return self * y.inv()

    def serial(self):
        return [str(a) for a in self.c]


def scalar(x):
    return K((F(x), F(0), F(0), F(0)))


R3 = K((F(0), F(1), F(0), F(0)))
R11 = K((F(0), F(0), F(1), F(0)))
ZERO = K()
ONE = scalar(1)


def minus(P, Q):
    return tuple(a - b for a, b in zip(P, Q))


def norm(P):
    return P[0] * P[0] + P[1] * P[1]


def distance(P, Q):
    return norm(minus(P, Q))


def graph(points):
    need(len(set(points)) == len(points), "point collision")
    return {(i, j) for i, j in combinations(range(len(points)), 2) if distance(points[i], points[j]) == ONE}


def lattice(a, b):
    return (scalar(F(2 * a + b, 2)), R3 * F(b, 2))


def partitions(n):
    P = [(0,)]
    for _ in range(1, n):
        P = [a + (j,) for a in P for j in range(max(a) + 2)]
    return P


def falling(q, k):
    return factorial(q) // factorial(q - k) if k <= q else 0


def adjacency(n, edges):
    A = [set() for _ in range(n)]
    for u, v in edges:
        A[u].add(v)
        A[v].add(u)
    return A


def extension_count(n, edges, pi, q, selected=None):
    """Exact list-colouring DP on the private forest, with full fixed ports."""
    m = len(pi)
    if max(pi) >= q or any(pi[u] == pi[v] for u, v in edges if v < m):
        return 0
    if selected is None:
        selected = set(range(m, n))
    A = adjacency(n, edges)
    lists = {u: set(range(q)) - {pi[v] for v in A[u] if v < m} for u in selected}
    seen = set()
    def dfs(u, parent):
        need(u not in seen, ("private graph not forest", u))
        seen.add(u)
        children = [v for v in A[u] if v in selected and v != parent]
        tables = [dfs(v, u) for v in children]
        return {c: prod(sum(value for d, value in table.items() if c != d) for table in tables)
                for c in lists[u]}
    total = 1
    for u in sorted(selected):
        if u not in seen:
            total *= sum(dfs(u, None).values())
    return total


def brute_extensions(n, edges, pi, q):
    if max(pi) >= q or any(pi[u] == pi[v] for u, v in edges if v < len(pi)):
        return 0
    return sum(all(colours[u] != colours[v] for u, v in edges)
               for tail in product(range(q), repeat=n - len(pi))
               for colours in [pi + tail])


def explicit_bank_count(edges, pi, q):
    if max(pi) >= q or any(pi[u] == pi[v] for u, v in edges if v < 6):
        return 0
    A = adjacency(14, edges)
    lists = {u: set(range(q)) - {pi[v] for v in A[u] if v < 6} for u in range(6, 14)}
    isolates = prod(len(lists[u]) for u in (6, 7, 12, 13))
    return isolates * sum(len(lists[8] - {x}) * len(lists[11] - {y})
                          for x in lists[9] for y in lists[10] if x != y)


def list_witness(A, lists):
    colours = [-1] * len(A)
    def dfs(left):
        if not left:
            return tuple(colours)
        possibilities = {u: lists[u] & ~sum(1 << c for c in set(colours[v] for v in A[u]) if c >= 0) for u in left}
        u = min(left, key=lambda v: (possibilities[v].bit_count(), v))
        available = possibilities[u]
        if not available:
            return None
        while available:
            bit = available & -available
            available -= bit
            colours[u] = bit.bit_length() - 1
            answer = dfs(left - {u})
            if answer is not None:
                return answer
        colours[u] = -1
        return None
    return dfs(set(range(len(A))))


def count_and_witness(n, edges, q):
    A = adjacency(n, edges)
    order = sorted(range(n), key=lambda u: -len(A[u]))
    colours = [-1] * n
    first = None
    def dfs(i):
        nonlocal first
        if i == n:
            if first is None:
                first = list(colours)
            return 1
        u = order[i]
        bad = {colours[v] for v in A[u] if colours[v] >= 0}
        total = 0
        for c in range(q):
            if c not in bad:
                colours[u] = c
                total += dfs(i + 1)
        colours[u] = -1
        return total
    return dfs(0), first


def circle_center(A, B, C):
    d, e = minus(B, A), minus(C, A)
    det = d[0] * e[1] - d[1] * e[0]
    if det == ZERO:
        return None
    b, c = norm(B) - norm(A), norm(C) - norm(A)
    center = ((b * e[1] - c * d[1]) / (2 * det),
              (d[0] * c - e[0] * b) / (2 * det))
    return center if distance(center, A) == ONE else None


def main():
    source = Path(__file__).read_bytes()
    root_path = Path(__file__).resolve().parents[4]
    pins = json.loads(Path("SOURCE_PINS.json").read_text(encoding="utf-8"))
    inherited_bytes = {}
    for pin in pins:
        data = (root_path / pin["path"]).read_bytes()
        need(hashlib.sha256(data).hexdigest() == pin["sha256"], ("source pin", pin["path"]))
        inherited_bytes[pin["path"]] = data
    # The actual triangular-lattice double hub, ports first.
    B = [(0, 1), (1, -1), (-1, 0), (-1, 1), (2, 0), (2, -1)]
    ports = [lattice(*z) for z in B]
    hubs = [lattice(0, 0), lattice(1, 0)]
    triangular = ports + hubs
    tri_edges = graph(triangular)
    need(len(tri_edges) == 13, "triangular double-hub incidence")
    target = (0, 1, 2, 3, 2, 3)
    alternate = (0, 1, 2, 3, 3, 2)
    need(extension_count(8, tri_edges, target, 5) == 0 and extension_count(8, tri_edges, target, 6) == 2,
         "interaction zero")

    # Same shared ports; rational free-angle leaves make the six ports independent.
    free_ports = [ports[0], ports[1], lattice(-1, 0), (scalar(F(-3, 5)), scalar(F(4, 5))),
                  lattice(2, 0), (scalar(F(8, 5)), scalar(F(-4, 5)))]
    free_points = free_ports + hubs
    free_edges = graph(free_points)
    expected = {(0, 6), (1, 6), (2, 6), (3, 6), (0, 7), (1, 7), (4, 7), (5, 7), (6, 7)}
    need(free_edges == expected, ("nonlattice exact graph", free_edges))

    # All external centres touching at least two triangular ports. Each pair's
    # computed list meets the Euclidean circle-intersection count, certifying
    # completeness among ALL real planar centres, not just lattice candidates.
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    candidates = sorted({(a + x, b + y) for a, b in B for x, y in steps})
    normL = lambda z: z[0] ** 2 + z[0] * z[1] + z[1] ** 2
    near = lambda z, w: normL((z[0] - w[0], z[1] - w[1])) == 1
    for i, j in combinations(range(6), 2):
        d2 = normL((B[i][0] - B[j][0], B[i][1] - B[j][1]))
        actual = [z for z in candidates if near(z, B[i]) and near(z, B[j])]
        need(len(actual) == (2 if d2 < 4 else 1 if d2 == 4 else 0), ("complete circle pair", i, j, d2))
    bank_coords = [z for z in candidates if z not in B and sum(near(z, b) for b in B) >= 2]
    need(bank_coords == [(-2, 1), (-1, 2), (0, -1), (0, 0), (1, 0), (1, 1), (2, -2), (3, -1)], "external bank")
    bank_points = ports + [lattice(*z) for z in bank_coords]
    bank_edges = graph(bank_points)
    hub_ids = {6 + bank_coords.index(z) for z in [(0, 0), (1, 0)]}
    need({e for e in bank_edges if e[0] >= 6} == {(8, 9), (9, 10), (10, 11)}, "bank private path")
    lattice_colouring = [(a + 2 * b) % 3 for a, b in B + bank_coords]
    need(all(lattice_colouring[u] != lattice_colouring[v] for u, v in bank_edges), "lattice 3-colouring")

    P = partitions(6)
    tables = []
    sums = {"triangular": {}, "free_angle": {}, "full_bank": {}}
    supports = {"triangular": {}, "free_angle": {}, "full_bank": {}}
    for q in (5, 6):
        for key in sums:
            sums[key][q] = 0
            supports[key][q] = []
    for pi in P:
        row = {"partition": list(pi)}
        for q in (5, 6):
            values = {}
            for key, n, edges in [("triangular", 8, tri_edges), ("free_angle", 8, free_edges),
                                  ("full_bank", 14, bank_edges)]:
                values[key] = extension_count(n, edges, pi, q)
                if key != "full_bank":
                    need(values[key] == brute_extensions(n, edges, pi, q), ("independent literal hubs", key, q, pi))
                else:
                    need(values[key] == explicit_bank_count(edges, pi, q), ("independent bank factorization", q, pi))
                sums[key][q] += falling(q, len(set(pi))) * values[key]
                if values[key] > 0:
                    supports[key][q].append(pi)
            need((values["triangular"] > 0) == (values["full_bank"] > 0), ("max bank support", q, pi))
            row[str(q)] = values
        tables.append(row)
    need(set(P) - set(supports["free_angle"][5]) == {target, alternate, (0, 1, 2, 3, 4, 5)}, "exact q5 free-angle zero orbits")
    need(len(supports["free_angle"][6]) == 203, "q6 free-angle support")
    need(sums["free_angle"] == {5: 46080, 6: 300000}, "free-angle polynomial counts")

    # All256 one-stage centre assemblies, retaining complete q5/q6 supports.
    bank_supports = []
    type_counts = {5: {}, 6: {}}
    for mask in range(256):
        selected = {6 + j for j in range(8) if mask >> j & 1}
        row = {"mask": mask}
        for q in (5, 6):
            support_mask = sum((1 << j) for j, pi in enumerate(P)
                               if extension_count(14, bank_edges, pi, q, selected) > 0)
            row[str(q)] = hex(support_mask)
            type_counts[q][support_mask] = type_counts[q].get(support_mask, 0) + 1
            expected_mask = sum((1 << j) for j, pi in enumerate(P)
                                if (extension_count(8, tri_edges, pi, q) if hub_ids <= selected
                                    else max(pi) < q and all(pi[u] != pi[v] for u, v in tri_edges if v < 6)))
            need(support_mask == expected_mask, ("complete subset support classification", mask, q))
        bank_supports.append(row)
    need(sorted(type_counts[5].values()) == [64, 192] and list(type_counts[6].values()) == [256], "subset support types")

    # Genuine nonlattice4-chromatic port carrier: exact Moser spindle.
    A = (ZERO, ZERO)
    B0 = (R3, ZERO)
    C = (R3 / 2, scalar(F(1, 2)))
    D = (R3 / 2, scalar(F(-1, 2)))
    def rotate(P):
        return ((5 * P[0] - R11 * P[1]) / 6, (R11 * P[0] + 5 * P[1]) / 6)
    spindle = [A, B0, C, D, rotate(B0), rotate(C), rotate(D)]
    spindle_edges = graph(spindle)
    need(spindle_edges == {(0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
                           (0, 5), (0, 6), (4, 5), (4, 6), (5, 6), (1, 4)}, "Moser complete incidence")
    c3, _ = count_and_witness(7, spindle_edges, 3)
    c4, witness4 = count_and_witness(7, spindle_edges, 4)
    c5, witness5 = count_and_witness(7, spindle_edges, 5)
    need(c3 == 0 and c4 > 0 and c5 > 0, "actual4-chromatic carrier")
    triple_records = []
    outside = set()
    for triple in combinations(range(7), 3):
        centre = circle_center(*(spindle[j] for j in triple))
        p0, p1, p2 = (spindle[j] for j in triple)
        d, e = minus(p1, p0), minus(p2, p0)
        determinant = d[0] * e[1] - d[1] * e[0]
        gram_residual = norm(d) * norm(e) * distance(p1, p2) - 4 * determinant * determinant
        radius_one = determinant != ZERO and gram_residual == ZERO
        need(radius_one == (centre is not None), ("independent unit circumradius", triple))
        old_centres = [j for j, point in enumerate(spindle) if all(distance(point, spindle[v]) == ONE for v in triple)]
        need(not radius_one or len(old_centres) == 1, ("old centre completeness", triple))
        is_external = centre is not None and centre not in spindle
        if is_external:
            outside.add(centre)
        triple_records.append(dict(triple=triple, unit_centre=centre is not None,
                                   external=is_external,
                                   circumradius_residual=gram_residual.serial(), old_centres=old_centres,
                                   centre=None if centre is None else [z.serial() for z in centre]))
    outside_contacts = [dict(centre=[z.serial() for z in centre],
                             contacts=[j for j, p in enumerate(spindle) if distance(p, centre) == ONE])
                        for centre in sorted(outside, key=repr)]
    maximum = max([2] + [len(item["contacts"]) for item in outside_contacts])
    need(not outside_contacts and maximum == 2, ("external spindle triple contact", outside_contacts))

    # A genuine unit triangular prism supplies a3-core; common-palette list
    # enumeration tests whether this first surviving private core can block.
    triangle = [(ZERO, ZERO), (ONE, ZERO), (scalar(F(1, 2)), R3 / 2)]
    shift = (scalar(F(3, 5)), scalar(F(4, 5)))
    prism = triangle + [tuple(a + b for a, b in zip(p, shift)) for p in triangle]
    prism_edges = graph(prism)
    need(prism_edges == {(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)},
         "complete unit prism incidence")
    AP = adjacency(6, prism_edges)
    triples = [sum(1 << c for c in block) for block in combinations(range(5), 3)]
    list_cases = 0
    witnesses_hash = hashlib.sha256()
    for tail in product(triples, repeat=5):
        lists = (7,) + tail
        witness = list_witness(AP, lists)
        need(witness is not None and all(lists[u] >> c & 1 for u, c in enumerate(witness))
             and all(witness[u] != witness[v] for u, v in prism_edges), ("prism list colouring", lists, witness))
        list_cases += 1
        witnesses_hash.update(bytes(lists + witness))
    need(list_cases == 100000, "complete normalized prism lists")

    coords = dict(field="Q(sqrt3,sqrt11); coefficients in basis1,sqrt3,sqrt11,sqrt33",
                  triangular=[[[str(x) for x in z.c] for z in p] for p in triangular],
                  free_angle=[[[str(x) for x in z.c] for z in p] for p in free_points],
                  full_bank=[[[str(x) for x in z.c] for z in p] for p in bank_points],
                  spindle=[[[str(x) for x in z.c] for z in p] for p in spindle],
                  prism=[[[str(x) for x in z.c] for z in p] for p in prism],
                  triangular_edges=sorted(map(list, tri_edges)), free_angle_edges=sorted(map(list, free_edges)),
                  full_bank_edges=sorted(map(list, bank_edges)), spindle_edges=sorted(map(list, spindle_edges)),
                  prism_edges=sorted(map(list, prism_edges)),
                  full_bank_colouring3=lattice_colouring, spindle_colouring4=witness4,
                  spindle_colouring5=witness5, blocked_partition=list(target))
    receipt = dict(boundary_partitions=203, palettes=[5, 6],
                   full_counts=sums,
                   support_counts={key: {q: len(values) for q, values in byq.items()} for key, byq in supports.items()},
                   free_angle_labelled_support5=5 ** 6 - 2 * falling(5, 4),
                   free_angle_labelled_support6=6 ** 6,
                   blocked_q5_extensions=0, same_state_q6_extensions=2,
                   complete_external_two_contact_centres=8, bank_assemblies=256,
                   q5_assembly_support_type_sizes=sorted(type_counts[5].values()),
                   q6_assembly_support_type_sizes=list(type_counts[6].values()),
                   spindle_counts={3: c3, 4: c4, 5: c5}, spindle_triples=len(triple_records),
                   spindle_external_triple_centres=len(outside),
                   spindle_maximum_external_contacts=maximum,
                   normalized_prism_3lists_in_5palette=list_cases,
                   represented_unrestricted_prism_list_assignments=10 ** 6,
                   prism_list_witness_sha256=witnesses_hash.hexdigest(),
                   inherited_source_pins=len(pins),
                   source_sha256=hashlib.sha256(source).hexdigest())
    need(Path(__file__).read_bytes() == source, "source changed during run")
    need(all((root_path / name).read_bytes() == data for name, data in inherited_bytes.items()),
         "inherited source changed during run")
    outputs = {"coordinates.json": coords, "boundary_tables.json": tables,
               "bank_supports.json": bank_supports,
               "spindle_contacts.json": dict(triples=triple_records, external_centres=outside_contacts)}
    for name, data in outputs.items():
        if "--write" in sys.argv:
            Path(name).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        else:
            need(json.loads(Path(name).read_text(encoding="utf-8")) == json.loads(json.dumps(data)), ("frozen artifact", name))
    print(json.dumps(receipt, indent=2, sort_keys=True))


def prod(items):
    total = 1
    for x in items:
        total *= x
    return total


if __name__ == "__main__":
    main()
