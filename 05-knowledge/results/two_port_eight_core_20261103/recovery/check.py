"""Independent endpoint-choice expansion and exact interpolation controls."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORES = HERE.parent / "cores"
CHECKS = 0


def need(condition, description):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(description)


def edges_of(n, mask):
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    return [edge for bit, edge in enumerate(pairs) if mask >> bit & 1]


def literal_coefficients(n, edges):
    """Enumerate every complete orientation, with no intermediate truncation."""
    coefficients = Counter()
    exponents = [0] * n
    leaves = 0

    def visit(k, sign):
        nonlocal leaves
        if k == len(edges):
            leaves += 1
            if all(e <= 2 for e in exponents):
                coefficients[tuple(exponents)] += sign
            return
        u, v = edges[k]
        exponents[u] += 1
        visit(k + 1, sign)
        exponents[u] -= 1
        exponents[v] += 1
        visit(k + 1, -sign)
        exponents[v] -= 1

    visit(0, 1)
    need(leaves == 2 ** len(edges), "complete edge-endpoint universe")
    need(exponents == [0] * n, "backtracking restored every exponent")
    return {e: c for e, c in coefficients.items() if c}, leaves


def denominator(values, a):
    return prod(a - b for b in values if b != a)


def coefficient_functional(sets, evaluate):
    return sum((Fraction(evaluate(a), prod(denominator(s, x)
               for s, x in zip(sets, a))) for a in product(*sets)), Fraction(0))


def interpolation_controls():
    monomials = 0
    for n in range(4):
        for t in product(range(3), repeat=n):
            total = sum(t)
            sets = [tuple(1 + 3 * i + 2 * j for j in range(e + 1))
                    for i, e in enumerate(t)]
            for u in product(range(total + 1), repeat=n):
                if sum(u) > total:
                    continue
                actual = coefficient_functional(
                    sets, lambda a, u=u: prod(x ** k for x, k in zip(a, u)))
                need(actual == int(u == t), "exact monomial interpolation")
                monomials += 1
    need(coefficient_functional([(1, 2)], lambda a: a[0] ** 2) == 3,
         "degree-bound hostile gives 3, not the zero x coefficient")

    graphs = 0
    certified = 0
    for n in range(5):
        for mask in range(1 << (n * (n - 1) // 2)):
            edges = edges_of(n, mask)
            coefficients, _ = literal_coefficients(n, edges)
            graphs += 1
            if not coefficients:
                continue
            t = min(coefficients, key=lambda e: sum(x * 3 ** i for i, x in enumerate(e)))
            # Overlapping lists, unequal selected sizes, and a singleton all occur.
            sets = [tuple((i % 2) + j for j in range(e + 1))
                    for i, e in enumerate(t)]
            value = coefficient_functional(sets, lambda a: prod(a[u] - a[v]
                                           for u, v in edges))
            need(value == coefficients[t], "graph interpolation equals literal coefficient")
            need(any(all(a[u] != a[v] for u, v in edges) for a in product(*sets)),
                 "certificate yields an actual proper list colouring")
            certified += 1

    triangle, _ = literal_coefficients(3, edges_of(3, 7))
    need(triangle.get((1, 1, 1), 0) == 0 and bool(triangle),
         "one zero coefficient does not refute three-choosability")
    k4, _ = literal_coefficients(4, edges_of(4, 63))
    need(not k4, "K4 has no capped polynomial support")
    need(not any(all(a[u] != a[v] for u, v in edges_of(4, 63))
                 for a in product(range(3), repeat=4)), "K4 common-three-list hostile")
    return {"monomials": monomials, "all_labelled_graphs_through4": graphs,
            "small_graph_certificates": certified}


def main():
    source_names = ("check.py", "LIST_CERTIFICATE.md", "INCOMING.md")
    pins = {name: sha256((HERE / name).read_bytes()).hexdigest() for name in source_names}
    census_raw = (CORES / "order8.json").read_bytes()
    census = json.loads(census_raw)
    need(census["n"] == 8 and census["classes"] == 38, "specified supplied bank")
    need(len(census["representatives"]) == 38, "all supplied representatives present")
    for name, expected in census["source_sha256"].items():
        need(sha256((CORES / name).read_bytes()).hexdigest() == expected,
             "producer receipt source pin: " + name)
    total_orientations = 0
    rows = []
    for item in census["representatives"]:
        mask = item["mask"]
        edges = edges_of(8, mask)
        need(len(edges) == item["edges"], "edge count " + str(mask))
        adjacency = [sum(1 << v for u, v in edges if u == i) +
                     sum(1 << u for u, v in edges if v == i) for i in range(8)]
        need(adjacency == item["rows"], "mask-to-adjacency convention " + str(mask))
        coefficients, orientations = literal_coefficients(8, edges)
        total_orientations += orientations
        need(len(coefficients) == item["polynomial_support"], "full capped support " + str(mask))
        chosen = item["list_certificate"]
        if not coefficients:
            need(chosen is None, "null certificate for zero support " + str(mask))
        else:
            t = min(coefficients, key=lambda e: sum(x * 3 ** i for i, x in enumerate(e)))
            need(chosen is not None, "nonempty support certificate " + str(mask))
            need(list(t) == chosen["exponents"], "first base-three exponent vector " + str(mask))
            need(coefficients[t] == chosen["coefficient"], "signed coefficient " + str(mask))
            need(sum(t) == len(edges) and max(t) <= 2, "homogeneous three-list bounds " + str(mask))
        rows.append({"mask": mask, "edges": len(edges), "orientations": orientations,
                     "support": len(coefficients), "certificate": chosen})
    residuals = [row["mask"] for row in rows if not row["support"]]
    need(residuals == [5675891, 5683827, 5691763, 5692274,
                       6215909, 6216436, 14587370, 48675698], "exact eight residuals")
    controls = interpolation_controls()
    need(pins == {name: sha256((HERE / name).read_bytes()).hexdigest() for name in source_names},
         "independent audit sources unchanged")
    need(census_raw == (CORES / "order8.json").read_bytes(), "input receipt unchanged")
    print(json.dumps({"status": "FINITE-EXACT INDEPENDENT COEFFICIENT AUDIT; THEOREM PROVED SEPARATELY",
                      "checks": CHECKS, "order8_certificates": 30, "order8_residuals": residuals,
                      "order8_complete_orientations": total_orientations,
                      "controls": controls, "representatives": rows,
                      "source_sha256": pins, "census_sha256": sha256(census_raw).hexdigest()},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
