"""Separate symbolic geometry/algebra path; requires SymPy (tested with 1.14.0)."""
import json
from itertools import combinations, product
from pathlib import Path
import sympy as s


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


data = json.loads(Path(__file__).with_name("forcing_check.json").read_text())
basis = (s.Integer(1), s.sqrt(3), s.sqrt(11), s.sqrt(33))
points = [tuple(sum(s.Rational(c) * b for c, b in zip(row, basis))
                for row in point) for point in data["points"]]
edges = []
for i, j in combinations(range(19), 2):
    squared = s.expand(sum((a - b) ** 2 for a, b in zip(points[i], points[j])))
    require(squared != 0, "coincident geometric points")
    if squared == 1:
        edges.append([i, j])
require(edges == data["graph19_edges"], "independent exact contact graph")
colouring = {int(i): c for i, c in data["global_four_colouring"].items()}
require(all(colouring[i] != colouring[j] for i, j in edges), "four-colouring")
core = [(i, j) for i, j in edges if j < 7]
counts = []
for central in range(6):
    count = 0
    for rest in product(range(3), repeat=6):
        colours = (central,) + rest
        count += all(colours[i] != colours[j] for i, j in core)
    counts.append(count)
require(counts == [0, 0, 0, 24, 24, 24], "independent literal palette count")
require([sum(v in e for e in edges) for v in range(7)] == [7, 6, 6, 6, 6, 6, 6],
        "sharp interior degree boundary")
print("SymPy", s.__version__, ": all171 contacts, four-colouring, degree boundary, palette counts PASS")

x, k = s.symbols("x k")
coefs = {
    2: (756, 792, 270, 30),
    3: (13380, 12974, 4140, 435),
    4: (142335, 130750, 39810, 4020),
    5: (893760, 788970, 231840, 22680),
    6: (2593080, 2222640, 635040, 60480),
}
p = 0
for j, c in coefs.items():
    polynomial = sum(v * x ** i for i, v in enumerate(c))
    falling = s.prod(x - i for i in range(j - 2))
    p += s.Rational(2 ** (j - 2), s.factorial(2 * j - 3)) * falling * polynomial
expected = (16*k**7 + 64*k**6 + 8*k**5 - 140*k**4 - 22*k**3 + 135*k**2 + 61*k - 110)/6
require(s.expand(p.subs(x, k - 2) - expected) == 0, "independent D3 polynomial expansion")
require(s.expand(sum(v*x**i for i, v in enumerate(coefs[5])) -
                 210*(2*x+7)*(6*x+19)*(9*x+32)) == 0, "C5 factorization")
require(s.expand(sum(v*x**i for i, v in enumerate(coefs[6])) -
                 7560*(2*x+7)**3) == 0, "C6 factorization")
print("Independent symbolic C2..C6 expansion, C5/C6 factors and N3 formula PASS")
print("PASS")
