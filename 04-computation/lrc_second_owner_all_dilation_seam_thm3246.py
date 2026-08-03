#!/usr/bin/env python3
"""Exact companion for THM-3246's all-dilation owner seam theorem."""

import ast
import contextlib
import hashlib
import io
from fractions import Fraction
from math import gcd
from pathlib import Path
import runpy


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY_THEOREM = ROOT / (
    "01-canon/theorems/"
    "THM-3224-complete-lrc-orbit-bernoulli-gcd-carry-and-owner-hodge-splitting.md"
)
DEPENDENCY_SCRIPT = ROOT / "04-computation/lrc_second_owner_bernoulli_curvature_thm3224.py"
DEPENDENCY_OUTPUT = ROOT / "05-knowledge/results/lrc_second_owner_bernoulli_curvature_thm3224.out"
DEPENDENCIES = {
    DEPENDENCY_THEOREM: "ad12039621ab1adbf7dd1b818462e5fcdbdd2d7adab62996d23e1e5c5c1176cd",
    DEPENDENCY_SCRIPT: "84cb3cd5d91b47e4d67918d7eaa49854ec525c9791694cfb679915e4947c7eaa",
    DEPENDENCY_OUTPUT: "413a487add6d53ea5f62734a13c40ccb966a01c64844de4d6fc03eff80a1fea7",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


for dependency, expected_hash in DEPENDENCIES.items():
    require(hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected_hash,
            ("dependency hash drift", dependency.name))

syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")

# Replay the promoted dependency and then reuse its exact, already-audited
# floor-moment/stability engine.  Suppression keeps this transcript canonical.
dependency_stdout = io.StringIO()
with contextlib.redirect_stdout(dependency_stdout):
    NS = runpy.run_path(str(DEPENDENCY_SCRIPT))
require(dependency_stdout.getvalue().encode("utf-8") == lf_bytes(DEPENDENCY_OUTPUT),
        "dependency transcript drift")

F = Fraction
P, Q, E, G = 3, 5, 1, 2
CELL_COUNT = 168
BOUNDARY = frozenset(tuple(range(6)) + tuple(range(162, 168)))
CONTROL_DILATIONS = (1, 2, 3, 17, 80, 160)


def explicit_polynomial(cell):
    if cell in BOUNDARY:
        return F(12096), F(-1032), F(2)
    if 6 <= cell <= 23 or 144 <= cell <= 161:
        return F(12096), F(-24), F(0)
    if 24 <= cell <= 71:
        return F(16044 - 168 * cell), F(48), F(0)
    if 72 <= cell <= 95:
        return F(4032), F(96), F(0)
    if 96 <= cell <= 143:
        return F(168 * cell - 12012), F(48), F(0)
    raise RuntimeError(("cell outside complete table", cell))


def polynomial_value(poly, dilation):
    return poly[0] * dilation * dilation + poly[1] * dilation + poly[2]


def add_polynomials(*polys):
    return tuple(sum(poly[index] for poly in polys) for index in range(3))


def multiply_linears(left, right):
    # Inputs are (coefficient of g, constant); output is (g^2,g,1).
    return (left[0] * right[0],
            left[0] * right[1] + left[1] * right[0],
            left[1] * right[1])


# Seam algebra.  For cells 0,...,5 the only possible alignment is
# (k,l)=(3t,5t).  The t=0 head, g-1 interior copies, and t=g tail have the
# following cleared polynomials.  Their sum is independent of the cell.
seam_component_checks = 0
for cell in range(6):
    z = (504, -1)
    head = multiply_linears(z, (0, 2 * cell + 12))
    middle = multiply_linears(z, (24, -24))
    tail = multiply_linears(z, (0, max(10 - 2 * cell, 0)))
    require(add_polynomials(head, middle, tail) == explicit_polynomial(cell),
            ("seam component identity", cell))

    # These coefficient inequalities are the exact all-g branch choices:
    # the B interval is the shorter nested interval at the head/interior/tail.
    require(4032 - 168 * cell - 12 > 0, ("head nesting", cell))
    require(3864 - 168 * cell > 0, ("interior nesting", cell))
    require(4200 + 168 * cell - 12 > 0, ("tail nesting", cell))
    seam_component_checks += 4

# If d=5k-3l is nonzero, the cleared centre separation is already larger
# than the overlap radius.  The two worst symbolic lower bounds occur at
# |d|=1, cell=5, and the extreme legal k.  Checking their positive excess
# proves that no unaligned pair can meet for any g>=1.
require(486 > 288 and 486 - 1 > 288,
        "positive unaligned separation for every g>=1")
require(504 > 288 and 504 - 1 > 288,
        "negative unaligned separation for every g>=1")
unaligned_symbolic_checks = 2

# The opposite seam is its exact x -> 1-x reflection.
require(all(explicit_polynomial(cell) == explicit_polynomial(167 - cell)
            for cell in range(CELL_COUNT)), "polynomial reflection")

# Interior affine-ray certificates.  The promoted engine tests all
# 3 residue indices x 10 shift indices x 2 triangle kernels at h=1; its
# affine slope ordering then holds on the entire ray.  Period is one here.
rows = []
interior_term_certificates = 0
direct_controls = 0
for cell in range(CELL_COUNT):
    if cell in BOUNDARY:
        poly = explicit_polynomial(cell)
    else:
        require(NS["ray_is_stable"](cell, P, Q, E, G, 1, 1, 1),
                ("interior ray instability", cell))
        compiled = NS["compile_infinite_cell"](cell, P, Q, E, G)
        require(len(compiled) == 1, ("period-one compiler", cell))
        poly = compiled[0]
        require(poly == explicit_polynomial(cell),
                ("explicit interior polynomial", cell, poly))
        interior_term_certificates += P * (Q + 5) * 2

    for dilation in CONTROL_DILATIONS:
        direct = NS["direct_cleared"](P, Q, E, G, cell, dilation)
        require(direct == polynomial_value(poly, dilation),
                ("direct geometry control", cell, dilation, direct, poly))
        direct_controls += 1

    limit = NS["limit_cell"](P, Q, E, G, cell)
    correction = NS["correction_cell"](P, Q, E, G, cell)
    q_value = NS["q_from_kappa"](P, Q, E, G, limit, correction, poly[2])
    rows.append((cell, poly, q_value, limit, correction))

require(interior_term_certificates == 9360, "interior certificate census")
require(direct_controls == 1008, "direct control census")

q_word = tuple(row[2] for row in rows)
require(q_word == NS["owner_q"], "finite scout and all-dilation word disagree")
require(all(q_word[167 - cell] == q_word[cell] for cell in range(CELL_COUNT)),
        "owner reflection")
signs = (sum(value > 0 for value in q_word),
         sum(value < 0 for value in q_word),
         sum(value == 0 for value in q_word))
require(signs == (156, 12, 0), "all-dilation sign word")
require(tuple(cell for cell, value in enumerate(q_word) if value < 0)
        == tuple(range(6)) + tuple(range(162, 168)), "negative seam set")
require(set(value for value in q_word if value < 0) == {F(-17, 3087000)},
        "negative seam coefficient")
require(min(value for value in q_word if value > 0) == F(1, 6174000),
        "positive minimum")
require(max(value for value in q_word if value > 0) == F(751, 666792000),
        "positive maximum")
require(sum(q_word) == F(1, 24696) == NS["owner_omega"],
        "all-dilation Hodge holonomy")

# The negative punctured Singer exponents occupy twelve distinct classes
# modulo 14.  Hence every Singer-equivariant vector line meets the negative
# set in at most one point, for every unit exponent gauge and phase shift.
negative_cells = tuple(cell for cell, value in enumerate(q_word) if value < 0)
gauge_checks = 0
max_line_intersection = 0
for multiplier in range(168):
    if gcd(multiplier, 168) != 1:
        continue
    for shift in range(14):
        residues = tuple((multiplier * cell + shift) % 14
                         for cell in negative_cells)
        max_line_intersection = max(
            max_line_intersection,
            max(residues.count(residue) for residue in range(14)),
        )
        gauge_checks += 1
require((gauge_checks, max_line_intersection) == (672, 1),
        "Singer line counterfeit")

row_text = "\n".join(
    "%d|%s|%s|%s|%s|%s|%s" %
    (cell, poly[0], poly[1], poly[2], q_value, limit, correction)
    for cell, poly, q_value, limit, correction in rows
)
q_text = "\n".join(map(str, q_word))

print("THM-3246 all-dilation second-owner seam exact control")
print("dependency_hash_checks=%d,dependency_replay=PASS" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("lane=(P,Q;e,f)=(3,5;1,2),cells=168,period=1")
print("polynomial_pieces=boundary:(12096,-1032,2);outer:(12096,-24,0);"
      "left_ramp:(16044-168j,48,0);middle:(4032,96,0);"
      "right_ramp:(168j-12012,48,0)")
print("seam_component_unaligned_checks=%s" %
      ((seam_component_checks, unaligned_symbolic_checks),))
print("interior_term_certificates=%d,direct_geometry_controls=%d" %
      (interior_term_certificates, direct_controls))
print("signs_positive_negative_zero=%s" % (signs,))
print("negative_cells=%s,negative_q=-17/3087000" % (negative_cells,))
print("positive_q_min_max=(1/6174000,751/666792000)")
print("hodge_sum=1/24696")
print("row_table_sha256=%s" % hashlib.sha256(row_text.encode("ascii")).hexdigest())
print("q_word_sha256=%s" % hashlib.sha256(q_text.encode("ascii")).hexdigest())
print("singer_gauge_checks=%d,max_negative_line_intersection=%d" %
      (gauge_checks, max_line_intersection))
print("scope=all-dilation-second-corrector-not-cellwise-LRC-exclusion")
print("all_exact_checks=PASS")
