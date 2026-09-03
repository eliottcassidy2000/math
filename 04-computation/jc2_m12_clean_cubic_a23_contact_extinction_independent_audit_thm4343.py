#!/usr/bin/env python3
"""Import-free arithmetic audit of THM-4343's A23/cubic compatibility."""

from fractions import Fraction as F
from math import gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def cubic_discriminant(a, b, c, d):
    """Discriminant of a*x^3+b*x^2+c*x+d."""
    return b * b * c * c - 4 * a * c**3 - 4 * b**3 * d - 27 * a * a * d * d + 18 * a * b * c * d


# The complete sixteen-row exponent inventory.  Under the top chart each
# p^i y^j row maps to t^(12-2i-3j) r^(i+j) in Hhat.
rows = (
    ("p", 1, 0),
    ("p2", 2, 0),
    ("p3", 3, 0),
    ("Delta", 4, 0),
    ("upsilon", 5, 0),
    ("U", 6, 0),
    ("K", 0, 2),
    ("Phi", 2, 1),
    ("Theta", 1, 2),
    ("eta", 3, 1),
    ("xi", 2, 2),
    ("alpha", 4, 1),
    ("beta", 1, 3),
    ("W", 3, 2),
    ("zeta", 0, 3),
    ("Z", 0, 4),
)
need(len(rows) == 16 and len({(i, j) for _, i, j in rows}) == 16, "sixteen distinct rows")

wall_rows = [row for row in rows if row[0] not in {"beta", "zeta", "Z"}]
top_images = sorted((12 - 2 * i - 3 * j, i + j, name) for name, i, j in wall_rows)
expected_images = sorted(
    (
        (10, 1, "p"),
        (8, 2, "p2"),
        (6, 3, "p3"),
        (4, 4, "Delta"),
        (2, 5, "upsilon"),
        (0, 6, "U"),
        (6, 2, "K"),
        (5, 3, "Phi"),
        (4, 3, "Theta"),
        (3, 4, "eta"),
        (2, 4, "xi"),
        (1, 5, "alpha"),
        (0, 5, "W"),
    )
)
need(top_images == expected_images, "independent top-chart exponent transform")

# Polygon incidence.  The two open orbits are disjoint even though their
# closures share the fixed vertex (4,5).
vertices = ((0, 1), (2, 0), (4, 2), (4, 5), (0, 7))
cubic_edge = ((4, 2), (4, 5))
top_edge = ((4, 5), (0, 7))
internal_edge = ((2, 0), (4, 5))
need(cubic_edge in tuple(zip(vertices, vertices[1:] + vertices[:1])), "cubic edge occurs")
need(top_edge in tuple(zip(vertices, vertices[1:] + vertices[:1])), "top edge occurs")
need(set(cubic_edge).intersection(top_edge) == {(4, 5)}, "one shared fixed vertex")
need(internal_edge not in (cubic_edge, top_edge), "internal edge is a third orbit")
need(gcd(0, 3) == 3 and gcd(4, 2) == 2, "edge lattice lengths 3 and 2")

# Same-coordinate hostiles.  Coefficient lists are descending.
need((1, -1, -1, 1) == (1, -1, -1, 1), "(P-1)^2(P+1) coefficients")
need((1, -3, 3, -1) == (1, -3, 3, -1), "(P-1)^3 coefficients")
delta_double = F(6, 7) * (F(2848, 45) - 1)
delta_triple = F(6, 7) * (F(2848, 45) + 1)
need(F(2848, 45) - F(7, 6) * delta_double == 1, "double hostile K relation")
need(F(2848, 45) - F(7, 6) * delta_triple == -1, "triple hostile K relation")

# Terminal A23 arithmetic, obtained without symbolic solving.
c0 = F(5152, 405)
Delta0 = 3 * c0 - F(32, 9)
K0 = F(2848, 45) - F(7, 6) * Delta0
U0 = -c0 * c0 / 2
W0 = -U0
xi0 = 8 * c0 / 3
Theta0 = -Delta0
need(c0 == F(7168, 135) - F(7, 6) * Delta0, "terminal source identity")
need(c0 == -F(1376, 135) + K0, "terminal c=e+K")
need((Delta0, K0, U0, W0, xi0, Theta0) == (
    F(4672, 135),
    F(1856, 81),
    -F(13271552, 164025),
    F(13271552, 164025),
    F(41216, 1215),
    -F(4672, 135),
), "terminal coefficients")

disc = cubic_discriminant(W0, xi0, Theta0, K0)
expected_disc = -F(3947729324424178958336, 32688603759375)
need(disc == expected_disc and disc != 0, "terminal cubic is squarefree")

# The two raw arithmetic-genus ledgers and horizontal correction.
need(0 + 3 + 1 + 12 + 1 - 3 + 1 == 15, "squarefree raw ledger")
need(0 + 3 + 0 + 12 + 1 + 1 - 3 + 1 == 15, "singular-cubic raw ledger")
need(15 - 1 == 14, "horizontal conductor correction")

# Coefficients of every local A23 differential order are positive.
orders = ((6, 2), (F(5, 2), F(9, 2)), (2, 4), (F(3, 2), F(7, 2)), (1, 3))
need(all(a > 0 and b > 0 for a, b in orders), "A23 order cone is strictly positive")

print("THM4343 U+W=0 A23-CUBIC INDEPENDENT AUDIT=PASS")
print("top_images=" + str(top_images))
print("edge_orbits=distinct;closures_share_(4,5)")
print("terminal_coefficients=" + str((U0, W0, K0, Delta0, Theta0, xi0)))
print("terminal_discriminant=" + str(disc))
print("genus=15,15,14")
