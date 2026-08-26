#!/usr/bin/env python3
"""Independent exact geometry audit for the THM-4147 repeated top-edge wall.

Read-only with respect to the repository.  This rebuilds the cleared source
fibre directly, rather than importing the producer's support routines.
"""

from itertools import permutations
from math import gcd
import hashlib
import json

import sympy as sp


s, p, Q, a, z, tau, W, q = sp.symbols("s p Q a z tau W q")
lam, alpha, eps, K, Phi, Delta, Theta, eta, C, B = sp.symbols(
    "lambda alpha epsilon K Phi Delta Theta eta C B"
)


def convex_hull(points):
    pts = sorted(set(points))

    def cross(o, x, y):
        return (x[0] - o[0]) * (y[1] - o[1]) - (x[1] - o[1]) * (y[0] - o[0])

    lo = []
    for pt in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], pt) <= 0:
            lo.pop()
        lo.append(pt)
    hi = []
    for pt in reversed(pts):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], pt) <= 0:
            hi.pop()
        hi.append(pt)
    return lo[:-1] + hi[:-1]


def area2(poly):
    return abs(sum(
        poly[i][0] * poly[(i + 1) % len(poly)][1]
        - poly[(i + 1) % len(poly)][0] * poly[i][1]
        for i in range(len(poly))
    ))


def boundary(poly):
    return sum(
        gcd(
            abs(poly[(i + 1) % len(poly)][0] - poly[i][0]),
            abs(poly[(i + 1) % len(poly)][1] - poly[i][1]),
        )
        for i in range(len(poly))
    )


def weighted_initial(poly, wa, wz):
    P = sp.Poly(sp.expand(poly), a, z)
    weights = [(wa * i + wz * j, c * a**i * z**j) for (i, j), c in P.terms()]
    m = min(v for v, _ in weights)
    return m, sp.factor(sum(t for v, t in weights if v == m))


def z_valuation(expr):
    P = sp.Poly(sp.expand(expr), z)
    return min(i[0] for i, c in P.terms() if c != 0)


def compose(x, y):
    return tuple(x[y[i]] for i in range(len(x)))


def inverse(x):
    out = [0] * len(x)
    for i, xi in enumerate(x):
        out[xi] = i
    return tuple(out)


def perm_index(x):
    seen = set()
    cycles = 0
    for i in range(len(x)):
        if i not in seen:
            cycles += 1
            j = i
            while j not in seen:
                seen.add(j)
                j = x[j]
    return len(x) - cycles


def support(x):
    return {i for i, xi in enumerate(x) if xi != i}


# Direct cleared fibre, with no producer helper imported.
H = (
    lam * p + alpha * p**2 + eps * p**3 + K * s**2 * p**2
    + Phi * s * p**3 + Delta * p**4 + Theta * s**2 * p**3
    + eta * s * p**3 * (p - s**2)
)
F = sp.expand((s**2 - p) * (1 - Q * H) - sp.Rational(1, 2) * Q * s**2)
FP = sp.Poly(F, s, p)
support_points = [m for m, coeff in FP.terms() if coeff != 0]
hull = convex_hull(support_points)
expected_hull = [(0, 1), (2, 0), (5, 3), (1, 5), (0, 5)]
assert hull == expected_hull, hull
A2 = area2(hull)
bdry = boundary(hull)
pick_genus = (A2 - bdry + 2) // 2
assert (A2, bdry, pick_genus) == (31, 11, 11)

# Exact face extraction by primitive inward inequalities.
faces = {}
edge_data = [
    ("e1", (1, 2), 2, 1, 1),
    ("cubic", (-1, 1), -2, 3, 2),
    ("repeated", (-1, -2), -11, 2, 8),
    ("delta", (0, -1), -5, 1, 4),
    ("vertical_affine", (1, 0), 0, 4, 1),
]
for name, (u, v), c0, lattice_len, residue_index in edge_data:
    terms = []
    for (i, j), coeff in FP.terms():
        if u * i + v * j == c0:
            terms.append(coeff * s**i * p**j)
    faces[name] = sp.factor(sum(terms))

assert sp.expand(faces["cubic"] - s**2 * ((1 - Q / 2) - K * Q * (s*p)**2 + eta * Q * (s*p)**3)) == 0
assert sp.expand(faces["repeated"] - Q * eta * s * p**3 * (s**2 - p)**2) == 0
assert sp.expand(faces["delta"] - Q * p**5 * (Delta + eta*s)) == 0

# Direct toric substitution and exact strict-transform identity.
r = 1 - a
L = sp.cancel(z**11 * F.subs({s: 1/z, p: r/z**2}))
L = sp.expand(L)
L_expected = sp.expand(
    a*z**9 - sp.Rational(1, 2)*Q*z**9 + Q*eta*a**2*r**3
    - Q*a*r**3*(Delta*r + Theta)*z - Q*Phi*a*r**3*z**2
    - Q*a*(eps*r**3 + K*r**2)*z**3 - Q*alpha*a*r**2*z**5
    - Q*lam*a*r*z**7
)
assert sp.expand(L - L_expected) == 0

# Four local Newton hierarchies.  Theta=C-Delta and eps+K=B.
L_C = sp.expand(L.subs(Theta, C - Delta))
in_A = weighted_initial(L_C, 1, 1)
assert in_A[0] == 2 and sp.expand(in_A[1] - Q*a*(eta*a - C*z)) == 0

L_B = sp.expand(L_C.subs(C, 0))
in_B1 = weighted_initial(L_B, 2, 1)
assert in_B1[0] == 4 and sp.expand(in_B1[1] - Q*a*(eta*a - Phi*z**2)) == 0
in_B2 = weighted_initial(L_B, 7, 1)
assert in_B2[0] == 9 and sp.expand(in_B2[1] + Q*z**2*(2*Phi*a + z**7)/2) == 0

L_Bsym = sp.expand(L_B.subs(K, B - eps))
L_Crow = sp.expand(L_Bsym.subs(Phi, 0))
in_C1 = weighted_initial(L_Crow, 3, 1)
assert in_C1[0] == 6 and sp.expand(in_C1[1] - Q*a*(eta*a - B*z**3)) == 0
in_C2 = weighted_initial(L_Crow, 6, 1)
assert in_C2[0] == 9 and sp.expand(in_C2[1] + Q*z**3*(2*B*a + z**6)/2) == 0

L_D = sp.expand(L_Crow.subs(B, 0))
in_D = weighted_initial(L_D, 9, 2)
assert in_D[0] == 18 and sp.expand(in_D[1] - Q*(2*a**2*eta - z**9)/2) == 0

# L_a valuations on the smooth branch leading terms.
La_A = sp.diff(L_C, a)
assert z_valuation(La_A.subs(a, C*z/eta)) == 1
assert z_valuation(La_A.subs(a, -z**8/(2*C))) == 1
La_B = sp.diff(L_B, a)
assert z_valuation(La_B.subs(a, Phi*z**2/eta)) == 2
assert z_valuation(La_B.subs(a, -z**7/(2*Phi))) == 2
La_C = sp.diff(L_Crow, a)
assert z_valuation(La_C.subs(a, B*z**3/eta)) == 3
assert z_valuation(La_C.subs(a, -z**6/(2*B))) == 3

# The cusp has z=tau^2, a=c*tau^9; L_a starts in degree 9.
cusp_c = sp.symbols("cusp_c", nonzero=True)
La_D_tau = sp.expand(sp.diff(L_D, a).subs({a: cusp_c*tau**9, z: tau**2}))
tau_val = min(i[0] for i, coeff in sp.Poly(La_D_tau, tau).terms() if coeff != 0)
assert tau_val == 9

rows = {
    "A": {"packet": (7, 7, 4, 2, 2, 2, 1), "delta": 1, "L": 23},
    "B": {"packet": (6, 6, 4, 2, 2, 2, 1), "delta": 2, "L": 21},
    "C": {"packet": (5, 5, 4, 2, 2, 2, 1), "delta": 3, "L": 19},
    "D": {"packet": (7, 4, 2, 2, 2, 1), "delta": 4, "L": 16},
}
for name, row in rows.items():
    row["genus"] = pick_genus - row["delta"]
    row["defect"] = sum(e - 1 for e in row["packet"])
    assert row["defect"] == 2 * row["genus"] - 2
    rational_packet = tuple(e for e in row["packet"] if e != 2)
    # Here the only e=2 entries are precisely the three conjugates of the cubic place.
    assert row["packet"].count(2) == 3
    row["finite_n"] = sum(rational_packet)
    row["beta"] = 3
    row["full_n"] = row["finite_n"] + 6
    assert row["L"] == row["finite_n"] + row["beta"] + 1
    assert row["full_n"] - row["L"] == 2
    assert row["defect"] > 4

# Carrier separability/Hurwitz polynomial.
carrier_poly = -eta*W**3 + K*W**2 - (q - sp.Rational(1, 2))
carrier_disc = sp.factor(sp.discriminant(carrier_poly, W))
carrier_disc_expected = sp.expand((q-sp.Rational(1, 2)) * (4*K**3 - 27*eta**2*(q-sp.Rational(1, 2))))
assert sp.expand(carrier_disc - carrier_disc_expected) == 0

# Independent exhaustive finite check of the general commutator-overlap inequality through S_5.
hostile_pairs = 0
for n in range(1, 6):
    perms = list(permutations(range(n)))
    for x in perms:
        ix = inverse(x)
        sx = support(x)
        for y in perms:
            hostile_pairs += 1
            comm = compose(compose(compose(x, y), ix), inverse(y))
            assert perm_index(comm) <= 2 * len(sx & support(y))

result = {
    "support": sorted(support_points),
    "hull": hull,
    "area2": A2,
    "boundary": bdry,
    "pick_genus": pick_genus,
    "faces": {k: str(v) for k, v in faces.items()},
    "initials": {
        "A": [in_A[0], str(in_A[1])],
        "B_first": [in_B1[0], str(in_B1[1])],
        "B_second": [in_B2[0], str(in_B2[1])],
        "C_first": [in_C1[0], str(in_C1[1])],
        "C_second": [in_C2[0], str(in_C2[1])],
        "D": [in_D[0], str(in_D[1])],
    },
    "rows": rows,
    "carrier_discriminant": str(carrier_disc),
    "commutator_hostile_pairs": hostile_pairs,
}
payload = json.dumps(result, sort_keys=True, separators=(",", ":"))
print(json.dumps(result, sort_keys=True, indent=2))
print("semantic_sha256=" + hashlib.sha256(payload.encode()).hexdigest())
