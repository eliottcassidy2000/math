"""Two-rung certificates for trinomial smaller endpoint degree min(a,c)<=8.

The parameter g is an indeterminate, not a bounded search variable. Every
shifted resultant coefficient is tested exactly over Q. A separate finite
exception list covers the opposite-endpoint orientation. Explicit gates
remain active under Python -O. Total exponent width a+c is unbounded;
the filename's width8 suffix is only a shorthand for min(a,c)<=8.
"""
from collections import defaultdict
from hashlib import sha256
import json
from math import factorial, gcd
from pathlib import Path

import sympy as sp

G, V = sp.symbols("g v")
DOMAIN = sp.QQ.poly_ring(G)
ROOT = Path(__file__).resolve().parents[1]
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def semigroup_contains(n, A, B):
    if n < 0:
        return False
    z0 = 0 if A == 1 else (n * pow(B, -1, A)) % A
    return B * z0 <= n


def symbolic_row(a, A, B, k):
    reps = [(y, z) for z in range(k * a // B + 1)
            if (k * a - B * z) % A == 0
            for y in [(k * a - B * z) // A]]
    require(bool(reps), "nonempty symbolic row")
    ymin = min(y for y, _ in reps)
    terms = {}
    for y, z in reps:
        require((y - ymin) % B == 0, "correct invariant exponent")
        terms[((y - ymin) // B,)] = sp.ff(k * G, y + z) / (
            factorial(y) * factorial(z))
    full = sp.Poly.from_dict(terms, V, domain=DOMAIN)
    content, primitive = full.primitive()
    return ymin, content, primitive


def direct_moment(a, b, c, m):
    """Independent literal Laurent convolution with middle-coefficient labels."""
    current = {(0, 0): 1}
    for _ in range(m):
        next_row = defaultdict(int)
        for (charge, exponent), coefficient in current.items():
            next_row[charge - a, exponent] += coefficient
            next_row[charge + b, exponent + 1] += coefficient
            next_row[charge + c, exponent] += coefficient
        current = next_row
    return sp.Poly.from_dict({(j,): value for (charge, j), value in current.items()
                             if charge == 0}, V, domain=sp.QQ)


def first_rows(a, b, c, m):
    """Direct unordered channels, independent of the semigroup normal form."""
    return [(x, y, m - x - y) for x in range(m + 1)
            for y in range(m - x + 1)
            if -a * x + b * y + c * (m - x - y) == 0]


def primitive_integer(poly):
    out = sp.Poly(poly, G, domain=sp.QQ).clear_denoms()[1].primitive()[1]
    return out if out.LC() > 0 else -out


def coefficients_ascending(poly):
    return [int(poly.nth(i)) for i in range(poly.degree() + 1)]


def main():
    manifest = {"status": "FINITE-EXACT symbolic certificates in Q[g]",
                "scope": "primitive trinomial min(a,c)<=8, arbitrary complex nonzero coefficients",
                "families": [], "exceptional_supports": []}
    family_count = 0
    shifted_coefficient_count = 0
    for a in range(1, 9):
        for A in range(1, a + 1):
            for B in range(A + 1, a + 1):
                if gcd(A, B) != 1 or not semigroup_contains(a - A * B, A, B):
                    continue
                gmin = a // A + 1
                y1, c1, P = symbolic_row(a, A, B, 1)
                y2, c2, Q = symbolic_row(a, A, B, 2)
                require(P.degree() >= 1, "collided first row")
                # These are polynomial resultants, not specialized numerical gcds.
                R = primitive_integer(P.resultant(Q))
                shifted = sp.Poly(R.as_expr().subs(G, G + gmin).expand(), G)
                coeffs = coefficients_ascending(shifted)
                require(all(v > 0 for v in coeffs), "positive shifted resultant coefficients")
                require(shifted.degree() == R.degree(), "no lost top coefficient")
                # Exact polynomial reconstruction checks the frozen certificate.
                require(sp.expand(sum(v * G ** i for i, v in enumerate(coeffs))
                                  - R.as_expr().subs(G, G + gmin)) == 0,
                        "shift certificate identity")
                record = {"a": a, "A": A, "B": B, "gmin": gmin,
                          "first_ymin": y1, "second_ymin": y2,
                          "first_content": str(sp.factor(c1)),
                          "second_content": str(sp.factor(c2)),
                          "first_primitive": str(P.as_expr()),
                          "second_primitive": str(Q.as_expr()),
                          "resultant_primitive_ascending": coefficients_ascending(R),
                          "resultant_factorization": str(sp.factor(R.as_expr())),
                          "shifted_ascending": coeffs}
                manifest["families"].append(record)
                family_count += 1
                shifted_coefficient_count += len(coeffs)
                print(f"family={(a,A,B)} g>={gmin} degrees={(P.degree(),Q.degree())} "
                      f"resultant_degree={R.degree()} positive_shift_coefficients={len(coeffs)}")

    require(family_count == 30, "complete symbolic family count")

    # If a>8 and c<=8, collision gives g|c-b and B<g, hence
    # a<=g(g-1)-c<=(c-1)(c-2)-c. This bounded rectangle is exhaustive.
    exceptions = []
    tested_rectangle = 0
    for c in range(2, 9):
        for b in range(1, c):
            for a in range(9, (c - 1) * (c - 2) - c + 1):
                if gcd(a, gcd(b, c)) != 1:
                    continue
                tested_rectangle += 1
                g = gcd(a + b, a + c)
                A, B = (a + b) // g, (a + c) // g
                if not semigroup_contains(a - A * B, A, B):
                    continue
                require(g > B, "collision width inequality")
                Praw, Qraw = direct_moment(a, b, c, g), direct_moment(a, b, c, 2 * g)
                common = sp.gcd(Praw, Qraw)
                require(len(common.terms()) == 1, "exceptional torus coprimality")
                rows1, rows2 = first_rows(a, b, c, g), first_rows(a, b, c, 2 * g)
                y1, c1, P = symbolic_row(a, A, B, 1)
                y2, c2, Q = symbolic_row(a, A, B, 2)
                require(sp.expand(V ** y1 * c1.subs(G, g) *
                                  P.as_expr().subs({G: g, V: V ** B}) - Praw.as_expr()) == 0,
                        "first exception literal/convolution identity")
                require(sp.expand(V ** y2 * c2.subs(G, g) *
                                  Q.as_expr().subs({G: g, V: V ** B}) - Qraw.as_expr()) == 0,
                        "second exception literal/convolution identity")
                pp, qq = P.subs(G, g), Q.subs(G, g)
                resultant = sp.resultant(pp.as_expr(), qq.as_expr(), V)
                require(resultant != 0, "exception exact resultant")
                record = {"support": [-a, b, c], "g": g, "A": A, "B": B,
                          "first_rows": rows1, "second_rows": rows2,
                          "first_raw_moment": str(Praw.as_expr()),
                          "second_raw_moment": str(Qraw.as_expr()),
                          "resultant": str(resultant)}
                manifest["exceptional_supports"].append(record)
                exceptions.append((a, b, c))
    require(exceptions == [(9, 1, 6), (11, 1, 7), (13, 1, 8), (20, 1, 8), (12, 3, 8)],
            "exhaustive opposite-endpoint exceptional list")

    # Raw coefficient-polynomial controls across different carry signatures.
    raw_controls = 0
    for a, A, B in [(3, 1, 2), (4, 1, 2), (6, 2, 3), (7, 1, 3), (8, 2, 3)]:
        gmin = a // A + 1
        for g in [gmin, gmin + 2, gmin + 5]:
            b, c = g * A - a, g * B - a
            for k in (1, 2):
                ymin, content, P = symbolic_row(a, A, B, k)
                expected = sp.expand(V ** ymin * content.subs(G, g) *
                                     P.as_expr().subs({G: g, V: V ** B}))
                actual = direct_moment(a, b, c, k * g).as_expr()
                require(expected == actual, "independent symbolic versus literal Laurent row")
                raw_controls += 1

    # Readable one-carry proof, including its exact first-return equality case.
    first = sp.ff(G, 3) * V ** 3 / 6 + sp.ff(G, 2) * V
    second = (sp.ff(2 * G, 6) * V ** 6 / 720 + sp.ff(2 * G, 5) * V ** 4 / 24
              + sp.ff(2 * G, 4) * V ** 2 / 4 + sp.ff(2 * G, 3) / 6)
    remainder = sp.rem(sp.Poly(second, V, domain=sp.QQ.frac_field(G)),
                       sp.Poly(1 + (G - 2) * V ** 2 / 6, V,
                               domain=sp.QQ.frac_field(G))).as_expr()
    closed = 2 * G * (G - 1) * (2 * G - 1) * (23 * G ** 2 - 47 * G + 20) / (15 * (G - 2) ** 2)
    require(sp.cancel(remainder - closed) == 0, "one-carry remainder identity")
    require(sp.Poly((23 * G ** 2 - 47 * G + 20).subs(G, G + 4).expand(), G).all_coeffs()
            == [23, 137, 200], "one-carry positivity proof")
    require(direct_moment(3, 1, 5, 4).as_expr() == 4 * V ** 3 + 12 * V,
            "sharp first-row control")
    require(direct_moment(3, 1, 5, 8).rem(sp.Poly(V ** 2 + 3, V)).as_expr() == 560,
            "sharp second-row control")

    manifest["symbolic_family_count"] = family_count
    manifest["shifted_coefficient_count"] = shifted_coefficient_count
    manifest["opposite_endpoint_primitive_rectangle"] = tested_rectangle
    manifest["literal_convolution_controls"] = raw_controls
    manifest["one_carry_remainder"] = str(closed)
    encoded = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    target = ROOT / "05-knowledge/results/overnight_20260906_moments_width8_certificates.json"
    target.write_bytes(encoded)
    print("symbolic_families=", family_count)
    print("positive_shift_coefficients=", shifted_coefficient_count)
    print("opposite_endpoint_primitive_rectangle=", tested_rectangle)
    print("opposite_endpoint_exceptions=", exceptions)
    print("literal_convolution_controls=", raw_controls, "+ 10 exceptional rows + 2 sharp controls")
    print("one_carry_remainder=", closed)
    print("checks=", CHECKS)
    print("certificate_sha256=", sha256(encoded).hexdigest())
    print("ALL EXACT GATES PASSED; general two-rung coprimality and arbitrary-support min(M,N)>=3 remain OPEN")


if __name__ == "__main__":
    main()
