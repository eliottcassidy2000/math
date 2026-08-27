#!/usr/bin/env python3
"""Independent hostile audit of THM-4251 hidden degree-24 exclusion.

This file deliberately imports no repository computation.  It rebuilds the
Eisenstein lattice, its legal symmetry group, the coefficient-form elliptic
addition/duplication calculation, four good finite-field specializations, and
the full-lattice projection rows used in the proposed subtraction.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from hashlib import sha256

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.fields import field


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Eisenstein arithmetic and the hidden Gram, rebuilt locally.
# ---------------------------------------------------------------------------


E = tuple[int, int]
V = tuple[E, E]


def e_add(x: E, y: E) -> E:
    return x[0] + y[0], x[1] + y[1]


def e_neg(x: E) -> E:
    return -x[0], -x[1]


def e_mul(x: E, y: E) -> E:
    # omega^2=-omega-1
    a, b = x
    c, d = y
    return a*c-b*d, a*d+b*c-b*d


def e_conj(x: E) -> E:
    # conjugate(omega)=omega^2=-1-omega
    return x[0]-x[1], -x[1]


def e_norm(x: E) -> int:
    return e_mul(x, e_conj(x))[0]


def e_trace(x: E) -> int:
    return 2*x[0]-x[1]


OMEGA: E = (0, 1)
OMEGA2: E = (-1, -1)
MINUS_OMEGA: E = (0, -1)
UNITS: tuple[E, ...] = (
    (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)
)
CROSS: E = (-4, -2)
GRAM = sp.Matrix(
    ((6, -3, -3, 0),
     (-3, 6, 3, -3),
     (-3, 3, 6, -3),
     (0, -3, -3, 6))
)


def degree(a: E, b: E) -> int:
    return 6*e_norm(a) + 6*e_norm(b) + e_trace(
        e_mul(e_mul(a, e_conj(b)), CROSS)
    )


def matrix_degree(vector: V) -> int:
    a, b = vector
    x = sp.Matrix((a[0], a[1], b[0], b[1]))
    return int((x.T*GRAM*x)[0])


def unit_action(vector: V, unit: E) -> V:
    a, b = vector
    return e_mul(unit, a), e_mul(unit, b)


def tau_action(vector: V) -> V:
    # T(a,b)=(-omega*b,a), so T^2 is scalar -omega.
    a, b = vector
    return e_mul(MINUS_OMEGA, b), a


def generated_orbit(seed: V) -> set[V]:
    orbit = {seed}
    queue = deque([seed])
    while queue:
        vector = queue.popleft()
        successors = [tau_action(vector)]
        successors.extend(unit_action(vector, unit) for unit in UNITS)
        for successor in successors:
            if successor not in orbit:
                orbit.add(successor)
                queue.append(successor)
    return orbit


def shell_audit() -> dict[str, object]:
    lam = sp.symbols("lambda")
    characteristic = sp.factor(GRAM.charpoly(lam).as_expr())
    # This is the ordinary real-coordinate Gram spectrum, not the relative
    # Hermitian spectrum.  Conflating the two would give a false charpoly.
    expected_characteristic = (lam-6)*(lam-3)*(lam**2-15*lam+18)
    require(sp.expand(characteristic-expected_characteristic) == 0,
            "integral Gram spectrum changed")

    # The Hermitian eigenvalues relative to N(a)+N(b) are 6+/-sqrt(12).
    # Hence q=24 gives N(a)+N(b)<12.  Since
    # N(m+n omega)>=(m^2+n^2)/2, every integer coordinate has abs <=4.
    vectors24: set[V] = set()
    vectors6: set[V] = set()
    for bound, sink24, sink6 in ((4, vectors24, vectors6),):
        for am in range(-bound, bound+1):
            for an in range(-bound, bound+1):
                for bm in range(-bound, bound+1):
                    for bn in range(-bound, bound+1):
                        vector = ((am, an), (bm, bn))
                        require(degree(*vector) == matrix_degree(vector),
                                "Hermitian and integral Grams disagree")
                        value = degree(*vector)
                        if value == 24:
                            sink24.add(vector)
                        if value == 6:
                            sink6.add(vector)

    # A widened hostile enumeration checks that the analytic box was not
    # implemented with an off-by-one error.
    widened24: set[V] = set()
    widened6: set[V] = set()
    for am in range(-6, 7):
        for an in range(-6, 7):
            for bm in range(-6, 7):
                for bn in range(-6, 7):
                    vector = ((am, an), (bm, bn))
                    value = degree(*vector)
                    if value == 24:
                        widened24.add(vector)
                    if value == 6:
                        widened6.add(vector)
    require(vectors24 == widened24 and vectors6 == widened6,
            "shell changed in widened box")
    require(len(vectors24) == 24 and len(vectors6) == 24,
            "degree-6/24 shell count is not 24")
    require(all(all(c % 2 == 0 for c in (*a, *b))
                for a, b in vectors24),
            "degree-24 shell has a non-even coordinate")
    halves = {
        ((a[0]//2, a[1]//2), (b[0]//2, b[1]//2))
        for a, b in vectors24
    }
    require(halves == vectors6, "halving is not a shell bijection")

    zero: V = ((0, 0), (0, 0))
    require(tau_action(tau_action(((1, 2), (3, 4))))
            == unit_action(((1, 2), (3, 4)), MINUS_OMEGA),
            "T^2=-omega failed")
    require(all(degree(*tau_action(v)) == degree(*v) for v in vectors24),
            "T failed to preserve the hidden degree")
    require(all(degree(*unit_action(v, u)) == degree(*v)
                for v in vectors24 for u in UNITS),
            "a target unit failed to preserve degree")

    seeds: tuple[V, V] = (((2, 0), (0, 0)), ((2, 0), (2, 0)))
    orbits = [generated_orbit(seed) for seed in seeds]
    require(all(len(orbit) == 12 for orbit in orbits),
            "a claimed legal orbit is not free of size 12")
    require(orbits[0].isdisjoint(orbits[1]), "the two legal orbits overlap")
    require(orbits[0] | orbits[1] == vectors24,
            "the two legal orbits do not exhaust degree 24")
    require(zero not in vectors24, "zero entered a positive shell")

    return {
        "gram_charpoly": str(characteristic),
        "hermitian_lambda_min": "6-sqrt(12)>2",
        "degree6": len(vectors6),
        "degree24": len(vectors24),
        "all_even": True,
        "orbit_sizes": tuple(sorted(map(len, orbits))),
        "seeds": seeds,
    }


# ---------------------------------------------------------------------------
# Characteristic-zero denominator and noncancellation audit.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SymbolicRow:
    name: str
    line: str
    line_norms: tuple[str, str, str, str]
    a_at_one_norm: str
    ab_resultant_norm: str
    numerator_at_zero_norm: str


def symbolic_audit() -> tuple[SymbolicRow, ...]:
    p, z, u = sp.symbols("p z u")
    phi = z**4-z**2+1
    p_relation = p**2-(1+2*z-z**3)*p+1
    groebner = sp.groebner([p_relation, phi], p, z, order="lex")

    def reduce_coefficients(expression: sp.Expr) -> sp.Expr:
        polynomial = sp.Poly(sp.expand(expression), u)
        return sp.factor(sum(
            groebner.reduce(sp.expand(coefficient))[1]*u**exponent[0]
            for exponent, coefficient in polynomial.terms()
        ))

    def reduced(expression: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.fraction(sp.cancel(expression))
        nr = reduce_coefficients(numerator)
        dr = reduce_coefficients(denominator)
        require(dr != 0, "a field scalar denominator vanished")
        return sp.cancel(nr/dr)

    def norm_nonzero(expression: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.fraction(reduced(expression))
        numerator = sp.factor(groebner.reduce(sp.expand(numerator))[1])
        denominator = sp.factor(groebner.reduce(sp.expand(denominator))[1])
        n_norm = sp.factor(sp.resultant(
            sp.resultant(numerator, p_relation, p), phi, z
        ))
        d_norm = sp.factor(sp.resultant(
            sp.resultant(denominator, p_relation, p), phi, z
        ))
        require(n_norm != 0 and d_norm != 0,
                "an asserted nonzero field element has zero norm")
        return sp.factor(n_norm/d_norm)

    a_f = u-p**2
    b_f = u+p**3
    a_g = z**2*(p**2*u-1)
    b_g = -z**3*(1+p**3*u)
    delta_a = sp.Poly(a_g-a_f, u)
    delta_b = sp.Poly(b_g-b_f, u)
    slope_num = delta_b.coeff_monomial(u)
    slope_den = delta_a.coeff_monomial(u)
    proportionality = reduced(
        slope_den*(b_g-b_f)-slope_num*(a_g-a_f)
    )
    require(proportionality == 0,
            "f+Tf addition slope is not constant")
    slope = sp.cancel(slope_num/slope_den)
    a_sum = reduced((u-1)*slope**2-a_f-a_g)
    b_sum = reduced(slope*(a_f-a_sum)-b_f)

    rows: list[SymbolicRow] = []
    for name, a_line, b_line in (
        ("f", a_f, b_f),
        ("f+Tf", a_sum, b_sum),
    ):
        a_num, a_den = sp.fraction(sp.cancel(a_line))
        b_num, b_den = sp.fraction(sp.cancel(b_line))
        a_num = reduce_coefficients(a_num)
        a_den = reduce_coefficients(a_den)
        b_num = reduce_coefficients(b_num)
        b_den = reduce_coefficients(b_den)
        require(a_den != 0 and b_den != 0, "line scalar vanished")
        require(sp.Poly(a_num, u).degree() == 1,
                f"{name}: A is not linear in u")
        b_poly = sp.Poly(b_num, u)
        require(b_poly.degree() == 1, f"{name}: B is not linear in u")
        b0 = b_poly.coeff_monomial(1)
        b1 = b_poly.coeff_monomial(u)

        # If A=A(u)/(2t), B=B(u)/(2t), then duplication has
        # numerator 9*A^4-8*A*(u-1)*B^2 and denominator
        # 8*t*(u-1)*B^2, up to nonzero field scalars.  We check each
        # possible t-polynomial cancellation independently.
        ab_resultant = sp.resultant(a_num, b_num, u)
        n_zero = reduced(
            9*a_line.subs(u, 0)**4
            + 8*a_line.subs(u, 0)*b_line.subs(u, 0)**2
        )
        line_norms = tuple(str(norm_nonzero(value)) for value in (
            b0, b1, b0+b1, b0-b1
        ))
        a_one_norm = norm_nonzero(a_line.subs(u, 1))
        ab_norm = norm_nonzero(ab_resultant)
        n_zero_norm = norm_nonzero(n_zero)

        # The reverse line is b1+b0*u.  Their resultant is the negative
        # product L(1)L(-1), so the only possible extra reciprocal roots
        # are u=+/-1.  The two displayed norms exclude both.
        reverse_resultant = sp.expand(b1*b1-b0*b0)
        require(sp.expand(reverse_resultant+(b0+b1)*(b0-b1)) == 0,
                "linear reciprocal resultant identity changed")
        require(norm_nonzero(reverse_resultant) != 0,
                f"{name}: an extra reciprocal root exists")

        rows.append(SymbolicRow(
            name=name,
            line=str(sp.factor(b_num)),
            line_norms=line_norms,
            a_at_one_norm=str(a_one_norm),
            ab_resultant_norm=str(ab_norm),
            numerator_at_zero_norm=str(n_zero_norm),
        ))
    return tuple(rows)


# ---------------------------------------------------------------------------
# Direct exact finite-field group law, independent of repository code.
# ---------------------------------------------------------------------------


def finite_field_rows() -> tuple[tuple[object, ...], ...]:
    embeddings = (
        (313, 29, 135, 21),
        (349, 24, 246, 28),
        (373, 69, 297, 33),
        (397, 157, 161, 27),
    )
    answer: list[tuple[object, ...]] = []
    for q, z, p, scale in embeddings:
        require((z**4-z**2+1) % q == 0, "bad cyclotomic embedding")
        require((p**2-(1+2*z-z**3)*p+1) % q == 0,
                "bad p embedding")
        scale_den = (2*p**3+3*p**2-1) % q
        require((scale**6*scale_den-4) % q == 0,
                "bad scale embedding")
        require((4*z**3-2*z) % q != 0, "ramified z embedding")
        require((2*p-(1+2*z-z**3)) % q != 0,
                "ramified p embedding")

        _, t = field("t", GF(q))
        u = t*t
        inv2 = pow(2, -1, q)
        a_scale = scale**2*inv2 % q
        b_scale = scale**3*inv2 % q
        af = a_scale*(u-p**2)/t
        bf = b_scale*(u+p**3)/t
        ag = a_scale*z**2*(p**2*u-1)/t
        bg = -b_scale*z**3*(1+p**3*u)/t
        c = (u-1)/(2*t)

        def add(left: tuple[object, object],
                right: tuple[object, object]) -> tuple[object, object]:
            ax, by = left
            cx, dy = right
            slope = (dy-by)/(cx-ax)
            sx = c*slope**2-ax-cx
            sy = slope*(ax-sx)-by
            return sx, sy

        def double(point: tuple[object, object]) -> tuple[object, object]:
            ax, by = point
            dx = 9*t*ax**4/(2*(u-1)*by**2)-2*ax
            dy = 3*t*ax**2*(ax-dx)/((u-1)*by)-by
            return dx, dy

        representatives = ((af, bf), add((af, bf), (ag, bg)))
        local_rows = []
        for name, representative in zip(("f", "f+Tf"), representatives):
            doubled = double(representative)
            denominator = doubled[0].denom
            numerator = doubled[0].numer
            require(denominator.degree() == 7,
                    f"{name}: reduced finite denominator is not degree 7")
            require(numerator.gcd(denominator).degree() == 0,
                    f"{name}: field fraction was not reduced")

            reciprocal = denominator.ring.zero
            for (exponent,), coefficient in denominator.to_dict().items():
                reciprocal[(7-exponent,)] = coefficient*((-1)**exponent)
            common = denominator.gcd(reciprocal)
            require(common.degree() == 2,
                    f"{name}: reciprocal gcd degree is not two")
            common = common.monic()
            gcd_coeffs = tuple(
                int(common.to_dict().get((k,), 0)) % q
                for k in range(2, -1, -1)
            )
            require(gcd_coeffs == (1, 0, q-1),
                    f"{name}: reciprocal gcd is not t^2-1")

            # Directly reconstruct the predicted shape from the B-line.
            b_num = representative[1].numer
            predicted = t*(u-1)*b_num**2
            require(predicted.denom.degree() == 0,
                    f"{name}: predicted shape acquired a t-denominator")
            require(predicted.numer.monic() == denominator.monic(),
                    f"{name}: denominator differs from t(u-1)L^2")
            den_monic = denominator.monic()
            coefficient_row = tuple(
                int(den_monic.to_dict().get((k,), 0)) % q
                for k in range(7, -1, -1)
            )
            local_rows.append((name, gcd_coeffs, coefficient_row))
        digest = sha256(repr(tuple(local_rows)).encode("ascii")).hexdigest()
        answer.append(((q, z, p, scale), tuple(local_rows), digest))
    return tuple(answer)


# ---------------------------------------------------------------------------
# Full-lattice visible/hidden projection rows, reconstructed from the gluing.
# ---------------------------------------------------------------------------


def e_residue(x: E) -> E:
    return x[0] % 2, x[1] % 2


def projection_histograms() -> dict[int, Counter[int]]:
    # q(v)=16N(a)+4N(d), ell=(A,B) with
    # A == omega^2*d and B == d modulo 2L.
    a_counts: Counter[int] = Counter()
    for am in range(-5, 6):
        for an in range(-5, 6):
            value = e_norm((am, an))
            if value <= 10:
                a_counts[value] += 1

    d_counts: Counter[tuple[E, E, int]] = Counter()
    for dm in range(-10, 11):
        for dn in range(-10, 11):
            d = (dm, dn)
            value = e_norm(d)
            if value <= 42:
                residue = (e_residue(e_mul(OMEGA2, d)), e_residue(d))
                d_counts[(residue[0], residue[1], value)] += 1

    hidden_counts: Counter[tuple[int, E, E]] = Counter()
    boundary_hit = False
    # The eigenvalue bound gives |coordinate|<=12 for q(ell)<=168.
    # Enumerate one layer wider as an implementation hostile.
    for am in range(-13, 14):
        for an in range(-13, 14):
            for bm in range(-13, 14):
                for bn in range(-13, 14):
                    ell = ((am, an), (bm, bn))
                    qell = degree(*ell)
                    if 0 <= qell <= 168:
                        if max(map(abs, (am, an, bm, bn))) == 13:
                            boundary_hit = True
                        hidden_counts[(qell, e_residue(ell[0]),
                                       e_residue(ell[1]))] += 1
    require(not boundary_hit, "hidden enumeration touched widened boundary")

    histograms: dict[int, Counter[int]] = {}
    for target in (34, 42):
        histogram: Counter[int] = Counter()
        for (qell, ra, rb), hidden_multiplicity in hidden_counts.items():
            qv = 4*target-qell
            if qell < 0 or qv < 0:
                continue
            for (dra, drb, nd), d_multiplicity in d_counts.items():
                if (dra, drb) != (ra, rb):
                    continue
                remainder = qv-4*nd
                if remainder < 0 or remainder % 16:
                    continue
                na = remainder//16
                histogram[qell] += (
                    hidden_multiplicity*d_multiplicity*a_counts[na]
                )
        histograms[target] = histogram

    expected = {
        34: Counter({12:1536, 24:2304, 36:5952, 60:5760,
                     72:8064, 84:5376, 108:4992, 120:1728, 132:576}),
        42: Counter({12:672, 24:288, 36:2304, 60:288, 72:3456,
                     84:3072, 108:3744, 120:1728, 132:576,
                     156:672, 168:192}),
    }
    require(histograms == expected, "independent projection histogram differs")
    require(sum(histograms[34].values()) == 36288,
            "degree-34 full theta total changed")
    require(sum(histograms[42].values()) == 16992,
            "degree-42 full theta total changed")
    return histograms


def main() -> None:
    shell = shell_audit()
    symbolic = symbolic_audit()
    finite = finite_field_rows()
    histograms = projection_histograms()

    print("INDEPENDENT HIDDEN DEGREE-24 REFEREE AUDIT")
    print("imports_repository_code=False")
    print("shell", shell)
    for row in symbolic:
        print("symbolic", row)
    for embedding, rows, digest in finite:
        print("finite", embedding, rows, digest)
    print("projection34", sorted(histograms[34].items()))
    print("projection42", sorted(histograms[42].items()))
    previous34 = sum(histograms[34].values())-histograms[34][12]-histograms[34][132]
    previous42 = sum(histograms[42].values())-histograms[42][12]-histograms[42][168]
    new34 = previous34-histograms[34][24]
    new42 = previous42-histograms[42][24]
    require((previous34, previous42) == (34176, 16128),
            "previous THM-4247 remainder changed")
    require((new34, new42) == (31872, 15840),
            "degree-24 row subtraction changed")
    print("thm4247_only_row_subtraction",
          {"previous": (previous34, previous42),
           "removed_qell24": (histograms[34][24], histograms[42][24]),
           "new": (new34, new42)})
    print("current_stronger_residual", "governed_by_THM4249")
    print("gate", "common reciprocal roots=t:+1,-1 => Z/U=0; gate UZ!=0 excludes both")
    print("verdict", "PASS relative to THM-4230/4241/4247 geometric inputs")


if __name__ == "__main__":
    main()
