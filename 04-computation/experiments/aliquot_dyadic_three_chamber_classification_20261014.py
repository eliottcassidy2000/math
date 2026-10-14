#!/usr/bin/env python3
"""Exact certificates for the complete dyadic three-chamber classification.

No third-party dependencies. Symbolic identities and positive-coefficient
certificates establish the displayed interval inequalities. The exponent
window 1..12 only corroborates the unbounded proof in the companion note.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class Poly:
    """Minimal exact Z[x,u,v] ring."""
    def __init__(self, coefficients):
        if isinstance(coefficients, int):
            coefficients = {(0, 0, 0): coefficients}
        self.coefficients = {m: c for m, c in coefficients.items() if c}

    @staticmethod
    def coerce(other):
        return other if isinstance(other, Poly) else Poly(other)

    def __add__(self, other):
        out = self.coefficients.copy()
        for m, c in self.coerce(other).coefficients.items():
            out[m] = out.get(m, 0) + c
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.coefficients.items()})

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        out = {}
        for a, c in self.coefficients.items():
            for b, d in self.coerce(other).coefficients.items():
                m = tuple(i+j for i, j in zip(a, b))
                out[m] = out.get(m, 0) + c*d
        return Poly(out)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        result = Poly(1)
        for _ in range(exponent):
            result *= self
        return result

    def __eq__(self, other):
        return self.coefficients == self.coerce(other).coefficients

    def substitute(self, variable, replacement):
        variables = [Poly({tuple(int(i == j) for i in range(3)): 1}) for j in range(3)]
        variables[variable] = self.coerce(replacement)
        result = Poly(0)
        for monomial, coefficient in self.coefficients.items():
            term = Poly(coefficient)
            for var, power in zip(variables, monomial):
                term *= var**power
            result += term
        return result


def numerator(x, y, z):
    return ((z-1)*(y-1)*(2*x-1)+(z-1)*y*(2*y-1)+y*z*(2*z-1))


def denominator(x, y, z):
    return x*y+y*z+z*x-x-y-z+1


def positivity_certificate(polynomial):
    x, u, v = [Poly({tuple(int(i == j) for i in range(3)): 1}) for j in range(3)]
    # Coordinates X=x-2, U=u-2, V=v-2u are nonnegative on the full domain.
    shifted = polynomial.substitute(2, v+2*u).substitute(1, u+2).substitute(0, x+2)
    require(all(c > 0 for c in shifted.coefficients.values()), "nonnegative coefficient certificate")
    require(shifted.coefficients.get((0, 0, 0), 0) > 0, "strict positive constant")
    return {"substitution": "x=X+2, u=U+2, v=2(U+2)+V; X,U,V>=0",
            "constant": shifted.coefficients[(0, 0, 0)],
            "coefficients_X_U_V": [[list(m), c] for m, c in sorted(shifted.coefficients.items())]}


def symbolic_checks():
    x, u, v = [Poly({tuple(int(i == j) for i in range(3)): 1}) for j in range(3)]
    s, t = 1+u+v, u*v+u+v
    d = x*x*t-x*s+1
    a = 2*x*x*u*s-x*(u+1)*(v+2)+v+1
    r = 2*x*x*u*s-x*(u*v+v+2*u*u+2*u)+u+v
    b = 2*u*u+u+(u-1)*v
    k = x*x*(u+3)+x*(u-2)-1
    f = 14*x*x-3*x-2+(u-2)*(7*x*x-2+x*(2*u-1))
    identities = {
        "denominator": d == denominator(x, x*u, x*v),
        "positive_denominator_product": d == x*(x*u)*(x*v)-(x-1)*(x*u-1)*(x*v-1),
        "ascending_H_formula": x*a == numerator(x*v, x, x*u)+d,
        "reverse_H_formula": x*r == numerator(x*v, x*u, x)+d,
        "lower_bound": r-d == x*(x-1)*b-x*(v-1)+u+v-1,
        "B_minus_v": b-v == 2*u*u+u+(u-2)*v,
        "order_gap": a-r == (u-1)*(2*x*(u+1)-1),
        "upper_bound": 3*d-a == (v-2*u)*k+f,
        "positive_F_baseline": 14*x*x-3*x-2 == (x-2)*(14*x+25)+48,
        "reverse_u2_contradiction": (r-2*d).substitute(1, 2)
            == (4-10*x)+(v-4)*(1-x-2*x*x),
        "reverse_next_edge": (x*v-1)*(2*x-1)+2*x*v-1
            == x*(v*(2*x+1)-2),
    }
    # Reuse u as an independent second core y for these polynomial identities.
    repeated_d = denominator(x, u, u)
    identities["repeated_maximum_difference"] = (
        numerator(u, u, x)-numerator(u, x, u) == (x-u)*(2*x+2*u-1))
    identities["repeated_maximum_interval"] = (
        2*repeated_d-(u-x)*(2*x+2*u-1) == 2*x*x+4*x*u-3*x-3*u+2)
    require(all(identities.values()), "a symbolic identity failed")
    positive = {"R-D": positivity_certificate(r-d),
                "A-R": positivity_certificate(a-r),
                "3D-A": positivity_certificate(3*d-a)}
    return {"ring": "Z[x,u,v]", "identities": identities,
            "identity_count": len(identities),
            "positive_coefficient_certificates": positive,
            "unbounded_interval": "1<R/D<A/D<3 for x>=2,u>=2,v>=2u"}


def fraction_closure(exponents):
    slope, intercept = Fraction(1), Fraction(0)
    for exponent in exponents:
        core = 1 << exponent
        alpha = Fraction(core-1, core)
        slope, intercept = alpha*slope, alpha*intercept+2*core-1
    n = intercept/(1-slope)
    out = []
    for exponent in exponents:
        core = 1 << exponent
        out.append(n/core)
        n = Fraction(core-1, core)*n+2*core-1
    return tuple(out)


def corroboration():
    bound = 12
    counts, integral = Counter(), []
    formula_checks = 0
    for exponents in product(range(1, bound+1), repeat=3):
        a, b, c = exponents
        if exponents != min(exponents, (b, c, a), (c, a, b)):
            continue
        multiplicity_min = exponents.count(min(exponents))
        if len(set(exponents)) == 1:
            kind = "constant"
        elif multiplicity_min == 2:
            kind = "repeated_minimum"
        elif len(set(exponents)) == 2:
            kind = "repeated_maximum"
        else:
            kind = "distinct_ascending" if b<c else "distinct_reverse"
        counts[kind] += 1
        values = fraction_closure(exponents)
        cores = tuple(1 << e for e in exponents)
        d = denominator(*cores)
        explicit = tuple(Fraction(numerator(cores[i], cores[(i+1)%3], cores[(i+2)%3]), d)
                         for i in range(3))
        require(values == explicit, "independent rational composition")
        formula_checks += 1
        closes = all(p.denominator == 1 for p in values)
        expected = len(set(exponents)) == 1 or exponents == (1, 1, 4)
        require(closes == expected, "full classification small-window check")
        if closes:
            integral.append({"exponents": exponents, "cofactors": tuple(int(p) for p in values)})
        if kind.startswith("distinct_"):
            i = exponents.index(max(exponents))
            h = (values[i]+1)/min(cores)
            require(1<h<3, "largest-core interval")
    require(dict(counts) == {"constant": 12, "repeated_minimum": 66,
                            "repeated_maximum": 66, "distinct_ascending": 220,
                            "distinct_reverse": 220}, "case partition")
    require(sum(counts.values()) == (bound**3+2*bound)//3, "cyclic representatives")
    return {"exponent_bound": bound, "ordered_templates_covered": bound**3,
            "cyclic_representatives": sum(counts.values()), "cases": dict(counts),
            "independent_fraction_formula_checks": formula_checks,
            "integral_closures": integral,
            "note": "This is the previously checked small window, not an enlarged census."}


def main():
    result = {
        "status": "complete unbounded dyadic affine three-chamber classification proved",
        "actual_aliquot_consequence": "No three-cycle of distinct even integers all with prime odd part",
        "scope_exclusions": ["odd members", "composite odd parts", "general sociable three-cycles"],
        "virtual_exception": {"exponents": [1, 1, 4], "cofactors": [23, 13, 1],
                              "formal_cycle": [46, 26, 16, 46], "actual_s_16": 15},
        "prior_theorem_used": "aliquot_repeated_minimum_classification_20261014.md",
        "symbolic": symbolic_checks(),
        "finite_corroboration": corroboration(),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
