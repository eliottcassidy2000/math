#!/usr/bin/env python3
"""Exact symbolic identities and bounded corroboration of an unbounded proof.

The mathematical proof is in the companion note. The finite window below
is unchanged from the preceding prime-chamber packet; it is not the proof.
"""
from fractions import Fraction
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class Poly:
    """Small exact Z[x,t] ring, with no computer-algebra dependency."""
    def __init__(self, coefficients):
        if isinstance(coefficients, int):
            coefficients = {(0, 0): coefficients}
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
        for (a, b), c in self.coefficients.items():
            for (d, e), f in self.coerce(other).coefficients.items():
                m = a + d, b + e
                out[m] = out.get(m, 0) + c * f
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
        replacement = self.coerce(replacement)
        x, t = Poly({(1, 0): 1}), Poly({(0, 1): 1})
        return sum((c * (replacement if variable == 0 else x)**a
                    * (replacement if variable == 1 else t)**b
                    for (a, b), c in self.coefficients.items()), Poly(0))


def symbolic_checks():
    x, t = Poly({(1, 0): 1}), Poly({(0, 1): 1})
    d = (2*x*x-x)*t + (x-1)**2
    q = 6*x**4 - 6*x**3 - x*x + 3*x - 1
    require(q == (2*x*x-1)*(3*x*x-3*x+1), "Q factorization")
    residual = d*(-1+(1-t)*x)-q
    require(all(a >= 2 for (a, b) in residual.coefficients), "quotient congruence mod x^2")
    large_boundary = ((3*x-1)*d-q).substitute(1, x)
    require(large_boundary == x*(4*x*x-5*x+2), "large-ratio H<3 boundary")
    e = 3*x*x*d-t*(q+d)
    bracket = -(2*x-1)*t*t+(3*x*x-1)*t+3*x*(x-1)**2
    require(e == x*bracket, "small-ratio H bound")
    lower = (x*x+x-1)*t+3*x*(x-1)**2
    require(bracket-lower == (2*x-1)*t*(x-t), "small-ratio nonnegative remainder")
    d2 = d.substitute(1, 2)
    require(q-d2*(x*x-x-1) == x*x*(x*x+3*x-1), "t=2 contradiction")
    require(q.substitute(0, 1) == 1, "Q mod x-1")
    return {"identities_checked": 7, "ring": "Z[x,t]",
            "quotient_congruence": "K = -1 + (1-t)x (mod x^2)",
            "sector_bounds": {"t>=x": "0<H<3", "0<t<=x": "0<H<3x/t"}}


def bounded_checks():
    checks = 0
    constants = 0
    exceptions = []
    sector_checks = 0
    for a in range(1, 65):
        x = 1 << a
        q = (2*x*x-1)*(3*x*x-3*x+1)
        for c in range(a, 3*a+3):
            z = 1 << c
            t = 1 << (c-a)
            d = (2*x-1)*z+(x-1)**2
            closes = q % d == 0
            expected = c == a or (a, c) == (1, 4)
            require(closes == expected, "classification corroboration")
            if closes:
                if c == a:
                    constants += 1
                else:
                    exceptions.append({"a": a, "c": c, "x": x, "t": t,
                                       "D": d, "Q": q, "quotient_K": q//d})
            if c > a:
                h = Fraction(q+d, x*d)
                if t >= x:
                    require(h < 3, "large-ratio sector")
                if t <= x:
                    require(h < Fraction(3*x, t), "small-ratio sector")
                if t < x:
                    require(x >= 2*t, "dyadic separation")
                    if t >= 8:
                        require(Fraction(x+1-t)-Fraction(3*x,t) >= t-5 > 0,
                                "t>=8 contradiction gap")
                sector_checks += 1
            checks += 1
    require(divmod(21463, 529) == (40, 303), "t=4 residual case")
    require(sum(d for d in range(1, 16) if 16 % d == 0) == 15, "virtual-loop gate")
    return {"a_range": [1, 64], "c_range": "a <= c <= 3a+2",
            "templates_checked": checks, "constant_closures": constants,
            "nonconstant_integral_closures": exceptions,
            "sector_checks": sector_checks,
            "t4_residual": {"x": 8, "D": 529, "Q": 21463, "remainder": 303},
            "virtual_loop_true_successor_of_16": 15}


def main():
    result = {
        "status": "proved repeated-minimum classification; full dyadic classification remains open",
        "theorem": "For a>=1,c>=a, integral closure occurs iff c=a or (a,c)=(1,4)",
        "symbolic": symbolic_checks(),
        "bounded_corroboration": bounded_checks(),
        "proof_scope": "No exponent bound in theorem; finite checks only corroborate the proof",
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
