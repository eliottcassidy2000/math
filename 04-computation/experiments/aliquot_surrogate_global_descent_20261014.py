#!/usr/bin/env python3
"""Exact certificates for global descent of the even prime-cofactor surrogate.

No third-party dependencies. The theorem is proved by the symbolic branch
identities and inequalities in the companion note, not by a cycle census.
All checks stay active with Python -O.
"""
from collections import Counter
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class Poly:
    """Minimal exact Z[r,y,q] ring."""
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


def symbolic_checks():
    r, y, q = [Poly({tuple(int(i == j) for i in range(3)): 1}) for j in range(3)]
    x = r*y
    n0 = x*(y*q-1)
    p1 = (x-1)*q+r
    n1 = y*p1
    n2 = (y-1)*p1+2*y-1
    delta = -(q-2)*(x+y-1)+1-r
    identities = {
        "first_edge": (x-1)*(y*q-1)+2*x-1 == n1,
        "two_step_difference": n2-n0 == delta,
        "q_one_second_iterate": n2.substitute(2, 1) == r*y*y-r+y,
    }
    start = r*y*(y-1)
    second = r*y*y-r+y
    # In the r<y branch, multiply by r to clear y/r exactly.
    r_third = (r-1)*second+r*(2*r-1)
    r_delta3 = r*(-y*y+r*y+y+r)-y
    identities["r_less_y_third_difference"] = r_third-r*start == r_delta3
    identities["r_less_y_upper_bound_gap"] = (
        -r*(y*y-3*y+4)-2*r_delta3 == (y-2*r)*(r*(y+1)+2))
    identities["r_less_y_strict_bound"] = y*y-3*y+4 == (y-4)*(y+1)+8
    # In the r>y branch, multiply by y to clear r/y exactly.
    y_third = (y-1)*second+y*(2*y-1)
    y_delta3 = -r*y+r+3*y*y-2*y
    identities["r_greater_y_third_difference"] = y_third-y*start == y_delta3
    identities["R_at_least_four_bound"] = (
        -q*(y-1)+3*y-2 == 2-y-(q-4)*(y-1))
    twice_start = 2*y*y*(y-1)
    twice_third = twice_start+y
    twice_fourth = (y-1)*(2*y*y-2*y+1)+2*y-1
    identities["R_two_third_iterate"] = (
        y_third.substitute(0, 2*y) == y*twice_third)
    identities["R_two_fourth_difference"] = (
        twice_fourth-twice_start == -(2*y-1)*(y-2))
    identities["equal_cores_power_hit"] = second.substitute(0, y) == y**3
    # r is an independent power M for these denominator-cleared identities.
    identities["power_first_edge"] = (r-1)+2*r-1 == 3*r-2
    identities["power_second_edge"] = (3*r-2)+6 == 3*r+4
    identities["power_third_edge"] = (3*r+4)+12 == 3*r+16
    identities["quadratic_escape_positive"] = (
        r*r-4*r-16 == (r-8)*(r+4)+16)
    identities["cubic_escape_positive"] = (
        y**3-4*y*y-16 == (y-8)*(y*y+4*y+32)+240)
    require(all(identities.values()), "symbolic identity failed")
    return {"ring": "Z[r,y,q]", "identities": identities,
            "identity_count": len(identities),
            "inequality_domains": {
                "two_step": "q>=3, x>=2, y>=2, r=x/y>=2",
                "r_less_y": "y>=4, 2<=r<=y/2",
                "R_at_least_four": "y>=4, R=r/y>=4",
                "R_two": "y>=4",
                "quadratic_escape": "x>=8",
                "cubic_escape": "y>=8"}}


def core(n):
    require(n > 0, "positive core argument")
    return n & -n


def f(n):
    require(n > 0 and n % 2 == 0, "positive even state")
    x = core(n)
    return (x-1)*(n//x)+2*x-1


def orbit(n, steps):
    out = [n]
    for _ in range(steps):
        out.append(f(out[-1]))
    return out


EXCEPTIONS = (2, 4, 8, 12, 16)
CYCLES = ((8, 22, 14, 10), (16, 46, 26))


def descent_certificate(n):
    x, p = core(n), n//core(n)
    if p == 2*x-1:
        require(f(n) == n, "fixed point")
        return {"n": n, "case": "fixed", "steps": 0, "endpoint": n}
    if n in EXCEPTIONS:
        return {"n": n, "case": "explicit_terminal_exception", "steps": 0, "endpoint": n}
    if p > 2*x-1:
        kind, steps, expected = "immediate_descent", 1, n-p+2*x-1
    elif p == 1:
        require(x >= 32, "small powers are explicit exceptions")
        kind, steps, expected = "power_descent", 3, 3*x//4+4
        values = orbit(n, 3)
        require(core(values[1]) == core(values[2]) == 2, "power successor cores")
    elif p == x-1:
        require(x >= 8, "small quadratic escape is explicit")
        kind, steps, expected = "quadratic_power_escape", 4, 3*x*x//4+4
        require(f(n) == x*x, "quadratic power hit")
    else:
        require(1 < p < 2*x-1, "remaining increasing range")
        y, q = core(p+1), (p+1)//core(p+1)
        require(y < x and x % y == 0 and q % 2 == 1, "first dyadic split")
        r = x//y
        values = orbit(n, 2)
        require(values[1] == y*((x-1)*q+r), "first edge formula")
        require(core(values[1]) == y, "first successor exact core")
        require(values[2]-n == -(q-2)*(x+y-1)+1-r, "two-step formula")
        if q >= 3:
            kind, steps, expected = "q_at_least_three", 2, n-(q-2)*(x+y-1)+1-r
        else:
            require(q == 1 and y >= 4, "Mersenne-shaped remaining branch")
            require(values[2] == r*y*y-r+y, "q-one second iterate")
            if r < y:
                require(y % r == 0 and r <= y//2, "smaller dyadic quotient")
                require(core(values[2]) == r, "smaller-quotient exact core")
                kind, steps = "r_less_than_y", 3
                expected = n-y*y+r*y+y-y//r+r
                require(2*(expected-n) <= -(y*y-3*y+4), "third-step upper bound")
            elif r > y:
                require(r % y == 0 and core(values[2]) == y, "larger-quotient exact core")
                R = r//y
                third = f(values[2])
                require(third-n == -R*(y-1)+3*y-2, "larger quotient formula")
                if R == 2:
                    require(third == n+y and core(third) == y, "twice-quotient exact core")
                    kind, steps, expected = "r_twice_y", 4, n-(2*y-1)*(y-2)
                else:
                    require(R >= 4, "dyadic larger-quotient partition")
                    kind, steps, expected = "r_at_least_four_y", 3, third
                    require(expected-n <= 2-y, "larger-quotient strict bound")
            else:
                require(values[2] == y**3, "cubic power hit")
                if y == 4:
                    require(n == 48, "small cubic escape")
                    kind, steps, expected = "cubic_48_escape", 6, 46
                else:
                    require(y >= 8, "remaining cubic power size")
                    kind, steps, expected = "cubic_power_escape", 5, 3*y**3//4+4
    values = orbit(n, steps)
    require(values[-1] == expected, "predicted descent endpoint")
    require(expected < n, "strict block descent")
    require(1 <= steps <= 6, "bounded block length")
    return {"n": n, "case": kind, "steps": steps, "endpoint": expected,
            "orbit_block": values}


def exceptional_checks():
    for cycle in CYCLES:
        require(all(f(n) == cycle[(i+1)%len(cycle)] for i, n in enumerate(cycle)),
                "explicit cycle")
        require(any(n//core(n) == 1 for n in cycle), "cofactor-one gate")
    chains = [[2,4,10], [4,10], [12,16], [48,76,64,190,98,52,46]]
    for chain in chains:
        require(all(f(n) == m for n, m in zip(chain, chain[1:])), "explicit exceptional chain")
    require(f(120) == 120 and 120//core(120) == 15, "composite-cofactor fixed point")
    return {"cycles": CYCLES, "exceptional_chains": chains,
            "surrogate_only_fixed_example": {"n": 120, "odd_part": 15}}


def bounded_branch_checks():
    bound = 4096
    counts = Counter()
    first_examples = {}
    for n in range(2, bound+1, 2):
        certificate = descent_certificate(n)
        kind = certificate["case"]
        counts[kind] += 1
        first_examples.setdefault(kind, certificate)
    expected_cases = {"fixed", "explicit_terminal_exception", "immediate_descent",
                      "power_descent", "quadratic_power_escape", "q_at_least_three",
                      "r_less_than_y", "r_twice_y", "r_at_least_four_y",
                      "cubic_48_escape", "cubic_power_escape"}
    require(set(counts) == expected_cases, "all proof branches covered")
    return {"input_bound": bound, "even_inputs": bound//2,
            "purpose": "Branch regression inside the predecessor window, not a cycle census",
            "counts": dict(counts), "first_examples": first_examples}


def main():
    result = {
        "status": "Global eventual periodicity and complete F-cycle classification proved",
        "F": "For positive even n=x*p with exact dyadic core x, F(n)=(x-1)*p+2*x-1",
        "descent_lemma": "Every positive even n is fixed, is in {2,4,8,12,16}, or has F^j(n)<n for 1<=j<=6",
        "fixed_points": "x*(2*x-1) for dyadic x>=2",
        "coarse_bounds": {"first_terminal_entry_steps_at_most": "3*n",
                          "every_orbit_value_at_most": "729*n"},
        "actual_aliquot_consequence": "Every finite cycle entirely in even prime-odd-part chambers is constant",
        "scope_exclusions": ["general aliquot convergence", "cycles with odd members", "cycles with composite odd parts"],
        "novelty": "Independent elementary derivation; global novelty unconfirmed",
        "symbolic_certificates": symbolic_checks(),
        "explicit_orbits": exceptional_checks(),
        "bounded_branch_checks": bounded_branch_checks(),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

