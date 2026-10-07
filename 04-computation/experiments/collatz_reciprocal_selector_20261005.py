"""Exact controls for reversal of shell selectors and the dyadic product.

No infinite product or divergent moment is numerically summed. Their proofs
are in the companion note; this checks finite identities and hostile cases.
"""
from fractions import Fraction as F
from math import comb, prod
import json

from collatz_refinement_energy_dual_20261005 import shell, multiply, evaluate
from collatz_floor_transport_deadlines_20261005 import step, word_weight


def need(ok, message):
    if not ok:
        raise ValueError(message)


def natural(n):
    need(type(n) is int and n >= 0, "exact natural required")


def exact(x):
    need(type(x) in (int, F), "exact rational required")
    return F(x)


def selector(n):
    natural(n)
    p = (F(1),)
    for d in range(1, n+1):
        x = shell(d)
        p = multiply(p, (-x/(1-x), 1/(1-x)))
    return p


def reversed_selector(n):
    natural(n)
    p = (F(1),)
    for d in range(1, n+1):
        x = shell(d)
        p = multiply(p, (1/(1-x), -x/(1-x)))
    return p


def norm(p):
    return sum(map(abs, p), F(0))


def trace_coordinate(alpha):
    alpha = exact(alpha)
    need(alpha != 0, "nonzero reciprocal coordinate required")
    return (alpha+1/alpha+2)/4


def dyadic_product(n, alpha):
    natural(n)
    alpha = exact(alpha)
    need(alpha != 0, "nonzero reciprocal coordinate required")
    return prod(((1-alpha/F(2**d))*(1-1/(F(2**d)*alpha)) /
                 (1-F(1, 2**d))**2 for d in range(1, n+1)), start=F(1))


def reversal_tail_bound(n):
    natural(n)
    need(n >= 5, "tail bound starts at degree five")
    return F(65536, 15*2**n)


def rooted_ray(t):
    natural(t)
    need(t >= 1, "positive ray index required")
    n = (2**(6*t)-1)//3
    return n, (n-3)//6, F(2, (3*t)*(3*t+1))


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        need(ok, message)
        checks += 1

    selectors = {n: selector(n) for n in range(0, 32)}
    reversed_polys = {n: reversed_selector(n) for n in range(0, 32)}
    for n in range(0, 32):
        q, r = selectors[n], reversed_polys[n]
        check(tuple(reversed(q)) == r, "independent coefficient reversal")
        check(evaluate(q, 1) == evaluate(r, 1) == 1, "normalization")
        check(norm(q) == norm(r) < 512, "uniform coefficient norm")
        for k in range(n+1):
            bound = 512*comb(n, k)*prod((shell(d) for d in range(1, n-k+1)), start=F(1))
            check(abs(q[k]) <= bound, "fixed forward coefficient escape bound")
        for d in range(1, n+1):
            check(evaluate(r, 1/shell(d)) == 0, "reciprocal shell zero")

    for n in range(5, 26):
        r = reversed_polys[n]
        for m in range(n+1, 32):
            s = reversed_polys[m]
            diff = tuple(s[k]-(r[k] if k < len(r) else 0) for k in range(len(s)))
            check(norm(diff) <= reversal_tail_bound(n), "exact finite tail bound")

    alphas = (F(-3), F(-1), F(-1, 2), F(1, 8), F(1, 2), F(1), F(2), F(3), F(8))
    for n in range(0, 19):
        for alpha in alphas:
            value = dyadic_product(n, alpha)
            check(value == evaluate(reversed_polys[n], trace_coordinate(alpha)), "factorized trace carrier")
            check(value == dyadic_product(n, 1/alpha), "reciprocal symmetry")
            lhs = (2*alpha-1)*(1-alpha/F(2**n))*dyadic_product(n, 2*alpha)
            rhs = 2*alpha*(1-alpha)*(1-1/(F(2**(n+1))*alpha))*value
            check(lhs == rhs, "finite shift identity including singular factors")
            z = trace_coordinate(alpha)
            check(trace_coordinate(alpha**2) == (2*z-1)**2, "trace-square recursion")

    for d in range(0, 80):
        h = shell(d)
        check(shell(2*d) == h*h/(2-h)**2, "doubled shell map")
    # At alpha=1 the infinite shifted identity only says F(2)=0. It cannot
    # recover F(1); even the identically zero function satisfies that identity.
    check(trace_coordinate(2) == 1/shell(1), "first dyadic zero")
    check(evaluate(reversed_polys[1], 1) == 1 and dyadic_product(1, 2) == 0,
          "normalization must be retained independently")

    ray = []
    for t in range(1, 21):
        n, j, mass = rooted_ray(t)
        check(n == 6*j+3 and n % 2 == 1, "known three-divisible rooted leaf")
        check(step(n) == (1, 6*t), "one-step root control")
        check(word_weight((6*t,)) == mass, "independent exact leaf weight")
        # mass >= 2^-b and h_0(j)^-1 >= 2^(j-2), without making 2^j.
        b = mass.denominator.bit_length()
        check(mass >= F(1, 2**b), "rational mass binary lower bound")
        log_lower = j-2-b
        if t >= 2:
            check(log_lower > t, "inverse-moment ray terms cannot tend to zero")
        if t <= 3:
            check(mass/shell(j) >= F(2)**log_lower, "direct finite inverse-moment lower bound")
        if t <= 5:
            ray.append({"t": t, "source": n, "index": j, "mass": str(mass),
                        "inverse_first_moment_term_log2_lower": log_lower})

    hostiles = (lambda: selector(True), lambda: selector(-1), lambda: reversed_selector(1.0),
                lambda: dyadic_product(2, 0), lambda: dyadic_product(2, 0.5),
                lambda: trace_coordinate(True), lambda: reversal_tail_bound(4),
                lambda: rooted_ray(0))
    for hostile in hostiles:
        try:
            hostile()
        except ValueError:
            check(True, "malformed input rejected")
        else:
            check(False, "malformed input accepted")

    print(json.dumps({"status": "PASS; finite identities, no positivity oracle", "checks": checks,
                      "degrees": "0..31", "shift_degrees": "0..18",
                      "coefficient_norm_bound": "<512",
                      "reversal_l1_tail_bound": "65536/(15*2^N), N>=5",
                      "rooted_ray_controls": ray, "type_hostiles": len(hostiles),
                      "universal_Collatz_coverage": "OPEN"}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
