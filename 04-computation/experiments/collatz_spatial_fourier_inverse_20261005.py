"""Rational controls for the degree-one spatial Fourier inverse.

The infinite spectral theorem is analytic, not inferred from these finite
matrices. No astronomical finite stencil from the coarse plan is expanded.
"""
from fractions import Fraction as F
from itertools import product
import json

C = F(1, 200000)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def natural(n):
    need(type(n) is int and n >= 0, "exact natural required")


def shell(d):
    natural(d)
    return F(4*(1 << d), ((1 << d)+1)**2)


def parameters(radius):
    natural(radius)
    delta = F(8, 1 << radius)
    gamma = C-delta
    need(gamma > 0, "band radius does not pay the spectral tail")
    return delta, gamma, 1-gamma/9


def ceil_fraction(x):
    return -(-x.numerator//x.denominator)


def plan(tolerance):
    need(type(tolerance) in (int, F) and 0 < tolerance <= 1,
         "rational tolerance in (0,1] required")
    eps = F(tolerance)
    radius = 22
    while F(8, 1 << radius) > C*eps/8:
        radius += 1
    delta, gamma, rho = parameters(radius)
    # rho^(K+1) <= 1/(1+(K+1)*(1-rho)); no enormous rational power made.
    depth = max(0, ceil_fraction((4/eps-1)/(1-rho))-1)
    bound = 1/(1+(depth+1)*(1-rho))+delta/gamma
    addresses = 2*radius*depth+1
    # The stencil l2 norm <=1/gamma, l1 norm <=addresses/gamma.
    noise = eps*gamma/(4*addresses)
    return {"radius": radius, "depth": depth, "address_count_bound": addresses,
            "approximation_bound": bound, "per_measurement_error": noise}


def ldlt_pivots(matrix):
    n = len(matrix)
    lower = [[F(0) for _ in range(n)] for _ in range(n)]
    pivots = []
    for i in range(n):
        lower[i][i] = F(1)
        for j in range(i):
            lower[i][j] = (matrix[i][j]-sum((lower[i][k]*lower[j][k]*pivots[k]
                                             for k in range(j)), F(0)))/pivots[j]
        pivot = matrix[i][i]-sum((lower[i][k]**2*pivots[k] for k in range(i)), F(0))
        need(pivot != 0, "zero exact LDL pivot")
        pivots.append(pivot)
    return tuple(pivots)


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        need(ok, message)
        checks += 1

    check(F(22, 7)**2/F(56, 81) < 15, "elementary upper exponent bound")
    check(F(72, 3**15) > C, "coarse Fourier coercivity constant")
    check(2*(F(1, 3)+F(1, 81)) > F(69, 100), "log two lower bound")
    for n in range(1, 13):
        a = [[shell(abs(i-j)) for j in range(n)] for i in range(n)]
        low = [[a[i][j]-(C if i == j else 0) for j in range(n)] for i in range(n)]
        high = [[(9 if i == j else 0)-a[i][j] for j in range(n)] for i in range(n)]
        for pivot in ldlt_pivots(low)+ldlt_pivots(high):
            check(pivot > 0, "independent exact finite spectral sandwich")

    for v in product((-1, 0, 1), repeat=5):
        q = sum((F(v[i]*v[j])*shell(abs(i-j)) for i in range(5) for j in range(5)), F(0))
        mass = sum(x*x for x in v)
        check(C*mass <= q <= 9*mass, "signed quadratic-form control")
    for radius in range(22, 45):
        delta, gamma, rho = parameters(radius)
        check(0 < gamma < C and 0 < rho < 1, "paid finite-band contraction")
        actual_partial_tail = 2*sum((shell(d) for d in range(radius+1, radius+41)), F(0))
        check(actual_partial_tail < delta, "independent rational tail partial control")
        for k in range(10):
            check(rho**(k+1) <= 1/(1+(k+1)*(1-rho)), "Bernoulli error bound")

    plans = []
    for exponent in (0, 1, 4, 8, 16, 32):
        eps = F(1, 2**exponent)
        p = plan(eps)
        delta, gamma, rho = parameters(p['radius'])
        noise_bill = p['per_measurement_error']*p['address_count_bound']/gamma
        check(p['approximation_bound'] < eps/2, "finite parameter choice without large expansion")
        check(noise_bill == eps/4, "safe address-dependent measurement bill")
        check(p['approximation_bound']+noise_bill < eps, "total requested tolerance")
        if exponent <= 8:
            plans.append({k: str(v) if type(v) is F else v for k, v in p.items()} | {"tolerance": str(eps)})

    # The finite observation map is injective in the displayed examples, yet
    # positivity of its entries cannot imply positivity of every input atom.
    check(all(shell(m) > 0 for m in range(101)), "delta-zero has everywhere positive first-moment field")
    for bad in (lambda: shell(True), lambda: shell(-1), lambda: parameters(20),
                lambda: parameters(22.0), lambda: plan(0), lambda: plan(0.5),
                lambda: plan(True), lambda: plan(F(2))):
        try:
            bad()
        except ValueError:
            check(True, "malformed input rejected")
        else:
            check(False, "malformed input accepted")
    print(json.dumps({"status": "PASS; finite spectral and parameter controls, no moment oracle",
                      "checks": checks, "coercivity_lower": str(C), "operator_upper": 9,
                      "finite_LDL_sizes": "1..12", "signed_vectors": 243,
                      "coarse_plans_not_expanded": plans,
                      "scope": "First-moment field identifies the law; source positivity remains OPEN"},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
