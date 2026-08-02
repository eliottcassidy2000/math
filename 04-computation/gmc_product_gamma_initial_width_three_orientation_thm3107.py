#!/usr/bin/env python3
"""Exact companion for THM-3107.

The script derives the two cleared divisibility numerators from Laurent
monomials in the five moment-response values, constructs every one-layer
coefficient vector, proves the unique dominant-mode tail inequalities, and
exhausts every histogram below the two exact tail thresholds.

Only Python integers and fractions are used.
"""

from fractions import Fraction
from itertools import product
from math import comb, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ZERO_EXP = (0, 0, 0, 0, 0)


def clean(poly):
    return {e: c for e, c in poly.items() if c}


def add(*polys):
    out = {}
    for poly in polys:
        for e, c in poly.items():
            out[e] = out.get(e, 0) + c
    return clean(out)


def scale(poly, scalar):
    return clean({e: scalar * c for e, c in poly.items()})


def multiply(left, right):
    out = {}
    for e, c in left.items():
        for f, d in right.items():
            g = tuple(x + y for x, y in zip(e, f))
            out[g] = out.get(g, 0) + c * d
    return clean(out)


def monomial(exponent, coefficient=1):
    return {tuple(exponent): coefficient}


def response_moment(indices):
    """Return L(prod f_i) as a Laurent monomial in r_1,...,r_5.

    We have scaled w_1 to one.  Thus w_n=prod_{j=1}^{n-1} r_j and
    L(prod f_i)=w_{sum i}/prod w_i.
    """

    total = sum(indices)
    exponent = []
    for j in range(1, 6):
        exponent.append(int(j < total) - sum(int(j < i) for i in indices))
    return monomial(exponent)


U = ((1, 1), (0, -1))
V = ((2, 1), (1, -1))


def contract(vectors):
    out = {}
    for choices in product(*vectors):
        indices = tuple(item[0] for item in choices)
        coefficient = prod(item[1] for item in choices)
        out = add(out, scale(response_moment(indices), coefficient))
    return out


def cleared_banks():
    g11 = contract((U, U))
    g12 = contract((U, V))
    g22 = contract((V, V))
    t111 = contract((U, U, U))
    t112 = contract((U, U, V))
    t122 = contract((U, V, V))
    t222 = contract((V, V, V))

    i1 = add(
        scale(multiply(multiply(t112, g11), g22), 3),
        scale(multiply(t222, multiply(g11, g11)), -1),
        scale(multiply(multiply(t111, g12), g22), -2),
    )
    i2 = add(
        scale(multiply(multiply(t122, g11), g22), 3),
        scale(multiply(multiply(t222, g12), g11), -2),
        scale(multiply(t111, multiply(g22, g22)), -1),
    )

    # F_j=-r_1^2 I_j.  Every Laurent exponent must now be nonnegative.
    shift = (2, 0, 0, 0, 0)
    banks = []
    for invariant in (i1, i2):
        bank = {
            tuple(a + b for a, b in zip(e, shift)): -c
            for e, c in invariant.items()
        }
        require(all(min(e) >= 0 for e in bank), "clearing left a Laurent pole")
        banks.append(clean(bank))
    return banks


def polynomial_multiply(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return out


def layer_vector(exponent):
    """Coefficients of prod_{j=1}^5 (1+j*x)^exponent_j."""

    out = [1]
    for j, power in enumerate(exponent, start=1):
        for _ in range(power):
            out = polynomial_multiply(out, [1, j])
    out += [0] * (7 - len(out))
    require(len(out) == 7, "one-layer degree exceeded six")
    return tuple(out)


def weak_compositions(total, parts, prefix=()):
    if parts == 1:
        yield prefix + (total,)
        return
    for first in range(total + 1):
        yield from weak_compositions(total - first, parts - 1, prefix + (first,))


EXPECTED_DOMINANTS = (
    ((2, 1, 1, 1, 1), 1, (1, 16, 100, 310, 499, 394, 120)),
    ((1, 2, 1, 1, 1), 2, (1, 17, 115, 395, 724, 668, 240)),
)
EXPECTED_THRESHOLDS = (14, 16)
EXPECTED_FINITE_COUNTS = (27132, 54264)
EXPECTED_MIN_POSITIVE = (16, 12)


def audit_bank(index, bank):
    terms = [(coefficient, exponent, layer_vector(exponent))
             for exponent, coefficient in sorted(bank.items(), reverse=True)]
    require(len(terms) == (20, 23)[index], "generic response-bank size drift")

    dominant = []
    for term in terms:
        coefficient, _, vector = term
        if coefficient <= 0:
            continue
        if all(all(vector[k] >= other[2][k] for k in range(1, 7))
               for other in terms):
            dominant.append(term)
    require(len(dominant) == 1, "dominant response mode is not unique")
    dom_coefficient, dom_exponent, dom_vector = dominant[0]
    expected_exponent, expected_coefficient, expected_vector = EXPECTED_DOMINANTS[index]
    require(dom_exponent == expected_exponent, "dominant exponent drift")
    require(dom_coefficient == expected_coefficient, "dominant coefficient drift")
    require(dom_vector == expected_vector, "dominant response vector drift")

    negative_ratios = []
    for coefficient, exponent, vector in terms:
        if coefficient >= 0:
            continue
        ratios = tuple(Fraction(vector[k], dom_vector[k]) for k in range(1, 7))
        maximum = max(ratios)
        require(ratios[0] == maximum, f"negative mode {exponent} not controlled at k=1")
        require(Fraction(0) <= maximum < 1,
                f"negative mode {exponent} has no contracting dominant ratio")
        negative_ratios.append((-coefficient, maximum))

    def tail_bound(active_layers):
        return sum(Fraction(multiplicity) * ratio ** active_layers
                   for multiplicity, ratio in negative_ratios)

    threshold = EXPECTED_THRESHOLDS[index]
    require(tail_bound(threshold) < dom_coefficient, "tail threshold does not pass")
    require(tail_bound(threshold - 1) >= dom_coefficient,
            "claimed tail threshold is not first for this certificate")

    checked = 0
    zero_count = 0
    minimum_positive = None
    for active_layers in range(threshold):
        for histogram in weak_compositions(active_layers, 6):
            coefficient = sum(
                c * prod(vector[k + 1] ** histogram[k] for k in range(6))
                for c, _, vector in terms
            )
            weighted_degree = sum((k + 1) * histogram[k] for k in range(6))
            require(coefficient >= 0, "negative finite histogram coefficient")
            require((coefficient == 0) == (weighted_degree < 5),
                    "order-five vanishing boundary drift")
            if coefficient == 0:
                zero_count += 1
            else:
                minimum_positive = coefficient if minimum_positive is None else min(
                    minimum_positive, coefficient
                )
            checked += 1

    require(checked == EXPECTED_FINITE_COUNTS[index], "histogram census drift")
    require(zero_count == 12, "finite zero histogram count drift")
    require(minimum_positive == EXPECTED_MIN_POSITIVE[index],
            "minimum positive finite coefficient drift")

    return {
        "terms": len(terms),
        "dominant": dom_vector,
        "threshold": threshold,
        "tail": tail_bound(threshold),
        "checked": checked,
        "zeros": zero_count,
        "minimum": minimum_positive,
    }


def evaluate_bank(bank, responses):
    return sum(
        coefficient * prod(responses[j] ** exponent[j] for j in range(5))
        for exponent, coefficient in bank.items()
    )


def direct_tensors_from_weights(weights):
    def moment(indices):
        return weights[sum(indices)] / prod(weights[i] for i in indices)

    def numeric_contract(vectors):
        return sum(
            prod(item[1] for item in choices) * moment(tuple(item[0] for item in choices))
            for choices in product(*vectors)
        )

    g11 = numeric_contract((U, U))
    g12 = numeric_contract((U, V))
    g22 = numeric_contract((V, V))
    t111 = numeric_contract((U, U, U))
    t112 = numeric_contract((U, U, V))
    t122 = numeric_contract((U, V, V))
    t222 = numeric_contract((V, V, V))
    i1 = 3 * t112 * g11 * g22 - t222 * g11 * g11 - 2 * t111 * g12 * g22
    i2 = 3 * t122 * g11 * g22 - 2 * t222 * g12 * g11 - t111 * g22 * g22
    atomic = 2 * t222 * g12 - 3 * t122 * g22
    return (g11, g12, g22, t111, t112, t122, t222, i1, i2, atomic)


def rising(theta, n):
    out = Fraction(1)
    for j in range(n):
        out *= theta + j
    return out


def boundary_controls(banks):
    # The shape-two atomic route fails although both exact invariants survive.
    theta = Fraction(2)
    weights = [rising(theta, n) for n in range(7)]
    values = direct_tensors_from_weights(weights)
    expected = (
        Fraction(1, 2), Fraction(1, 2), Fraction(5, 6),
        Fraction(1, 2), Fraction(1), Fraction(13, 6), Fraction(16, 3),
        Fraction(-1, 2), Fraction(-11, 36), Fraction(-1, 12),
    )
    require(values == expected, "shape-two atomic hostile drift")

    # Positive coefficients without a negative-real-root factorization do not
    # preserve the orientation.  P(n)=n^3+n^2+6n+1000 has discriminant <0.
    def response_polynomial(n):
        return n ** 3 + n ** 2 + 6 * n + 1000

    responses = tuple(Fraction(response_polynomial(n), response_polynomial(0))
                      for n in range(1, 6))
    hostile = tuple(evaluate_bank(bank, responses) for bank in banks)
    require(hostile == (
        Fraction(-1229376, 3814697265625),
        Fraction(1539056, 3814697265625),
    ), "positive-coefficient response hostile drift")
    discriminant = -26896828
    require(discriminant < 0, "cubic response unexpectedly real-rooted")

    # Independent direct-tensor versus response-bank fixtures guard both the
    # sign convention and the clearing F_j=-r_1^2 I_j.
    shape_fixtures = (
        (Fraction(1, 2),),
        (Fraction(2),),
        (Fraction(7, 3),),
        (Fraction(1, 2), Fraction(4, 3)),
        (Fraction(2), Fraction(3)),
        (Fraction(5, 4), Fraction(7, 2)),
        (Fraction(1, 3), Fraction(5, 4), Fraction(7, 2)),
        (Fraction(2, 3), Fraction(3, 2), Fraction(11, 4)),
        (Fraction(1), Fraction(2), Fraction(5)),
    )
    direct_fixtures = 0
    for shapes in shape_fixtures:
        weights = [prod(rising(theta, n) for theta in shapes) for n in range(7)]
        direct = direct_tensors_from_weights(weights)
        responses = tuple(prod(1 + Fraction(j, 1) / theta for theta in shapes)
                          for j in range(1, 6))
        for index, bank in enumerate(banks):
            invariant = direct[7 + index]
            require(evaluate_bank(bank, responses) == -responses[0] ** 2 * invariant,
                    "direct tensor versus cleared response bank mismatch")
            require(invariant < 0, "rational product-Gamma fixture lost orientation")
        direct_fixtures += 1

    # Exact Gauss multiplication on rational controls.
    gauss_checks = 0
    for theta in (Fraction(1, 3), Fraction(2), Fraction(7, 4)):
        for dilation in range(1, 7):
            shapes = tuple((theta + r) / dilation for r in range(dilation))
            for n in range(7):
                left = rising(theta, dilation * n)
                right = dilation ** (dilation * n) * prod(rising(shape, n) for shape in shapes)
                require(left == right, "Gauss multiplication identity drift")
                gauss_checks += 1
    return hostile, discriminant, gauss_checks, direct_fixtures


def main():
    banks = cleared_banks()
    audits = [audit_bank(index, bank) for index, bank in enumerate(banks)]
    hostile, discriminant, gauss_checks, direct_fixtures = boundary_controls(banks)

    print("THM3107 exact product-Gamma initial width-three orientation")
    print(
        "generic_banks="
        f"I1:{audits[0]['terms']};I2:{audits[1]['terms']};"
        "clearing=-r1^2"
    )
    print(
        "dominant_vectors="
        f"I1:{','.join(map(str, audits[0]['dominant']))};"
        f"I2:{','.join(map(str, audits[1]['dominant']))}"
    )
    print(
        "tail_thresholds="
        f"I1:{audits[0]['threshold']}:{audits[0]['tail']};"
        f"I2:{audits[1]['threshold']}:{audits[1]['tail']}"
    )
    print(
        "finite_histograms="
        f"I1:{audits[0]['checked']}:zeros{audits[0]['zeros']}:min{audits[0]['minimum']};"
        f"I2:{audits[1]['checked']}:zeros{audits[1]['zeros']}:min{audits[1]['minimum']}"
    )
    print("coefficient_boundary=box_0_to_6;zero_below_degree_5;positive_in_box_from_degree_5")
    print("shape2=I1:-1/2;I2:-11/36;atomic_D:-1/12")
    print(
        "nonreal_response_hostile="
        f"P:n3+n2+6n+1000;disc:{discriminant};F1:{hostile[0]};F2:{hostile[1]}"
    )
    print(
        f"direct_bank_fixtures={direct_fixtures};"
        f"gauss_multiplication_checks={gauss_checks};anchored_AP=0,d,2d"
    )
    print("translated_AP_boundary=quadratic_tilt_2a_vs_cubic_tilt_3a")


if __name__ == "__main__":
    main()
