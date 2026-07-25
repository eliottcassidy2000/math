#!/usr/bin/env python3
"""Exact referee for the conditional same-core exclusion in THM-2232.

The two high-first profiles left by THM-2229 are (4,6,8) and (5,7,9).
This companion assumes that their normalized blocker cores coincide:

    c_j = 13^(a+2j) u,  j=0,1,2,  a in {4,5},  13 does not divide u.

It enumerates the exact stationary two-state danger process

    X_t(x) = 1_{D_u}(13^t x)

and evaluates an explicit signed-transfer dual certificate.  It also checks
the pointwise positive-part majorant at every endpoint of the deliberately
enlarged residual box, checks the transfer centering exactly, and records
hostile controls showing that the unsigned three-checkpoint carrier itself
survives above the target.
"""

from fractions import Fraction
from itertools import product


P = 13
TARGET = Fraction(961, 6930)
STATIONARY_NUMERATORS = (6, 1)
STATIONARY = tuple(Fraction(value, 7) for value in STATIONARY_NUMERATORS)
TRANSITION_NUMERATORS = ((11, 2), (12, 1))
TRANSITION = tuple(
    tuple(Fraction(value, P) for value in row)
    for row in TRANSITION_NUMERATORS
)
WEIGHTS = (1, 2, 2, 1)
DUAL_TIMES = {
    4: (1, 3, 5, 7),
    5: (2, 4, 6, 8),
}
EXPECTED_CAPACITY = {
    4: Fraction(2303649491556761, 51185893014090757),
    5: Fraction(2711005672365842843, 60552911435669365531),
}
EXPECTED_CHECKPOINT_MASS = Fraction(916159, 4826809)
EXPECTED_CORRECTION = {
    4: Fraction(-4884270, 62748517),
    5: Fraction(4884270, 815730721),
}
EXPECTED_DUAL_RANGE = {
    4: (Fraction(-255878338, 62748517), Fraction(2)),
    5: (Fraction(-4), Fraction(1636345712, 815730721)),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def word_probability(word):
    """Stationary Markov probability of one finite bit word."""
    value = STATIONARY[word[0]]
    for left, right in zip(word, word[1:]):
        value *= TRANSITION[left][right]
    return value


def word_integer_weight(word):
    """Numerator of the same probability over 7*13^(len(word)-1)."""
    value = STATIONARY_NUMERATORS[word[0]]
    for left, right in zip(word, word[1:]):
        value *= TRANSITION_NUMERATORS[left][right]
    return value


def checkpoint_support(word, first_depth, parity):
    """Intersection of all divided-blocker windows of the given parity."""
    checkpoints = range(parity, first_depth + 1, 2)
    return all(
        any(word[first_depth - checkpoint + 2 * j] for j in range(3))
        for checkpoint in checkpoints
    )


def centered_dual(word, first_depth):
    """The explicit q_a = 2-sum w_t(X_t-(-1/13)^t X_0)."""
    value = Fraction(2)
    for weight, time in zip(WEIGHTS, DUAL_TIMES[first_depth]):
        transfer_eigenvalue = Fraction((-1) ** time, P**time)
        value -= weight * (
            word[time] - transfer_eigenvalue * word[0]
        )
    return value


def positive_part(value):
    return max(Fraction(0), value)


def pointwise_majorant(positive_support, negative_support, dual_value):
    """Dual upper envelope for 0<=R_+<=1_P and 0<=R_-<=5*1_N."""
    return (
        int(positive_support) * positive_part(1 - dual_value)
        + 5 * int(negative_support) * positive_part(dual_value)
    )


def audit_markov_kernel():
    require(sum(STATIONARY) == 1, "stationary law does not sum to one")
    require(
        all(sum(row) == 1 for row in TRANSITION),
        "transition row does not sum to one",
    )
    for state in (0, 1):
        pushed = sum(
            STATIONARY[previous] * TRANSITION[previous][state]
            for previous in (0, 1)
        )
        require(pushed == STATIONARY[state], "stationarity drift")
    require(
        STATIONARY[0] * TRANSITION[0][1]
        == STATIONARY[1] * TRANSITION[1][0],
        "detailed balance drift",
    )
    for terminal_bit in (0, 1):
        backward_one = (
            STATIONARY[1]
            * TRANSITION[1][terminal_bit]
            / STATIONARY[terminal_bit]
        )
        require(
            backward_one == Fraction(2 - terminal_bit, P),
            "root-count conditional law drift",
        )


def audit_case(first_depth):
    last_time = first_depth + 4
    words = tuple(product((0, 1), repeat=last_time + 1))
    direct_total_probability = Fraction(0)
    integer_total_weight = 0
    direct_capacity = Fraction(0)
    integer_capacity_numerator = Fraction(0)
    positive_support_mass = Fraction(0)
    negative_support_mass = Fraction(0)
    majorant_checks = 0
    naive_failure = None
    dual_minimum = None
    dual_maximum = None

    correction = sum(
        weight * Fraction((-1) ** time, P**time)
        for weight, time in zip(WEIGHTS, DUAL_TIMES[first_depth])
    )
    require(
        correction == EXPECTED_CORRECTION[first_depth],
        "centered X_0 correction drift",
    )

    # If the centering is omitted, integral R*q_raw has the generally
    # nonzero coefficient -correction multiplying integral R*X_0.
    raw_correlation_defect = -correction
    centered_correlation_defect = raw_correlation_defect + correction
    require(
        raw_correlation_defect != 0,
        "hostile uncentered control unexpectedly cancels",
    )
    require(
        centered_correlation_defect == 0,
        "signed-transfer centering does not cancel",
    )

    for word in words:
        probability = word_probability(word)
        integer_weight = word_integer_weight(word)
        direct_total_probability += probability
        integer_total_weight += integer_weight

        positive_support = checkpoint_support(word, first_depth, 0)
        negative_support = checkpoint_support(word, first_depth, 1)
        dual_value = centered_dual(word, first_depth)
        dual_minimum = (
            dual_value
            if dual_minimum is None
            else min(dual_minimum, dual_value)
        )
        dual_maximum = (
            dual_value
            if dual_maximum is None
            else max(dual_maximum, dual_value)
        )
        envelope = pointwise_majorant(
            positive_support,
            negative_support,
            dual_value,
        )
        direct_capacity += probability * envelope
        integer_capacity_numerator += integer_weight * envelope
        if positive_support:
            positive_support_mass += probability
        if negative_support:
            negative_support_mass += probability

        # Hostile pointwise audit.  The actual residual has disjoint positive
        # and negative parts with integral values.  We enlarge it to every
        # endpoint of the product box, including simultaneous positive and
        # negative masses, so passing this check is strictly stronger than
        # needed.
        positive_values = (0, 1) if positive_support else (0,)
        negative_values = range(6) if negative_support else (0,)
        for residual_positive in positive_values:
            for residual_negative in negative_values:
                left = (
                    residual_positive * (1 - dual_value)
                    + residual_negative * dual_value
                )
                require(
                    left <= envelope,
                    "pointwise positive-part majorant failed",
                )
                majorant_checks += 1

                # Omitting positive parts is unsafe.  Freeze the first exact
                # witness rather than merely asserting that one exists.
                naive = (
                    int(positive_support) * (1 - dual_value)
                    + 5 * int(negative_support) * dual_value
                )
                if naive_failure is None and left > naive:
                    naive_failure = (
                        "".join(map(str, word)),
                        dual_value,
                        residual_positive,
                        residual_negative,
                        left,
                        naive,
                    )

    common_denominator = 7 * P**last_time
    require(direct_total_probability == 1, "word probabilities do not sum to one")
    require(
        integer_total_weight == common_denominator,
        "independent integer Markov normalization drift",
    )
    integer_capacity = integer_capacity_numerator / common_denominator
    require(
        direct_capacity == integer_capacity,
        "direct/integer-weight capacity mismatch",
    )
    require(
        direct_capacity == EXPECTED_CAPACITY[first_depth],
        "frozen signed-dual capacity drift",
    )
    require(
        direct_capacity < TARGET,
        "signed-dual capacity does not beat target",
    )
    require(
        positive_support_mass == EXPECTED_CHECKPOINT_MASS,
        "geometric checkpoint mass drift",
    )
    require(
        positive_support_mass > TARGET,
        "hostile unsigned checkpoint control no longer survives",
    )
    require(naive_failure is not None, "missing hostile truncation witness")
    require(
        (dual_minimum, dual_maximum) == EXPECTED_DUAL_RANGE[first_depth],
        "frozen dual range drift",
    )

    positive_lipschitz = positive_part(1 - dual_minimum)
    negative_lipschitz = positive_part(dual_maximum)
    positive_checkpoint_count = len(range(0, first_depth + 1, 2))
    negative_checkpoint_count = len(range(1, first_depth + 1, 2))
    atom_stability_coefficient = (
        positive_checkpoint_count * positive_lipschitz
        + 5 * negative_checkpoint_count * negative_lipschitz
    )
    atom_stability_threshold = (
        TARGET - direct_capacity
    ) / atom_stability_coefficient
    require(
        atom_stability_coefficient > 0 and atom_stability_threshold > 0,
        "stability constants are not positive",
    )

    return {
        "capacity": direct_capacity,
        "positive_support_mass": positive_support_mass,
        "negative_support_mass": negative_support_mass,
        "correction": correction,
        "raw_correlation_defect": raw_correlation_defect,
        "majorant_checks": majorant_checks,
        "naive_failure": naive_failure,
        "dual_range": (dual_minimum, dual_maximum),
        "positive_lipschitz": positive_lipschitz,
        "negative_lipschitz": negative_lipschitz,
        "atom_stability_coefficient": atom_stability_coefficient,
        "atom_stability_threshold": atom_stability_threshold,
    }


def format_naive_failure(witness):
    word, dual_value, residual_positive, residual_negative, left, naive = witness
    return (
        f"word={word} q={dual_value} "
        f"Rplus={residual_positive} Rminus={residual_negative} "
        f"left={left} naive={naive}"
    )


def main():
    audit_markov_kernel()
    print("stationary=(6/7,1/7)")
    print("transition=((11/13,2/13),(12/13,1/13))")
    print("reverse_root_law=P(previous=1|terminal=b)=(2-b)/13")
    print(f"target={TARGET}")

    for first_depth in (4, 5):
        result = audit_case(first_depth)
        print(
            f"a={first_depth} profile="
            f"({first_depth},{first_depth + 2},{first_depth + 4})"
        )
        print(
            f"a={first_depth} centered_X0_correction="
            f"{result['correction']}"
        )
        print(
            f"a={first_depth} raw_correlation_defect_coefficient="
            f"{result['raw_correlation_defect']} "
            "centered_defect_coefficient=0"
        )
        print(
            f"a={first_depth} positive_checkpoint_mass="
            f"{result['positive_support_mass']} "
            f"unsigned_margin={result['positive_support_mass'] - TARGET}"
        )
        print(
            f"a={first_depth} negative_checkpoint_mass="
            f"{result['negative_support_mass']}"
        )
        print(
            f"a={first_depth} signed_dual_capacity={result['capacity']} "
            f"target_margin={TARGET - result['capacity']}"
        )
        print(
            f"a={first_depth} dual_range="
            f"({result['dual_range'][0]},{result['dual_range'][1]}) "
            f"positive_lipschitz={result['positive_lipschitz']} "
            f"negative_lipschitz={result['negative_lipschitz']}"
        )
        print(
            f"a={first_depth} atom_stability_coefficient="
            f"{result['atom_stability_coefficient']} "
            f"required_total_core_distance="
            f"{result['atom_stability_threshold']}"
        )
        print(
            f"a={first_depth} hostile_majorant_endpoint_checks="
            f"{result['majorant_checks']}"
        )
        print(
            f"a={first_depth} naive_untruncated_majorant_failure="
            f"{format_naive_failure(result['naive_failure'])}"
        )

    print("direct_vs_integer_weight_markov_evaluation=PASS")
    print("pointwise_positive_part_majorant=PASS")
    print("same_core_profiles_excluded=2")
    print("unrestricted_profiles_excluded=0")
    print("full_profile_closure_requires_inverse_or_stability_step=YES")


if __name__ == "__main__":
    main()
