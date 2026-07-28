#!/usr/bin/env python3
"""Exact referee for THM-2644.

All logical checks use ``require`` so optimized Python executes the same
certificate as ordinary Python.
"""

from itertools import combinations, product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def convolution(left, right):
    n = len(left)
    require(len(right) == n, "convolution size mismatch")
    return tuple(
        sum(left[a] * right[(x - a) % n] for a in range(n))
        for x in range(n)
    )


def convolution_power(weight, exponent):
    require(exponent >= 1, "positive convolution exponent required")
    n = len(weight)
    out = tuple(1 if a == 0 else 0 for a in range(n))
    base = tuple(weight)
    k = exponent
    while k:
        if k & 1:
            out = convolution(out, base)
        base = convolution(base, base)
        k //= 2
    return out


def statistics(weight):
    n = len(weight)
    mass = sum(weight)
    energy = sum(value * value for value in weight)
    defect = mass * mass - energy
    same_return = sum(weight[a] * weight[(-a) % n] for a in range(n))
    involution_energy = sum(
        weight[a] * weight[a] for a in range(n) if (2 * a) % n == 0
    )
    reverse = tuple(weight[(-a) % n] for a in range(n))
    require(convolution(reverse, weight)[0] == energy,
            "reverse/Gram identity coefficient mismatch")
    require(convolution(weight, weight)[0] == same_return,
            "same-orientation identity coefficient mismatch")
    return mass, energy, defect, same_return, involution_energy


def equality_shape(weight):
    n = len(weight)
    support = [a for a, value in enumerate(weight) if value > 0]
    return (
        len(support) <= 1
        or (len(support) == 2 and (support[0] + support[1]) % n == 0)
    )


def audit_nonnegative(weight):
    require(all(value >= 0 for value in weight), "nonnegative audit received a signed weight")
    n = len(weight)
    mass, energy, defect, same_return, involution_energy = statistics(weight)
    require(defect >= 0, "purity defect became negative")
    require(involution_energy >= same_return - defect,
            "involution-energy inequality failure")
    require((involution_energy == same_return - defect) == equality_shape(weight),
            "equality classification failure")
    support_size = sum(value > 0 for value in weight)
    require((defect == 0) == (support_size <= 1), "purity classification failure")
    if same_return > defect:
        require(any(weight[a] > 0 for a in range(n) if (2 * a) % n == 0),
                "strict gate missed the involution locus")
        if n % 2 == 1:
            require(weight[0] > 0, "odd cyclic strict gate missed the identity")
    if mass > 0 and energy == mass * mass and same_return > 0 and n % 2 == 1:
        require(support_size == 1 and weight[0] == mass,
                "pure odd return did not decode the identity")


def sparse_positive_weights(n, max_support, max_height):
    yield (0,) * n
    for size in range(1, max_support + 1):
        for support in combinations(range(n), size):
            for values in product(range(1, max_height + 1), repeat=size):
                weight = [0] * n
                for index, value in zip(support, values):
                    weight[index] = value
                yield tuple(weight)


def main():
    small_profiles = 0
    for n in (3, 5, 7):
        for weight in product(range(3), repeat=n):
            audit_nonnegative(weight)
            small_profiles += 1

    boolean_profiles = 0
    boolean_strict = 0
    for weight in product((0, 1), repeat=13):
        audit_nonnegative(weight)
        _, _, defect, same_return, _ = statistics(weight)
        boolean_strict += same_return > defect
        boolean_profiles += 1
    require(boolean_profiles == 8192 and boolean_strict == 1,
            "C13 Boolean strict-gate census")

    sparse_profiles = 0
    sparse_strict = 0
    for weight in sparse_positive_weights(13, 3, 3):
        audit_nonnegative(weight)
        _, _, defect, same_return, _ = statistics(weight)
        sparse_strict += same_return > defect
        sparse_profiles += 1
    require(sparse_profiles == 8464, "C13 sparse weighted profile count")

    identity = tuple(1 if a == 0 else 0 for a in range(13))
    translated = tuple(1 if a == 5 else 0 for a in range(13))
    inverse_pair = tuple(1 if a in (5, 8) else 0 for a in range(13))
    robust = tuple(3 if a == 0 else 1 if a == 1 else 0 for a in range(13))
    uniform = (1,) * 13
    require(statistics(identity) == (1, 1, 0, 1, 1), "identity control")
    require(statistics(translated) == (1, 1, 0, 0, 0), "translated singleton hostile")
    require(statistics(inverse_pair) == (2, 2, 2, 2, 0), "inverse-pair boundary")
    require(statistics(robust) == (4, 10, 6, 9, 9), "robust non-pure control")
    require(statistics(uniform) == (13, 13, 156, 13, 1), "dense control")

    signed = tuple(1 if a == 1 else -1 if a == 2 else 0 for a in range(13))
    signed_stats = statistics(signed)
    require(signed_stats[:4] == (0, 2, -2, 0) and signed[0] == 0,
            "signed nonnegativity hostile")
    even = (0, 1)
    require(statistics(even) == (1, 1, 0, 1, 1) and even[0] == 0,
            "even involution hostile")

    kfold_checks = 0
    coprime_positive = 0
    divisible_hostiles = 0
    for a in range(13):
        singleton = tuple(1 if x == a else 0 for x in range(13))
        for k in range(2, 27):
            observed = convolution_power(singleton, k)[0]
            expected = int((k * a) % 13 == 0)
            require(observed == expected, "cyclic k-return formula failure")
            if gcd(k, 13) == 1 and observed:
                require(a == 0, "coprime return failed to decode identity")
                coprime_positive += 1
            if k % 13 == 0 and a != 0:
                require(observed == 1, "divisible-length hostile failed")
                divisible_hostiles += 1
            kfold_checks += 1
    require(kfold_checks == 325 and divisible_hostiles == 24,
            "k-fold control census")

    print("THM-2644 ODD-TORSOR PURITY/RETURN EXACT REFEREE")
    print(f"small_odd_weight_profiles={small_profiles}")
    print(f"C13_boolean_profiles={boolean_profiles} strict_gate_profiles={boolean_strict}")
    print(f"C13_sparse_height3_profiles={sparse_profiles} strict_gate_profiles={sparse_strict}")
    print("controls identity=(1,1,0,1) translated=(1,1,0,0) inverse_pair=(2,2,2,2)")
    print("robust_nonpure=(M4,E10,delta6,R9) dense=(M13,E13,delta156,R13)")
    print("hostiles signed=(M0,E2,delta-2,R0) even_C2=(M1,E1,delta0,R1)")
    print(f"kfold_C13_checks={kfold_checks} coprime_positive={coprime_positive} p_divisible_nonzero_hostiles={divisible_hostiles}")
    print("reverse_gram_and_same_orientation_compositions=EXACT")
    print("PASS")


if __name__ == "__main__":
    main()
