#!/usr/bin/env python3
"""Dependency-free exact referee for candidate THM-2410."""

from fractions import Fraction
from itertools import combinations, product


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cyclotomic_reduce(coefficients, prime):
    """Reduce sum_j coefficients[j] zeta^j modulo Phi_prime."""
    require(len(coefficients) == prime, "wrong coefficient-vector length")
    top = coefficients[-1]
    return tuple(coefficients[j] - top for j in range(prime - 1))


def cyclic_convolution(left, right):
    require(len(left) == P and len(right) == P, "wrong convolution length")
    out = [0] * P
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in enumerate(right):
            if b:
                out[(i + j) % P] += a * b
    return out


def fourier_vector(mask, character, safe=False):
    coefficients = [0] * P
    chosen = set(mask)
    for residue in range(P):
        value = int((residue not in chosen) if safe else (residue in chosen))
        coefficients[(character * residue) % P] += value
    return coefficients


def fourier_nonzero(mask, character, safe=False):
    return any(cyclotomic_reduce(fourier_vector(mask, character, safe), P))


def cyclic_consecutive(mask):
    chosen = set(mask)
    size = len(chosen)
    return any(
        chosen == {(start + offset) % P for offset in range(size)}
        for start in range(P)
    )


def all_masks_prime_cyclotomic_check():
    checked = 0
    residues = tuple(range(P))
    for size in range(1, P):
        for mask in combinations(residues, size):
            for character in range(1, P):
                require(
                    fourier_nonzero(mask, character, safe=False),
                    "proper danger mask lost a nontrivial character",
                )
                require(
                    fourier_nonzero(mask, character, safe=True),
                    "proper safe mask lost a nontrivial character",
                )
            checked += 1
    require(checked == 2**P - 2, "proper-mask census mismatch")
    return checked


def relation_character_census():
    weights = (1, 2, 3, 4, 5, 11)
    inverse_last = pow(weights[-1], -1, P)
    count = 0
    for prefix in product(range(1, P), repeat=5):
        last = (-sum(a * b for a, b in zip(weights[:-1], prefix))) % P
        last = (last * inverse_last) % P
        if last:
            count += 1
    formula = ((P - 1) ** 6 + (P - 1)) // P
    require(count == formula == 229_692, "unit-annihilator count mismatch")
    endpoint = count * (P - 1) ** 3
    eligible = endpoint * (P - 2)
    require(endpoint == 396_907_776, "endpoint census mismatch")
    require(eligible == 4_365_985_536, "endpoint/deep census mismatch")
    return count, endpoint, eligible


def quotient_normalization_check():
    full_states = P**10
    quotient_states = P**9
    require(full_states == P * quotient_states, "gauge orbit census mismatch")
    require(
        Fraction(P, full_states) == Fraction(1, quotient_states),
        "full and quotient normalizations differ",
    )
    return full_states, quotient_states


def hostile_check():
    unit_masks = (
        (9, 10, 11, 12),
        (0, 1),
        (2, 3),
        (3, 4),
        (5, 6),
        (7, 8),
    )
    blocker_masks = ((1, 2), (3, 4), (5, 6))
    masks = unit_masks + blocker_masks
    weights = (1, 2, 3, 4, 5, 11, 0, 0, 0)
    character = (1,) * 9
    trajectory = (1, 1, 1, 1, 1, 1, 0, 0, 0)

    union = set().union(*(set(mask) for mask in unit_masks))
    incidence = sum(len(mask) for mask in unit_masks)
    require(
        all(cyclic_consecutive(mask) for mask in masks),
        "hostile contains a nonconsecutive physical mask",
    )
    require(union == set(range(P)), "unit danger masks do not cover F_13")
    require(incidence == 14, "unit incidence is not fourteen")
    multiplicities = {
        residue: sum(residue in mask for mask in unit_masks) for residue in range(P)
    }
    require(
        [residue for residue, count in multiplicities.items() if count == 2] == [3],
        "hostile does not have its unique double root at three",
    )
    require(
        all(count in (1, 2) for count in multiplicities.values()),
        "hostile unit incidence profile is not physical",
    )
    require(
        sum(a * b for a, b in zip(character, weights)) % P == 0,
        "hostile character is not in the relation annihilator",
    )
    phase = sum(a * b for a, b in zip(character, trajectory)) % P
    require(phase == 6, "hostile trajectory phase mismatch")

    product_vector = [1] + [0] * (P - 1)
    for mask in masks:
        factor = fourier_vector(mask, 1, safe=True)
        require(any(cyclotomic_reduce(factor, P)), "local hostile factor vanished")
        product_vector = cyclic_convolution(product_vector, factor)
    require(
        any(cyclotomic_reduce(product_vector, P)),
        "rank-one hostile product vanished locally",
    )

    anchored = 0
    for t in range(P):
        covered = any(((-t * trajectory[i]) % P) in unit_masks[i] for i in range(6))
        require(covered, "translated hostile is not anchored at zero")
        anchored += int(covered)
    require(anchored == P, "not every translated hostile is anchored")

    rotation_sum = [0] * P
    for t in range(P):
        rotation_sum[(phase * t) % P] += 1
    require(
        not any(cyclotomic_reduce(rotation_sum, P)),
        "translated local coefficients did not cancel",
    )
    return incidence, phase


def septimal_convolution_hostile():
    prime = 7
    f = [1] + [0] * (prime - 1)
    g = [0] + [1] * (prime - 1)
    for character in range(1, prime):
        fhat = [0] * prime
        ghat = [0] * prime
        for residue in range(prime):
            fhat[(character * residue) % prime] += f[residue]
            ghat[(character * residue) % prime] += g[residue]
        require(any(cyclotomic_reduce(fhat, prime)), "delta lost septimal support")
        require(
            any(cyclotomic_reduce(ghat, prime)),
            "complement lost septimal support",
        )
    require(all(a * b == 0 for a, b in zip(f, g)), "septimal hostile product live")
    return prime - 1


def main():
    mask_count = all_masks_prime_cyclotomic_check()
    relation_count, endpoint_count, eligible_count = relation_character_census()
    full_states, quotient_states = quotient_normalization_check()
    incidence, hostile_phase = hostile_check()
    septimal_colours = septimal_convolution_hostile()

    endpoint_floor = Fraction(2, 169) ** 18
    deep_floor = Fraction(2, 169) ** 20

    print("THM-2410 FULL-COORDINATE PROJECTOR -- exact audit")
    print(f"proper nonempty F13 masks checked: {mask_count}")
    print(f"unit annihilator characters: {relation_count}")
    print(f"all-nonzero endpoint characters: {endpoint_count}")
    print(f"eligible endpoint/deep pairs: {eligible_count}")
    print(f"full states / quotient states: {full_states} / {quotient_states}")
    print("normalization identity: 13/13^10 = 1/13^9")
    print(f"hostile unit incidence / trajectory phase: {incidence} / {hostile_phase}")
    print("hostile physical masks: cyclic-consecutive / unique double root 3")
    print("hostile local factors: all nonzero")
    print("hostile anchored slices: 13 / 13")
    print("hostile integrated coefficient: 0 in Q(zeta_13)")
    print(f"endpoint Gram rational floor: {endpoint_floor}")
    print(f"deep Gram rational floor: {deep_floor}")
    print(f"septimal nontrivial colours checked per factor: {septimal_colours}")
    print("septimal full-support product hostile: pointwise zero")
    print("PASS")


if __name__ == "__main__":
    main()
