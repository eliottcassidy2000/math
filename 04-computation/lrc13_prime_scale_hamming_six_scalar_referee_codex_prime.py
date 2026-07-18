#!/usr/bin/env python3
"""Exact referee for the prime common-scale H6 scalar obstruction.

For effective order p, scalar mask cardinality is a residue-class count in
(-p,p].  This script derives the complete mod-13 bonus table, proves the
large-q inequality, and scans the only exceptional six-support banks p=23,29.
"""

from hashlib import sha256
from itertools import combinations


P = 13
LABELS = tuple(range(1, P))
SUPPORTS = tuple(combinations(LABELS, 6))
EXPECTED_B6 = (1, 3, 5, 6, 6, 6, 7, 9, 11, 12, 12, 12)
EXPECTED_B5 = (1, 3, 5, 5, 5, 5, 6, 8, 10, 10, 10, 10)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_prime(value: int) -> bool:
    return value >= 2 and all(
        value % divisor for divisor in range(2, int(value**0.5) + 1)
    )


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


def scalar_cardinality(scale: int, provider_owner_ratio: int) -> int:
    """Count the fixed nonzero residue class in the half-open band."""
    target = scale * provider_owner_ratio % P
    return sum(
        value % P == target for value in range(-scale + 1, scale + 1)
    )


def bonus_row(residue: int) -> tuple[int, ...]:
    # At scale 13q+s every residue class gains exactly two representatives
    # when q increases by one, so q=0 is the complete bonus row.
    return tuple(scalar_cardinality(residue, ratio_) for ratio_ in LABELS)


BONUS = tuple(bonus_row(residue) for residue in LABELS)
B6 = tuple(sum(sorted(row, reverse=True)[:6]) for row in BONUS)
B5 = tuple(sum(sorted(row, reverse=True)[:5]) for row in BONUS)


def capacity(scale: int, support: tuple[int, ...], owner: int) -> int:
    return sum(
        scalar_cardinality(scale, ratio(provider, owner))
        for provider in support
    )


def scalar_supports(scale: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        support
        for support in SUPPORTS
        if all(capacity(scale, support, owner) >= scale for owner in support)
    )


def tournament_fingerprint(scale: int) -> tuple:
    """Orient pair-capacity comparisons; label order is the tie path."""
    out = [0] * 12
    ties = 0
    flips = 0
    for left, right in combinations(LABELS, 2):
        left_observable = scalar_cardinality(scale, ratio(right, left))
        right_observable = scalar_cardinality(scale, ratio(left, right))
        winner = left
        if left_observable == right_observable:
            ties += 1
        elif right_observable > left_observable:
            winner = right
            flips += 1
        loser = left + right - winner
        out[winner - 1] |= 1 << (loser - 1)

    scores = tuple(sorted(mask.bit_count() for mask in out))
    triangles = sum(
        bool(
            ((out[i] >> j) & 1 and (out[j] >> k) & 1 and (out[k] >> i) & 1)
            or ((out[i] >> k) & 1 and (out[k] >> j) & 1 and (out[j] >> i) & 1)
        )
        for i, j, k in combinations(range(12), 3)
    )

    reach = [
        [i == j or bool((out[i] >> j) & 1) for j in range(12)]
        for i in range(12)
    ]
    for middle in range(12):
        for source in range(12):
            for target in range(12):
                reach[source][target] |= (
                    reach[source][middle] and reach[middle][target]
                )
    assigned = set()
    scc_sizes = []
    for vertex in range(12):
        if vertex in assigned:
            continue
        component = {
            other
            for other in range(12)
            if reach[vertex][other] and reach[other][vertex]
        }
        assigned |= component
        scc_sizes.append(len(component))

    paths = [[0] * 12 for _mask in range(1 << 12)]
    for vertex in range(12):
        paths[1 << vertex][vertex] = 1
    for mask in range(1, 1 << 12):
        if mask & (mask - 1) == 0:
            continue
        for last in range(12):
            if not (mask >> last) & 1:
                continue
            previous_mask = mask ^ (1 << last)
            paths[mask][last] = sum(
                paths[previous_mask][previous]
                for previous in range(12)
                if (previous_mask >> previous) & 1
                and (out[previous] >> last) & 1
            )
    return (
        scores,
        triangles,
        tuple(sorted(scc_sizes, reverse=True)),
        sum(paths[-1]),
        ties,
        flips,
    )


def main() -> None:
    require(len(SUPPORTS) == 924, "support census mismatch")
    require(B6 == EXPECTED_B6, "six-largest bonus table mismatch")
    require(B5 == EXPECTED_B5, "five-largest bonus table mismatch")
    require(max(B6[index] - (index + 1) for index in range(12)) == 2,
            "large-q bonus bound mismatch")
    require(max(B5[index] - (index + 1) for index in range(12)) == 2,
            "order-one five-provider bound mismatch")

    # Exact recurrence behind the structural formula.  The finite range covers
    # THM-860's p<=840 bank; the proof is the translation by the two new
    # endpoints in the same residue class when scale increases by thirteen.
    for q in range(65):
        for residue in LABELS:
            scale = 13 * q + residue
            for ratio_ in LABELS:
                require(
                    scalar_cardinality(scale, ratio_)
                    == 2 * q + BONUS[residue - 1][ratio_ - 1],
                    "residue-cardinality recurrence mismatch",
                )

    # A D1 provider contributes zero at any distinct owner.  Five Dp
    # providers are uniformly short for every q>=1 and residue s.
    for q in range(1, 65):
        for residue in LABELS:
            require(
                10 * q + B5[residue - 1] < 13 * q + residue,
                "order-one exclusion inequality failed",
            )

    # Six Dp providers are uniformly short once q>=3.
    for q in range(3, 65):
        for residue in LABELS:
            require(
                12 * q + B6[residue - 1] < 13 * q + residue,
                "large-prime six-provider inequality failed",
            )

    require(
        scalar_supports(17)
        == ((1, 3, 4, 9, 10, 12), (2, 5, 6, 7, 8, 11)),
        "scale-seventeen anchor mismatch",
    )
    require(scalar_supports(23) == (), "p=23 exceptional support survived")
    require(scalar_supports(29) == (), "p=29 exceptional support survived")

    small_prime_margins = {}
    for prime in (19, 23, 29, 31, 37):
        q, residue = divmod(prime, 13)
        small_prime_margins[prime] = prime - (12 * q + B6[residue - 1])
    require(
        small_prime_margins == {19: 1, 23: -1, 29: 0, 31: 1, 37: 1},
        "small-prime margin table mismatch",
    )

    primes_to_840 = tuple(value for value in range(19, 841) if is_prime(value))
    for prime in primes_to_840:
        q, residue = divmod(prime, 13)
        if prime in (23, 29):
            require(scalar_supports(prime) == (), "exceptional prime survived")
        else:
            require(
                12 * q + B6[residue - 1] < prime,
                "prime not excluded by the six-largest bound",
            )

    fingerprint23 = tournament_fingerprint(23)
    fingerprint29 = tournament_fingerprint(29)
    expected_fingerprint = (
        (2, 3, 3, 3, 5, 5, 6, 6, 8, 8, 8, 9),
        40,
        (12,),
        124_961,
        42,
        12,
    )
    require(fingerprint23 == expected_fingerprint,
            "p=23 tournament fingerprint mismatch")
    require(fingerprint29 == expected_fingerprint,
            "p=29 tournament fingerprint mismatch")

    table_digest = sha256(
        bytes(value for row in BONUS for value in row)
        + bytes(B6)
        + bytes(B5)
    ).hexdigest()
    exceptional_digest = sha256(
        repr((scalar_supports(17), scalar_supports(23), scalar_supports(29))).encode()
    ).hexdigest()
    require(
        table_digest
        == "5fa2bb54d6ebdc665e62fdd88054ecc976ec5811bdbc8f366f87b7b41e93f6ae",
        "bonus-table digest mismatch",
    )
    require(
        exceptional_digest
        == "1f089b4656bf5896974e11beeaa0cfa37ac0e321007f61947f63da1f49b87f47",
        "exceptional-support digest mismatch",
    )

    print("prime common-scale Hamming-six scalar referee")
    print("order-p cardinality #{x:-p<x<=p, x=p*r mod13}; recurrence a_(p+13)=a_p+2")
    print("residues s 1,2,3,4,5,6,7,8,9,10,11,12")
    print("six-largest bonuses " + ",".join(map(str, B6)))
    print("five-largest bonuses " + ",".join(map(str, B5)))
    print("D1 exclusion max(B5-s)=2<3q; all-Dp exclusion max(B6-s)=2<q for q>=3")
    print("small prime six-capacity margins p-top6 "
          + " ".join(f"{prime}:{margin}" for prime, margin in small_prime_margins.items()))
    print("exceptional support scans p23 0/924; p29 0/924; p17 anchor QR,NQR")
    print(f"primes 19..839 excluded {len(primes_to_840)}; prime scales with scalar support 0")
    print("ratio tournament p23=p29 scores 2,3,3,3,5,5,6,6,8,8,8,9; "
          "triangles 40; SCC 12; Hamiltonian paths 124961; ties 42; flips 12")
    print("faithful carrier is the labelled induced ratio-capacity digraph with absolute owner sums; tournament orientation is lossy")
    print(f"bonus-table SHA256 {table_digest}")
    print(f"exceptional-support SHA256 {exceptional_digest}")


if __name__ == "__main__":
    main()
