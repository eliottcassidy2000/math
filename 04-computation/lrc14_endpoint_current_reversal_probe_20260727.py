#!/usr/bin/env python3
"""Exact pointed-origin and reversal audit for THM-2625.

THM-2625 proves that its allocated endpoint current

    J(L,R) = L*(L) R*(R)

is nonzero on all 13^4 endpoint pairs.  This companion reconstructs the same
two exact finite-field specializations and asks two finer questions.

1. Does choosing an origin on a nondegenerate determinant fibre ever lose the
   current?  (It cannot: every one of the 13 possible pointed edges survives.)
2. Does row reversal (L,R)->(R,L) preserve the current, even up to a scalar
   depending only on the fibre (q,Delta)?

The second question is tested without division: two reversal ratios differ iff

    J(R1,L1) J(L2,R2) - J(R2,L2) J(L1,R1)

is nonzero.  A nonzero specialization proves the corresponding cyclotomic
integer is nonzero in characteristic zero.  Thus this is a typed obstruction,
not a comparison of unrelated marginal statistics.
"""

from collections import Counter
from hashlib import sha256

import lrc14_canonical_endpoint_current_thm2625 as current


P = 13
EXPECTED_RATIO_DIGESTS = (
    "3aaeebf1ad57458e3e345626fea34629f519b4bbee9a94e89c98d5a8fdb6c372",
    "d1953a4518963fbd0d7f07e9ad91ac9e499bddd54615a5eabec206de8f6c510f",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def endpoint_banks():
    """Reconstruct THM-2625's left/right endpoint banks in both fields."""
    e_zero = current.build_set(current.PAT_E, current.ZERO_ELL)
    q_intervals = current.build_set(current.PAT_QA, current.ZERO_ELL)
    q_starts = [start for start, _ in q_intervals]
    tabs = current.make_tabs(q_intervals, current.X, current.MODS)
    z13 = [pow(root, current.NN // P, prime)
           for prime, root in current.MODS]
    p_banks = [dict(), dict()]
    q_banks = [dict(), dict()]
    for address, ell in current.REPS.items():
        e_intervals = current.build_set(current.PAT_E, ell)
        a_values, _ = current.x_sweep(
            e_intervals, q_intervals, q_starts,
            current.X, current.MODS, tabs,
        )
        b_values = current.endpoint_sum(e_intervals, -current.Y, current.MODS)
        for field_index, (prime, _) in enumerate(current.MODS):
            phase = pow(
                z13[field_index],
                current.M_DEEP * ell[current.TB] % P,
                prime,
            )
            p_banks[field_index][address] = phase * a_values[field_index] % prime
            q_banks[field_index][address] = b_values[field_index] % prime

    points = [(x, y) for x in range(P) for y in range(P)]
    banks = []
    for field_index, (prime, _) in enumerate(current.MODS):
        powers = [pow(z13[field_index], exponent, prime)
                  for exponent in range(P)]
        left = {}
        right = {}
        for point in points:
            lx = 0
            rx = 0
            for address in points:
                tau0, tau1 = current.TAU[address]
                pairing = (tau0 * point[0] + tau1 * point[1]) % P
                lx = (lx + p_banks[field_index][address]
                      * powers[-pairing % P]) % prime
                rx = (rx + q_banks[field_index][address]
                      * powers[pairing]) % prime
            left[point] = lx
            right[point] = rx
        require(all(left.values()) and all(right.values()),
                "THM-2625 endpoint support changed")
        banks.append((prime, left, right))
    return banks


def det(left, right):
    return (left[0] * right[1] - left[1] * right[0]) % P


def add(point, vector, scalar=1):
    return ((point[0] + scalar * vector[0]) % P,
            (point[1] + scalar * vector[1]) % P)


def canonical_origin(q, delta):
    """A coordinate-pivot origin R0 with det(q,R0)=Delta.

    This chart is deliberately only a reproducible enumeration device.  The
    theorem-level statement is independent of it because all 13 cells survive.
    """
    q0, q1 = q
    if q0:
        return (0, delta * pow(q0, -1, P) % P)
    require(q1, "nonzero q required")
    return (-delta * pow(q1, -1, P) % P, 0)


def field_rank(matrix, prime):
    matrix = [list(row) for row in matrix]
    rank = 0
    columns = len(matrix[0]) if matrix else 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(matrix))
                      if matrix[row][column] % prime), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [value * inverse % prime for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column] % prime:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (matrix[row][j] - factor * matrix[rank][j]) % prime
                for j in range(columns)
            ]
        rank += 1
    return rank


def cross_ratio(a, b, c, d, prime):
    denominator = (a - d) * (b - c) % prime
    require(denominator, "cross-ratio denominator vanished")
    return (a - c) * (b - d) * pow(denominator, -1, prime) % prime


def main():
    banks = endpoint_banks()
    points = [(x, y) for x in range(P) for y in range(P)]
    sectors = [
        ((q0, q1), delta)
        for q0, q1 in points if (q0, q1) != (0, 0)
        for delta in range(1, P)
    ]
    require(len(sectors) == 2016, "nondegenerate sector universe changed")

    # The h-chart merely permutes the 13 cells of each fibre.  This makes the
    # origin-free corollary of THM-2625 explicit, including h=3.
    pointed = []
    for q, delta in sectors:
        r0 = canonical_origin(q, delta)
        fibre = []
        for h in range(P):
            right = add(r0, q, h)
            left = add(right, q)
            require(det(left, right) == delta, "pointed fibre escaped sector")
            fibre.append((left, right))
        require(len(set(fibre)) == P, "pointed chart did not enumerate fibre")
        pointed.extend(fibre)
    require(len(pointed) == 26208, "pointed-edge census changed")

    first_field_ratios = {}
    first_field_pairs = {}
    field_digests = []
    ratio_histograms = []
    spectral_histograms = []
    mobius_label_counts = []
    mobius_character_counts = []
    cross_ratio_profile_counts = []
    for field_index, (prime, left_bank, right_bank) in enumerate(banks):
        ratio_hist = Counter()
        spectral_hist = Counter()
        digest = sha256()
        constant_sectors = []
        mobius_label_sectors = []
        mobius_character_sectors = []
        cross_ratio_profiles = set()
        z13 = pow(current.MODS[field_index][1], current.NN // P, prime)
        z13_powers = [pow(z13, exponent, prime) for exponent in range(P)]

        # The gauge-free quadratic survivor A(L)A(R), A=Lstar*Rstar, is
        # exactly symmetric and remains nonzero on every endpoint pair.
        vertex_product = {
            point: left_bank[point] * right_bank[point] % prime
            for point in points
        }
        psi = {
            point: left_bank[point] * pow(right_bank[point], -1, prime) % prime
            for point in points
        }
        require(all(vertex_product.values()), "quadratic vertex factor vanished")
        symmetric_support = 0
        for left in points:
            for right in points:
                symmetric = vertex_product[left] * vertex_product[right] % prime
                reverse_symmetric = (
                    vertex_product[right] * vertex_product[left] % prime
                )
                require(symmetric == reverse_symmetric and symmetric,
                        "quadratic two-copy current lost reversal invariance")
                symmetric_support += 1
        require(symmetric_support == P**4, "quadratic survivor support changed")

        for q, delta in sectors:
            r0 = canonical_origin(q, delta)
            ratios = []
            cells = []
            for h in range(P):
                right = add(r0, q, h)
                left = add(right, q)
                forward = left_bank[left] * right_bank[right] % prime
                reverse = left_bank[right] * right_bank[left] % prime
                require(forward and reverse,
                        "pointed current or its row reversal vanished")
                ratio = reverse * pow(forward, -1, prime) % prime
                ratios.append(ratio)
                cells.append((left, right, forward, reverse))
                digest.update(ratio.to_bytes((prime.bit_length() + 7) // 8, "big"))
            ratio_hist[len(set(ratios))] += 1
            product = 1
            for ratio in ratios:
                product = product * ratio % prime
            require(product == 1, "full-cycle reversal coboundary has holonomy")
            for h, (left, right, _forward, _reverse) in enumerate(cells):
                require(
                    ratios[h] == psi[right] * pow(psi[left], -1, prime) % prime,
                    "rank-one reversal ratio is not the vertex coboundary",
                )
            punctured_product = 1
            for h in range(1, 12):
                punctured_product = punctured_product * ratios[h] % prime
            r1 = cells[1][1]
            r12 = cells[11][0]
            require(
                punctured_product == psi[r1] * pow(psi[r12], -1, prime) % prime,
                "punctured path did not retain only its endpoint ratio",
            )

            # Natural cyclic Fourier complexity.  One character would be a
            # geometric H-drift; two characters would be affine-character.
            spectral_support = 0
            for frequency in range(P):
                coefficient = sum(
                    ratios[h] * z13_powers[-frequency * h % P]
                    for h in range(P)
                ) % prime
                spectral_support += int(coefficient != 0)
            spectral_hist[spectral_support] += 1

            # Coordinate-label Mobius law rho_h=(a h+b)/(c h+d).
            label_matrix = [
                (h, 1, -ratios[h] * h, -ratios[h])
                for h in range(P)
            ]
            if field_rank(label_matrix, prime) < 4:
                mobius_label_sectors.append((q, delta))

            # More intrinsic cyclic version: a degree-one Mobius transform of
            # one nontrivial character zeta^(a h).
            character_laws = []
            for frequency in range(1, P):
                character_matrix = []
                for h in range(P):
                    z = z13_powers[frequency * h % P]
                    character_matrix.append(
                        (z, 1, -ratios[h] * z, -ratios[h])
                    )
                if field_rank(character_matrix, prime) < 4:
                    character_laws.append(frequency)
            if character_laws:
                mobius_character_sectors.append((q, delta, tuple(character_laws)))

            cross_ratio_profiles.add(tuple(
                cross_ratio(ratios[0], ratios[1], ratios[2], ratios[h], prime)
                for h in range(3, P)
            ))
            if len(set(ratios)) == 1:
                constant_sectors.append((q, delta))
            if field_index == 0:
                first_field_ratios[(q, delta)] = tuple(ratios)
                first_field_pairs[(q, delta)] = tuple(cells)
        require(not constant_sectors,
                "a nondegenerate sector unexpectedly has scalar reversal")
        require(ratio_hist == Counter({P: 2016}),
                "reversal-ratio fibre census changed")
        require(spectral_hist == Counter({P: 2016}),
                "reversal drift lost full cyclic Fourier support")
        require(not mobius_label_sectors,
                "a coordinate-label Mobius reversal law appeared")
        require(not mobius_character_sectors,
                "a degree-one character-Mobius reversal law appeared")
        require(len(cross_ratio_profiles) == 2016,
                "two sectors acquired the same pointed cross-ratio profile")
        ratio_histograms.append(tuple(sorted(ratio_hist.items())))
        spectral_histograms.append(tuple(sorted(spectral_hist.items())))
        mobius_label_counts.append(len(mobius_label_sectors))
        mobius_character_counts.append(len(mobius_character_sectors))
        cross_ratio_profile_counts.append(len(cross_ratio_profiles))
        field_digests.append(digest.hexdigest())
        require(digest.hexdigest() == EXPECTED_RATIO_DIGESTS[field_index],
                "reversal-ratio digest changed")
        print(
            f"field={field_index + 1} p={prime} "
            f"pointed_support=26208 reversed_support=26208 "
            f"constant_reversal_sectors={len(constant_sectors)} "
            f"distinct_ratio_hist={tuple(sorted(ratio_hist.items()))}"
        )
        print(f"field={field_index + 1} ratio_digest={digest.hexdigest()}")
        print(
            f"field={field_index + 1} spectral_support_hist="
            f"{tuple(sorted(spectral_hist.items()))} "
            f"mobius_label_sectors={len(mobius_label_sectors)} "
            f"mobius_character_sectors={len(mobius_character_sectors)} "
            f"cross_ratio_profiles={len(cross_ratio_profiles)}"
        )

    # Pick the lexicographically first field-1 witness and verify the same
    # cross-multiplied inequality in both embeddings.
    witness_sector = sectors[0]
    ratios = first_field_ratios[witness_sector]
    h0 = 0
    h1 = next(h for h in range(1, P) if ratios[h] != ratios[h0])
    witness = []
    for prime, left_bank, right_bank in banks:
        q, delta = witness_sector
        r0 = canonical_origin(q, delta)
        row = []
        for h in (h0, h1):
            right = add(r0, q, h)
            left = add(right, q)
            forward = left_bank[left] * right_bank[right] % prime
            reverse = left_bank[right] * right_bank[left] % prime
            row.append((left, right, forward, reverse))
        cross_difference = (
            row[0][3] * row[1][2] - row[1][3] * row[0][2]
        ) % prime
        require(cross_difference,
                "selected reversal-ratio witness collapsed in an embedding")
        witness.append(cross_difference)

    q, delta = witness_sector
    r0 = canonical_origin(q, delta)
    print(
        f"witness_sector=(q={q},Delta={delta}), R0={r0}, "
        f"h_pair=({h0},{h1}), cross_differences={tuple(witness)}"
    )
    print(f"ratio_histograms={tuple(ratio_histograms)}")
    print(f"spectral_histograms={tuple(spectral_histograms)}")
    print(f"mobius_label_counts={tuple(mobius_label_counts)}")
    print(f"mobius_character_counts={tuple(mobius_character_counts)}")
    print(f"cross_ratio_profile_counts={tuple(cross_ratio_profile_counts)}")
    print(f"ratio_digests={tuple(field_digests)}")
    require(tuple(witness)
            == (98876734886111000, 382444748747993057),
            "minimal reversal witness changed")
    print("pointed-origin conclusion: every h, including h=3, survives")
    print("reversal conclusion: support survives, current does not descend to a sector scalar")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
