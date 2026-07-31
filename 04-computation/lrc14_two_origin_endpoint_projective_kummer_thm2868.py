#!/usr/bin/env python3
"""Exact full-current Prony splitter on the THM-2847 42-cell horn.

This is the exact companion for the THM-2868 proved candidate.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


PINNED = {
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    COMP / "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
}
for path, expected in PINNED.items():
    actual = sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, f"pinned dependency changed: {path.name}")


import lrc14_q3_q11_transverse_endpoint_horn_thm2847 as horn


allocation = horn.allocation
endpoint_base = allocation.endpoint_base
endpoint = endpoint_base.endpoint
P = 13
ORIGIN = allocation.RIGHT_ORIGIN
TARGET_STEP = allocation.TARGET_STEP
STEPPED = allocation.add(ORIGIN, TARGET_STEP)
COMMON_SOURCE = {
    352341050142921841: 254455016269350867,
    956354278959359281: 318932490657369324,
}
PRONY_WINDOW = (1, 2, 3, 4)
FREQUENCY_SECTIONS = tuple(1 + 42 * r for r in range(P))
SECTION_OFFSETS = tuple(
    1 if r == 4 else 2 if r == 8 else 0 for r in range(P)
)
UNIT_SECTIONS = tuple(
    m + offset
    for m, offset in zip(FREQUENCY_SECTIONS, SECTION_OFFSETS)
)
FREQUENCY_MEASUREMENTS = tuple(sorted(
    {
        sample
        for m in UNIT_SECTIONS
        for sample in (m, m + 1)
    }
    | set(PRONY_WINDOW)
))


def signature(interval, s, t, clock, full_module, e3, clocks):
    return (
        allocation.contained(interval, e3),
        allocation.contained(interval, clocks[clock]),
        horn.safe_comb_contains(
            interval, full_module, full_module.W[1], 182,
            -14 * s - 13, -14 * s + 13,
        ),
        horn.safe_comb_contains(
            interval, full_module, full_module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        ),
        horn.safe_comb_contains(
            interval, full_module, full_module.C2, 182,
            14 * s - 13, 14 * s + 13,
        ),
        horn.safe_comb_contains(
            interval, full_module, full_module.C3, 182,
            14 * t - 13, 14 * t + 13,
        ),
    )


def full_cells(q, interval, target, unit, period, full_module, e3, clocks):
    shifted = horn.circular_shift_interval(target, q * unit, period)
    return tuple(
        (s, t, clock)
        for s in allocation.COMMON_S
        for t in allocation.COMMON_T
        for clock in range(7)
        if all(signature(interval, s, t, clock, full_module, e3, clocks))
        and all(signature(target, s, t, clock, full_module, e3, clocks))
        and all(signature(shifted, s, t, clock, full_module, e3, clocks))
    )


def indexed_restriction(carrier, address):
    ell = endpoint_base.REPS[address]
    present = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
    starts = tuple(left for left, _right in present)
    return allocation.indexed_weighted_intersection(
        carrier, present, starts
    )


def shifted_chart(address, carrier):
    ell = endpoint_base.REPS[address]
    present = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
    shifted_ell = tuple(
        (ell[index] + endpoint_base.WMOD[index]) % P
        for index in range(9)
    )
    shifted_present = allocation.physical.overlap.shift_union(
        present, -(endpoint_base.T // P)
    )
    starts = tuple(left for left, _right in shifted_present)
    restricted = allocation.indexed_weighted_intersection(
        carrier, shifted_present, starts
    )
    return shifted_ell, restricted


def weighted_endpoint_sum(restricted, frequency, embedding):
    return allocation.endpoint_sum(restricted, frequency, embedding)


def main() -> None:
    projective_exponents = tuple(
        (1111 - 156 * m) % 2366 for m in PRONY_WINDOW
    )
    projective_norm_exponents = tuple(
        13 * exponent % 2366 for exponent in projective_exponents
    )
    require(
        projective_exponents == (955, 799, 643, 487)
        and projective_norm_exponents == (585, 923, 1261, 1599)
        and all(gcd(exponent, 2366) == 1
                for exponent in projective_exponents)
        and all((1275 * exponent - exponent) % 2366 == 546
                for exponent in projective_exponents),
        "projective Kummer torsor arithmetic failed",
    )
    require(
        SECTION_OFFSETS
        == (0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0)
        and tuple(
            r for r, m in enumerate(FREQUENCY_SECTIONS)
            if m % P == 0
        ) == (4,)
        and tuple(
            r for r, m in enumerate(FREQUENCY_SECTIONS)
            if (m + 1) % P == 0
        ) == (8,)
        and all(
            gcd(sample, 91) == 1
            for m in UNIT_SECTIONS for sample in (m, m + 1)
        )
        and all(
            gcd(12 + 26 * sample, 91) == 1
            for m in UNIT_SECTIONS for sample in (m, m + 1)
        )
        and tuple(m % 7 for m in UNIT_SECTIONS)
        == (1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 1)
        and tuple((m + 1) % 7 for m in UNIT_SECTIONS)
        == (2, 2, 2, 2, 3, 2, 2, 2, 4, 2, 2, 2, 2),
        "variable-offset 91-unit atlas arithmetic failed",
    )
    physical_projective_exponents = tuple(
        (1111 - 156 * m) % 2366 for m in FREQUENCY_SECTIONS
    )
    require(
        physical_projective_exponents
        == tuple((955 + 546 * r) % 2366 for r in range(P))
        and len(set(physical_projective_exponents)) == P
        and all(
            gcd(exponent, 2366) == 1
            for exponent in physical_projective_exponents
        )
        and all(
            (
                physical_projective_exponents[(r + 1) % P]
                - physical_projective_exponents[r]
            ) % 2366 == 546
            for r in range(P)
        ),
        "physical projective C13 orbit arithmetic failed",
    )
    require(
        (0, 4, 8) == (0, 7 - 3, 11 - 3)
        and (7 - 11) % P == 9
        and FREQUENCY_SECTIONS[4] == 169
        and FREQUENCY_SECTIONS[8] + 1 == 338,
        "q3/q7/q11 seam-triangle alignment changed",
    )
    (
        _module, full_module, _details, e3, clocks, _q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    unit = period // P
    interval = allocation.ATOM_INTERVAL
    target = tuple(x + allocation.physical.SHIFT for x in interval)
    source_atom = ((*interval, 1),)
    target_atom = ((*target, 1),)

    cells3 = full_cells(
        3, interval, target, unit, period, full_module, e3, clocks
    )
    cells11 = full_cells(
        11, interval, target, unit, period, full_module, e3, clocks
    )
    common = tuple(sorted(set(cells3) & set(cells11)))
    expected_common = tuple(
        (s, t, 1)
        for s in (0, 3, 8, 9, 10, 11, 12)
        for t in range(5, 11)
    )
    require(common == expected_common and len(common) == 42,
            "42-cell common horn changed")

    source_restricted = indexed_restriction(source_atom, ORIGIN)
    terminal = tuple(
        endpoint.build_set(endpoint_base.PAT_Q12, endpoint_base.ZERO)
    )
    q_starts = tuple(left for left, _right in terminal)
    direct_restricted = {}
    gauge_restricted = {}
    for address in (ORIGIN, STEPPED):
        for q in (3, 7, 11):
            direct_restricted[address, q] = indexed_restriction(
                allocation.physical.overlap.shift_weighted(
                    target_atom, q * unit
                ),
                address,
            )
            _shifted_ell, gauge_restricted[address, q] = shifted_chart(
                address,
                allocation.physical.overlap.shift_weighted(
                    target_atom, ((q - 1) % P) * unit
                ),
            )

    field_rows = []
    for embedding in endpoint.MODS:
        prime, root = embedding
        z13 = pow(root, endpoint_base.NRED // P, prime)
        xi2366 = pow(root, endpoint.NN // 2366, prime)
        require(
            pow(xi2366, 2366, prime) == 1
            and all(
                pow(xi2366, 2366 // divisor, prime) != 1
                for divisor in (2, 7, 13)
            )
            and z13 == pow(xi2366, 182, prime),
            "local primitive 2366th root certificate failed",
        )

        # These are the literal source values proved and exactly replayed
        # in THM-2806, equation (28).  Reusing them avoids rerunning its
        # expensive x-sweep; the present probe pins that immutable script.
        source_value = COMMON_SOURCE[prime]
        require(source_value != 0 and len(source_restricted) == 1,
                "common source coefficient vanished")
        galois_u = 547
        require(
            gcd(galois_u, endpoint.NN) == 1
            and galois_u % 91 == 1,
            "chosen full-field automorphism does not fix zeta_91",
        )
        conjugate_embedding = (
            prime, pow(root, galois_u, prime)
        )
        conjugate_tabs = endpoint.make_tabs(
            terminal, endpoint_base.X0, (conjugate_embedding,)
        )
        conjugate_source = allocation.x_sweep(
            source_restricted,
            terminal,
            q_starts,
            conjugate_embedding,
            conjugate_tabs,
        )[0]
        require(
            conjugate_source != source_value,
            "common source coefficient unexpectedly became F-rational",
        )

        # THM-2806 already proves the marked-sheet source-gauge identity
        # for this exact x-sweep.  Varying m changes only the separate
        # deepest phase; W has zero deepest coordinate, so the proved
        # identity is literally the same for the whole multiplier bank.
        require(
            endpoint_base.WMOD[endpoint_base.DEEP] == 0,
            "source gauge acquired a deepest-coordinate phase",
        )

        raw = {}
        restricted_bank = {}
        for m in FREQUENCY_MEASUREMENTS:
            frequency = -(12 + 26 * m)
            for address in (ORIGIN, STEPPED):
                for q in (3, 7, 11):
                    restricted = direct_restricted[address, q]
                    restricted_bank[(m, address, q)] = restricted
                    raw[(m, address, q)] = weighted_endpoint_sum(
                        restricted, frequency, embedding
                    )

                    # Exact representative-gauge replay.  The target
                    # carrier index shifts q -> q-1 when the present chart
                    # shifts by -unit.
                    shifted_value = weighted_endpoint_sum(
                        gauge_restricted[address, q],
                        frequency,
                        embedding,
                    )
                    require(
                        shifted_value == raw[(m, address, q)],
                        "target representative gauge failed",
                    )

        require(
            all(
                len(restricted_bank[(m, ORIGIN, q)]) == 1
                for m in FREQUENCY_MEASUREMENTS for q in (3, 11)
            )
            and all(
                len(restricted_bank[(m, STEPPED, 11)]) == 1
                for m in FREQUENCY_MEASUREMENTS
            )
            and all(
                not restricted_bank[(m, STEPPED, 3)]
                and not restricted_bank[(m, ORIGIN, 7)]
                and not restricted_bank[(m, STEPPED, 7)]
                for m in FREQUENCY_MEASUREMENTS
            ),
            "E3 support pattern changed",
        )
        require(
            all(
                raw[(m, ORIGIN, 3)]
                == raw[(m, ORIGIN, 11)]
                == raw[(m, STEPPED, 11)]
                != 0
                and raw[(m, STEPPED, 3)] == 0
                and raw[(m, ORIGIN, 7)] == 0
                and raw[(m, STEPPED, 7)] == 0
                for m in FREQUENCY_MEASUREMENTS
            ),
            "multiplier-bank endpoint support/value pattern changed",
        )

        endpoint_interval = restricted_bank[(1, ORIGIN, 3)][0]
        left, right, weight = endpoint_interval
        require(weight == 1, "endpoint interval weight changed")
        alpha_left = pow(
            root, (12 * endpoint.RDIL * left) % endpoint.NN, prime
        )
        alpha_right = pow(
            root, (12 * endpoint.RDIL * right) % endpoint.NN, prime
        )
        lambda_left = pow(
            root, (26 * endpoint.RDIL * left) % endpoint.NN, prime
        )
        lambda_right = pow(
            root, (26 * endpoint.RDIL * right) % endpoint.NN, prime
        )
        require(
            lambda_left == pow(xi2366, 13, prime)
            and lambda_right == pow(xi2366, 169, prime)
            and lambda_left != lambda_right,
            "Prony nodes or their primitive-2366 normalization changed",
        )

        signed_endpoint = tuple(
            (
                raw[(m, ORIGIN, 3)]
                - raw[(m, STEPPED, 3)]
            ) % prime
            for m in PRONY_WINDOW
        )
        signed_full = tuple(
            source_value * value % prime for value in signed_endpoint
        )
        signed_q11 = tuple(
            source_value
            * (
                raw[(m, ORIGIN, 11)]
                - raw[(m, STEPPED, 11)]
            )
            % prime
            for m in PRONY_WINDOW
        )
        require(all(signed_full) and not any(signed_q11),
                "signed full-current selector failed")

        expected_endpoint = tuple(
            (
                alpha_left * pow(lambda_left, m, prime)
                - alpha_right * pow(lambda_right, m, prime)
            ) % prime
            for m in PRONY_WINDOW
        )
        require(signed_endpoint == expected_endpoint,
                "literal endpoint sequence is not two-node Prony")
        node_sum = (lambda_left + lambda_right) % prime
        node_product = lambda_left * lambda_right % prime
        require(
            signed_full[2]
            == (
                node_sum * signed_full[1]
                - node_product * signed_full[0]
            ) % prime
            and signed_full[3]
            == (
                node_sum * signed_full[2]
                - node_product * signed_full[1]
            ) % prime,
            "full-current Prony recurrence failed",
        )

        inverse_difference = pow(
            (lambda_left - lambda_right) % prime, -1, prime
        )
        split_left = (
            signed_full[1] - lambda_right * signed_full[0]
        ) * inverse_difference % prime
        split_minus_right = (
            lambda_left * signed_full[0] - signed_full[1]
        ) * inverse_difference % prime
        require(
            split_left
            == source_value * alpha_left * lambda_left % prime
            and split_minus_right
            == -source_value * alpha_right * lambda_right % prime
            and split_left != 0
            and split_minus_right != 0,
            "oriented full-current endpoint split failed",
        )

        # The endpoint numerator projector is already charged before the
        # common source is restored.
        h_endpoint = (
            signed_endpoint[1]
            - lambda_right * signed_endpoint[0]
        ) % prime
        require(
            h_endpoint
            == (
                (lambda_left - lambda_right)
                * alpha_left * lambda_left
            ) % prime
            and h_endpoint != 0,
            "two-sample charged endpoint projector failed",
        )

        # Build the physical frequency-dual C13 atlas from actual 91-unit
        # multiplier pairs.  At the two index-zero sections r=4,8, shift
        # the local Prony window by e_r=1,2 and transport the isolated
        # nodes back by their known powers.  Thus no raw multiplier is
        # divisible by 7 or 13, even though the reconstructed formal
        # sections are the desired m_r=1+42r.
        require(
            pow(lambda_left, 42, prime) == pow(z13, 3, prime)
            and pow(lambda_right, 42, prime) == 1,
            "endpoint nodes do not have chi_3 plus trivial monodromy",
        )
        physical_left = []
        physical_minus_right = []
        physical_ratio = []
        for r, (m, offset, n) in enumerate(zip(
            FREQUENCY_SECTIONS, SECTION_OFFSETS, UNIT_SECTIONS
        )):
            current_n = (
                source_value
                * (
                    raw[(n, ORIGIN, 3)]
                    - raw[(n, STEPPED, 3)]
                )
            ) % prime
            current_next = (
                source_value
                * (
                    raw[(n + 1, ORIGIN, 3)]
                    - raw[(n + 1, STEPPED, 3)]
                )
            ) % prime
            split_left_n = (
                current_next - lambda_right * current_n
            ) * inverse_difference % prime
            split_minus_right_n = (
                lambda_left * current_n - current_next
            ) * inverse_difference % prime
            transported_left = (
                split_left_n
                * pow(pow(lambda_left, offset, prime), -1, prime)
            ) % prime
            transported_minus_right = (
                split_minus_right_n
                * pow(pow(lambda_right, offset, prime), -1, prime)
            ) % prime
            require(
                transported_left
                == source_value
                * alpha_left
                * pow(lambda_left, m, prime)
                % prime
                and transported_minus_right
                == -source_value
                * alpha_right
                * pow(lambda_right, m, prime)
                % prime
                and transported_left != 0
                and transported_minus_right != 0,
                f"unit-chart Prony transport failed at section {r}",
            )
            physical_left.append(transported_left)
            physical_minus_right.append(transported_minus_right)
            physical_ratio.append(
                transported_left
                * pow(transported_minus_right, -1, prime)
                % prime
            )
        physical_left = tuple(physical_left)
        physical_minus_right = tuple(physical_minus_right)
        physical_ratio = tuple(physical_ratio)
        require(
            all(
                physical_left[(r + 1) % P]
                == pow(z13, 3, prime) * physical_left[r] % prime
                and physical_minus_right[(r + 1) % P]
                == physical_minus_right[r]
                and physical_ratio[(r + 1) % P]
                == pow(z13, 3, prime) * physical_ratio[r] % prime
                for r in range(P)
            )
            and len(set(physical_ratio)) == P
            and all(
                physical_ratio[r]
                == pow(xi2366, physical_projective_exponents[r], prime)
                for r in range(P)
            ),
            "transported physical frequency orbit failed",
        )

        # Put the transported split-left values over the q3 target row:
        # H(r,q)=U_r delta_3(q).  Its normalized two-axis DFT has one
        # multiplier-frequency row and all thirteen target columns.
        inverse_169 = pow(P * P, -1, prime)
        multiplier_target_dft = {}
        for multiplier_character in range(P):
            for target_frequency in range(P):
                total = 0
                for frequency_section in range(P):
                    value = physical_left[frequency_section]
                    total += (
                        value
                        * pow(
                            z13,
                            (
                                -multiplier_character * frequency_section
                                - target_frequency * 3
                            ) % P,
                            prime,
                        )
                    )
                multiplier_target_dft[
                    multiplier_character, target_frequency
                ] = total * inverse_169 % prime
        dft_support = tuple(
            key for key, value in multiplier_target_dft.items() if value
        )
        expected_support = tuple((3, b) for b in range(P))
        require(
            dft_support == expected_support
            and all(
                multiplier_target_dft[3, b]
                == physical_left[0]
                * pow(P, -1, prime)
                * pow(z13, -3 * b, prime)
                % prime
                for b in range(P)
            ),
            "physical multiplier-target DFT support or normalization failed",
        )
        diagonal_channels = tuple(
            key for key in dft_support if sum(key) % P == 0
        )
        require(
            diagonal_channels == ((3, 10),)
            and multiplier_target_dft[3, 10] != 0,
            "unique diagonal-invariant multiplier channel failed",
        )

        # The literal translation increments vanish for every multiplier:
        # no hidden endpoint phase is introduced by representative gauge.
        phase_rows = tuple(
            (
                (-12 * endpoint.RDIL * unit) % endpoint.NN,
                ((12 + 26 * m) * endpoint.RDIL * unit) % endpoint.NN,
            )
            for m in FREQUENCY_MEASUREMENTS
        )
        require(
            phase_rows == ((0, 0),) * len(FREQUENCY_MEASUREMENTS),
            "multiplier bank acquired a hidden gauge phase",
        )
        field_rows.append((
            prime,
            source_value,
            signed_full,
            split_left,
            split_minus_right,
            physical_left[0],
            physical_minus_right[0],
            physical_ratio[0],
            multiplier_target_dft[3, 10],
            conjugate_source,
        ))

    print("THM-2868 COMMON SIGNED FULL-CURRENT PRONY PROBE")
    print(
        "horn=q3/q11 common 42 cells; "
        f"first={common[0]}; last={common[-1]}"
    )
    print("Prony_window=(1,2,3,4); right_frequencies=(38,64,90,116)")
    print(
        f"frequency_sections={FREQUENCY_SECTIONS}; "
        f"unit_offsets={SECTION_OFFSETS}; "
        f"unit_pairs={tuple((n, n + 1) for n in UNIT_SECTIONS)}"
    )
    print(
        "unit_right_frequency_pairs="
        f"{tuple((12 + 26 * n, 12 + 26 * (n + 1)) for n in UNIT_SECTIONS)}"
    )
    for row in field_rows:
        print(
            f"field={row[0]}; common_source={row[1]}; "
            f"signed_full={row[2]}; "
            f"split_left={row[3]}; split_minus_right={row[4]}; "
            f"physical_U0={row[5]}; physical_V0={row[6]}; "
            f"physical_t0={row[7]}; "
            f"diagonal_clutch_3_10={row[8]}; "
            f"sigma547_common_source={row[9]}"
        )
    print(
        "support=q3 signed selector nonzero; q11 cancels; "
        "q7 zero at both origins"
    )
    print(
        "gauge=marked-sheet representative transport exact for the full "
        "28-sample replay bank in both fields"
    )
    print(
        "physical_multiplier_target_DFT_support={3}xF13; "
        "U_(r+1)=omega^3*U_r; V_(r+1)=V_r; "
        "unique_diagonal_invariant_channel=(3,10); "
        "normalized_value=U0/13*omega^9"
    )
    print(
        f"projective_split_exponents={projective_exponents}; "
        f"13th_power_exponents={projective_norm_exponents}; "
        "each_primitive2366=1; sigma_ratio=omega^3; "
        "minimal_polynomial=X^13-t^13_over_Q(zeta91)"
    )
    print(
        f"physical_projective_exponents={physical_projective_exponents}; "
        "13_distinct_primitive2366=1; physical_ratio_rotation=omega^3"
    )
    print(
        "index_zero_repairs=r4:+1,r8:+2; all_26_raw_multipliers_and_"
        "right_frequencies_are_91_units=1; "
        "seam_triangle_sections=(0,4,8)<->(q3,q7,q11)"
    )
    print(
        "common_source_not_Q(zeta91)=1; witness automorphism u=547 "
        "fixes zeta91 and moves P in both certified fields"
    )
    print(
        "scope=frequency-cleared signed coefficient measurements with "
        "variable local Prony offsets, not one raw affine pair or a positive "
        "packet; common E3 ancestry, same marked triangle, and "
        "q11-to-q7 support transport remain absent"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
