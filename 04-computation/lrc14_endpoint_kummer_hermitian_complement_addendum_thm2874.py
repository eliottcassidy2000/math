#!/usr/bin/env python3
"""Exact Hermitian and full-pattern-complement addendum for THM-2874.

The script proves two coefficient-level statements:

* recombining the THM-2868 Prony summands gives exactly the normalized
  THM-2861 Hermitian edge; and
* the full endpoint-pattern Boolean complement supplies the q7 chart in
  THM-2874's cheapest test, but the resulting three-chart seam is flat.

It does not construct a typed E3/complement contraction or an LRC current.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_endpoint_kummer_galois_bockstein_thm2874.py":
        "3f15c44dc5f66c660ac3605cc25814adc39594bf193aa130a0f5353d6a6178b0",
    RESULTS / "lrc14_endpoint_kummer_galois_bockstein_thm2874.out":
        "90b993b56508ef3603f94104596b899ed9ec7084a2b58ead1604882873ef5455",
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    RESULTS / "lrc14_two_origin_endpoint_projective_kummer_thm2868.out":
        "ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9",
    COMP / "lrc14_endpoint_hermitian_edge_holonomy_thm2861.py":
        "57bad76968ec9c61d2202331e007860f2817d15c606d8ba558ab8b8d3c41f20c",
    RESULTS / "lrc14_endpoint_hermitian_edge_holonomy_thm2861.out":
        "9bc846b6269b6ca967d32b5b4091ec506b3ede632c58a249c78211e1ecc8b43d",
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    RESULTS / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.out":
        "155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas


allocation = atlas.allocation
endpoint_base = atlas.endpoint_base
endpoint = atlas.endpoint
horn = atlas.horn
P = 13
SEAM = (
    # frequency section, target residue, truth block, local Prony offset
    (0, 3, "E3", 0),
    (4, 7, "not-E3", 1),
    (8, 11, "E3", 2),
)


def complement(intervals, modulus):
    result = []
    cursor = 0
    for left, right in intervals:
        require(0 <= left < right <= modulus, "bad endpoint interval")
        require(cursor <= left, "endpoint intervals overlap")
        if cursor < left:
            result.append((cursor, left))
        cursor = right
    if cursor < modulus:
        result.append((cursor, modulus))
    return tuple(result)


def certify_partition(left, right, modulus):
    cursor = 0
    for start, stop in sorted((*left, *right)):
        require(start == cursor and start < stop, "Boolean partition has a gap")
        cursor = stop
    require(cursor == modulus, "Boolean partition misses the circle endpoint")


def dft(values, omega, prime):
    return tuple(
        sum(
            value * pow(omega, (-character * index) % P, prime)
            for index, value in enumerate(values)
        )
        % prime
        for character in range(P)
    )


def natural_lift(state, step):
    ancestry, residue = state
    return (
        ancestry + (residue + step) // P,
        (residue + step) % P,
    )


def main() -> None:
    (
        _module, full_module, _details, e3, clocks, _q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    unit = period // P
    interval = allocation.ATOM_INTERVAL
    target = tuple(value + allocation.physical.SHIFT for value in interval)
    target_atom = ((*target, 1),)

    ell = endpoint_base.REPS[allocation.RIGHT_ORIGIN]
    present = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
    absent = complement(present, endpoint_base.T)
    require(
        ell == (0,) * 9
        and endpoint_base.T == period
        and present == tuple(e3),
        "endpoint full pattern is not the physical macro-E3 locus",
    )
    certify_partition(present, absent, period)

    def literal_bits(atom_interval):
        bits = []
        for index, speed in enumerate(endpoint.W):
            if index == 0:
                dangerous = full_module.make_comb(
                    speed, 91, -13 - 7 * ell[index], 13 - 7 * ell[index]
                )
            else:
                dangerous = full_module.make_comb(
                    speed, 182,
                    -13 - 14 * ell[index], 13 - 14 * ell[index],
                )
            hit = allocation.intersect_sorted((atom_interval,), dangerous)
            if allocation.contained(atom_interval, dangerous):
                bits.append("D")
            elif not hit:
                bits.append("S")
            else:
                raise RuntimeError(f"factor {index} is partial on the atom")
        return "".join(bits)

    orbit_patterns = tuple(
        (
            residue,
            literal_bits(horn.circular_shift_interval(
                target, residue * unit, period
            )),
        )
        for residue in range(P)
    )
    expected_orbit_patterns = (
        (0, "SSSSSSSSD"),
        (1, "SSDSSSSSD"),
        (2, "SSDSSSSSD"),
        (3, "SSSSSSSSD"),
        (4, "SDSSSSSSD"),
        (5, "DDSSSSSSD"),
        (6, "DSSSSDSSD"),
        (7, "DSSSSDSSD"),
        (8, "DSSSSSSSD"),
        (9, "SSSSDSSSD"),
        (10, "SSSSDSSSD"),
        (11, "SSSSSSSSD"),
        (12, "SSSDSSSSD"),
    )
    guard_q5_pairs = tuple(
        (bits[0], bits[5]) for _residue, bits in orbit_patterns
    )
    require(
        orbit_patterns == expected_orbit_patterns
        and set(guard_q5_pairs) == {("S", "S"), ("D", "S"), ("D", "D")}
        and ("S", "D") not in guard_q5_pairs
        and all(
            (representative[0], representative[5]) == (0, 0)
            for representative in endpoint_base.REPS.values()
        ),
        "uniform guard/q5 three-corner horn changed",
    )

    truth_blocks = {"E3": present, "not-E3": absent}

    def restricted_piece(residue, block):
        truth = truth_blocks[block]
        shifted = allocation.physical.overlap.shift_weighted(
            target_atom, residue * unit
        )
        return allocation.indexed_weighted_intersection(
            shifted, truth, tuple(left for left, _right in truth)
        )

    cells3 = atlas.full_cells(
        3, interval, target, unit, period, full_module, e3, clocks
    )
    cells11 = atlas.full_cells(
        11, interval, target, unit, period, full_module, e3, clocks
    )
    common = tuple(sorted(set(cells3) & set(cells11)))
    target7 = horn.circular_shift_interval(target, 7 * unit, period)
    q7_rows = tuple(
        (
            cell,
            atlas.signature(target7, *cell, full_module, e3, clocks),
        )
        for cell in common
    )
    horn20 = tuple(
        cell
        for cell, bits in q7_rows
        if bits == (False, True, True, True, True, True)
    )
    expected20 = tuple(
        (s, t, 1)
        for s in (0, 3, 8, 9, 12)
        for t in (5, 6, 9, 10)
    )
    require(horn20 == expected20, "20-cell E3-only horn changed")

    seam_pieces = {}
    for _r, residue, block, _offset in SEAM:
        piece = restricted_piece(residue, block)
        opposite = "not-E3" if block == "E3" else "E3"
        expected_piece = allocation.physical.overlap.shift_weighted(
            target_atom, residue * unit
        )
        require(
            piece == expected_piece
            and len(piece) == 1
            and piece[0][2] == 1
            and piece[0][1] - piece[0][0] == 26444880
            and not restricted_piece(residue, opposite),
            f"truth-polarized endpoint atom changed at q={residue}",
        )
        require(
            allocation.physical.overlap.shift_weighted(
                piece, -residue * unit
            ) == target_atom,
            f"q={residue} atom is not the full cyclic translate",
        )
        seam_pieces[residue, block] = piece

    require(
        (12 * endpoint.RDIL * unit) % endpoint.NN == 0
        and (26 * endpoint.RDIL * unit) % endpoint.NN == 0,
        "target translation acquired an endpoint phase",
    )

    raw_multipliers = tuple(
        sample
        for r, _q, _block, offset in SEAM
        for sample in (
            1 + 42 * r + offset,
            2 + 42 * r + offset,
        )
    )
    raw_frequencies = tuple(12 + 26 * sample for sample in raw_multipliers)
    require(
        raw_multipliers == (1, 2, 170, 171, 339, 340)
        and all(gcd(value, 91) == 1 for value in raw_multipliers)
        and all(gcd(value, 91) == 1 for value in raw_frequencies),
        "six seam samples ceased to be 91-units",
    )

    start = (0, 3)
    require(
        natural_lift(start, 8) == (0, 11)
        and natural_lift(natural_lift(start, 8), 9) == (1, 7)
        and natural_lift(start, 4) == (0, 7),
        "THM-2851 Bockstein triangle changed",
    )

    field_rows = []
    q3_piece = seam_pieces[3, "E3"]
    for embedding in endpoint.MODS:
        prime, root = embedding
        xi = pow(root, endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        zeta = pow(xi, 2, prime)
        source_value = atlas.COMMON_SOURCE[prime]
        require(
            pow(xi, 2366, prime) == 1
            and pow(xi, 1183, prime) != 1
            and pow(omega, P, prime) == 1
            and omega != 1,
            "endpoint root normalization changed",
        )

        def endpoint_parameters(piece):
            left, right, weight = piece[0]
            require(weight == 1, "seam atom weight changed")
            return (
                pow(root, 12 * endpoint.RDIL * left % endpoint.NN, prime),
                pow(root, 12 * endpoint.RDIL * right % endpoint.NN, prime),
                pow(root, 26 * endpoint.RDIL * left % endpoint.NN, prime),
                pow(root, 26 * endpoint.RDIL * right % endpoint.NN, prime),
            )

        alpha_left, alpha_right, lambda_left, lambda_right = (
            endpoint_parameters(q3_piece)
        )
        require(
            all(
                endpoint_parameters(piece)
                == (alpha_left, alpha_right, lambda_left, lambda_right)
                for piece in seam_pieces.values()
            )
            and lambda_left == pow(xi, 13, prime)
            and lambda_right == pow(xi, 169, prime)
            and pow(lambda_left, 42, prime) == pow(omega, 3, prime)
            and pow(lambda_right, 42, prime) == 1,
            "translated seam atoms do not share the two Prony nodes",
        )
        inverse_difference = pow(
            (lambda_left - lambda_right) % prime, -1, prime
        )

        def split(piece, r, offset):
            formal = 1 + 42 * r
            actual = formal + offset
            raw = []
            for sample in (actual, actual + 1):
                frequency = -(12 + 26 * sample)
                endpoint_value = allocation.endpoint_sum(
                    piece, frequency, embedding
                )
                require(endpoint_value != 0, "endpoint sample vanished")
                raw.append(source_value * endpoint_value % prime)
            current, current_next = raw
            left_at_actual = (
                current_next - lambda_right * current
            ) * inverse_difference % prime
            right_at_actual = (
                lambda_left * current - current_next
            ) * inverse_difference % prime
            transported_left = (
                left_at_actual
                * pow(pow(lambda_left, offset, prime), -1, prime)
            ) % prime
            transported_right = (
                right_at_actual
                * pow(pow(lambda_right, offset, prime), -1, prime)
            ) % prime
            require(
                transported_left
                == (
                    source_value
                    * alpha_left
                    * pow(lambda_left, formal, prime)
                )
                % prime
                and transported_right
                == (
                    -source_value
                    * alpha_right
                    * pow(lambda_right, formal, prime)
                )
                % prime
                and transported_left != 0
                and transported_right != 0,
                f"Prony split formula failed at r={r}",
            )
            return tuple(raw), transported_left, transported_right

        full_u = []
        full_v = []
        for r, offset in enumerate(atlas.SECTION_OFFSETS):
            _raw, left, right = split(q3_piece, r, offset)
            full_u.append(left)
            full_v.append(right)
        full_u = tuple(full_u)
        full_v = tuple(full_v)
        require(
            all(
                full_u[(r + 1) % P]
                == pow(omega, 3, prime) * full_u[r] % prime
                and full_v[(r + 1) % P] == full_v[r]
                for r in range(P)
            ),
            "full frequency split cycle changed",
        )

        canonical_a = pow(zeta, 624, prime)
        canonical_b = pow(zeta, 510, prime)
        canonical_c = tuple(
            (
                canonical_a
                - canonical_b * pow(omega, 3 * r, prime)
            )
            % prime
            for r in range(P)
        )
        canonical_bar_c = tuple(
            (
                pow(canonical_a, -1, prime)
                - pow(canonical_b, -1, prime)
                * pow(omega, -3 * r, prime)
            )
            % prime
            for r in range(P)
        )
        full_s = tuple(
            (left + right) % prime for left, right in zip(full_u, full_v)
        )
        ratio = tuple(
            left * pow(right, -1, prime) % prime
            for left, right in zip(full_u, full_v)
        )
        require(
            full_v[0] == source_value * canonical_a % prime
            and all(value == full_v[0] for value in full_v)
            and all(
                full_u[r]
                == (
                    -source_value
                    * canonical_b
                    * pow(omega, 3 * r, prime)
                )
                % prime
                and full_s[r] == source_value * canonical_c[r] % prime
                and ratio[r]
                == (
                    -canonical_b
                    * pow(canonical_a, -1, prime)
                    * pow(omega, 3 * r, prime)
                )
                % prime
                and (1 + ratio[r]) % prime
                == canonical_c[r] * pow(canonical_a, -1, prime) % prime
                for r in range(P)
            ),
            "Kummer/Galois full-current recombination changed",
        )

        normalized_s = tuple((1 + value) % prime for value in ratio)
        hermitian_edge = tuple(
            normalized_s[(r + 1) % P]
            * (1 + pow(ratio[r], -1, prime))
            % prime
            for r in range(P)
        )
        canonical_edge = tuple(
            canonical_c[(r + 1) % P] * canonical_bar_c[r] % prime
            for r in range(P)
        )
        conjugate_edge = tuple(
            canonical_bar_c[(r + 1) % P] * canonical_c[r] % prime
            for r in range(P)
        )
        transform = dft(hermitian_edge, omega, prime)
        edge_support = tuple(
            character for character, value in enumerate(transform) if value
        )
        require(
            canonical_a * pow(canonical_a, -1, prime) % prime == 1
            and hermitian_edge == canonical_edge
            and all(
                hermitian_edge[r]
                == pow(omega, 3, prime) * conjugate_edge[r] % prime
                for r in range(P)
            )
            and edge_support == (0, 3, 10)
            and len(set(hermitian_edge)) == P
            and all(hermitian_edge),
            "recombined THM-2861 Hermitian edge changed",
        )
        trace = tuple(
            (hermitian_edge[r] + conjugate_edge[r]) % prime
            for r in range(P)
        )
        trace_recovery = (
            pow(omega, 3, prime)
            * pow((1 + pow(omega, 3, prime)) % prime, -1, prime)
        ) % prime
        require(
            all(trace)
            and all(
                hermitian_edge[r] == trace_recovery * trace[r] % prime
                for r in range(P)
            )
            and tuple(
                (conjugate_edge[r] + hermitian_edge[r]) % prime
                for r in range(P)
            ) == trace,
            "Hermitian symmetric-trace/reversal hostile changed",
        )

        seam_u = []
        seam_v = []
        seam_raw = []
        for r, residue, block, offset in SEAM:
            raw, left, right = split(
                seam_pieces[residue, block], r, offset
            )
            require(
                left == full_u[r] and right == full_v[r],
                "truth-polarized chart changed the split coefficient",
            )
            seam_raw.append(raw)
            seam_u.append(left)
            seam_v.append(right)
        seam_u = tuple(seam_u)
        seam_v = tuple(seam_v)
        seam_rq = tuple((r, residue) for r, residue, *_ in SEAM)
        invariant_characters = tuple(
            target_character
            for target_character in range(P)
            if len({
                (3 * r + target_character * residue) % P
                for r, residue in seam_rq
            }) == 1
        )
        compensated_u = tuple(
            seam_u[index] * pow(omega, 10 * residue, prime) % prime
            for index, (_r, residue) in enumerate(seam_rq)
        )
        compensated_v = tuple(
            seam_v[index] * pow(omega, 10 * residue, prime) % prime
            for index, (_r, residue) in enumerate(seam_rq)
        )
        require(
            invariant_characters == (10,)
            and tuple(
                (3 * r + 10 * residue) % P for r, residue in seam_rq
            ) == (4, 4, 4)
            and all(residue - r == 3 for r, residue in seam_rq)
            and len(set(compensated_u)) == 1
            and len(set(compensated_v)) == 3,
            "unique (3,10) full-complement seam invariant changed",
        )

        def rho(step):
            return pow(omega, 3 * step, prime)

        frequency_holonomy = (
            rho(9) * rho(8) * pow(rho(4), -1, prime)
        ) % prime
        carry_character_three = pow(omega, 3, prime)
        require(
            rho(9) == omega
            and frequency_holonomy == 1
            and carry_character_three != 1
            and all(
                rho(left) * rho(right) % prime
                == rho((left + right) % P)
                for left in range(P) for right in range(P)
            ),
            "flat frequency versus character-three Bockstein check failed",
        )
        field_rows.append((
            prime,
            full_u[0],
            full_v[0],
            edge_support,
            trace_recovery,
            compensated_u[0],
            rho(9),
            frequency_holonomy,
            carry_character_three,
            tuple(seam_raw),
        ))

    print("THM-2874 HERMITIAN / FULL-PATTERN-COMPLEMENT ADDENDUM")
    print(
        f"horn20={len(horn20)}; first={horn20[0]}; last={horn20[-1]}; "
        "q7_signature=(0,1,1,1,1,1)"
    )
    print(
        "truth_locus=endpoint_PAT_E3_equals_macro_E3_at_zero_chart; "
        "full_pattern_complement_is_disjoint_and_exhaustive"
    )
    print(
        f"literal_q_patterns={orbit_patterns}; "
        "guard_q5_corners={(S,S),(D,S),(D,D)}; "
        "missing_corner=(S,D); "
        "guard_q5_projection_uniform_over_169_addresses=1"
    )
    print(
        f"seam={SEAM}; raw_m={raw_multipliers}; raw_Y={raw_frequencies}; "
        "all_are_91_units=1"
    )
    for row in field_rows:
        print(
            f"field={row[0]}; U0={row[1]}; V0={row[2]}; "
            f"Hermitian_support={row[3]}; trace_recovery={row[4]}; "
            f"seam_invariant={row[5]}; q11_to_q7_ratio={row[6]}; "
            f"frequency_holonomy={row[7]}; chi3_T={row[8]}; "
            f"raw_pairs={row[9]}"
        )
    print(
        "recombination=S_r=P*c_r and V0=P*A; "
        "S_(r+1)bar(S_r)=P*bar(P)*c_(r+1)bar(c_r)"
    )
    print(
        "Hermitian=edge_is_exact_THM2861; support={0,3,10}; "
        "13_distinct=1; symmetric_trace_recovers_only_with_forward_"
        "frequency_orientation=1"
    )
    print(
        "complement_test=q3/E3,q7/not-E3,q11/E3 are full translated atoms; "
        "U_r*omega^(10q) is constant and b=10 is unique"
    )
    print(
        "Bockstein_verdict=q11_to_q7_projective_ratio_is_omega, but "
        "formal_coefficient_frequency_triangle_holonomy=1 while "
        "aligned_chi3(T)=omega^3"
    )
    print(
        "scope=exact graded coefficient charts and decisive flat hostile; "
        "no typed E3/complement contraction, QA/QAB alignment, positive "
        "current, ancestry action, row exclusion, or LRC14 proof"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
