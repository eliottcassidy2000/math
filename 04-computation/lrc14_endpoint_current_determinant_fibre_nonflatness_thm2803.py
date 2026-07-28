#!/usr/bin/env python3
"""Exact companion for THM-2803 determinant-fibre projective nonflatness.

For the canonical THM-2625 endpoint factors reconstructed by the immutable
THM-2779 companion, this script forms

    J_s(R) = P(R+s) Q(R)

on every determinant fibre ``det(s,R)=delta``.  A normalized transverse
vector ``t`` with ``det(s,t)=1`` gives the cycle profile

    j_(s,delta)(w) = J_s(w*s + delta*t),  w in F_13.

It checks all pairs of the thirteen profiles for each of the 168 nonzero
directions.  The comparison permits an independent cyclic origin change and
an arbitrary nonzero scalar on either profile; a stronger diagnostic also
permits reversal.  Projective comparison is performed by exact 2-by-2
minors, never floating-point division.

The script also checks all normalized transverse-frame changes and the
first-forgotten-digit pointwise-power separator used in the proof.

Script:
  04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_thm2803.py
Output:
  05-knowledge/results/lrc14_endpoint_current_determinant_fibre_nonflatness_thm2803.out
"""

import hashlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_bockstein_symplectic_endpoint_gate_thm2779.py"
)
SPEC = importlib.util.spec_from_file_location("thm2779_endpoint_gate_2803", SOURCE)
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)

P = 13
POINTS = GATE.POINTS
NONZERO = GATE.NONZERO
EXPECTED_PRIMES = (352341050142921841, 956354278959359281)
EXPECTED_DIGESTS = (
    (
        "d5b84c7702cead0b4fa13842d6207a0ba3acc7d68c9d0ac334518df9e6cc5cfe",
        "39dbf7d36cfd37a637499591bf45418007a9b26edd1e15b469ab7111a7bba9c2",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
        "0fa02865c34d02f6de3c453f82b5563b1c15c4cbfc1ab6de28bd4d09e5f995e8",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
        "117b842476e1666fb685063ea73e717cd9ba114c27df582af8bd58ada85becdf",
    ),
    (
        "2b0fe21a1f03eccf176398892d47924c764a9dc6850e4e48704452bc4517f31b",
        "671906fa806d5cf4653945cd85d3f2ef0fd15580c6a916e5272bfe0b14a265c7",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
        "9c2289ab8f0e4e44fd88d2623284d7dfadfdf353a78cf9605b05d0ff0cdc5d08",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
        "1af1aed658585a06a71fdc1e2d357da0ffa0b6ce226d3ab8ab5898d64ad6786b",
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(values):
    payload = b"".join(int(value).to_bytes(8, "big") for value in values)
    return hashlib.sha256(payload).hexdigest()


def add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def scale(a, x):
    return (a * x[0] % P, a * x[1] % P)


def transverse(step):
    """Choose a deterministic t with det(step,t)=1."""
    a, b = step
    if a:
        return (0, pow(a, -1, P))
    return (-pow(b, -1, P) % P, 0)


def rotate(values, shift):
    return tuple(values[(index + shift) % P] for index in range(P))


def reverse(values):
    return tuple(values[-index % P] for index in range(P))


def current(field, step, origin):
    prime, left, right = field
    return left[add(origin, step)] * right[origin] % prime


def make_profile(field, step, across, delta):
    return tuple(
        current(field, step, add(scale(w, step), scale(delta, across)))
        for w in range(P)
    )


def projective_minor_data(left, right, prime):
    """Return the number of nonzero Pluecker minors and the first witness."""
    nonzero = 0
    first = None
    for i in range(P):
        for j in range(i + 1, P):
            value = (left[i] * right[j] - left[j] * right[i]) % prime
            if value:
                nonzero += 1
                if first is None:
                    first = value
    return nonzero, 0 if first is None else first


def verify_synthetic_controls(profile, prime):
    """Positive controls for equality, scale, translation, and reversal."""
    require(profile == tuple(profile), "identity equality control")
    count, witness = projective_minor_data(profile, profile, prime)
    require(count == 0 and witness == 0, "identity projective control")

    scalar = 2
    shift = 3
    synthetic = tuple(
        scalar * profile[(index + shift) % P] % prime
        for index in range(P)
    )
    matches = []
    for origin_shift in range(P):
        count, _witness = projective_minor_data(
            profile,
            rotate(synthetic, origin_shift),
            prime,
        )
        if count == 0:
            matches.append(origin_shift)
    require((-shift) % P in matches, "translation-scale positive control")

    synthetic_reflected = reverse(profile)
    reflected = reverse(synthetic_reflected)
    reflected_matches = []
    for origin_shift in range(P):
        count, _witness = projective_minor_data(
            profile,
            rotate(reflected, origin_shift),
            prime,
        )
        if count == 0:
            reflected_matches.append(origin_shift)
    require(0 in reflected_matches, "reversal positive control")
    return len(matches), len(reflected_matches)


def verify_first_forgotten_digit_separator(prime):
    """Finite controls for the all-m pointwise algebra identity.

    The proof itself is formal: a coarse pullback has equal values at r and
    r+13, while H(r)=0 and H(r+13)=1.  We replay the THM-2791 normalized
    two-point chord and several powers to guard the typing and sign.
    """
    chord = tuple(1 if residue in (6, 7) else 0 for residue in range(P))
    coarse = tuple(chord[index % P] for index in range(P * P))
    high_digit = tuple(index // P for index in range(P * P))
    coarse_zero_checks = 0
    marked_formula_checks = 0
    marked_nonzero_checks = 0
    for exponent in range(1, 2 * P + 1):
        coarse_power = tuple(pow(value, exponent, prime) for value in coarse)
        marked_power = tuple(
            high_digit[index] * coarse_power[index] % prime
            for index in range(P * P)
        )
        for residue in range(P):
            coarse_difference = (
                coarse_power[residue] - coarse_power[residue + P]
            ) % prime
            marked_difference = (
                marked_power[residue] - marked_power[residue + P]
            ) % prime
            expected = -pow(chord[residue], exponent, prime) % prime
            require(coarse_difference == 0, "coarse power escaped fibre difference")
            require(
                marked_difference == expected,
                "first-forgotten-digit separator sign or value drift",
            )
            coarse_zero_checks += 1
            marked_formula_checks += 1
            marked_nonzero_checks += marked_difference != 0
    require(marked_nonzero_checks == 2 * (2 * P), "chord separator support drift")
    return coarse_zero_checks, marked_formula_checks, marked_nonzero_checks


def analyze_field(field, field_index):
    prime = field[0]
    require(prime == EXPECTED_PRIMES[field_index], "certified field drift")
    primitive_root = GATE.BASE.MODS[field_index][1]
    zeta = pow(primitive_root, GATE.BASE.NN // P, prime)

    profile_bank = []
    mode_bank = []
    minor_count_bank = []
    witness_bank = []
    reflected_minor_count_bank = []
    reflected_witness_bank = []

    canonical_raw_matches = 0
    canonical_projective_matches = 0
    cyclic_projective_matches = 0
    dihedral_projective_matches = 0
    pair_count = 0
    cyclic_comparisons = 0
    dihedral_comparisons = 0
    frame_entry_checks = 0
    nontrivial_mode_support = 0
    zero_mode_support = 0
    per_step_distinct = []

    first_profile = None
    for step in NONZERO:
        across = transverse(step)
        require(GATE.det(step, across) == 1, "bad normalized transverse frame")
        profiles = tuple(
            make_profile(field, step, across, delta)
            for delta in range(P)
        )
        require(all(all(profile) for profile in profiles), "current support hole")
        if first_profile is None:
            first_profile = profiles[0]
        profile_bank.extend(value for profile in profiles for value in profile)

        for frame_shear in range(P):
            alternate = add(across, scale(frame_shear, step))
            require(
                GATE.det(step, alternate) == 1,
                "normalized frame shear changed determinant",
            )
            for delta in range(P):
                alternate_profile = make_profile(field, step, alternate, delta)
                expected = rotate(profiles[delta], frame_shear * delta % P)
                require(
                    alternate_profile == expected,
                    "transverse-frame shear is not the predicted origin change",
                )
                frame_entry_checks += P

        for profile in profiles:
            zero_mode = sum(profile) % prime
            zero_mode_support += zero_mode != 0
            for character in range(1, P):
                mode = sum(
                    profile[w] * pow(zeta, -character * w, prime)
                    for w in range(P)
                ) % prime
                require(mode != 0, "inherited nontrivial mode support hole")
                nontrivial_mode_support += 1
                mode_bank.append(mode)

        step_raw_distinct = True
        for left_delta in range(P):
            for right_delta in range(left_delta + 1, P):
                left_profile = profiles[left_delta]
                right_profile = profiles[right_delta]
                pair_count += 1
                if left_profile == right_profile:
                    canonical_raw_matches += 1
                    step_raw_distinct = False

                count, _witness = projective_minor_data(
                    left_profile,
                    right_profile,
                    prime,
                )
                if count == 0:
                    canonical_projective_matches += 1

                cyclic_match = False
                dihedral_match = False
                reflected = reverse(right_profile)
                for shift in range(P):
                    count, witness = projective_minor_data(
                        left_profile,
                        rotate(right_profile, shift),
                        prime,
                    )
                    minor_count_bank.append(count)
                    witness_bank.append(witness)
                    cyclic_comparisons += 1
                    cyclic_match = cyclic_match or count == 0

                    reflected_count, reflected_witness = projective_minor_data(
                        left_profile,
                        rotate(reflected, shift),
                        prime,
                    )
                    reflected_minor_count_bank.append(reflected_count)
                    reflected_witness_bank.append(reflected_witness)
                    dihedral_comparisons += 1
                    dihedral_match = dihedral_match or reflected_count == 0
                cyclic_projective_matches += cyclic_match
                dihedral_projective_matches += cyclic_match or dihedral_match
        per_step_distinct.append(step_raw_distinct)

    require(first_profile is not None, "empty profile universe")
    synthetic_matches, reflection_controls = verify_synthetic_controls(
        first_profile,
        prime,
    )
    separator = verify_first_forgotten_digit_separator(prime)

    expected_profiles = len(NONZERO) * P
    expected_pairs = len(NONZERO) * P * (P - 1) // 2
    expected_comparisons = expected_pairs * P
    expected_frame_entries = len(NONZERO) * P * P * P
    expected_modes = expected_profiles * (P - 1)
    require(len(profile_bank) == expected_profiles * P, "profile universe drift")
    require(pair_count == expected_pairs, "profile-pair universe drift")
    require(cyclic_comparisons == expected_comparisons, "cyclic comparison drift")
    require(dihedral_comparisons == expected_comparisons, "dihedral comparison drift")
    require(frame_entry_checks == expected_frame_entries, "frame audit drift")
    require(nontrivial_mode_support == expected_modes, "mode census drift")
    require(zero_mode_support == expected_profiles, "zero Fourier mode support hole")
    require(all(per_step_distinct), "raw determinant-profile collision")
    require(canonical_raw_matches == 0, "raw profile collision")
    require(canonical_projective_matches == 0, "canonical projective collision")
    require(cyclic_projective_matches == 0, "cyclic-projective profile collision")
    require(dihedral_projective_matches == 0, "dihedral-projective profile collision")
    require(
        all(count == P * (P - 1) // 2 for count in minor_count_bank),
        "cyclic pointwise-ratio collision",
    )
    require(all(witness_bank), "cyclic first-minor witness vanished")
    require(
        all(
            count == P * (P - 1) // 2
            for count in reflected_minor_count_bank
        ),
        "dihedral pointwise-ratio collision",
    )
    require(all(reflected_witness_bank), "dihedral first-minor witness vanished")

    result = {
        "prime": prime,
        "profiles": expected_profiles,
        "profile_pairs": expected_pairs,
        "cyclic_comparisons": expected_comparisons,
        "dihedral_comparisons": expected_comparisons,
        "raw_matches": canonical_raw_matches,
        "projective_matches": canonical_projective_matches,
        "cyclic_matches": cyclic_projective_matches,
        "dihedral_matches": dihedral_projective_matches,
        "minor_min": min(minor_count_bank),
        "minor_max": max(minor_count_bank),
        "reflected_minor_min": min(reflected_minor_count_bank),
        "reflected_minor_max": max(reflected_minor_count_bank),
        "frame_entry_checks": frame_entry_checks,
        "nontrivial_modes": nontrivial_mode_support,
        "zero_modes": zero_mode_support,
        "synthetic_matches": synthetic_matches,
        "reflection_controls": reflection_controls,
        "separator": separator,
        "profile_digest": digest(profile_bank),
        "mode_digest": digest(mode_bank),
        "minor_count_digest": digest(minor_count_bank),
        "witness_digest": digest(witness_bank),
        "reflected_minor_count_digest": digest(reflected_minor_count_bank),
        "reflected_witness_digest": digest(reflected_witness_bank),
    }
    require(
        (
            result["profile_digest"],
            result["mode_digest"],
            result["minor_count_digest"],
            result["witness_digest"],
            result["reflected_minor_count_digest"],
            result["reflected_witness_digest"],
        )
        == EXPECTED_DIGESTS[field_index],
        "deterministic field digest drift",
    )
    return result


def main():
    fields = GATE.build_endpoint_factors()
    require(
        tuple(field[0] for field in fields) == EXPECTED_PRIMES,
        "certified dual-field universe drift",
    )
    results = tuple(
        analyze_field(field, field_index)
        for field_index, field in enumerate(fields)
    )

    print("THM-2803 ENDPOINT-CURRENT DETERMINANT-FIBRE PROJECTIVE NONFLATNESS")
    print("status=VERIFIED-EXACT dual-field coefficient-side companion")
    print("current=J_s(R)=P(R+s)Q(R); fibres=det(s,R)=delta")
    print("equivalence=independent_nonzero_scalar_times_cyclic_origin_shift")
    print("projective_test=all_2_by_2_minors; no division")
    print("pair_basis=lex_s_then_0<=left_delta<right_delta<13_then_shift")
    print("frame_gauge=t_to_t+c*s gives profile_delta(w)_to_profile_delta(w+c*delta)")
    print("separator_operation=pointwise_product_on_high_sheet_function_algebra")
    print("separator_identity=Lambda_r((pi*a)^m)=0; Lambda_r(H*(pi*a)^m)=-a(r)^m for every m>=1")
    for result in results:
        print(
            f"p={result['prime']} profiles={result['profiles']} "
            f"pairs={result['profile_pairs']} "
            f"raw_matches={result['raw_matches']} "
            f"projective_matches={result['projective_matches']} "
            f"cyclic_projective_matches={result['cyclic_matches']} "
            f"dihedral_projective_matches={result['dihedral_matches']}"
        )
        print(
            f"p={result['prime']} cyclic_comparisons="
            f"{result['cyclic_comparisons']} "
            f"nonzero_minor_count={result['minor_min']}.."
            f"{result['minor_max']} "
            f"reflected_comparisons={result['dihedral_comparisons']} "
            f"reflected_nonzero_minor_count={result['reflected_minor_min']}.."
            f"{result['reflected_minor_max']}"
        )
        print(
            f"p={result['prime']} frame_entry_checks="
            f"{result['frame_entry_checks']} "
            f"nontrivial_modes={result['nontrivial_modes']} "
            f"zero_mode_support={result['zero_modes']}/{result['profiles']} "
            f"synthetic_translation_matches={result['synthetic_matches']} "
            f"synthetic_reflection_matches={result['reflection_controls']}"
        )
        coarse, marked, marked_nonzero = result["separator"]
        print(
            f"p={result['prime']} separator_coarse_zero={coarse}/{coarse} "
            f"separator_marked_formula={marked}/{marked} "
            f"separator_chord_nonzero={marked_nonzero}/52"
        )
        print(
            f"p={result['prime']} profile_sha256={result['profile_digest']} "
            f"mode_sha256={result['mode_digest']}"
        )
        print(
            f"p={result['prime']} minor_count_sha256="
            f"{result['minor_count_digest']} "
            f"witness_sha256={result['witness_digest']}"
        )
        print(
            f"p={result['prime']} reflected_minor_count_sha256="
            f"{result['reflected_minor_count_digest']} "
            f"reflected_witness_sha256={result['reflected_witness_digest']}"
        )
    print("scope=canonical_coefficient_current_only; no high-sheet physical descent, semantic ancestry, positive floor, row exclusion, or LRC14")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
