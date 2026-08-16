#!/usr/bin/env python3
"""Exact folded-C7/K4 drift-transporter probe.

This companion compares two already typed finite objects:

* the raw 7 x 13 common-base response table in proved THM-2594; and
* the four K4-Walsh drift rows in proved THM-3514.

The comparison is representation-theoretic.  It does not construct the
missing U_full common-ancestry relation.  In particular, a channel-dependent
circulant transporter is only a marked split-field map, not a physical
current or a character-independent Boolean transplant.
"""

from __future__ import annotations

import ast
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_r5_folded_c7_to_ufull_k4_drift_transporter_probe_20260816.py"
UFULL_PATH = ROOT / "04-computation/lrc_ufull_guard_bucket_all_role_spectral_probe_20260816.py"
R5_OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_stage2_theta_contraction_opus_20260728.out"
UFULL_SHA256 = "98a4cf5c82ca10027302baf2c7fb59acb0f305143e22453d0a8660fef8d90cf0"
R5_OUTPUT_SHA256 = "bef4ee9a18ff3e2f455bad66a95252dd9989b2f60953e26e8ea0c2dc6ae7f5df"
EXPECTED_GAMMA_SHA256 = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_SEMANTIC_SHA256 = "bae3749e01361ebbaf6c9cc0e7160d2e4d22d6df420d4512ba94f778718af3a1"
P = 13


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    body = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    require(lf_sha256(path) == expected_hash, (name, lf_sha256(path), expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


U = load_module(UFULL_PATH, "lrc_ufull_all_role_parent", UFULL_SHA256)


def worker(alpha: int) -> tuple[object, ...]:
    return U.worker(alpha)


def parse_r5_table() -> tuple[tuple[tuple[int, ...], ...], int]:
    require(lf_sha256(R5_OUTPUT_PATH) == R5_OUTPUT_SHA256, "r5 output hash drift")
    text = R5_OUTPUT_PATH.read_text(encoding="utf-8")
    match = re.search(r"entry denominators: DENC = .*? = ([0-9]+)", text)
    require(match is not None, "missing DENC")
    denominator = int(match.group(1))
    section = text.split("[8] the response array", 1)[1].split(
        "fibrewise ANOVA", 1
    )[0]
    rows: dict[int, tuple[int, ...]] = {}
    pattern = re.compile(r"ell=(\d+): theta=0\.\.2: (\[[^\]]+\])")
    for ell_text, values_text in pattern.findall(section):
        values = ast.literal_eval(values_text)
        require(isinstance(values, list) and len(values) == 3, values)
        rows[int(ell_text)] = tuple(int(value) for value in values) + (0,) * 10
    require(tuple(sorted(rows)) == tuple(range(7)), sorted(rows))
    return tuple(rows[index] for index in range(7)), denominator


def interaction_numerators(table: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    row_sums = tuple(sum(row) for row in table)
    column_sums = tuple(sum(table[ell][theta] for ell in range(7)) for theta in range(P))
    total = sum(row_sums)
    return tuple(
        tuple(
            91 * table[ell][theta]
            - 7 * row_sums[ell]
            - 13 * column_sums[theta]
            + total
            for theta in range(P)
        )
        for ell in range(7)
    )


def dft(values: tuple[int, ...], frequency: int, root: int, prime: int) -> int:
    return sum(
        value * pow(root, -frequency * index % P, prime)
        for index, value in enumerate(values)
    ) % prime


def spectrum(values: tuple[int, ...], root: int, prime: int) -> tuple[int, ...]:
    return tuple(dft(values, frequency, root, prime) for frequency in range(P))


def inverse_spectrum(values: tuple[int, ...], root: int, prime: int) -> tuple[int, ...]:
    inverse = pow(P, -1, prime)
    return tuple(
        inverse
        * sum(
            values[frequency] * pow(root, frequency * index % P, prime)
            for frequency in range(P)
        )
        % prime
        for index in range(P)
    )


def convolution(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> tuple[int, ...]:
    return tuple(
        sum(left[shift] * right[(index - shift) % P] for shift in range(P)) % prime
        for index in range(P)
    )


def centered(values: tuple[int, ...], prime: int) -> tuple[int, ...]:
    mean = sum(values) * pow(P, -1, prime) % prime
    return tuple((value - mean) % prime for value in values)


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(row) for row in matrix]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column] % prime),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column] % prime, -1, prime)
        rows[rank] = [entry * inverse % prime for entry in rows[rank]]
        for row in range(len(rows)):
            if row == rank:
                continue
            factor = rows[row][column] % prime
            if factor:
                rows[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(rows[row], rows[rank])
                ]
        rank += 1
    return rank


def ratio_equivalent(
    left: tuple[int, ...], right: tuple[int, ...], root: int, prime: int
) -> bool:
    """Equal up to a nonzero scalar and a cyclic gauge shift."""
    require(all(left) and all(right), "ratio equivalence needs units")
    scalar = left[0] * pow(right[0], -1, prime) % prime
    for shift in range(P):
        if all(
            left[frequency]
            == scalar
            * pow(root, frequency * shift % P, prime)
            % prime
            * right[frequency]
            % prime
            for frequency in range(P)
        ):
            return True
    return False


def equivalence_class_count(
    rows: tuple[tuple[int, ...], ...], root: int, prime: int
) -> int:
    representatives: list[tuple[int, ...]] = []
    for row in rows:
        if not any(ratio_equivalent(row, representative, root, prime)
                   for representative in representatives):
            representatives.append(row)
    return len(representatives)


def projective_channel_mixing_rank(
    source_rows: tuple[tuple[int, ...], ...],
    target_rows: tuple[tuple[int, ...], ...],
    prime: int,
) -> int:
    """Rank of the linear system Y_k proportional to M X_k.

    The sixteen unknowns are the entries of one frequency-independent 4 x 4
    channel matrix M.  A nonzero solution leaves a scalar multiplier at each
    frequency, which is exactly one common circulant convolution in Fourier
    coordinates.
    """
    require(len(source_rows) == len(target_rows) == 4, "four channels required")
    equations = []
    for frequency in range(P):
        source = tuple(row[frequency] for row in source_rows)
        target = tuple(row[frequency] for row in target_rows)
        for left in range(4):
            for right in range(left + 1, 4):
                equation = [0] * 16
                for column in range(4):
                    equation[right * 4 + column] = (
                        equation[right * 4 + column]
                        + target[left] * source[column]
                    ) % prime
                    equation[left * 4 + column] = (
                        equation[left * 4 + column]
                        - target[right] * source[column]
                    ) % prime
                equations.append(tuple(equation))
    return rank_mod(tuple(equations), prime)


def paley_fold() -> tuple[object, ...]:
    residues = {1, 2, 4}
    arcs = tuple(
        (left, right)
        for left in range(7)
        for right in range(7)
        if left != right and (right - left) % 7 in residues
    )
    require(len(arcs) == 21, len(arcs))
    orbits = ((0,), (1, 6), (2, 5), (3, 4))
    counts = tuple(
        tuple(
            sum((left, right) in arcs for left in source for right in target)
            for target in orbits
        )
        for source in orbits
    )
    require(
        all(counts[i][j] == counts[j][i] for i in range(4) for j in range(4)),
        counts,
    )
    doubled = tuple(
        next(
            index
            for index, target_orbit in enumerate(orbits)
            if (2 * source_orbit[0]) % 7 in target_orbit
        )
        for source_orbit in orbits[1:]
    )
    require(doubled == (2, 3, 1), doubled)

    vertices = ((0, 0), (0, 1), (1, 0), (1, 1))
    characters = ((0, 0), (1, 0), (0, 1), (1, 1))
    matchings = []
    for character in characters[1:]:
        matching = tuple(
            sorted(
                tuple(sorted((i, j)))
                for bit in (0, 1)
                for i in range(4)
                for j in range(i + 1, 4)
                if (
                    (character[0] * vertices[i][0] + character[1] * vertices[i][1]) % 2
                    == bit
                    == (character[0] * vertices[j][0] + character[1] * vertices[j][1]) % 2
                )
            )
        )
        require(len(matching) == 2, matching)
        matchings.append(matching)
    require(len(set(matchings)) == 3, matchings)
    require(
        sorted(edge for matching in matchings for edge in matching)
        == [(i, j) for i in range(4) for j in range(i + 1, 4)],
        matchings,
    )

    # x^3+x^2-2x-1 is the minimal polynomial of zeta_7+zeta_7^-1.
    a, b, c, d = 1, 1, -2, -1
    discriminant = b * b * c * c - 4 * a * c**3 - 4 * b**3 * d - 27 * a * a * d * d + 18 * a * b * c * d
    require(discriminant == 49, discriminant)
    require(all(value**3 + value**2 - 2 * value - 1 for value in (-1, 1)),
            "real septic cubic unexpectedly reducible")
    return arcs, orbits, counts, doubled, tuple(matchings), discriminant


def main() -> None:
    raw_integer, denominator = parse_r5_table()
    interaction_integer = interaction_numerators(raw_integer)

    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(value for chunk in chunks for value in chunk[5])
    require(U.digest_integers(gamma) == EXPECTED_GAMMA_SHA256, "gamma digest")

    (
        _word,
        _t_den,
        nn,
        prime,
        root,
        zeta13,
        *_rest,
    ) = U.B.context()
    require(nn % 7 == 0, nn)
    zeta7 = pow(root, nn // 7, prime)
    require(pow(zeta7, 7, prime) == 1 and zeta7 != 1, "zeta7 order")
    require(pow(zeta13, P, prime) == 1 and zeta13 != 1, "zeta13 order")
    require(denominator % prime and (91 * denominator) % prime, "bad denominator")

    normalizer = pow(P**3, -1, prime)
    q_buckets = {
        q: tuple(
            sum(chunk[6][q_index][bucket] for chunk in chunks)
            % prime
            * normalizer
            % prime
            for bucket in range(len(U.B.BUCKETS))
        )
        for q_index, q in enumerate(U.Q_CLASSES)
    }
    q_h = (1, 0, 1)
    q_q5 = (1, 0, 0)
    corners = tuple(
        tuple(
            (
                q_buckets[q_h][U.B.BUCKET_INDEX[(left, right, drift)]]
                - q_buckets[q_q5][U.B.BUCKET_INDEX[(left, right, drift)]]
            )
            % prime
            for drift in range(P)
        )
        for left, right in U.CORNER_PAIRS
    )
    walsh = {
        "constant": tuple(sum(corners[row][d] for row in range(4)) % prime for d in range(P)),
        "left": tuple((corners[0][d] + corners[1][d] - corners[2][d] - corners[3][d]) % prime for d in range(P)),
        "right": tuple((corners[0][d] - corners[1][d] + corners[2][d] - corners[3][d]) % prime for d in range(P)),
        "mixed": tuple((corners[0][d] - corners[1][d] - corners[2][d] + corners[3][d]) % prime for d in range(P)),
    }
    names = ("constant", "left", "right", "mixed")
    target_spectra = {name: spectrum(walsh[name], zeta13, prime) for name in names}
    require(all(all(walsh[name]) and all(target_spectra[name]) for name in names),
            "U_full Walsh support regression")
    require(sum(walsh["constant"]) % prime == 389266878372286537904,
            "global bridge regression")
    require(rank_mod(tuple(walsh[name] for name in names), prime) == 4,
            "Walsh rank regression")

    inverse_denominator = pow(denominator, -1, prime)
    inverse_interaction_denominator = pow(91 * denominator, -1, prime)
    raw = tuple(
        tuple(value % prime * inverse_denominator % prime for value in row)
        for row in raw_integer
    )
    interaction = tuple(
        tuple(value % prime * inverse_interaction_denominator % prime for value in row)
        for row in interaction_integer
    )
    source_profiles = {}
    interaction_profiles = {}
    source_spectra = {}
    for septimal in range(7):
        source_profiles[septimal] = tuple(
            sum(
                raw[ell][theta] * pow(zeta7, -septimal * ell % 7, prime)
                for ell in range(7)
            )
            % prime
            for theta in range(P)
        )
        interaction_profiles[septimal] = tuple(
            sum(
                interaction[ell][theta] * pow(zeta7, -septimal * ell % 7, prime)
                for ell in range(7)
            )
            % prime
            for theta in range(P)
        )
        source_spectra[septimal] = spectrum(source_profiles[septimal], zeta13, prime)
    require(all(all(source_spectra[septimal]) for septimal in range(7)),
            "a raw septimal profile is not cyclic")
    full_source_rank = rank_mod(
        tuple(source_profiles[septimal] for septimal in range(7)), prime
    )
    require(full_source_rank == 3, ("three-window source rank", full_source_rank))
    require(all(value == 0 for value in interaction_profiles[0]),
            "double centring retained the trivial septimal mode")
    for septimal in range(1, 7):
        interaction_hat = spectrum(interaction_profiles[septimal], zeta13, prime)
        require(interaction_hat[0] == 0 and all(interaction_hat[1:]),
                ("augmentation spectrum", septimal))
        require(centered(source_profiles[septimal], prime) == interaction_profiles[septimal],
                ("raw/interaction margin split", septimal))

    folded_pairs = ((1, 6), (2, 5), (3, 4))
    nontrivial_names = names[1:]

    def allocation_rows(allocation: dict[str, int]) -> tuple[tuple[int, ...], ...]:
        return tuple(
            tuple(
                target_spectra[name][frequency]
                * pow(source_spectra[allocation[name]][frequency], -1, prime)
                % prime
                for frequency in range(P)
            )
            for name in names
        )

    trivial_preserving = []
    for pair_order in itertools.permutations(folded_pairs):
        for signs in itertools.product((0, 1), repeat=3):
            allocation = {"constant": 0}
            allocation.update(
                {
                    name: pair_order[index][signs[index]]
                    for index, name in enumerate(nontrivial_names)
                }
            )
            trivial_preserving.append(allocation)
    require(len(trivial_preserving) == 48, len(trivial_preserving))

    all_allocations = []
    folded_classes = ((0,),) + folded_pairs
    for class_order in itertools.permutations(folded_classes):
        nonzero_positions = tuple(index for index, cls in enumerate(class_order) if len(cls) == 2)
        for signs in itertools.product((0, 1), repeat=3):
            selected = []
            sign_index = 0
            for cls in class_order:
                if len(cls) == 1:
                    selected.append(0)
                else:
                    selected.append(cls[signs[sign_index]])
                    sign_index += 1
            require(len(nonzero_positions) == 3, nonzero_positions)
            all_allocations.append(dict(zip(names, selected)))
    require(len(all_allocations) == 192, len(all_allocations))

    def allocation_census(allocations: list[dict[str, int]]) -> tuple[int, int, tuple[int, ...]]:
        exact_common = 0
        gauge_common = 0
        class_counts = []
        for allocation in allocations:
            rows = allocation_rows(allocation)
            if len(set(rows)) == 1:
                exact_common += 1
            classes = equivalence_class_count(rows, zeta13, prime)
            class_counts.append(classes)
            if classes == 1:
                gauge_common += 1
        return exact_common, gauge_common, tuple(sorted(set(class_counts)))

    preserving_census = allocation_census(trivial_preserving)
    arbitrary_census = allocation_census(all_allocations)
    unrestricted_mode_allocations = [
        dict(zip(names, selected))
        for selected in itertools.product(range(7), repeat=4)
    ]
    require(len(unrestricted_mode_allocations) == 7**4, len(unrestricted_mode_allocations))
    unrestricted_census = allocation_census(unrestricted_mode_allocations)

    target_rows = tuple(target_spectra[name] for name in names)
    projective_ranks = []
    selected_source_ranks = []
    for signs in itertools.product((0, 1), repeat=3):
        selected = (0,) + tuple(folded_pairs[index][signs[index]] for index in range(3))
        source_rows = tuple(source_spectra[index] for index in selected)
        selected_source_ranks.append(rank_mod(source_rows, prime))
        projective_ranks.append(projective_channel_mixing_rank(source_rows, target_rows, prime))
    require(len(projective_ranks) == 8, projective_ranks)
    require(tuple(selected_source_ranks) == (3,) * 8, selected_source_ranks)
    require(tuple(projective_ranks) == (12,) * 8, projective_ranks)
    # Each 4-channel source has left nullity one, so the space of 4 x 4
    # matrices annihilating it has dimension 4.  The projective system also
    # has nullity 16-12=4.  Hence every formal solution is an annihilator;
    # none sends a source frequency vector to the nonzero target vector.
    projective_solution_space = "annihilator-only"

    calibration = {"constant": 0, "left": 1, "right": 2, "mixed": 3}
    ratio_rows = allocation_rows(calibration)
    full_kernels = tuple(inverse_spectrum(row, zeta13, prime) for row in ratio_rows)
    for name, septimal, kernel in zip(names, (0, 1, 2, 3), full_kernels):
        require(
            convolution(kernel, source_profiles[septimal], prime) == walsh[name],
            ("full transporter", name),
        )

    source_augmentations = tuple(centered(source_profiles[index], prime) for index in (0, 1, 2, 3))
    target_augmentations = tuple(centered(walsh[name], prime) for name in names)
    augmentation_kernels = []
    margin_multipliers = []
    for source_vector, target_vector, full_ratio in zip(
        source_augmentations, target_augmentations, ratio_rows
    ):
        source_hat = spectrum(source_vector, zeta13, prime)
        target_hat = spectrum(target_vector, zeta13, prime)
        require(source_hat[0] == target_hat[0] == 0, "centred zero mode")
        require(all(source_hat[1:]) and all(target_hat[1:]), "augmentation cyclicity")
        multiplier = (0,) + tuple(
            target_hat[frequency] * pow(source_hat[frequency], -1, prime) % prime
            for frequency in range(1, P)
        )
        require(multiplier[1:] == full_ratio[1:], "mixed modes changed under centring")
        kernel = inverse_spectrum(multiplier, zeta13, prime)
        require(convolution(kernel, source_vector, prime) == target_vector,
                "augmentation transporter")
        augmentation_kernels.append(kernel)
        margin_multipliers.append(full_ratio[0])
        difference = tuple((a - b) % prime for a, b in zip(
            inverse_spectrum(full_ratio, zeta13, prime), kernel
        ))
        require(len(set(difference)) == 1, "margin correction is not constant")

    paley = paley_fold()
    kernel_supports = tuple(sum(bool(value) for value in row) for row in full_kernels)
    augmentation_supports = tuple(sum(bool(value) for value in row) for row in augmentation_kernels)
    semantic = (
        UFULL_SHA256,
        R5_OUTPUT_SHA256,
        (prime, root, nn, zeta7, zeta13),
        raw_integer,
        denominator,
        tuple(source_profiles[index] for index in range(7)),
        tuple(source_spectra[index] for index in range(7)),
        tuple((name, walsh[name], target_spectra[name]) for name in names),
        preserving_census,
        arbitrary_census,
        unrestricted_census,
        tuple(projective_ranks),
        tuple(selected_source_ranks),
        projective_solution_space,
        calibration,
        full_kernels,
        tuple(augmentation_kernels),
        tuple(margin_multipliers),
        paley,
        "split-field channel transport only; no common ancestry/current/LRC14",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("R5 FOLDED-C7 TO U_FULL K4xF13 DRIFT TRANSPORTER PROBE")
    print("status=FINITE-EXACT marked split-field comparison; LRC(14) OPEN")
    print(f"dependencies=(ufull={UFULL_SHA256},r5_output={R5_OUTPUT_SHA256})")
    print(f"field=(prime={prime},zeta7={zeta7},zeta13={zeta13})")
    print("source_raw_spectrum=(7 septimal profiles x 13 root frequencies; nonzero=91/91)")
    print("source_interaction_spectrum=(6 primitive septimal profiles x 12 primitive root frequencies; nonzero=72/72; double centring deletes the 13-mode septimal-trivial sector and six remaining root-zero margins)")
    print("target_spectrum=(4 Walsh channels x 13 drift frequencies; nonzero=52/52; rank=4)")
    print(f"paley_negation_fold=(orbits={paley[1]},arc_counts={paley[2]},times2_cycle={paley[3]})")
    print(f"k4_nontrivial_walsh_matchings={paley[4]}")
    print("rational_algebra_boundary=(Q[C7]^inversion = Q plus cyclic cubic x^3+x^2-2x-1 of discriminant 49; Q[F2^2]=Q^4; no Q-algebra isomorphism; both split after marked base extension)")
    print(f"allocation_census=(trivial_preserving_48={preserving_census},all_folded_192={arbitrary_census},all_mode_choices_2401={unrestricted_census})")
    print("allocation_fields=(exact_common_kernel,gauge/amplitude_common_kernel,observed_numbers_of_kernel_classes)")
    print(f"fixed_channel_mixing_system=(all_7_source_profiles_rank={full_source_rank},selected_4_source_ranks={tuple(selected_source_ranks)},target_rank=4,projective_equation_ranks={tuple(projective_ranks)},solution_space={projective_solution_space})")
    print(f"calibration={calibration}; full_kernel_supports={kernel_supports}; full_kernel_rank={rank_mod(full_kernels, prime)}")
    print(f"full_kernel_sha256={digest(full_kernels)}")
    print(f"augmentation_kernel_supports={augmentation_supports}; augmentation_kernel_rank={rank_mod(tuple(augmentation_kernels), prime)}")
    print(f"augmentation_kernel_sha256={digest(tuple(augmentation_kernels))}")
    print(f"margin_multipliers={tuple(margin_multipliers)}")
    print("margin_split=for each selected raw C13 profile, full minus augmentation transport is a constant kernel; for septimal modes 1,2,3 the centred source is the THM-2594 interaction, while the septimal-trivial profile is wholly absent after double centring")
    print("hostile=no folded-mode allocation yields one scalar drift convolution after independent channel amplitude/cyclic-gauge calibration; all fixed channel-mixing projective solutions annihilate the rank-three, three-window source, so one common convolution cannot reach the rank-four Walsh target")
    print("survivor=four marked channel-dependent invertible circulants transport the raw r5 profiles to the four U_full Walsh rows; this is formal spectral equivalence only")
    print("next_test=a shared linear drift transplant needs a fourth independent source channel (new deep window/sidecar); a lawful relation could instead be nonlinear, nonconvolutional, channel-dependent, or frequency-dependent")
    print(f"semantic_sha256={semantic_hash}")
    print("nonconsequence=no Boolean coupling, physical current, H1 flux, grouped coefficient, row exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
