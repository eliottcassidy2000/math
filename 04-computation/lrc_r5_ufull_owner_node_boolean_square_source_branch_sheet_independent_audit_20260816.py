#!/usr/bin/env python3
"""Independent audit of the source-time P_169 high-branch sidecar.

The submitted sidecar implementation is not imported.  The audit rebuilds the
thirteen high-digit profiles directly from inverse images of E, applies Q and
P_(13^5), and couples them to the independently audited Boolean-square
endpoint integrand.

The retained digit b=floor(n/13) is a source-time P_169 branch sheet.  It is
not the current-leg r_owner=a mod 13, an arrival atom, a deep root, an exact
THM-2334 address, or a chronology.
"""

from __future__ import annotations

from bisect import bisect_right
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SQUARE_PATH = (
    ROOT
    / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_refiner_"
    "independent_audit_20260816.py"
)
SQUARE_SHA256 = "a75300a81efebef83683c41ac073ffa4d3268da83e96071d7b1b576b36e5bbc7"
SQUARE_SEMANTIC = "af0d543232869e82ee8d0191478ba7a833954cb19b387dedb6fb6f44a6fa272c"

P = 13
V = 4
ROWS = P * V
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (12, 12, 0))
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_NONEMPTY = (0,) + (1,) * 11 + (0,)
EXPECTED_PARENT_GAMMA = "8b8f2a2785b084e1578ba0512e4577ab79fd674b84b588fbdb8186f2009242c2"
EXPECTED_PARENT_TABLE = "5f4b9609faaa5f148d112a7cde5cfba0ab2c1385b4c53ea9c4bcfc6e93d106fc"
EXPECTED_FULL_SHAPE = (
    ("dc", 1), ("branch", 12), ("state", 3), ("residue", 12),
    ("branch_state", 36), ("branch_residue", 144),
    ("state_residue", 36), ("triple", 432),
)
EXPECTED_THREE_WAY_SHAPE = (
    ("dc", 0), ("branch", 0), ("state", 0), ("residue", 0),
    ("branch_state", 0), ("branch_residue", 0),
    ("state_residue", 0), ("triple", 432),
)
EXPECTED_DIGESTS = (
    "ba0143302119c93b7b7b4f778478c14c265efb8c0bcb9c2945e9eaae9009b493",
    "63d9304d28689ffa04b1568d538c67f3da59ecdc7e0d1c6535ad323c800ac0e5",
    "504e91431cf4d55f6e3cb5aa7c6e6f2c34db538cb692c67857a5c7627b252261",
    "61cf3999a1ef902db958851e99b8ea923c20099d81464c159613c0cf0054efdf",
    "8874b55a71e975e20c7b2b29efd9eec0ea93dcaecc0209d43c61712bbc9c0acf",
    "5d02353fd1cfa642be290264353e5e02441a2b69fa4144be72a7ae5ed9776396",
    "d8d388655fe34b27122c92da84c3cfb29d96303451334fbb4fe9d8c8fc13f2b8",
    "a70c2937fa755de450cc913d031b7e5f1eaf4b756900e79832db247582f76cba",
)
EXPECTED_SEMANTIC_SHA256 = "6cac7c90f3ebe3f33a27c1979b74006c545260bb6d41a68f4f76442456114fe6"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_square():
    require(lf_sha256(SQUARE_PATH) == SQUARE_SHA256, "square parent hash drift")
    spec = importlib.util.spec_from_file_location("branch_sheet_square", SQUARE_PATH)
    require(spec is not None and spec.loader is not None, "square parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.EXPECTED_SEMANTIC_SHA256 == SQUARE_SEMANTIC,
            "square parent semantic drift")
    return module


SQ = load_square()
R = SQ.R
PRIME = R.JOINT_PRIME


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def temporal_sheet_certificate():
    cells = 91 * P**2
    counts = [0] * P
    for cell in range(cells):
        midpoint = 2 * cell + 1
        source_guard_sheet = midpoint // (2 * P**2 * 7)
        inverse_branch = midpoint // (2 * 91)
        high_digit = inverse_branch // P
        require(0 <= inverse_branch < P**2, ("inverse branch", cell, inverse_branch))
        require(source_guard_sheet == high_digit,
                ("source-sheet certificate", cell, source_guard_sheet, high_digit))
        counts[high_digit] += 1
    require(tuple(counts) == (91 * P,) * P, ("sheet census", counts))
    return cells, tuple(counts)


def high_digit_packet_profiles(e_intervals, grid: int):
    """Build sum_{n=13b}^{13b+12} 1_E((y+n)/169) by direct inverse images."""
    require(grid % P**2 == 0, "source grid lacks P_169 cells")
    unit = grid // P**2
    groups = [[] for _branch in range(P)]
    input_mass = sum(right - left for left, right in e_intervals)
    lifted_mass = 0
    for inverse_branch in range(P**2):
        cell_left = inverse_branch * unit
        cell_right = cell_left + unit
        for left, right in e_intervals:
            overlap_left = max(left, cell_left)
            overlap_right = min(right, cell_right)
            if overlap_left >= overlap_right:
                continue
            mapped_left = P**2 * overlap_left - inverse_branch * grid
            mapped_right = P**2 * overlap_right - inverse_branch * grid
            require(0 <= mapped_left < mapped_right <= grid,
                    ("inverse-image interval", inverse_branch, mapped_left, mapped_right))
            groups[inverse_branch // P].append((mapped_left, mapped_right, 1))
            lifted_mass += mapped_right - mapped_left
    require(lifted_mass == P**2 * input_mass,
            ("inverse-image mass", lifted_mass, P**2 * input_mass))
    return tuple(R.fold_weighted(group, 1, grid) for group in groups)


def branch_profiles():
    grid = R.SRC.T_DEN
    e_intervals = R.SRC.build_set(R.SRC.PAT_E, R.SRC.ZELL)
    q_intervals = R.SRC.build_set(R.SRC.PAT_QB, R.SRC.ZELL)
    packet = high_digit_packet_profiles(e_intervals, grid)
    restricted = tuple(
        R.intersect_profile_with_set(profile, q_intervals, grid) for profile in packet
    )
    folded = tuple(R.fold_weighted(list(pieces), R.SRC.DCOLL, grid)
                   for pieces in restricted)
    roots = tuple(tuple(R.pull_profile_to_root(profile, root, grid)
                        for root in range(P)) for profile in folded)
    scale = R.JOINT_COORDINATE // grid
    scaled = tuple(tuple((tuple(position * scale for position in profile[0]), profile[1])
                         for profile in branch) for branch in roots)
    boundaries = tuple(sorted(
        {0, R.JOINT_COORDINATE}
        | {position for branch in scaled for profile in branch for position in profile[0]}
    ))
    nonempty = tuple(int(any(value for profile in branch for value in profile[1]))
                     for branch in scaled)
    require(nonempty == EXPECTED_NONEMPTY, ("nonempty branches", nonempty))
    return scaled, boundaries, digest_json(scaled), nonempty


def restore_profiles(branches, source_u):
    for root in range(P):
        points = sorted(
            {0, R.JOINT_COORDINATE}
            | set(source_u[root][0])
            | {point for branch in branches for point in branch[root][0]}
        )
        for point in points[:-1]:
            require(sum(profile_value(branch[root], point) for branch in branches)
                    == profile_value(source_u[root], point),
                    ("rootwise branch restoration", root, point))


def geometry(branches, boundaries, source_u):
    support = [[False] * V for _branch in range(P)]
    components = [[0] * V for _branch in range(P)]
    previous = [[False] * V for _branch in range(P)]
    merged = sorted(set(boundaries) | {point for profile in source_u for point in profile[0]})
    for left, right in zip(merged, merged[1:]):
        midpoint_twice = left + right
        owner = (
            7 * midpoint_twice < R.JOINT_COORDINATE
            or 7 * midpoint_twice > 13 * R.JOINT_COORDINATE
        )
        current = [[False] * V for _branch in range(P)]
        if owner:
            state = SQ.state_index(source_u, left)
            for branch in range(P):
                present = any(profile_value(branches[branch][root], left)
                              for root in range(P))
                if present:
                    support[branch][state] = True
                    current[branch][state] = True
                    if not previous[branch][state]:
                        components[branch][state] += 1
        previous = current
    matrix = tuple(tuple(int(value) for value in row) for row in support)
    component_matrix = tuple(tuple(row) for row in components)
    counts_state = tuple(sum(matrix[branch][state] for branch in range(P))
                         for state in range(V))
    counts_branch = tuple(sum(row) for row in matrix)
    expected_matrix = ((0, 0, 0, 0),) + ((1, 1, 1, 1),) * 11 + ((0, 0, 0, 0),)
    require(matrix == expected_matrix and component_matrix == expected_matrix,
            ("branch-state geometry", matrix, component_matrix))
    return matrix, component_matrix, counts_state, counts_branch


def integrate_pair(
    alpha, beta, word, endpoint_grid, branches, branch_boundaries,
    source_u, source_v, harmonic, danger, literal_tau=None,
):
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    rows = [[0] * ROWS for _tau in tau_values]
    positions = sorted(set(events) | set(branch_boundaries)
                       | {point for profile in source_u + source_v for point in profile[0]})
    mask = 0
    active = q_active = weighted = 0
    primitive_left = harmonic.value(positions[0])
    for left, right in zip(positions, positions[1:]):
        mask ^= events.get(left, 0)
        primitive_right = harmonic.value(right)
        jump = (primitive_right - primitive_left) % PRIME
        primitive_left = primitive_right
        if left == right or not mask:
            continue
        active += 1
        require(R.cell_index(left, right) == 0,
                ("branch sidecar escaped owner cell", alpha, beta, left, right))
        state = SQ.state_index(source_u, left)
        if not jump:
            continue
        q_active += 1
        whole_u = tuple(profile_value(profile, left) for profile in source_u)
        v_values = tuple(profile_value(profile, left) for profile in source_v)
        root_support = tuple(root for root, value in enumerate(whole_u) if value)
        branch_by_root = {
            root: tuple(profile_value(branches[branch][root], left) for branch in range(P))
            for root in root_support
        }
        require(all(sum(branch_by_root[root]) == whole_u[root] for root in root_support),
                ("local branch restoration", left, root_support))
        weighted += 1
        chamber = R.chamber(left, right)
        for row_index, tau in enumerate(tau_values):
            selected_mask = mask if literal_tau is not None else (
                mask & R.guard_mask(chamber, tau, danger)
            )
            selected = tuple(root for root in range(P) if (selected_mask >> root) & 1)
            selected_u = tuple(root for root in root_support if root in selected)
            right_value = sum(v_values[root] for root in selected)
            for root in selected_u:
                require(v_values[root] == 0, ("branch same-root", left, root))
            left_by_branch = tuple(sum(branch_by_root[root][branch] for root in selected_u)
                                   for branch in range(P))
            require(sum(left_by_branch) == sum(whole_u[root] for root in selected_u),
                    ("selected branch marginal", alpha, beta, tau, left))
            for branch, left_value in enumerate(left_by_branch):
                coordinate = branch * V + state
                rows[row_index][coordinate] = (
                    rows[row_index][coordinate] + left_value * right_value * jump
                ) % PRIME
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("branch endpoint mask", alpha, beta, literal_tau))
    return tuple(tuple(row) for row in rows), (
        interval_count, mapped_count, active, q_active, weighted
    )


def build_gamma(word, endpoint_grid, branches, boundaries, source_u, source_v,
                harmonic, danger):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    gamma = []
    totals = [0] * 5
    for alpha in range(P):
        for beta in range(P):
            rows, counts = integrate_pair(
                alpha, beta, word, endpoint_grid, branches, boundaries,
                source_u, source_v, harmonic, danger,
            )
            phase = pow(zeta, beta, PRIME)
            gamma.extend(tuple(phase * value % PRIME for value in row) for row in rows)
            totals = [left + right for left, right in zip(totals, counts)]
    return tuple(gamma), tuple(totals)


def inverse_table(gamma, zeta):
    table = [[0] * P for _coordinate in range(ROWS)]
    normalizer = pow(P**3, -1, PRIME)
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                for residue in range(P):
                    phase = pow(zeta, -(alpha + tau * residue) % P, PRIME)
                    for coordinate in range(ROWS):
                        table[coordinate][residue] = (
                            table[coordinate][residue] + row[coordinate] * phase
                        ) % PRIME
    return tuple(tuple(value * normalizer % PRIME for value in row) for row in table)


def triple_fourier(table, zeta):
    roots = tuple(tuple(pow(zeta, -frequency * value % P, PRIME)
                        for value in range(P)) for frequency in range(P))
    walsh = tuple(tuple(
        -1 if (character[0] * label[0] + character[1] * label[1]) % 2 else 1
        for label in SQ.STATE_LABELS
    ) for character in SQ.STATE_LABELS)
    return tuple(tuple(tuple(sum(
        table[branch * V + state][residue]
        * roots[branch_frequency][branch]
        * walsh[channel][state]
        * roots[residue_frequency][residue]
        for branch in range(P) for state in range(V) for residue in range(P)
    ) % PRIME for residue_frequency in range(P))
        for channel in range(V)) for branch_frequency in range(P))


def spectrum_shape(spectrum):
    names = (
        "dc", "branch", "state", "residue", "branch_state",
        "branch_residue", "state_residue", "triple",
    )
    counts = {name: 0 for name in names}
    for branch_frequency in range(P):
        for channel in range(V):
            for residue_frequency in range(P):
                if not spectrum[branch_frequency][channel][residue_frequency]:
                    continue
                flags = (branch_frequency != 0, channel != 0, residue_frequency != 0)
                name = {
                    (False, False, False): "dc",
                    (True, False, False): "branch",
                    (False, True, False): "state",
                    (False, False, True): "residue",
                    (True, True, False): "branch_state",
                    (True, False, True): "branch_residue",
                    (False, True, True): "state_residue",
                    (True, True, True): "triple",
                }[flags]
                counts[name] += 1
    return tuple((name, counts[name]) for name in names)


def three_way_interaction(table):
    tensor = [[[table[branch * V + state][residue] for residue in range(P)]
               for state in range(V)] for branch in range(P)]
    inv_p = pow(P, -1, PRIME)
    inv_v = pow(V, -1, PRIME)
    for state in range(V):
        for residue in range(P):
            mean = sum(tensor[branch][state][residue] for branch in range(P))
            mean = mean * inv_p % PRIME
            for branch in range(P):
                tensor[branch][state][residue] = (
                    tensor[branch][state][residue] - mean
                ) % PRIME
    for branch in range(P):
        for residue in range(P):
            mean = sum(tensor[branch][state][residue] for state in range(V))
            mean = mean * inv_v % PRIME
            for state in range(V):
                tensor[branch][state][residue] = (
                    tensor[branch][state][residue] - mean
                ) % PRIME
    for branch in range(P):
        for state in range(V):
            mean = sum(tensor[branch][state]) * inv_p % PRIME
            for residue in range(P):
                tensor[branch][state][residue] = (
                    tensor[branch][state][residue] - mean
                ) % PRIME
    return tuple(tuple(tuple(row) for row in branch) for branch in tensor)


def main() -> None:
    certificate = temporal_sheet_certificate()
    source_u, source_v, source_boundaries, source_digest, _total, _types = R.source_profiles()
    source_grid = R.SRC.T_DEN
    source_u = tuple(R.scale_profile(profile, R.JOINT_COORDINATE, source_grid)
                     for profile in source_u)
    source_v = tuple(R.scale_profile(profile, R.JOINT_COORDINATE, source_grid)
                     for profile in source_v)
    branches, branch_boundaries, profile_digest, nonempty = branch_profiles()
    restore_profiles(branches, source_u)
    geometry_record = geometry(branches, branch_boundaries, source_u)
    literal_boundaries = (
        set(branch_boundaries) | R.fixed_boundaries(source_boundaries, source_grid)
    )
    word, endpoint_grid = R.endpoint_word_and_grid()
    harmonic = R.HarmonicPrimitive(word, endpoint_grid)
    danger = R.danger_arcs()
    gamma, work_counts = build_gamma(
        word, endpoint_grid, branches, branch_boundaries,
        source_u, source_v, harmonic, danger,
    )
    require(work_counts == EXPECTED_WORK_COUNTS, ("branch work counts", work_counts))
    parent_gamma = tuple(tuple(sum(row[branch * V + state] for branch in range(P))
                                     % PRIME for state in range(V)) for row in gamma)
    require(digest_json(parent_gamma) == EXPECTED_PARENT_GAMMA,
            "branch gamma misses square")

    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    guard_records = []
    for alpha, beta, tau in CONTROL_TRIPLES:
        R.literal_guard_restoration(word, endpoint_grid, alpha, beta, tau,
                                    literal_boundaries, danger)
        direct, _counts = integrate_pair(
            alpha, beta, word, endpoint_grid, branches, branch_boundaries,
            source_u, source_v, harmonic, danger, literal_tau=tau,
        )
        phase = pow(zeta, beta, PRIME)
        index = (alpha * P + beta) * P + tau
        require(tuple(phase * value % PRIME for value in direct[0]) == gamma[index],
                ("literal branch guard", alpha, beta, tau))
        guard_records.append((alpha, beta, tau))

    table = inverse_table(gamma, zeta)
    parent_table = tuple(tuple(sum(table[branch * V + state][residue]
                                   for branch in range(P)) % PRIME
                               for residue in range(P)) for state in range(V))
    require(digest_json(parent_table) == EXPECTED_PARENT_TABLE,
            "branch table misses square")
    table_rank = R.rank_mod(table)
    spectrum = triple_fourier(table, zeta)
    shape = spectrum_shape(spectrum)
    three_way = three_way_interaction(table)
    three_way_flat = tuple(three_way[branch][state]
                           for branch in range(P) for state in range(V))
    three_way_rank = R.rank_mod(three_way_flat)
    three_way_spectrum = triple_fourier(three_way_flat, zeta)
    three_way_shape = spectrum_shape(three_way_spectrum)
    fixed = tuple(tuple(table[branch * V + state][6] for state in range(V))
                  for branch in range(P))
    branch_marginal = tuple(tuple(sum(table[branch * V + state][residue]
                                      for state in range(V)) % PRIME
                                  for residue in range(P)) for branch in range(P))
    fixed_record = (R.rank_mod(fixed), sum(value != 0 for row in fixed for value in row))
    rank_record = (R.rank_mod(parent_table), table_rank, R.rank_mod(branch_marginal))
    require(rank_record == (4, 6, 3), ("branch rank jump", rank_record))
    require(shape == EXPECTED_FULL_SHAPE, ("branch full shape", shape))
    require((three_way_rank, three_way_shape) == (6, EXPECTED_THREE_WAY_SHAPE),
            ("branch three-way", three_way_rank, three_way_shape))
    require(fixed_record == (3, 44), ("branch fixed relation", fixed_record))

    digests = (
        profile_digest, digest_json(gamma), digest_json(table), digest_json(spectrum),
        digest_json(three_way), digest_json(three_way_spectrum), digest_json(fixed),
        digest_json(branch_marginal),
    )
    require(digests == EXPECTED_DIGESTS, ("branch candidate digests", digests))
    record = (
        source_digest, certificate, nonempty, geometry_record, tuple(guard_records),
        work_counts, rank_record, shape, three_way_rank, three_way_shape,
        fixed_record, digests,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("branch audit semantic", semantic))

    print("== independent hostile audit: source-time P169 high-branch x Boolean square x relation ==")
    print(f"parent=(square_sha256={SQUARE_SHA256},square_semantic={SQUARE_SEMANTIC})")
    print(f"sidecar=b=floor(n/13);temporal_source_sheet_certificate={certificate}: PASS")
    print(f"branch_profile=(sha256={profile_digest},nonempty={nonempty},rootwise_restoration=True)")
    print(f"geometric_branch_x_state_support={geometry_record}")
    print(f"source_profile_sha256={source_digest};work_counts={work_counts};literal_controls={tuple(guard_records)}: PASS")
    print("marginals=(branch->audited_square gamma and table): PASS;same_root_pointwise_zero=PASS")
    print(f"rank_jump_(square,branch_x_state,branch_marginal)={rank_record}")
    print(f"full_spectral_support={shape}")
    print(f"three_way_interaction=(rank={three_way_rank},support={three_way_shape})")
    print(f"fixed_(1,0,6)=(branch_x_state_rank={fixed_record[0]},nonzero_entries={fixed_record[1]}/52,digest={digests[6]})")
    print(f"digests_(profiles,gamma,table,spectrum,three_way,three_way_spectrum,fixed,branch_marginal)={digests}")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT scoped nonredundant source-time inverse-branch sidecar")
    print("typing=source P169 high digit;not current-leg r_owner,not arrival atom,not exact address,not chronology")
    print("scope=one owner base;root difference and pointed tail marginalized;no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
