#!/usr/bin/env python3
"""Refine the Boolean square to its six realized pointed cut states.

The four owner-visible support states contain exactly six lawful pairs
``(state,u)`` with ``u`` an active THM-2471 source root.  In geometric order,

    (0,0), (1,0), (1,6), (3,6), (3,12), (2,12).

This is a pointed six-state path, not the abstract six-cut completion and not
a six-vertex tournament.  The script retains this marked tail together with
the source root difference ``s=u-q`` before integrating the actual U_full
endpoint factors.  Summing marked tails must recover the pinned
Boolean-square x root-difference tensor exactly.

There is no canonical cyclic group on the six pointed states.  Consequently
the decisive invariant is tensor flattening rank after three-axis centering.
For a secondary coordinate-aware test, the five geometrically ordered path
edge contrasts are Fourier transformed only in the two genuine F_13 axes.

This is a source-root ancestry refinement.  It is not a THM-2334 exact
address, an inverse-ancestry sheet, a U_clock chronology, or a physical
current.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = (
    ROOT
    / "04-computation/lrc_r5_ufull_owner_node_boolean_square_root_difference_probe_20260816.py"
)
PARENT_SHA256 = "dddeea995e9ab7e8abfd010e00c798c93a28ff924da8172dd114756184942a57"

EXPECTED_DIGESTS = (
    (
        "416f03a90894e526bc767a34839c2489aa4cac7051ac04687e0feacfa36d58d2",
        "a24603fe10909d5fec3f4a3668ae4a2eaffb06ba36c7d6a7121665fc8a996c38",
    ),
    (
        "9c5e227d9c142373973a562a54c6a67cac60a82da1121a028b9658920d155a19",
        "4231a2033312b13ebc3070c5084f55dff421186120831142dd8a57722f0d0c7b",
        "4958e3b37ff473644a610671763404d5aa2fe850c4e57de3c2c94eb0a89bb184",
    ),
    (
        "ae9744c987a3b2b6824feab0b64b234a239f7efa250193304e9aa5beb18e6fb5",
        "34356f952b40ee3f8ed423d732114c2bb2c1107dddf89bbe5eb3d956215d753f",
        "154bae6efd632d281fa96acafd386cde4d75ae38f460e706d02f0746371e3310",
    ),
    (
        "9d3733c61ce0dacc056d382db70ad0a10ad5aff62a8c57c3968229e21fdfad58",
        "4d9857743743203f86380b3a162c0e4314de119b95b3cbd137d0e3f180993453",
        "e59b7082fe2e25718b01c07d7d29a843d4e695075a23b85e83da960d70a8bd8d",
    ),
)
EXPECTED_RANKS = (
    ((6, 12, 13), (6, 12, 13), (4, 12, 13)),
    ((5, 12, 12), (5, 12, 12), (3, 12, 12)),
)
EXPECTED_SUPPORT_SHA256 = "4e7bd717552c0976344f55ee878eb6ae89d5c74b9584214007942f734aed94af"
EXPECTED_FIXED = (
    (6, 5, (12, 12),
     "37daad2d2da203cc2ac57342765b379096aedcbf17b903bd286210a34cdfa356"),
    (6, 5, (12, 12),
     "43c2035d3635fa8309aa870a6f8e0e42e625da32e2573091879737a758fe1308"),
    (4, 3, (0, 0),
     "5cbc6f2b579968bf49df54873fc7ad55436c9137fd7bc9f5a703f5bcd8850310"),
)
EXPECTED_SEMANTIC_SHA256 = "5001bed534a9b8f953529101f0b7a51cf6994c7dcb29b7aed576aec239078384"

P = 13
W = 6
POINTED_STATES = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))
PATH_EDGE_TYPES = (
    "enter_doubleton_at_0",
    "retail_0_to_6",
    "owner_gap_bridge_at_6",
    "retail_6_to_12",
    "leave_doubleton_at_12",
)
POINT_INDEX = {pair: index for index, pair in enumerate(POINTED_STATES)}
STATE_POINTS = tuple(
    tuple(index for index, (point_state, _root) in enumerate(POINTED_STATES)
          if point_state == state)
    for state in range(4)
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_parent():
    observed = lf_sha256(PARENT_PATH)
    require(observed == PARENT_SHA256,
            ("root-difference parent drift", observed, PARENT_SHA256))
    spec = importlib.util.spec_from_file_location("pointed_six_parent", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "pointed-six parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load_parent()
B = R.B
C = R.C
PRIME = C.JOINT_PRIME


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def integrate_pointed(alpha: int, beta: int, literal_tau: int | None = None):
    ctx = C.context()
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    weighted = [
        [[0 for _difference in range(P)] for _point in range(W)]
        for _tau in tau_values
    ]
    support_only = [
        [[0 for _difference in range(P)] for _point in range(W)]
        for _tau in tau_values
    ]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    point_segments = [0] * W
    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        require(C.cell_of_segment(left, right) == 0,
                ("pointed owner escaped cell zero", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("pointed active middle chamber", alpha, beta, left, right))
        state = B.state_of_segment(left, right)
        active_segments += 1
        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])
        u_support = tuple(root for root, value in enumerate(u_values) if value)
        v_support = tuple(root for root, value in enumerate(v_values) if value)
        if not u_support or not v_support:
            continue
        weighted_segments += 1
        require(set(u_support).isdisjoint(v_support),
                ("pointed source supports overlap", left, right, u_support, v_support))
        for u in u_support:
            require((state, u) in POINT_INDEX,
                    ("unrealized pointed state", left, right, state, u, u_support))
            point_segments[POINT_INDEX[(state, u)]] += 1
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet
                    for sheet in range(P)
                    if (mask >> sheet) & 1 and C.T.safe(chamber, sheet + tau)
                )
            else:
                selected = tuple(sheet for sheet in range(P) if (mask >> sheet) & 1)
            selected_set = set(selected)
            selected_u = tuple(root for root in u_support if root in selected_set)
            selected_v = tuple(root for root in v_support if root in selected_set)
            for u in selected_u:
                point = POINT_INDEX[(state, u)]
                for q in selected_v:
                    difference = (u - q) % P
                    require(difference != 0, ("pointed same-root pair", left, right, u, q))
                    weighted[row_index][point][difference] = (
                        weighted[row_index][point][difference]
                        + u_values[u] * v_values[q] * jump
                    ) % PRIME
                    support_only[row_index][point][difference] = (
                        support_only[row_index][point][difference] + jump
                    ) % PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("pointed endpoint mask", alpha, beta, literal_tau, mask))
    counts = (
        interval_count,
        mapped,
        active_segments,
        q_active_segments,
        weighted_segments,
        tuple(point_segments),
    )
    freeze = lambda bank: tuple(
        tuple(tuple(row) for row in point_rows) for point_rows in bank
    )
    return freeze(weighted), freeze(support_only), counts


def worker(alpha: int):
    zeta = C.context()["zeta"]
    weighted_rows = []
    support_rows = []
    scalar_counts = [0] * 5
    point_counts = [0] * W
    for beta in range(P):
        weighted, support_only, record = integrate_pointed(alpha, beta)
        phase = pow(zeta, beta, PRIME)

        def twist(bank):
            return tuple(
                tuple(
                    tuple(phase * value % PRIME for value in difference_row)
                    for difference_row in point_rows
                )
                for point_rows in bank
            )

        weighted_rows.append(twist(weighted))
        support_rows.append(twist(support_only))
        scalar_counts = [left + right for left, right in zip(scalar_counts, record[:5])]
        point_counts = [left + right for left, right in zip(point_counts, record[5])]
    return alpha, tuple(weighted_rows), tuple(support_rows), (tuple(scalar_counts), tuple(point_counts))


def marginal_gamma(gamma_rows):
    return tuple(
        tuple(
            tuple(
                sum(row[point][difference] for point in STATE_POINTS[state]) % PRIME
                for difference in range(P)
            )
            for state in range(4)
        )
        for row in gamma_rows
    )


def inverse_tensor(gamma_rows, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [[[0 for _relation in range(P)] for _difference in range(P)] for _point in range(W)]
    phases = tuple(
        tuple(pow(zeta, -(alpha + tau * relation) % P, PRIME) for relation in range(P))
        for alpha in range(P) for tau in range(P)
    )
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_rows[index]
                phase_row = phases[alpha * P + tau]
                index += 1
                for point in range(W):
                    for difference in range(P):
                        value = row[point][difference]
                        if not value:
                            continue
                        for relation, phase in enumerate(phase_row):
                            tensor[point][difference][relation] = (
                                tensor[point][difference][relation] + value * phase
                            ) % PRIME
    require(index == P**3, ("pointed gamma size", index))
    return tuple(
        tuple(
            tuple(value * normalizer % PRIME for value in relation_row)
            for relation_row in difference_rows
        )
        for difference_rows in tensor
    )


def marginal_tensor(tensor):
    return tuple(
        tuple(
            tuple(
                sum(tensor[point][difference][relation] for point in STATE_POINTS[state]) % PRIME
                for relation in range(P)
            )
            for difference in range(P)
        )
        for state in range(4)
    )


def flat_tail_lift(parent_tensor):
    return tuple(
        tuple(
            tuple(
                parent_tensor[state][difference][relation]
                * pow(len(STATE_POINTS[state]), -1, PRIME) % PRIME
                for relation in range(P)
            )
            for difference in range(P)
        )
        for state, _root in POINTED_STATES
    )


def centre_axis(tensor, axis: int):
    dimensions = (W, P, P)
    inverse = pow(dimensions[axis], -1, PRIME)
    mutable = [[[tensor[point][difference][relation] for relation in range(P)]
                for difference in range(P)] for point in range(W)]
    if axis == 0:
        for difference in range(P):
            for relation in range(P):
                mean = sum(mutable[point][difference][relation] for point in range(W)) * inverse % PRIME
                for point in range(W):
                    mutable[point][difference][relation] = (mutable[point][difference][relation] - mean) % PRIME
    elif axis == 1:
        for point in range(W):
            for relation in range(P):
                mean = sum(mutable[point][difference][relation] for difference in range(P)) * inverse % PRIME
                for difference in range(P):
                    mutable[point][difference][relation] = (mutable[point][difference][relation] - mean) % PRIME
    else:
        for point in range(W):
            for difference in range(P):
                mean = sum(mutable[point][difference][relation] for relation in range(P)) * inverse % PRIME
                for relation in range(P):
                    mutable[point][difference][relation] = (mutable[point][difference][relation] - mean) % PRIME
    return tuple(tuple(tuple(row) for row in plane) for plane in mutable)


def three_way_interaction(tensor):
    return centre_axis(centre_axis(centre_axis(tensor, 0), 1), 2)


def flattening_ranks(tensor):
    point_flat = tuple(
        tuple(tensor[point][difference][relation] for difference in range(P) for relation in range(P))
        for point in range(W)
    )
    difference_flat = tuple(
        tuple(tensor[point][difference][relation] for point in range(W) for relation in range(P))
        for difference in range(P)
    )
    relation_flat = tuple(
        tuple(tensor[point][difference][relation] for point in range(W) for difference in range(P))
        for relation in range(P)
    )
    return C.rank_mod(point_flat), C.rank_mod(difference_flat), C.rank_mod(relation_flat)


def fourier_2d(matrix, zeta: int):
    roots = tuple(
        tuple(pow(zeta, -frequency * value % P, PRIME) for value in range(P))
        for frequency in range(P)
    )
    first = tuple(
        tuple(
            sum(matrix[difference][relation] * roots[difference_frequency][difference]
                for difference in range(P)) % PRIME
            for relation in range(P)
        )
        for difference_frequency in range(P)
    )
    return tuple(
        tuple(
            sum(first[difference_frequency][relation] * roots[relation_frequency][relation]
                for relation in range(P)) % PRIME
            for relation_frequency in range(P)
        )
        for difference_frequency in range(P)
    )


def path_edge_spectra(tensor, zeta: int):
    mixed = centre_axis(centre_axis(tensor, 1), 2)
    return tuple(
        fourier_2d(
            tuple(
                tuple((mixed[right][difference][relation]
                       - mixed[left][difference][relation]) % PRIME
                      for relation in range(P))
                for difference in range(P)
            ),
            zeta,
        )
        for left, right in zip(range(W - 1), range(1, W))
    )


def edge_support_records(edge_spectra):
    records = []
    for spectrum in edge_spectra:
        zero_pairs = {
            (difference, relation)
            for difference in range(1, P) for relation in range(1, P)
            if spectrum[difference][relation] == 0
        }
        zero_relation_lines = tuple(
            relation
            for relation in range(1, P)
            if all((difference, relation) in zero_pairs for difference in range(1, P))
        )
        covered = {
            (difference, relation)
            for difference in range(1, P) for relation in zero_relation_lines
        }
        residual_zeros = tuple(sorted(zero_pairs - covered))
        require(not residual_zeros,
                ("path-edge zeros are not complete relation lines", residual_zeros))
        records.append(
            (
                sum(value != 0 for row in spectrum for value in row),
                sum(spectrum[difference][relation] != 0
                    for difference in range(1, P) for relation in range(1, P)),
                zero_relation_lines,
                digest_json(spectrum),
            )
        )
    return tuple(records)


def fixed_relation_record(tensor, relation: int):
    matrix = tuple(
        tuple(tensor[point][difference][relation] for difference in range(P))
        for point in range(W)
    )
    inverse_w = pow(W, -1, PRIME)
    inverse_p = pow(P, -1, PRIME)
    inverse_total = pow(W * P, -1, PRIME)
    row_sums = tuple(sum(row) % PRIME for row in matrix)
    column_sums = tuple(sum(matrix[point][difference] for point in range(W)) % PRIME
                        for difference in range(P))
    grand = sum(row_sums) % PRIME
    interaction = tuple(
        tuple((matrix[point][difference]
               - row_sums[point] * inverse_p
               - column_sums[difference] * inverse_w
               + grand * inverse_total) % PRIME
              for difference in range(P))
        for point in range(W)
    )
    tail_contrasts = tuple(
        tuple((matrix[right][difference] - matrix[left][difference]) % PRIME
              for difference in range(P))
        for left, right in ((1, 2), (3, 4))
    )
    return (
        C.rank_mod(matrix),
        C.rank_mod(interaction),
        tuple(sum(value != 0 for value in contrast) for contrast in tail_contrasts),
        digest_json(matrix),
    )


def main() -> None:
    require(STATE_POINTS == ((0,), (1, 2), (5,), (3, 4)),
            ("pointed-state fibres", STATE_POINTS))
    ctx = C.context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "pointed worker order")

    def flatten_bank(position: int):
        return tuple(
            row
            for chunk in chunks
            for beta_rows in chunk[position]
            for row in beta_rows
        )

    weighted_gamma = flatten_bank(1)
    support_gamma = flatten_bank(2)
    require(len(weighted_gamma) == len(support_gamma) == P**3,
            "pointed bank size")
    require(digest_json(marginal_gamma(weighted_gamma)) == R.EXPECTED_DIGESTS[0][0],
            "pointed weighted gamma marginal")
    require(digest_json(marginal_gamma(support_gamma)) == R.EXPECTED_DIGESTS[0][1],
            "pointed support gamma marginal")
    require(all(row[point][0] == 0 for row in weighted_gamma for point in range(W)),
            "pointed weighted same-root")
    require(all(row[point][0] == 0 for row in support_gamma for point in range(W)),
            "pointed support same-root")

    weighted_tensor = inverse_tensor(weighted_gamma, ctx["zeta"])
    support_tensor = inverse_tensor(support_gamma, ctx["zeta"])
    weighted_parent = marginal_tensor(weighted_tensor)
    support_parent = marginal_tensor(support_tensor)
    require(digest_json(weighted_parent) == R.EXPECTED_DIGESTS[1][0],
            "pointed weighted tensor marginal")
    require(digest_json(support_parent) == R.EXPECTED_DIGESTS[1][1],
            "pointed support tensor marginal")
    flat_tensor = flat_tail_lift(weighted_parent)
    require(marginal_tensor(flat_tensor) == weighted_parent,
            "flat tail lift misses parent")

    tensors = (weighted_tensor, support_tensor, flat_tensor)
    interactions = tuple(three_way_interaction(tensor) for tensor in tensors)
    ranks = tuple(flattening_ranks(tensor) for tensor in tensors)
    interaction_ranks = tuple(flattening_ranks(tensor) for tensor in interactions)
    edge_spectra = tuple(path_edge_spectra(tensor, ctx["zeta"]) for tensor in tensors)
    edge_records = tuple(edge_support_records(spectra) for spectra in edge_spectra)
    relation = 6
    fixed_records = tuple(fixed_relation_record(tensor, relation) for tensor in tensors)
    difference_support = tuple(
        tuple(difference for difference in range(P)
              if any(tensor[point][difference][relation_value] != 0
                     for point in range(W) for relation_value in range(P)))
        for tensor in tensors
    )

    scalar_counts = tuple(sum(chunk[3][0][index] for chunk in chunks) for index in range(5))
    point_counts = tuple(sum(chunk[3][1][point] for chunk in chunks) for point in range(W))
    require(scalar_counts == C.EXPECTED_WORK_COUNTS,
            ("pointed work counts", scalar_counts))
    require(all(count > 0 for count in point_counts),
            ("empty pointed state", point_counts))

    direct_controls = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
    direct_record = []
    for alpha, beta, tau in direct_controls:
        direct = integrate_pointed(alpha, beta, tau)
        index = (alpha * P + beta) * P + tau
        phase = pow(ctx["zeta"], beta, PRIME)
        for bank_index, gamma_bank in enumerate((weighted_gamma, support_gamma)):
            direct_row = tuple(
                tuple(phase * value % PRIME for value in difference_row)
                for difference_row in direct[bank_index][0]
            )
            require(direct_row == gamma_bank[index],
                    ("pointed literal guard", alpha, beta, tau, bank_index))
        direct_record.append(((alpha, beta, tau), direct[2]))

    digests = (
        tuple(digest_json(bank) for bank in (weighted_gamma, support_gamma)),
        tuple(digest_json(tensor) for tensor in tensors),
        tuple(digest_json(tensor) for tensor in interactions),
        tuple(digest_json(spectra) for spectra in edge_spectra),
    )
    support_record = (difference_support, edge_records)
    support_digest = digest_json(support_record)
    rank_record = (ranks, interaction_ranks)
    if EXPECTED_DIGESTS != "TO_BE_PINNED":
        require(digests == EXPECTED_DIGESTS, ("pointed digests", digests))
    if EXPECTED_RANKS != "TO_BE_PINNED":
        require(rank_record == EXPECTED_RANKS, ("pointed ranks", rank_record))
    if EXPECTED_SUPPORT_SHA256 != "TO_BE_PINNED":
        require(support_digest == EXPECTED_SUPPORT_SHA256,
                ("pointed support record", support_digest, EXPECTED_SUPPORT_SHA256))
    if EXPECTED_FIXED != "TO_BE_PINNED":
        require(fixed_records == EXPECTED_FIXED,
                ("pointed fixed relation", fixed_records))

    record = (
        PARENT_SHA256,
        R.EXPECTED_SEMANTIC_SHA256,
        POINTED_STATES,
        PATH_EDGE_TYPES,
        STATE_POINTS,
        C.JOINT_PRIME,
        C.JOINT_ROOT,
        ctx["zeta"],
        scalar_counts,
        point_counts,
        tuple(direct_record),
        digests,
        rank_record,
        support_record,
        relation,
        fixed_records,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("pointed semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 owner-node pointed-six x root-difference x relation probe ==")
    print(f"parent=(sha256={PARENT_SHA256},semantic={R.EXPECTED_SEMANTIC_SHA256})")
    print(f"pointed_path={POINTED_STATES};edge_types={PATH_EDGE_TYPES};state_fibres={STATE_POINTS}")
    print(f"field=(prime={C.JOINT_PRIME},root={C.JOINT_ROOT},zeta13={ctx['zeta']})")
    print(f"work_counts={scalar_counts};point_segment_counts={point_counts};literal_controls={direct_controls}: PASS")
    print("marginals=(weighted,support_only)->four-state_root-difference_parent: PASS;flat_tail_same_parent: PASS")
    print("same_root=(weighted=0,support_only=0): PASS")
    print(f"difference_support_(weighted,support_only,flat_tail)={difference_support}")
    print(f"flattening_ranks_(pointed,difference,relation)={ranks}")
    print(f"three_way_ANOVA_flattening_ranks={interaction_ranks}")
    print(f"path_edge_support_(all,mixed,zero_relation_lines,digest)={edge_records}")
    print(f"support_record_sha256={support_digest}")
    print(f"fixed_relation=(1,0,{relation});records_(rank,ANOVA_rank,tail_contrast_support,digest)={fixed_records}")
    print(f"digests=(gamma,tensor,ANOVA,path_edge_spectra)={digests}")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT six-pointed-state source-root refinement candidate on one owner base")
    print("scope=marked source tail is lawful;no cyclic C6 claim,not exact address,not inverse ancestry,not chronology,not physical current,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
