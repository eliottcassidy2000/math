#!/usr/bin/env python3
"""Cross current-leg inverse ancestry with the ordered source cut colour.

The two pinned parents retain, separately,

    r_owner = a mod 13                    (current-leg inverse digit),
    s = u-q mod 13                        (ordered source-root difference).

This probe retains both before endpoint multiplication and integration.  Its
exact table has axes

    V_4(state) x F_13(r_owner) x F_13(s) x F_13(relation).

Summing ``r_owner`` must recover the root-difference parent, while summing
``s`` must recover the inverse-owner parent.  Four-axis ANOVA and all three
2+2 flattenings decide whether the two refinements have a joint component.
Branch-flat and difference-flat lifts preserve one parent apiece and are the
primary hostiles.

This is a finite common-base source/current statistic.  It is not an exact
THM-2334 address, a U_clock chronology, a physical current, or LRC(14).
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OWNER_PATH = (
    ROOT
    / "04-computation/lrc_r5_ufull_owner_node_boolean_square_inverse_owner_branch_probe_20260816.py"
)
DIFFERENCE_PATH = (
    ROOT
    / "04-computation/lrc_r5_ufull_owner_node_boolean_square_root_difference_probe_20260816.py"
)
POINTED_PATH = (
    ROOT
    / "04-computation/lrc_r5_ufull_owner_node_pointed_six_state_root_difference_probe_20260816.py"
)
OWNER_SHA256 = "ae1cf021ea23f325eeded42ff1dea8df903d837b1d7ab551289d62f7ab7a0348"
DIFFERENCE_SHA256 = "dddeea995e9ab7e8abfd010e00c798c93a28ff924da8172dd114756184942a57"
POINTED_SHA256 = "52b5af635b394e8f6dda59746d369b4b62a73da5ee89c6ca8e758426a7e81b76"
POINTED_GAMMA_SHA256 = "416f03a90894e526bc767a34839c2489aa4cac7051ac04687e0feacfa36d58d2"
POINTED_TENSOR_SHA256 = "9c5e227d9c142373973a562a54c6a67cac60a82da1121a028b9658920d155a19"

EXPECTED_DIGESTS = (
    (
        "843fe780896a23471f9e844383cc03cbdc70836a138a6addc218e87e319649de",
        "f5b3e6144659e701c4bb1229d102cca1a03e6fe3f89477a97d4bfc71ec669c1b",
        "416f03a90894e526bc767a34839c2489aa4cac7051ac04687e0feacfa36d58d2",
    ),
    (
        "ccfdb2373578190156dca9e9012c0a98ac50cdbacea84f1ff21d5ec0ac94db5b",
        "e52847bde897843610721db4d1d42478f5e308250b15f6459b82ab0c36374f61",
        "3f73e4273bf04332e4187d73d4ead9862f5dbf7c084818e5b637d4614bcbc989",
        "0ba2b6072c22ed3ef8058f54c428923b633892849b4cc3ddad0d137479cceca8",
    ),
    "9c5e227d9c142373973a562a54c6a67cac60a82da1121a028b9658920d155a19",
    (
        "ebd2587a1d30db95c20c5b891c67949127e87ea31aa66876d3390434e37f6489",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
    ),
    (
        "b971c6881f2be41ed30bb4400d398073290ab0bce8069be3b702b91dba3ddb6e",
        "0499eed3a2e3899bd84105a5b0dc5a72b83399743542e7c4f695c06ccb0b9449",
        "18568a89c220bb23f51d2e1a889c3cc91781797bad44e0dd086690f46c32e203",
        "eb4ad01f374e55b0bca39a47290aa5a9bf96953f1c5aff4693995d0a3b1421e5",
    ),
    (
        "383a35df28f6a6ef305ddc69d027b5cfd4134c09aaa47166842cafe73a9b9232",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
    ),
)
EXPECTED_RANKS = (
    ((4, 4, 12, 13), (4, 1, 12, 13), (4, 1, 12, 13), (4, 4, 1, 6)),
    ((6, 40, 42), (4, 13, 12), (4, 13, 12), (6, 4, 4)),
    ((1, 2, 1, 2), (1, 1, 1, 1), (1, 1, 1, 1), (1, 2, 1, 2)),
    ((1, 1, 1, 1), (2, 2, 2, 2), (1, 1, 1, 1), (2, 2, 2, 2)),
    ((3, 4, 12, 12), (0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)),
    ((6, 36, 34), (0, 0, 0), (0, 0, 0), (0, 0, 0)),
    ((4, 4, 4, 4), (0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)),
)
EXPECTED_CENSUSES = (
    (
        (8788, 1, 12, 12, 144, 12, 144, 144, 1728, 3, 36, 36, 432, 36, 432, 432, 5184),
        (676, 1, 12, 12, 144, 0, 0, 0, 0, 3, 36, 36, 432, 0, 0, 0, 0),
        (676, 1, 12, 12, 144, 0, 0, 0, 0, 3, 36, 36, 432, 0, 0, 0, 0),
        (676, 1, 12, 0, 0, 12, 144, 0, 0, 3, 36, 0, 0, 36, 432, 0, 0),
    ),
    (
        (5184, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5184),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    ),
)
EXPECTED_FIXED = (
    (
        (4, 4, 6), (3, 4, 6),
        (676, 1, 12, 12, 144, 3, 36, 36, 432),
        (432, 0, 0, 0, 0, 0, 0, 0, 432), 572,
        "748db7f48fce3756f5d0c8d4178b9da0a91020dc262e84beceb4e0c26b3ceef1",
        "1583ce224a363d748156df52ac8286d09862898784fbc9c4417dc12c70703582",
    ),
    (
        (4, 1, 4), (0, 0, 0),
        (52, 1, 12, 0, 0, 3, 36, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0), 572,
        "87f1fc0dd8846903deb3ac7bbbf9f2c9daccc34059ee01f9f3b6503222edf74a",
        "197447345f85eb9bcafa5ca23a56870cb07d8ba0635314ff498abf7f9ab9d5bc",
    ),
    (
        (4, 1, 4), (0, 0, 0),
        (52, 1, 12, 0, 0, 3, 36, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0), 572,
        "d722bcb22b6d2ddef2c392252c07aa26d05034a2fedc4147a463db90f47eaad1",
        "197447345f85eb9bcafa5ca23a56870cb07d8ba0635314ff498abf7f9ab9d5bc",
    ),
    (
        (4, 4, 1), (0, 0, 0),
        (52, 1, 0, 12, 0, 3, 0, 36, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0), 676,
        "31d83bc96f7398e0a4ea30f549cf2f8e3364b5afebb5277602a0acbd27574976",
        "197447345f85eb9bcafa5ca23a56870cb07d8ba0635314ff498abf7f9ab9d5bc",
    ),
)
EXPECTED_SEMANTIC_SHA256 = (
    "0e527913f22fc7d165f795df42e25692fba9289ccaf85aee4a3a2402613bdb41"
)

P = 13
V = 4
POINTED_STATES = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))
POINT_INDEX = {pair: index for index, pair in enumerate(POINTED_STATES)}
STATE_POINTS = tuple(
    tuple(index for index, (point_state, _root) in enumerate(POINTED_STATES)
          if point_state == state)
    for state in range(V)
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


def load_module(path: Path, expected: str, name: str):
    observed = lf_sha256(path)
    require(observed == expected, (name, "source drift", observed, expected))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load_module(OWNER_PATH, OWNER_SHA256, "joint_inverse_owner_parent")
D = load_module(DIFFERENCE_PATH, DIFFERENCE_SHA256, "joint_root_difference_parent")
B = R.B
C = R.C
PRIME = C.JOINT_PRIME

require(D.C.JOINT_PRIME == PRIME, "parent field mismatch")
require(D.B.EXPECTED_SEMANTIC_SHA256 == B.EXPECTED_SEMANTIC_SHA256,
        "parent Boolean-square mismatch")
require(lf_sha256(POINTED_PATH) == POINTED_SHA256, "pointed parent source drift")


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def integrate_joint(alpha: int, beta: int, literal_tau: int | None = None):
    ctx = C.context()
    branch_profiles, branch_boundaries, _source_record = R.owner_branch_context()
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    for boundary in branch_boundaries:
        events.setdefault(boundary, 0)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    weighted = [
        [
            [[0 for _difference in range(P)] for _branch in range(P)]
            for _state in range(V)
        ]
        for _tau in tau_values
    ]
    support_normalized = [
        [
            [[0 for _difference in range(P)] for _branch in range(P)]
            for _state in range(V)
        ]
        for _tau in tau_values
    ]
    pointed = [
        [[0 for _difference in range(P)] for _point in POINTED_STATES]
        for _tau in tau_values
    ]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    state_segments = [0] * V
    branch_segments = [0] * P
    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        require(C.cell_of_segment(left, right) == 0,
                ("joint branch/difference escaped cell zero", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("joint branch/difference middle chamber", alpha, beta, left, right))
        state = B.state_of_segment(left, right)
        active_segments += 1
        state_segments[state] += 1
        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])
        branch_values = tuple(
            tuple(profile_value(branch_profiles[root][branch], left)
                  for branch in range(P))
            for root in range(P)
        )
        require(all(sum(branch_values[root]) == u_values[root] for root in range(P)),
                ("joint local branch partition", left, right))
        u_support = tuple(root for root, value in enumerate(u_values) if value)
        v_support = tuple(root for root, value in enumerate(v_values) if value)
        if not u_support or not v_support:
            continue
        weighted_segments += 1
        require(set(u_support).isdisjoint(v_support),
                ("joint source supports overlap", left, right, u_support, v_support))
        for branch in range(P):
            if any(branch_values[root][branch] for root in u_support):
                branch_segments[branch] += 1
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
                require((state, u) in POINT_INDEX,
                        ("joint unrealized pointed state", left, right, state, u))
                point = POINT_INDEX[(state, u)]
                active_branches = tuple(
                    branch for branch in range(P) if branch_values[u][branch]
                )
                require(active_branches, ("joint empty branch fibre", left, right, u))
                support_weight = pow(len(active_branches), -1, PRIME)
                for q in selected_v:
                    difference = (u - q) % P
                    require(difference != 0,
                            ("joint same-root pair", left, right, u, q))
                    right_weight = v_values[q] * jump % PRIME
                    pointed[row_index][point][difference] = (
                        pointed[row_index][point][difference]
                        + u_values[u] * right_weight
                    ) % PRIME
                    for branch in active_branches:
                        weighted[row_index][state][branch][difference] = (
                            weighted[row_index][state][branch][difference]
                            + branch_values[u][branch] * right_weight
                        ) % PRIME
                        support_normalized[row_index][state][branch][difference] = (
                            support_normalized[row_index][state][branch][difference]
                            + support_weight * jump
                        ) % PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("joint endpoint mask", alpha, beta, literal_tau, mask))
    counts = (
        interval_count,
        mapped,
        active_segments,
        q_active_segments,
        weighted_segments,
        tuple(state_segments),
        tuple(branch_segments),
    )
    freeze = lambda bank: tuple(
        tuple(tuple(tuple(row) for row in state_rows) for state_rows in tau_rows)
        for tau_rows in bank
    )
    frozen_pointed = tuple(
        tuple(tuple(row) for row in point_rows) for point_rows in pointed
    )
    return freeze(weighted), freeze(support_normalized), frozen_pointed, counts


def worker(alpha: int):
    zeta = C.context()["zeta"]
    weighted_rows = []
    support_rows = []
    pointed_rows = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    branch_counts = [0] * P
    for beta in range(P):
        weighted, support_normalized, pointed, record = integrate_joint(alpha, beta)
        phase = pow(zeta, beta, PRIME)

        def twist(bank):
            return tuple(
                tuple(
                    tuple(
                        tuple(phase * value % PRIME for value in difference_row)
                        for difference_row in branch_rows
                    )
                    for branch_rows in state_rows
                )
                for state_rows in bank
            )

        weighted_rows.append(twist(weighted))
        support_rows.append(twist(support_normalized))
        pointed_rows.append(tuple(
            tuple(
                tuple(phase * value % PRIME for value in difference_row)
                for difference_row in point_rows
            )
            for point_rows in pointed
        ))
        scalar_counts = [left + right for left, right in zip(scalar_counts, record[:5])]
        state_counts = [left + right for left, right in zip(state_counts, record[5])]
        branch_counts = [left + right for left, right in zip(branch_counts, record[6])]
    return (
        alpha,
        tuple(weighted_rows),
        tuple(support_rows),
        tuple(pointed_rows),
        (tuple(scalar_counts), tuple(state_counts), tuple(branch_counts)),
    )


def inverse_tensor(gamma_rows, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [
        [
            [[0 for _relation in range(P)] for _difference in range(P)]
            for _branch in range(P)
        ]
        for _state in range(V)
    ]
    phases = tuple(
        tuple(pow(zeta, -(alpha + tau * relation) % P, PRIME)
              for relation in range(P))
        for alpha in range(P) for tau in range(P)
    )
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_rows[index]
                phase_row = phases[alpha * P + tau]
                index += 1
                for state in range(V):
                    for branch in range(P):
                        for difference in range(P):
                            value = row[state][branch][difference]
                            if not value:
                                continue
                            for relation, phase in enumerate(phase_row):
                                tensor[state][branch][difference][relation] = (
                                    tensor[state][branch][difference][relation]
                                    + value * phase
                                ) % PRIME
    require(index == P**3, ("joint gamma size", index))
    return tuple(
        tuple(
            tuple(
                tuple(value * normalizer % PRIME for value in relation_row)
                for relation_row in difference_rows
            )
            for difference_rows in branch_rows
        )
        for branch_rows in tensor
    )


def inverse_pointed_tensor(gamma_rows, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [
        [[0 for _relation in range(P)] for _difference in range(P)]
        for _point in POINTED_STATES
    ]
    phases = tuple(
        tuple(pow(zeta, -(alpha + tau * relation) % P, PRIME)
              for relation in range(P))
        for alpha in range(P) for tau in range(P)
    )
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_rows[index]
                phase_row = phases[alpha * P + tau]
                index += 1
                for point in range(len(POINTED_STATES)):
                    for difference in range(P):
                        value = row[point][difference]
                        if not value:
                            continue
                        for relation, phase in enumerate(phase_row):
                            tensor[point][difference][relation] = (
                                tensor[point][difference][relation] + value * phase
                            ) % PRIME
    require(index == P**3, ("joint pointed gamma size", index))
    return tuple(
        tuple(
            tuple(value * normalizer % PRIME for value in relation_row)
            for relation_row in difference_rows
        )
        for difference_rows in tensor
    )


def gamma_branch_marginal(bank):
    return tuple(
        tuple(
            tuple(
                sum(row[state][branch][difference] for branch in range(P)) % PRIME
                for difference in range(P)
            )
            for state in range(V)
        )
        for row in bank
    )


def gamma_difference_marginal(bank):
    return tuple(
        tuple(
            tuple(
                sum(row[state][branch][difference] for difference in range(P)) % PRIME
                for branch in range(P)
            )
            for state in range(V)
        )
        for row in bank
    )


def tensor_branch_marginal(tensor):
    return tuple(
        tuple(
            tuple(
                sum(tensor[state][branch][difference][relation]
                    for branch in range(P)) % PRIME
                for relation in range(P)
            )
            for difference in range(P)
        )
        for state in range(V)
    )


def tensor_difference_marginal(tensor):
    return tuple(
        tuple(
            tuple(
                sum(tensor[state][branch][difference][relation]
                    for difference in range(P)) % PRIME
                for relation in range(P)
            )
            for branch in range(P)
        )
        for state in range(V)
    )


def branch_flat_lift(tensor):
    marginal = tensor_branch_marginal(tensor)
    inverse = pow(P, -1, PRIME)
    return tuple(
        tuple(
            tuple(
                tuple(marginal[state][difference][relation] * inverse % PRIME
                      for relation in range(P))
                for difference in range(P)
            )
            for _branch in range(P)
        )
        for state in range(V)
    )


def difference_flat_lift(tensor):
    marginal = tensor_difference_marginal(tensor)
    inverse = pow(P, -1, PRIME)
    return tuple(
        tuple(
            tuple(
                tuple(marginal[state][branch][relation] * inverse % PRIME
                      for relation in range(P))
                for _difference in range(P)
            )
            for branch in range(P)
        )
        for state in range(V)
    )


def centre_axis(tensor, axis: int):
    dimensions = (V, P, P, P)
    inverse = pow(dimensions[axis], -1, PRIME)
    mutable = [[[[tensor[s][b][d][t] for t in range(P)] for d in range(P)]
                for b in range(P)] for s in range(V)]
    if axis == 0:
        for b in range(P):
            for d in range(P):
                for t in range(P):
                    mean = sum(mutable[s][b][d][t] for s in range(V)) * inverse % PRIME
                    for s in range(V):
                        mutable[s][b][d][t] = (mutable[s][b][d][t] - mean) % PRIME
    elif axis == 1:
        for s in range(V):
            for d in range(P):
                for t in range(P):
                    mean = sum(mutable[s][b][d][t] for b in range(P)) * inverse % PRIME
                    for b in range(P):
                        mutable[s][b][d][t] = (mutable[s][b][d][t] - mean) % PRIME
    elif axis == 2:
        for s in range(V):
            for b in range(P):
                for t in range(P):
                    mean = sum(mutable[s][b][d][t] for d in range(P)) * inverse % PRIME
                    for d in range(P):
                        mutable[s][b][d][t] = (mutable[s][b][d][t] - mean) % PRIME
    else:
        for s in range(V):
            for b in range(P):
                for d in range(P):
                    mean = sum(mutable[s][b][d][t] for t in range(P)) * inverse % PRIME
                    for t in range(P):
                        mutable[s][b][d][t] = (mutable[s][b][d][t] - mean) % PRIME
    return tuple(
        tuple(tuple(tuple(row) for row in plane) for plane in block)
        for block in mutable
    )


def four_way_interaction(tensor):
    for axis in range(4):
        tensor = centre_axis(tensor, axis)
    return tensor


def axis_ranks(tensor):
    state = tuple(
        tuple(tensor[s][b][d][t] for b in range(P) for d in range(P) for t in range(P))
        for s in range(V)
    )
    branch = tuple(
        tuple(tensor[s][b][d][t] for s in range(V) for d in range(P) for t in range(P))
        for b in range(P)
    )
    difference = tuple(
        tuple(tensor[s][b][d][t] for s in range(V) for b in range(P) for t in range(P))
        for d in range(P)
    )
    relation = tuple(
        tuple(tensor[s][b][d][t] for s in range(V) for b in range(P) for d in range(P))
        for t in range(P)
    )
    return tuple(C.rank_mod(matrix) for matrix in (state, branch, difference, relation))


def bipartition_ranks(tensor):
    state_branch = tuple(
        tuple(tensor[s][b][d][t] for d in range(P) for t in range(P))
        for s in range(V) for b in range(P)
    )
    state_difference = tuple(
        tuple(tensor[s][b][d][t] for b in range(P) for t in range(P))
        for s in range(V) for d in range(P)
    )
    state_relation = tuple(
        tuple(tensor[s][b][d][t] for b in range(P) for d in range(P))
        for s in range(V) for t in range(P)
    )
    return tuple(C.rank_mod(matrix)
                 for matrix in (state_branch, state_difference, state_relation))


def state_branch_block_ranks(tensor):
    """Rank contributed by the thirteen branch rows inside each fixed state."""

    return tuple(
        C.rank_mod(tuple(
            tuple(tensor[state][branch][difference][relation]
                  for difference in range(P) for relation in range(P))
            for branch in range(P)
        ))
        for state in range(V)
    )


def pointed_carrier_record(tensor, pointed_tensor):
    """Certify each state-branch row space equals its pointed-tail row space."""

    record = []
    for state in range(V):
        branch_rows = tuple(
            tuple(tensor[state][branch][difference][relation]
                  for difference in range(P) for relation in range(P))
            for branch in range(P)
        )
        point_rows = tuple(
            tuple(pointed_tensor[point][difference][relation]
                  for difference in range(P) for relation in range(P))
            for point in STATE_POINTS[state]
        )
        ranks = (
            len(STATE_POINTS[state]),
            C.rank_mod(branch_rows),
            C.rank_mod(point_rows),
            C.rank_mod(branch_rows + point_rows),
        )
        require(ranks == (len(STATE_POINTS[state]),) * 4,
                ("pointed carrier row-space mismatch", state, ranks))
        record.append(ranks)
    return tuple(record)


def fourier_tensor(tensor, zeta: int):
    roots = tuple(
        tuple(pow(zeta, -frequency * value % P, PRIME) for value in range(P))
        for frequency in range(P)
    )
    walsh = [[[[
        sum(tensor[state][branch][difference][relation]
            * B.WALSH_SIGNS[character][state] for state in range(V)) % PRIME
        for relation in range(P)] for difference in range(P)] for branch in range(P)]
        for character in range(V)]
    branch_stage = [[[[
        sum(walsh[character][branch][difference][relation]
            * roots[branch_frequency][branch] for branch in range(P)) % PRIME
        for relation in range(P)] for difference in range(P)]
        for branch_frequency in range(P)] for character in range(V)]
    difference_stage = [[[[
        sum(branch_stage[character][branch_frequency][difference][relation]
            * roots[difference_frequency][difference] for difference in range(P)) % PRIME
        for relation in range(P)] for difference_frequency in range(P)]
        for branch_frequency in range(P)] for character in range(V)]
    return tuple(
        tuple(
            tuple(
                tuple(
                    sum(difference_stage[character][branch_frequency]
                        [difference_frequency][relation]
                        * roots[relation_frequency][relation]
                        for relation in range(P)) % PRIME
                    for relation_frequency in range(P)
                )
                for difference_frequency in range(P)
            )
            for branch_frequency in range(P)
        )
        for character in range(V)
    )


def support_census(spectrum):
    bins = [0] * 16
    for character in range(V):
        for branch_frequency in range(P):
            for difference_frequency in range(P):
                for relation_frequency in range(P):
                    if spectrum[character][branch_frequency][difference_frequency][relation_frequency]:
                        mask = (
                            ((character != 0) << 3)
                            | ((branch_frequency != 0) << 2)
                            | ((difference_frequency != 0) << 1)
                            | (relation_frequency != 0)
                        )
                        bins[mask] += 1
    return (sum(bins),) + tuple(bins)


def fixed_relation_record(tensor, relation: int):
    fixed = tuple(
        tuple(
            tuple(tensor[state][branch][difference][relation]
                  for difference in range(P))
            for branch in range(P)
        )
        for state in range(V)
    )
    interaction = D.three_way_interaction(fixed)
    spectrum = D.fourier_tensor(fixed, C.context()["zeta"])
    interaction_spectrum = D.fourier_tensor(interaction, C.context()["zeta"])
    return (
        D.flattening_ranks(fixed),
        D.flattening_ranks(interaction),
        D.support_census(spectrum),
        D.support_census(interaction_spectrum),
        sum(value != 0 for state in fixed for branch in state for value in branch),
        digest_json(fixed),
        digest_json(interaction),
    )


def main() -> None:
    ctx = C.context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "joint worker order")

    def flatten_bank(position: int):
        return tuple(
            row
            for chunk in chunks
            for beta_rows in chunk[position]
            for row in beta_rows
        )

    weighted_gamma = flatten_bank(1)
    support_gamma = flatten_bank(2)
    pointed_gamma = flatten_bank(3)
    require(len(weighted_gamma) == len(support_gamma) == len(pointed_gamma) == P**3,
            "joint bank size")
    require(digest_json(pointed_gamma) == POINTED_GAMMA_SHA256,
            "joint pointed gamma misses audited candidate")

    weighted_branch_marginal = gamma_branch_marginal(weighted_gamma)
    support_branch_marginal = gamma_branch_marginal(support_gamma)
    weighted_difference_marginal = gamma_difference_marginal(weighted_gamma)
    require(digest_json(weighted_branch_marginal) == D.EXPECTED_DIGESTS[0][0],
            "joint branch marginal misses weighted root-difference parent")
    require(digest_json(support_branch_marginal) == D.EXPECTED_DIGESTS[0][1],
            "joint branch marginal misses support root-difference parent")
    require(digest_json(weighted_difference_marginal) == R.EXPECTED_DIGESTS[0],
            "joint difference marginal misses inverse-owner parent")
    require(all(row[state][branch][0] == 0
                for row in weighted_gamma for state in range(V) for branch in range(P)),
            "joint weighted same-root slice")
    require(all(row[state][branch][0] == 0
                for row in support_gamma for state in range(V) for branch in range(P)),
            "joint support same-root slice")

    actual = inverse_tensor(weighted_gamma, ctx["zeta"])
    support = inverse_tensor(support_gamma, ctx["zeta"])
    pointed_tensor = inverse_pointed_tensor(pointed_gamma, ctx["zeta"])
    require(digest_json(pointed_tensor) == POINTED_TENSOR_SHA256,
            "joint pointed tensor misses audited candidate")
    require(digest_json(tensor_branch_marginal(actual)) == D.EXPECTED_DIGESTS[1][0],
            "joint inverse branch marginal misses root-difference tensor")
    require(digest_json(tensor_branch_marginal(support)) == D.EXPECTED_DIGESTS[1][1],
            "joint inverse support marginal misses root-difference tensor")
    require(digest_json(tensor_difference_marginal(actual)) == R.EXPECTED_DIGESTS[1],
            "joint inverse difference marginal misses inverse-owner tensor")

    branch_flat = branch_flat_lift(actual)
    difference_flat = difference_flat_lift(actual)
    require(tensor_branch_marginal(branch_flat) == tensor_branch_marginal(actual),
            "branch-flat parent")
    require(tensor_difference_marginal(difference_flat)
            == tensor_difference_marginal(actual), "difference-flat parent")
    tensors = (actual, support, branch_flat, difference_flat)
    interactions = tuple(four_way_interaction(tensor) for tensor in tensors)
    axis_rank_records = tuple(axis_ranks(tensor) for tensor in tensors)
    interaction_axis_ranks = tuple(axis_ranks(tensor) for tensor in interactions)
    bipartition_records = tuple(bipartition_ranks(tensor) for tensor in tensors)
    interaction_bipartition_ranks = tuple(
        bipartition_ranks(tensor) for tensor in interactions
    )
    state_block_records = tuple(state_branch_block_ranks(tensor) for tensor in tensors)
    carrier_record = pointed_carrier_record(actual, pointed_tensor)
    interaction_state_block_ranks = tuple(
        state_branch_block_ranks(tensor) for tensor in interactions
    )
    spectra = tuple(fourier_tensor(tensor, ctx["zeta"]) for tensor in tensors)
    interaction_spectra = tuple(
        fourier_tensor(tensor, ctx["zeta"]) for tensor in interactions
    )
    censuses = tuple(support_census(spectrum) for spectrum in spectra)
    interaction_censuses = tuple(
        support_census(spectrum) for spectrum in interaction_spectra
    )
    relation = 6
    fixed_records = tuple(fixed_relation_record(tensor, relation) for tensor in tensors)

    scalar_counts = tuple(sum(chunk[4][0][index] for chunk in chunks) for index in range(5))
    state_counts = tuple(sum(chunk[4][1][state] for chunk in chunks) for state in range(V))
    branch_counts = tuple(sum(chunk[4][2][branch] for chunk in chunks) for branch in range(P))
    require(all(count > 0 for count in branch_counts),
            ("joint empty branch", branch_counts))

    direct_controls = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
    direct_record = []
    for alpha, beta, tau in direct_controls:
        direct = integrate_joint(alpha, beta, tau)
        phase = pow(ctx["zeta"], beta, PRIME)
        index = (alpha * P + beta) * P + tau
        for bank_index, gamma_bank in enumerate((weighted_gamma, support_gamma)):
            direct_row = tuple(
                tuple(
                    tuple(phase * value % PRIME for value in difference_row)
                    for difference_row in branch_rows
                )
                for branch_rows in direct[bank_index][0]
            )
            require(direct_row == gamma_bank[index],
                    ("joint literal guard", alpha, beta, tau, bank_index))
        direct_pointed = tuple(
            tuple(phase * value % PRIME for value in difference_row)
            for difference_row in direct[2][0]
        )
        require(direct_pointed == pointed_gamma[index],
                ("joint pointed literal guard", alpha, beta, tau))
        direct_record.append(((alpha, beta, tau), direct[3]))

    digests = (
        tuple(digest_json(bank) for bank in (weighted_gamma, support_gamma, pointed_gamma)),
        tuple(digest_json(tensor) for tensor in tensors),
        digest_json(pointed_tensor),
        tuple(digest_json(tensor) for tensor in interactions),
        tuple(digest_json(spectrum) for spectrum in spectra),
        tuple(digest_json(spectrum) for spectrum in interaction_spectra),
    )
    ranks = (
        axis_rank_records,
        bipartition_records,
        state_block_records,
        carrier_record,
        interaction_axis_ranks,
        interaction_bipartition_ranks,
        interaction_state_block_ranks,
    )
    census_record = (censuses, interaction_censuses)
    if EXPECTED_DIGESTS != "TO_BE_PINNED":
        require(digests == EXPECTED_DIGESTS, ("joint digests", digests))
    if EXPECTED_RANKS != "TO_BE_PINNED":
        require(ranks == EXPECTED_RANKS, ("joint ranks", ranks))
    if EXPECTED_CENSUSES != "TO_BE_PINNED":
        require(census_record == EXPECTED_CENSUSES,
                ("joint censuses", census_record))
    if EXPECTED_FIXED != "TO_BE_PINNED":
        require(fixed_records == EXPECTED_FIXED,
                ("joint fixed relation", fixed_records))

    record = (
        OWNER_SHA256,
        R.EXPECTED_SEMANTIC_SHA256,
        DIFFERENCE_SHA256,
        D.EXPECTED_SEMANTIC_SHA256,
        POINTED_SHA256,
        POINTED_GAMMA_SHA256,
        POINTED_TENSOR_SHA256,
        C.JOINT_PRIME,
        C.JOINT_ROOT,
        ctx["zeta"],
        scalar_counts,
        state_counts,
        branch_counts,
        tuple(direct_record),
        digests,
        ranks,
        census_record,
        relation,
        fixed_records,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("joint semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 current-leg inverse branch x source-root difference four-way probe ==")
    print(f"parents=(owner_sha={OWNER_SHA256},owner_semantic={R.EXPECTED_SEMANTIC_SHA256},difference_sha={DIFFERENCE_SHA256},difference_semantic={D.EXPECTED_SEMANTIC_SHA256},pointed_sha={POINTED_SHA256})")
    print("coordinates=(state in V4,r_owner=a mod13,source_difference=u-q mod13,relation=(1,0,t))")
    print(f"field=(prime={C.JOINT_PRIME},root={C.JOINT_ROOT},zeta13={ctx['zeta']})")
    print(f"work_counts={scalar_counts};state_counts={state_counts};branch_counts={branch_counts};literal_controls={direct_controls}: PASS")
    print("marginals=(sum_r->weighted/support_root_difference,sum_s->inverse_owner): PASS;same_root=0")
    print("hostiles=(branch_flat_same_root_difference_parent,difference_flat_same_inverse_owner_parent): PASS")
    print(f"axis_ranks_(state,branch,difference,relation)_(actual,support,branch_flat,difference_flat)={axis_rank_records}")
    print(f"bipartition_ranks_(state_branch|difference_relation,state_difference|branch_relation,state_relation|branch_difference)={bipartition_records}")
    print(f"fixed_state_branch_block_ranks={state_block_records}")
    print(f"pointed_carrier_rowspace_(point_count,branch_rank,point_rank,union_rank)={carrier_record}: PASS")
    print(f"four_way_ANOVA_axis_ranks={interaction_axis_ranks}")
    print(f"four_way_ANOVA_bipartition_ranks={interaction_bipartition_ranks}")
    print(f"four_way_ANOVA_fixed_state_branch_block_ranks={interaction_state_block_ranks}")
    print("spectrum_mask_bits=(state,branch,difference,relation);census=(total,mask0,...,mask15)")
    print(f"spectral_censuses={censuses}")
    print(f"four_way_ANOVA_spectral_censuses={interaction_censuses}")
    print(f"fixed_relation=(1,0,{relation});records_(axis_ranks,ANOVA_ranks,spectral_census,ANOVA_census,nonzero,digests)={fixed_records}")
    print(f"digests=(gamma_with_pointed,tensor,pointed_tensor,ANOVA,spectrum,ANOVA_spectrum)={digests}")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT four-coordinate current-ancestry/source-cut candidate on one owner base")
    print("scope=joint finite statistic;not exact address,not chronology,not physical current,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
