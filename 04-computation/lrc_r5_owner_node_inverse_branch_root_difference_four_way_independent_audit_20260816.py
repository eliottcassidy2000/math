#!/usr/bin/env python3
"""Independent hostile audit of current digit x source-root difference.

The submitted four-way implementation is never imported.  This audit starts
from two hash-pinned independent parents: the current-leg digit audit and the
six-pointed source-root-difference audit.  It reconstructs the stronger

    point=(state,u) x r_owner x (u-q) x relation

tensor directly from source profiles, literal endpoint events, and ordered
root pairs.  The submitted state x digit x difference x relation tensor is
then only a pointed-state marginal.

The statewise row-space equality is tested by union ranks and explicit exact
coordinates in the pointed basis.  Branch amplitude rank, the six relation
carrier lines, and whole-tensor bipartition ranks are kept distinct.  No
complete address, chronology, physical current, row exclusion, or LRC(14)
conclusion is asserted.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CURRENT_PARENT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_inverse_owner_branch_"
    "independent_audit_20260816.py"
)
POINTED_PARENT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_ufull_owner_node_pointed_six_state_root_difference_"
    "independent_audit_20260816.py"
)
TWO_CURRENT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_two_current_digits_pointed_root_difference_carrier_transition_"
    "probe_20260816.py"
)
SOURCE_CURRENT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_nested_ancestry_digits_"
    "probe_20260816.py"
)

CURRENT_PARENT_SHA256 = (
    "4b0bd05ffa6195ff484433329e334d771bc27e7cd380136b50b45e7248bb98ba"
)
CURRENT_PARENT_SEMANTIC = (
    "7063720d0e0e4847ce752102de83274ea47d7740fc435a64bae425dbd7100121"
)
POINTED_PARENT_SHA256 = (
    "8c7cb5f98b15a768d4f4d6060074e0815a8f089f857ec4f3c55a0e7d877e1fec"
)
POINTED_PARENT_SEMANTIC = (
    "66db2301f88db1ced7784868095e198e3e12f1fb79175ebd902d0f569a5decef"
)
TWO_CURRENT_SHA256 = (
    "9d1671e0f823fdbaa9ab79915ba05dbb4dda4c6eabb97fe4484baf4c2e3205f2"
)
TWO_CURRENT_SEMANTIC = (
    "38725dc1d7129b326634c99bd70e1eb414590dc24fb83bd9522e2095e41f204c"
)
SOURCE_CURRENT_SHA256 = (
    "1188df8aa2a7a84c1e8ada5fc3cc8d3b839ece70298b94f1d94c9d440caa88f3"
)
SOURCE_CURRENT_SEMANTIC = (
    "6e5605f58b7a94ea5ea4e8f62cfa7ee135b0d52512225f4aaee248ad6e21a9ae"
)

P = 13
V = 4
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_STATE_COUNTS = (300115, 300116, 300115, 300116)
EXPECTED_BRANCH_COUNTS = (186244,) * P
EXPECTED_RANKS = (
    ((4, 4, 12, 13), (4, 1, 12, 13), (4, 1, 12, 13), (4, 4, 1, 6)),
    ((6, 40, 42), (4, 13, 12), (4, 13, 12), (6, 4, 4)),
    ((1, 2, 1, 2), (1, 1, 1, 1), (1, 1, 1, 1), (1, 2, 1, 2)),
    ((3, 4, 12, 12), (0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)),
    ((6, 36, 34), (0, 0, 0), (0, 0, 0), (0, 0, 0)),
)
EXPECTED_CARRIER = (
    (1, 1, 1, 1),
    (2, 2, 2, 2),
    (1, 1, 1, 1),
    (2, 2, 2, 2),
)
EXPECTED_CENSUSES = (
    (
        (8788, 1, 12, 12, 144, 12, 144, 144, 1728,
         3, 36, 36, 432, 36, 432, 432, 5184),
        (676, 1, 12, 12, 144, 0, 0, 0, 0,
         3, 36, 36, 432, 0, 0, 0, 0),
        (676, 1, 12, 12, 144, 0, 0, 0, 0,
         3, 36, 36, 432, 0, 0, 0, 0),
        (676, 1, 12, 0, 0, 12, 144, 0, 0,
         3, 36, 0, 0, 36, 432, 0, 0),
    ),
    (
        (5184, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 5184),
        (0,) * 17,
        (0,) * 17,
        (0,) * 17,
    ),
)
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
EXPECTED_SEMANTIC_SHA256 = "57bb908c381608415091feb14035354e7cef6dd3b11a4e9ce78c363fc703f2a0"


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


D = load_module(CURRENT_PARENT_PATH, CURRENT_PARENT_SHA256, "joint_current_parent")
PS = load_module(POINTED_PARENT_PATH, POINTED_PARENT_SHA256, "joint_pointed_parent")
require(D.EXPECTED_SEMANTIC_SHA256 == CURRENT_PARENT_SEMANTIC,
        "current parent semantic drift")
require(PS.EXPECTED_SEMANTIC_SHA256 == POINTED_PARENT_SEMANTIC,
        "pointed parent semantic drift")
require(lf_sha256(TWO_CURRENT_PATH) == TWO_CURRENT_SHA256,
        "two-current comparison source drift")
require(lf_sha256(SOURCE_CURRENT_PATH) == SOURCE_CURRENT_SHA256,
        "source-current comparison source drift")

RD = PS.RD
SQ = D.SQ
R = D.R
PRIME = D.PRIME
POINTS = PS.POINTS
POINT_INDEX = PS.POINT_INDEX
STATE_FIBRES = PS.STATE_FIBRES
require(POINTS == ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12)),
        ("pointed order", POINTS))
require(PRIME == PS.PRIME, "split-field mismatch")


def freeze4(value):
    return tuple(tuple(tuple(tuple(line) for line in block) for block in slab)
                 for slab in value)


def freeze3(value):
    return tuple(tuple(tuple(line) for line in slab) for slab in value)


def integrate_pair(alpha: int, beta: int, literal_tau: int | None = None):
    (
        branches, branch_boundaries, source_u, source_v,
        _source_boundaries, _source_digest, _source_record, _source_sha,
        _reflection_record,
    ) = D.owner_digit_profiles()
    word, endpoint_grid, harmonic, danger = D.endpoint_context()
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    actual = [[[[[0] * P for _r in range(P)] for _point in POINTS]
               for _tau in tau_values]][0]
    support = [[[[[0] * P for _r in range(P)] for _state in range(V)]
                for _tau in tau_values]][0]
    positions = sorted(
        set(events) | set(branch_boundaries)
        | {point for profile in source_v for point in profile[0]}
    )
    inv_p = pow(P, -1, PRIME)
    mask = 0
    active = q_active = weighted_intervals = 0
    state_counts = [0] * V
    branch_counts = [0] * P
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
                ("joint escaped cell zero", alpha, beta, left, right))
        chamber = R.chamber(left, right)
        require(chamber in ("left", "right"),
                ("joint middle chamber", alpha, beta, left, right))
        state = SQ.state_index(source_u, left)
        state_counts[state] += 1
        if not jump:
            continue
        q_active += 1
        u_values = tuple(R.profile_value(profile, left) for profile in source_u)
        v_values = tuple(R.profile_value(profile, left) for profile in source_v)
        digit_values = tuple(tuple(D.profile_value(branches[root][digit], left)
                                   for digit in range(P))
                             for root in range(P))
        require(all(sum(digit_values[root]) == u_values[root]
                    for root in range(P)), ("pointwise digit marginal", left))
        require(all(u_values[root] * v_values[root] == 0 for root in range(P)),
                ("pointwise same-root", left))
        for digit in range(P):
            if any(digit_values[root][digit] for root in range(P)):
                branch_counts[digit] += 1
        if any(u_values) and any(v_values):
            weighted_intervals += 1

        for tau_index, tau in enumerate(tau_values):
            selected_mask = mask if literal_tau is not None else (
                mask & R.guard_mask(chamber, tau, danger)
            )
            selected = tuple(root for root in range(P)
                             if (selected_mask >> root) & 1)
            direct_weight = 0
            direct_support = 0
            for u in selected:
                if not u_values[u]:
                    continue
                point = POINT_INDEX.get((state, u))
                require(point is not None, ("unrealized pointed tail", state, u))
                require(all(digit_values[u][digit] != 0 for digit in range(P)),
                        ("branch support not full", left, state, u))
                for q in selected:
                    if not v_values[q]:
                        continue
                    require(u != q, ("selected same-root", alpha, beta, tau, left, u))
                    difference = (u - q) % P
                    direct_weight += u_values[u] * v_values[q]
                    direct_support += 1
                    for digit in range(P):
                        weighted = digit_values[u][digit] * v_values[q]
                        actual[tau_index][point][digit][difference] = (
                            actual[tau_index][point][digit][difference]
                            + weighted * jump
                        ) % PRIME
                        support[tau_index][state][digit][difference] = (
                            support[tau_index][state][digit][difference]
                            + inv_p * jump
                        ) % PRIME
            require(direct_weight == (
                sum(u_values[root] for root in selected)
                * sum(v_values[root] for root in selected)
            ), ("weighted pair factor", left, tau))
            require(direct_support == (
                sum(u_values[root] != 0 for root in selected)
                * sum(v_values[root] != 0 for root in selected)
            ), ("support pair factor", left, tau))
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("joint endpoint mask closure", alpha, beta, literal_tau))
    return (
        tuple(tuple(tuple(tuple(line) for line in point) for point in tau)
              for tau in actual),
        tuple(tuple(tuple(tuple(line) for line in state) for state in tau)
              for tau in support),
        (interval_count, mapped_count, active, q_active, weighted_intervals),
        tuple(state_counts),
        tuple(branch_counts),
    )


def state_row(point_row):
    return tuple(tuple(tuple(
        sum(point_row[point][digit][difference]
            for point in STATE_FIBRES[state]) % PRIME
        for difference in range(P)) for digit in range(P))
        for state in range(V))


def pointed_row(point_row):
    return tuple(tuple(
        sum(point_row[point][digit][difference] for digit in range(P)) % PRIME
        for difference in range(P)) for point in range(len(POINTS)))


def worker(alpha: int):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    gamma_actual = []
    gamma_support = []
    gamma_pointed = []
    point_core = [[[[0] * P for _digit in range(P)] for _point in POINTS]
                  for _tau in range(P)]
    support_core = [[[[0] * P for _digit in range(P)] for _state in range(V)]
                    for _tau in range(P)]
    totals = [0] * 5
    state_totals = [0] * V
    branch_totals = [0] * P
    for beta in range(P):
        actual, support, counts, state_counts, branch_counts = integrate_pair(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        for tau in range(P):
            twisted_point = tuple(tuple(tuple(
                phase * actual[tau][point][digit][difference] % PRIME
                for difference in range(P)) for digit in range(P))
                for point in range(len(POINTS)))
            twisted_support = tuple(tuple(tuple(
                phase * support[tau][state][digit][difference] % PRIME
                for difference in range(P)) for digit in range(P))
                for state in range(V))
            gamma_actual.append(state_row(twisted_point))
            gamma_support.append(twisted_support)
            gamma_pointed.append(pointed_row(twisted_point))
            for point in range(len(POINTS)):
                for digit in range(P):
                    for difference in range(P):
                        point_core[tau][point][digit][difference] = (
                            point_core[tau][point][digit][difference]
                            + twisted_point[point][digit][difference]
                        ) % PRIME
            for state in range(V):
                for digit in range(P):
                    for difference in range(P):
                        support_core[tau][state][digit][difference] = (
                            support_core[tau][state][digit][difference]
                            + twisted_support[state][digit][difference]
                        ) % PRIME
        totals = [left + right for left, right in zip(totals, counts)]
        state_totals = [left + right for left, right in zip(state_totals, state_counts)]
        branch_totals = [left + right for left, right in zip(branch_totals, branch_counts)]
    return (
        alpha,
        tuple(gamma_actual), tuple(gamma_support), tuple(gamma_pointed),
        tuple(tuple(tuple(tuple(line) for line in point) for point in tau)
              for tau in point_core),
        tuple(tuple(tuple(tuple(line) for line in state) for state in tau)
              for tau in support_core),
        tuple(totals), tuple(state_totals), tuple(branch_totals),
    )


def build_banks():
    with ProcessPoolExecutor(max_workers=4) as executor:
        chunks = list(executor.map(worker, range(P)))
    chunks.sort(key=lambda item: item[0])
    gamma_actual = tuple(row for chunk in chunks for row in chunk[1])
    gamma_support = tuple(row for chunk in chunks for row in chunk[2])
    gamma_pointed = tuple(row for chunk in chunks for row in chunk[3])
    point_cores = tuple(chunk[4] for chunk in chunks)
    support_cores = tuple(chunk[5] for chunk in chunks)
    counts = tuple(sum(chunk[6][index] for chunk in chunks) for index in range(5))
    state_counts = tuple(sum(chunk[7][index] for chunk in chunks) for index in range(V))
    branch_counts = tuple(sum(chunk[8][index] for chunk in chunks) for index in range(P))
    return (
        gamma_actual, gamma_support, gamma_pointed,
        point_cores, support_cores, counts, state_counts, branch_counts,
    )


def invert_point_cores(cores, zeta: int):
    answer = [[[[0] * P for _difference in range(P)] for _digit in range(P)]
              for _point in POINTS]
    normalizer = pow(P**3, -1, PRIME)
    for alpha in range(P):
        for tau in range(P):
            row = cores[alpha][tau]
            for relation in range(P):
                phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                for point in range(len(POINTS)):
                    for digit in range(P):
                        for difference in range(P):
                            answer[point][digit][difference][relation] = (
                                answer[point][digit][difference][relation]
                                + row[point][digit][difference] * phase
                            ) % PRIME
    return tuple(tuple(tuple(tuple(value * normalizer % PRIME for value in line)
                             for line in block) for block in slab)
                 for slab in answer)


def invert_support_cores(cores, zeta: int):
    answer = [[[[0] * P for _difference in range(P)] for _digit in range(P)]
              for _state in range(V)]
    normalizer = pow(P**3, -1, PRIME)
    for alpha in range(P):
        for tau in range(P):
            row = cores[alpha][tau]
            for relation in range(P):
                phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                for state in range(V):
                    for digit in range(P):
                        for difference in range(P):
                            answer[state][digit][difference][relation] = (
                                answer[state][digit][difference][relation]
                                + row[state][digit][difference] * phase
                            ) % PRIME
    return freeze4(tuple(tuple(tuple(tuple(value * normalizer % PRIME
                                           for value in line)
                                   for line in block) for block in slab)
                         for slab in answer))


def state_marginal(point_tensor):
    return tuple(tuple(tuple(tuple(
        sum(point_tensor[point][digit][difference][relation]
            for point in STATE_FIBRES[state]) % PRIME
        for relation in range(P)) for difference in range(P))
        for digit in range(P)) for state in range(V))


def pointed_marginal(point_tensor):
    return tuple(tuple(tuple(
        sum(point_tensor[point][digit][difference][relation] for digit in range(P))
        % PRIME for relation in range(P)) for difference in range(P))
        for point in range(len(POINTS)))


def branch_marginal(tensor):
    return tuple(tuple(tuple(
        sum(tensor[state][digit][difference][relation] for digit in range(P))
        % PRIME for relation in range(P)) for difference in range(P))
        for state in range(V))


def difference_marginal(tensor):
    return tuple(tuple(tuple(
        sum(tensor[state][digit][difference][relation] for difference in range(P))
        % PRIME for relation in range(P)) for digit in range(P))
        for state in range(V))


def branch_flat(parent):
    inv_p = pow(P, -1, PRIME)
    return tuple(tuple(tuple(tuple(
        parent[state][difference][relation] * inv_p % PRIME
        for relation in range(P)) for difference in range(P))
        for _digit in range(P)) for state in range(V))


def difference_flat(parent):
    inv_p = pow(P, -1, PRIME)
    return tuple(tuple(tuple(tuple(
        parent[state][digit][relation] * inv_p % PRIME
        for relation in range(P)) for _difference in range(P))
        for digit in range(P)) for state in range(V))


def centre4(tensor):
    dims = (V, P, P, P)
    data = [[[[tensor[a][r][s][t] for t in range(P)] for s in range(P)]
             for r in range(P)] for a in range(V)]
    for axis, size in enumerate(dims):
        inverse = pow(size, -1, PRIME)
        old = freeze4(data)
        for a in range(V):
            for r in range(P):
                for s in range(P):
                    for t in range(P):
                        if axis == 0:
                            total = sum(old[x][r][s][t] for x in range(V))
                        elif axis == 1:
                            total = sum(old[a][x][s][t] for x in range(P))
                        elif axis == 2:
                            total = sum(old[a][r][x][t] for x in range(P))
                        else:
                            total = sum(old[a][r][s][x] for x in range(P))
                        data[a][r][s][t] = (
                            old[a][r][s][t] - total * inverse
                        ) % PRIME
    answer = freeze4(data)
    require(all(sum(answer[a][r][s][t] for a in range(V)) % PRIME == 0
                for r in range(P) for s in range(P) for t in range(P)),
            "state centering")
    require(all(sum(answer[a][r][s][t] for r in range(P)) % PRIME == 0
                for a in range(V) for s in range(P) for t in range(P)),
            "branch centering")
    require(all(sum(answer[a][r][s][t] for s in range(P)) % PRIME == 0
                for a in range(V) for r in range(P) for t in range(P)),
            "difference centering")
    require(all(sum(answer[a][r][s][t] for t in range(P)) % PRIME == 0
                for a in range(V) for r in range(P) for s in range(P)),
            "relation centering")
    return answer


def axis_ranks(tensor):
    rows = (
        tuple(tuple(tensor[a][r][s][t] for r in range(P) for s in range(P)
                    for t in range(P)) for a in range(V)),
        tuple(tuple(tensor[a][r][s][t] for a in range(V) for s in range(P)
                    for t in range(P)) for r in range(P)),
        tuple(tuple(tensor[a][r][s][t] for a in range(V) for r in range(P)
                    for t in range(P)) for s in range(P)),
        tuple(tuple(tensor[a][r][s][t] for a in range(V) for r in range(P)
                    for s in range(P)) for t in range(P)),
    )
    return tuple(R.rank_mod(matrix) for matrix in rows)


def bipartition_ranks(tensor):
    matrices = (
        tuple(tuple(tensor[a][r][s][t] for s in range(P) for t in range(P))
              for a in range(V) for r in range(P)),
        tuple(tuple(tensor[a][r][s][t] for r in range(P) for t in range(P))
              for a in range(V) for s in range(P)),
        tuple(tuple(tensor[a][r][s][t] for r in range(P) for s in range(P))
              for a in range(V) for t in range(P)),
    )
    return tuple(R.rank_mod(matrix) for matrix in matrices)


def state_branch_ranks(tensor):
    return tuple(R.rank_mod(tuple(tuple(
        tensor[state][digit][difference][relation]
        for difference in range(P) for relation in range(P)
    ) for digit in range(P))) for state in range(V))


def fourier4(tensor, zeta: int):
    walsh = tuple(tuple(
        -1 if (character[0] * label[0] + character[1] * label[1]) % 2 else 1
        for label in SQ.STATE_LABELS
    ) for character in SQ.STATE_LABELS)
    roots = tuple(tuple(pow(zeta, -frequency * value % P, PRIME)
                        for value in range(P)) for frequency in range(P))
    stage0 = [[[[sum(tensor[a][r][s][t] * walsh[c][a] for a in range(V))
                 % PRIME for t in range(P)] for s in range(P)]
               for r in range(P)] for c in range(V)]
    stage1 = [[[[sum(stage0[c][r][s][t] * roots[k][r] for r in range(P))
                 % PRIME for t in range(P)] for s in range(P)]
               for k in range(P)] for c in range(V)]
    stage2 = [[[[sum(stage1[c][k][s][t] * roots[l][s] for s in range(P))
                 % PRIME for t in range(P)] for l in range(P)]
               for k in range(P)] for c in range(V)]
    stage3 = [[[[sum(stage2[c][k][l][t] * roots[m][t] for t in range(P))
                 % PRIME for m in range(P)] for l in range(P)]
               for k in range(P)] for c in range(V)]
    return freeze4(stage3)


def census4(spectrum):
    bins = [0] * 16
    for c in range(V):
        for k in range(P):
            for l in range(P):
                for m in range(P):
                    if spectrum[c][k][l][m]:
                        mask = ((c != 0) << 3) | ((k != 0) << 2) | ((l != 0) << 1) | (m != 0)
                        bins[mask] += 1
    return (sum(bins),) + tuple(bins)


def centre3(tensor):
    data = [[[tensor[a][r][s] for s in range(P)] for r in range(P)]
            for a in range(V)]
    dims = (V, P, P)
    for axis, size in enumerate(dims):
        inverse = pow(size, -1, PRIME)
        old = freeze3(data)
        for a in range(V):
            for r in range(P):
                for s in range(P):
                    if axis == 0:
                        total = sum(old[x][r][s] for x in range(V))
                    elif axis == 1:
                        total = sum(old[a][x][s] for x in range(P))
                    else:
                        total = sum(old[a][r][x] for x in range(P))
                    data[a][r][s] = (old[a][r][s] - total * inverse) % PRIME
    return freeze3(data)


def ranks3(tensor):
    return (
        R.rank_mod(tuple(tuple(tensor[a][r][s] for r in range(P) for s in range(P))
                         for a in range(V))),
        R.rank_mod(tuple(tuple(tensor[a][r][s] for a in range(V) for s in range(P))
                         for r in range(P))),
        R.rank_mod(tuple(tuple(tensor[a][r][s] for a in range(V) for r in range(P))
                         for s in range(P))),
    )


def fourier3(tensor, zeta):
    walsh = tuple(tuple(
        -1 if (character[0] * label[0] + character[1] * label[1]) % 2 else 1
        for label in SQ.STATE_LABELS
    ) for character in SQ.STATE_LABELS)
    roots = tuple(tuple(pow(zeta, -frequency * value % P, PRIME)
                        for value in range(P)) for frequency in range(P))
    return tuple(tuple(tuple(sum(
        tensor[a][r][s] * walsh[c][a] * roots[k][r] * roots[l][s]
        for a in range(V) for r in range(P) for s in range(P)
    ) % PRIME for l in range(P)) for k in range(P)) for c in range(V))


def census3(spectrum):
    bins = [0] * 8
    for c in range(V):
        for k in range(P):
            for l in range(P):
                if spectrum[c][k][l]:
                    bins[((c != 0) << 2) | ((k != 0) << 1) | (l != 0)] += 1
    return (sum(bins),) + tuple(bins)


def fixed_record(tensor):
    matrix = tuple(tuple(tuple(tensor[a][r][s][6] for s in range(P))
                         for r in range(P)) for a in range(V))
    interaction = centre3(matrix)
    spectrum = fourier3(matrix, pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME))
    interaction_spectrum = fourier3(
        interaction, pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    )
    return (
        ranks3(matrix), ranks3(interaction), census3(spectrum),
        census3(interaction_spectrum),
        sum(value != 0 for slab in matrix for row in slab for value in row),
        digest_json(matrix), digest_json(interaction),
    )


def solve_square(matrix, right):
    n = len(matrix)
    work = [list(matrix[row]) + [right[row] % PRIME] for row in range(n)]
    for column in range(n):
        pivot = next((row for row in range(column, n)
                      if work[row][column] % PRIME), None)
        require(pivot is not None, ("singular coordinate minor", matrix))
        work[column], work[pivot] = work[pivot], work[column]
        inverse = pow(work[column][column], -1, PRIME)
        work[column] = [value * inverse % PRIME for value in work[column]]
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [(left - factor * right_value) % PRIME
                             for left, right_value in zip(work[row], work[column])]
    return tuple(work[row][-1] for row in range(n))


def basis_coordinates(basis, rows):
    rank = len(basis)
    width = len(basis[0])
    pivots = None
    for columns in itertools.combinations(range(width), rank):
        minor = tuple(tuple(basis[row][column] for row in range(rank))
                      for column in columns)
        if R.rank_mod(minor) == rank:
            pivots = columns
            break
    require(pivots is not None, "no coordinate pivot")
    matrix = tuple(tuple(basis[row][column] for row in range(rank))
                   for column in pivots)
    coordinates = []
    for row in rows:
        coeffs = solve_square(matrix, tuple(row[column] for column in pivots))
        reconstructed = tuple(sum(coeffs[index] * basis[index][column]
                                  for index in range(rank)) % PRIME
                              for column in range(width))
        require(reconstructed == row, ("coordinate reconstruction", pivots))
        coordinates.append(coeffs)
    return tuple(pivots), tuple(coordinates)


def carrier_records(actual, pointed, point_branch):
    state_records = []
    pointed_digit_records = []
    coordinate_records = []
    all_base = []
    all_digit = []
    for state in range(V):
        branch_rows = tuple(tuple(actual[state][digit][s][t]
                                  for s in range(P) for t in range(P))
                            for digit in range(P))
        base_rows = tuple(tuple(pointed[point][s][t]
                                for s in range(P) for t in range(P))
                          for point in STATE_FIBRES[state])
        digit_rows = tuple(tuple(point_branch[point][digit][s][t]
                                 for s in range(P) for t in range(P))
                           for point in STATE_FIBRES[state] for digit in range(P))
        state_record = (
            len(base_rows), R.rank_mod(branch_rows), R.rank_mod(base_rows),
            R.rank_mod(branch_rows + base_rows),
        )
        digit_record = (
            len(base_rows), R.rank_mod(base_rows), R.rank_mod(digit_rows),
            R.rank_mod(base_rows + digit_rows),
        )
        state_records.append(state_record)
        pointed_digit_records.append(digit_record)
        coordinate_records.append(basis_coordinates(base_rows, branch_rows))
        all_base.extend(base_rows)
        all_digit.extend(digit_rows)
    global_record = (
        R.rank_mod(tuple(all_base)), R.rank_mod(tuple(all_digit)),
        R.rank_mod(tuple(all_base + all_digit)),
    )
    return (
        tuple(state_records), tuple(pointed_digit_records), global_record,
        digest_json(tuple(coordinate_records)), tuple(coordinate_records),
    )


def main() -> None:
    R.split_field_certificate()
    (
        gamma_actual, gamma_support, gamma_pointed,
        point_cores, support_cores, work_counts, state_counts, branch_counts,
    ) = build_banks()
    require(len(gamma_actual) == len(gamma_support) == len(gamma_pointed) == P**3,
            "gamma size")
    require(work_counts == EXPECTED_WORK_COUNTS, ("work counts", work_counts))
    require(state_counts == EXPECTED_STATE_COUNTS, ("state counts", state_counts))
    require(branch_counts == EXPECTED_BRANCH_COUNTS, ("branch counts", branch_counts))

    gamma_branch_actual = tuple(tuple(tuple(
        sum(row[state][digit][difference] for digit in range(P)) % PRIME
        for difference in range(P)) for state in range(V)) for row in gamma_actual)
    gamma_branch_support = tuple(tuple(tuple(
        sum(row[state][digit][difference] for digit in range(P)) % PRIME
        for difference in range(P)) for state in range(V)) for row in gamma_support)
    gamma_difference_actual = tuple(tuple(tuple(
        sum(row[state][digit][difference] for difference in range(P)) % PRIME
        for digit in range(P)) for state in range(V)) for row in gamma_actual)
    require(digest_json(gamma_branch_actual) == RD.EXPECTED_DIGESTS[0][0],
            "gamma branch marginal misses weighted root difference")
    require(digest_json(gamma_branch_support) == RD.EXPECTED_DIGESTS[0][1],
            "gamma branch marginal misses support root difference")
    require(digest_json(gamma_difference_actual) == D.EXPECTED_DIGESTS[0],
            "gamma difference marginal misses current digit")
    require(digest_json(gamma_pointed) == PS.EXPECTED_DIGESTS[0][0],
            "gamma pointed marginal")

    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    point_branch = invert_point_cores(point_cores, zeta)
    actual = state_marginal(point_branch)
    support = invert_support_cores(support_cores, zeta)
    pointed = pointed_marginal(point_branch)
    weighted_root = branch_marginal(actual)
    support_root = branch_marginal(support)
    current = difference_marginal(actual)
    require(digest_json(weighted_root) == RD.EXPECTED_DIGESTS[1][0],
            "tensor branch marginal misses weighted root difference")
    require(digest_json(support_root) == RD.EXPECTED_DIGESTS[1][1],
            "tensor branch marginal misses support root difference")
    require(digest_json(current) == D.EXPECTED_DIGESTS[1],
            "tensor difference marginal misses current digit")
    require(digest_json(pointed) == PS.EXPECTED_DIGESTS[1][0],
            "tensor pointed marginal")
    require(all(actual[state][digit][0][relation] == 0
                and support[state][digit][0][relation] == 0
                for state in range(V) for digit in range(P)
                for relation in range(P)), "same-root slice")

    flat_branch = branch_flat(weighted_root)
    flat_difference = difference_flat(current)
    require(branch_marginal(flat_branch) == weighted_root, "branch-flat parent")
    require(difference_marginal(flat_difference) == current, "difference-flat parent")
    tensors = (actual, support, flat_branch, flat_difference)

    guard_records = []
    for alpha, beta, tau in CONTROL_TRIPLES:
        word, endpoint_grid, _harmonic, danger = D.endpoint_context()
        branches = D.owner_digit_profiles()
        R.literal_guard_restoration(
            word, endpoint_grid, alpha, beta, tau, set(branches[1]), danger
        )
        local_actual, local_support, _counts, _states, _branches = integrate_pair(
            alpha, beta, literal_tau=tau
        )
        phase = pow(zeta, beta, PRIME)
        index = (alpha * P + beta) * P + tau
        twisted_point = tuple(tuple(tuple(
            phase * local_actual[0][point][digit][difference] % PRIME
            for difference in range(P)) for digit in range(P))
            for point in range(len(POINTS)))
        twisted_support = tuple(tuple(tuple(
            phase * local_support[0][state][digit][difference] % PRIME
            for difference in range(P)) for digit in range(P))
            for state in range(V))
        require(state_row(twisted_point) == gamma_actual[index],
                ("literal actual guard", alpha, beta, tau))
        require(twisted_support == gamma_support[index],
                ("literal support guard", alpha, beta, tau))
        require(pointed_row(twisted_point) == gamma_pointed[index],
                ("literal pointed guard", alpha, beta, tau))
        guard_records.append((alpha, beta, tau))

    axis = tuple(axis_ranks(tensor) for tensor in tensors)
    bipartitions = tuple(bipartition_ranks(tensor) for tensor in tensors)
    state_blocks = tuple(state_branch_ranks(tensor) for tensor in tensors)
    interactions = tuple(centre4(tensor) for tensor in tensors)
    interaction_axis = tuple(axis_ranks(tensor) for tensor in interactions)
    interaction_bipartitions = tuple(bipartition_ranks(tensor) for tensor in interactions)
    require((axis, bipartitions, state_blocks, interaction_axis,
             interaction_bipartitions) == EXPECTED_RANKS,
            ("rank record", axis, bipartitions, state_blocks,
             interaction_axis, interaction_bipartitions))

    carrier = carrier_records(actual, pointed, point_branch)
    require(carrier[0] == EXPECTED_CARRIER, ("state carrier", carrier[0]))
    require(carrier[1] == EXPECTED_CARRIER, ("pointed digit carrier", carrier[1]))
    require(carrier[2] == (6, 6, 6), ("global pointed digit carrier", carrier[2]))

    spectra = tuple(fourier4(tensor, zeta) for tensor in tensors)
    interaction_spectra = tuple(fourier4(tensor, zeta) for tensor in interactions)
    censuses = tuple(census4(spectrum) for spectrum in spectra)
    interaction_censuses = tuple(census4(spectrum) for spectrum in interaction_spectra)
    require((censuses, interaction_censuses) == EXPECTED_CENSUSES,
            ("spectral censuses", censuses, interaction_censuses))

    fixed = tuple(fixed_record(tensor) for tensor in tensors)
    digests = (
        (digest_json(gamma_actual), digest_json(gamma_support), digest_json(gamma_pointed)),
        tuple(digest_json(tensor) for tensor in tensors),
        digest_json(pointed),
        tuple(digest_json(tensor) for tensor in interactions),
        tuple(digest_json(spectrum) for spectrum in spectra),
        tuple(digest_json(spectrum) for spectrum in interaction_spectra),
    )
    require(digests == EXPECTED_DIGESTS, ("candidate digest comparison", digests))

    # The two newer packages distinguish source/profile amplitude from the
    # relation carrier.  These are comparison records, not dependencies for
    # the reconstruction above.
    amplitude_carrier = (
        ("source_current", 17, 4, 6, 6, SOURCE_CURRENT_SEMANTIC),
        ("one_current", P, axis[0][1], carrier[2][0], carrier[2][2]),
        ("two_current", 13, 12, 6, 6, TWO_CURRENT_SEMANTIC),
    )
    record = (
        CURRENT_PARENT_SHA256, CURRENT_PARENT_SEMANTIC,
        POINTED_PARENT_SHA256, POINTED_PARENT_SEMANTIC,
        TWO_CURRENT_SHA256, SOURCE_CURRENT_SHA256,
        work_counts, state_counts, branch_counts, tuple(guard_records),
        axis, bipartitions, state_blocks, interaction_axis,
        interaction_bipartitions, carrier[:4], censuses,
        interaction_censuses, fixed, digests, amplitude_carrier,
        "static finite response; no address/clock/current/row/LRC",
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== independent hostile audit: current digit x source-root difference x pointed carrier ==")
    print(f"parents=(current_sha={CURRENT_PARENT_SHA256},semantic={CURRENT_PARENT_SEMANTIC},pointed_sha={POINTED_PARENT_SHA256},semantic={POINTED_PARENT_SEMANTIC})")
    print(f"comparison_packages=(source_current_sha={SOURCE_CURRENT_SHA256},semantic={SOURCE_CURRENT_SEMANTIC},two_current_sha={TWO_CURRENT_SHA256},semantic={TWO_CURRENT_SEMANTIC})")
    print(f"coordinates=(point={POINTS},r=a_mod13,s=u-q,relation=(1,0,t));field=(prime={PRIME},zeta13={zeta})")
    print(f"work_counts={work_counts};state_counts={state_counts};branch_counts={branch_counts};literal_controls={tuple(guard_records)}: PASS")
    print("marginals=(sum_r->weighted/support_root_difference,sum_s->independent_current_digit,sum_r->independent_pointed_root): SEGMENTWISE+GAMMA+INVERSE PASS")
    print("same_root_s0=(actual,support)=0;hostiles=(support_normalized,branch_flat,difference_flat): PASS")
    print(f"axis_ranks_(state,branch,difference,relation)={axis}")
    print(f"bipartition_ranks_(state_branch|difference_relation,state_difference|branch_relation,state_relation|branch_difference)={bipartitions}")
    print(f"fixed_state_branch_block_ranks={state_blocks}")
    print(f"four_way_ANOVA_axis_ranks={interaction_axis}")
    print(f"four_way_ANOVA_bipartition_ranks={interaction_bipartitions}")
    print(f"statewise_pointed_equality_(point_count,branch_rank,point_rank,union_rank)={carrier[0]}: PASS")
    print(f"stronger_pointed_digit_carrier_(point_count,base_rank,digit_rank,union_rank)={carrier[1]};global={carrier[2]}: PASS")
    print(f"pointed_coordinate_digest={carrier[3]}")
    print(f"spectral_censuses_(actual,support,branch_flat,difference_flat)={censuses}")
    print(f"four_way_ANOVA_spectral_censuses={interaction_censuses}")
    print(f"fixed_(1,0,6)_records={fixed}")
    print(f"digests=(gamma,tensors,pointed,ANOVA,spectra,ANOVA_spectra)={digests}")
    print(f"amplitude_vs_carrier={amplitude_carrier}")
    print("interpretation=branch_amplitude_rank4 lives inside pointed_relation_carrier_rank6;whole_tensor_bipartition_ranks40/42 remain distinct")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT scoped four-way interaction and exact statewise pointed-six factorization")
    print("scope=not complete address,not chronology,not physical current,no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
