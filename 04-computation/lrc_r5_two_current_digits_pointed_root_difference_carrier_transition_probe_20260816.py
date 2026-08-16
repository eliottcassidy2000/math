#!/usr/bin/env python3
"""Retain two inverse-cylinder digits on the pointed-six root-difference bundle.

The exact pre-integration axes are

    point=(Boolean state, absolute U-tail root) in the audited six points,
    r0=a mod 13,
    r1=floor(a/13) mod 13,
    s=u-q mod 13,
    relation=(1,0,t).

The endpoint sweep is rebuilt from the two-digit source profiles.  It reopens
each selected ordered root pair (u,q) before multiplying its endpoint jump;
the one-digit current/root-difference candidate is not used as an aggregation
oracle.  Two mandatory quotient maps are checked before inversion and after
inversion:

  * sum over s and over pointed tails in a fixed state gives the audited
    state x r0 x r1 two-digit tensor;
  * sum over r1 and over pointed tails in a fixed state gives the existing
    state x r0 x s one-digit/root-difference tensor.

The stronger sum over r0,r1 is also checked against the independently audited
six-pointed root-difference parent.  The decisive linear-algebra test asks for
one matrix K[r1] in Mat_6(F) satisfying, simultaneously for every r0,s,t,

    child[:,r0,r1,s,t] = K[r1] parent[:,r0,s,t].

State-block preservation is tested separately from unrestricted six-carrier
closure.  A support-normalized hostile preserves every one-digit pointed row,
while an r1-flat hostile has K[r1]=13^-1 I.  None of these finite response
coordinates is interpreted as chronology, a complete address, an arrival
class, a physical current, a row exclusion, or LRC(14).
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TWO_DIGIT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_probe_20260816.py"
)
TWO_DIGIT_AUDIT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_independent_audit_20260816.py"
)
ONE_DIGIT_ROOT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_owner_node_inverse_branch_root_difference_four_way_probe_20260816.py"
)
POINTED_AUDIT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_ufull_owner_node_pointed_six_state_root_difference_independent_audit_20260816.py"
)
NESTED_SOURCE_CURRENT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_nested_ancestry_digits_probe_20260816.py"
)
TWO_DIGIT_SHA256 = "3dab580e479e4ba7ac8801c1e5d8523018e0b3dc1c2176c072e7c609033eb6c8"
TWO_DIGIT_AUDIT_SHA256 = (
    "126be106a34a990f22e66b26d79ea3568ee4f394419936ed47f1bfbe0656788f"
)
TWO_DIGIT_AUDIT_SEMANTIC = (
    "f2af726e9d5abd1487e841623ce2f62ca647c86e5a1a68e41eda4d9dda6c81ac"
)
ONE_DIGIT_ROOT_SHA256 = "d66378cc8db99c4de087fa78413721de6aa7be3960f51ee39f29feae9313eeba"
POINTED_AUDIT_SHA256 = "8c7cb5f98b15a768d4f4d6060074e0815a8f089f857ec4f3c55a0e7d877e1fec"
POINTED_AUDIT_SEMANTIC = "66db2301f88db1ced7784868095e198e3e12f1fb79175ebd902d0f569a5decef"
NESTED_SOURCE_CURRENT_SHA256 = (
    "1188df8aa2a7a84c1e8ada5fc3cc8d3b839ece70298b94f1d94c9d440caa88f3"
)
NESTED_SOURCE_CURRENT_SEMANTIC = (
    "6e5605f58b7a94ea5ea4e8f62cfa7ee135b0d52512225f4aaee248ad6e21a9ae"
)
# (source-profile rank, output-amplitude rank, pointed carrier rank, union rank)
NESTED_AMPLITUDE_CARRIER_CONTROL = (17, 4, 6, 6)

P = 13
V = 4
DIRECT_CONTROLS = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
INVERSION_PROBES = (
    (0, 0, 1, 1),
    (1, 6, 7, 6),
    (3, 4, 9, 7),
    (5, 12, 11, 12),
)
PROBE_RELATIONS = (0, 6, 12)

# Pinned after the first normal exact replay.
EXPECTED_JOINT_ALPHA_MERKLE_SHA256 = (
    "b902b40e2ef44a7e75b107bb28b446aca3aaae5442dd82a00381d1962e4d77a8"
)
EXPECTED_TENSOR_DIGESTS = (
    (
        "f08ced17b1f727c8032a692e59df174de14ed06a9bc821e72788c8f347b28986",
        "6aa28175baa458ae8646e03b4b23234febb3a5ef60714ec7e35d623e4d2267a0",
        "e4f96fe95854d1f78649c9bbe71c91a5276ab29813934c90e11c6cd616bc7964",
        "e36cedbce618771be7ac5e94d74fa6be948ad8236c3145222877ee274a6b5412",
    ),
    (
        "10a70010c5d305df4f2e24e265f8d811d2080c25b12d8693c6ce5086bb670dab",
        "218b0a8e967649c238963bf4d8d3bfe076b74dc31d45d96f15fbc0e7a969915e",
    ),
)
EXPECTED_RANK_RECORD_SHA256 = (
    "05e3f52e72b0e81fcc7352cbf7dc2b2ad99b4159e1e30afaeefc9a594efe48f0"
)
EXPECTED_TRANSITION_RECORD_SHA256 = (
    "26f55dfca08659f36d6c13729847348e56f082c217443f645562076015dbf19d"
)
EXPECTED_WORK_RECORD_SHA256 = (
    "8f1703d00e79cad05a8c5e3809053597aed193f9de2f72ed472f560051ce69da"
)
EXPECTED_SEMANTIC_SHA256 = (
    "38725dc1d7129b326634c99bd70e1eb414590dc24fb83bd9522e2095e41f204c"
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


T2 = load_module(TWO_DIGIT_PATH, TWO_DIGIT_SHA256, "pointed_two_digit_parent")
J1 = load_module(ONE_DIGIT_ROOT_PATH, ONE_DIGIT_ROOT_SHA256, "pointed_one_digit_root_parent")
require(lf_sha256(TWO_DIGIT_AUDIT_PATH) == TWO_DIGIT_AUDIT_SHA256,
        "independent two-digit audit source drift")
require(lf_sha256(POINTED_AUDIT_PATH) == POINTED_AUDIT_SHA256,
        "independent pointed audit source drift")
require(lf_sha256(NESTED_SOURCE_CURRENT_PATH) == NESTED_SOURCE_CURRENT_SHA256,
        "nested source/current amplitude-carrier control drift")

C = T2.C
B = T2.B
PRIME = T2.PRIME
POINTS = J1.POINTED_STATES
POINT_INDEX = J1.POINT_INDEX
STATE_POINTS = J1.STATE_POINTS

require(POINTS == ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12)),
        ("pointed carrier order", POINTS))
require(PRIME == J1.PRIME, "finite-field mismatch")
require(T2.EXPECTED_SEMANTIC_SHA256
        == "61743457afc0cff984c87affa7f2e67bf3a21e08a401ea69b679319f2f51e826",
        "two-digit semantic drift")
require(J1.POINTED_GAMMA_SHA256
        == "416f03a90894e526bc767a34839c2489aa4cac7051ac04687e0feacfa36d58d2",
        "pointed gamma namespace drift")
require(J1.POINTED_TENSOR_SHA256
        == "9c5e227d9c142373973a562a54c6a67cac60a82da1121a028b9658920d155a19",
        "pointed tensor namespace drift")


def freeze5(bank):
    return tuple(
        tuple(
            tuple(tuple(tuple(line) for line in sheet) for sheet in block)
            for block in slab
        )
        for slab in bank
    )


def freeze_joint_row(row):
    return tuple(
        tuple(tuple(tuple(line) for line in block) for block in slab)
        for slab in row
    )


def source_support_record(cells):
    counts = set()
    patterns = set()
    for state, u_values, _v_values, _mask, address_values in cells:
        if state is None:
            continue
        for root, value in enumerate(u_values):
            if not value:
                continue
            pattern = tuple(
                sum(address_values[root][r0][r1] != 0 for r1 in range(P))
                for r0 in range(P)
            )
            counts.update(pattern)
            patterns.add(pattern)
    return tuple(sorted(counts)), tuple(sorted(patterns))


def source_ratio_kernels(cells):
    """Test source-cell proportionality before endpoint integration.

    A kernel entry is retained only when child/parent is constant over every
    source cell carrying the pointed tail.  The support-normalized profile is
    expected to pass everywhere; the actual amplitudes are allowed to expose
    exceptional boundary fibres, which are recorded rather than suppressed.
    """

    kernels = []
    records = []
    for normalized in (False, True):
        kernel = [[[None for _point in POINTS] for _r1 in range(P)]
                  for _r0 in range(P)]
        size_histogram = {}
        exceptions = []
        for point, (state, root) in enumerate(POINTS):
            for r0 in range(P):
                for r1 in range(P):
                    ratios = set()
                    for cell_state, u_values, _v_values, _mask, addresses in cells:
                        if cell_state != state or not u_values[root]:
                            continue
                        row = addresses[root][r0]
                        total = sum(row) % PRIME
                        if not total:
                            continue
                        if normalized:
                            active = tuple(index for index, value in enumerate(row)
                                           if value)
                            child = (total * pow(len(active), -1, PRIME) % PRIME
                                     if r1 in active else 0)
                        else:
                            child = row[r1]
                        ratios.add(child * pow(total, -1, PRIME) % PRIME)
                    size_histogram[len(ratios)] = size_histogram.get(len(ratios), 0) + 1
                    if len(ratios) == 1:
                        kernel[r0][r1][point] = next(iter(ratios))
                    else:
                        exceptions.append((point, POINTS[point], r0, r1,
                                           len(ratios), digest_json(tuple(sorted(ratios)))))
        kernels.append(tuple(tuple(tuple(row) for row in plane) for plane in kernel))
        grouped = {}
        for point, pair, r0, r1, size, _ratio_digest in exceptions:
            grouped.setdefault((point, pair, r0), []).append((r1, size))
        group_record = tuple(
            (key, tuple(values)) for key, values in sorted(grouped.items())
        )
        records.append((
            tuple(sorted(size_histogram.items())),
            len(exceptions),
            group_record,
            digest_json(tuple(exceptions)),
            digest_json(kernels[-1]),
        ))
    return tuple(kernels), tuple(records)


def expand_coefficients(coefficients, cells):
    actual = [[[[[0 for _s in range(P)] for _r1 in range(P)]
                for _r0 in range(P)] for _point in POINTS]
              for _tau in coefficients]
    support = [[[[[0 for _s in range(P)] for _r1 in range(P)]
                 for _r0 in range(P)] for _point in POINTS]
               for _tau in coefficients]

    for tau_index, coefficient_map in enumerate(coefficients):
        for (cell_index, root, difference), scalar in coefficient_map.items():
            if scalar == 0:
                continue
            state, _u_values, _v_values, _mask, address_values = cells[cell_index]
            require(state is not None, ("untyped joint source cell", cell_index))
            point = POINT_INDEX.get((state, root))
            require(point is not None,
                    ("unsupported pointed carrier", cell_index, state, root))
            require(difference != 0, ("same-root coefficient", cell_index, root))
            for r0 in range(P):
                address_row = address_values[root][r0]
                total = sum(address_row) % PRIME
                active = tuple(r1 for r1, value in enumerate(address_row) if value)
                if not total:
                    require(not active, ("zero address total with support", cell_index,
                                         root, r0, active))
                    continue
                require(active, ("empty r1 support", cell_index, root, r0))
                normalized = total * pow(len(active), -1, PRIME) % PRIME
                for r1 in active:
                    actual[tau_index][point][r0][r1][difference] = (
                        actual[tau_index][point][r0][r1][difference]
                        + address_row[r1] * scalar
                    ) % PRIME
                    support[tau_index][point][r0][r1][difference] = (
                        support[tau_index][point][r0][r1][difference]
                        + normalized * scalar
                    ) % PRIME

    frozen = (freeze5(actual), freeze5(support))
    for tau in range(len(coefficients)):
        for point in range(len(POINTS)):
            for r0 in range(P):
                for difference in range(P):
                    left = sum(frozen[0][tau][point][r0][r1][difference]
                               for r1 in range(P)) % PRIME
                    right = sum(frozen[1][tau][point][r0][r1][difference]
                                for r1 in range(P)) % PRIME
                    require(left == right,
                            ("support-normalized parent", tau, point, r0,
                             difference, left, right))
    return frozen


def integrate_joint(alpha: int, beta: int, literal_tau: int | None = None):
    _profiles, boundaries, cells, _source_record = T2.two_digit_context()
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    for boundary in boundaries:
        events.setdefault(boundary, 0)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coefficients = [dict() for _tau in tau_values]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    state_segments = [0] * V
    coefficient_updates = 0

    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        require(C.cell_of_segment(left, right) == 0,
                ("pointed two-digit escaped cell zero", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("pointed two-digit middle chamber", alpha, beta, left, right))
        cell_index = bisect_right(boundaries, left) - 1
        require(0 <= cell_index < len(cells)
                and right <= boundaries[cell_index + 1],
                ("endpoint crosses source cell", left, right, cell_index))
        state, u_values, v_values, _u_mask, _addresses = cells[cell_index]
        require(state is not None, ("active untyped cell", left, right))
        active_segments += 1
        state_segments[state] += 1

        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_support = tuple(root for root, value in enumerate(u_values) if value)
        v_support = tuple(root for root, value in enumerate(v_values) if value)
        if not u_support or not v_support:
            continue
        weighted_segments += 1
        require(set(u_support).isdisjoint(v_support),
                ("source supports overlap", left, right, u_support, v_support))

        for tau_index, tau in enumerate(tau_values):
            selected_mask = mask if literal_tau is not None else (
                mask & T2.SAFE_MASKS[chamber][tau]
            )
            selected_u = tuple(root for root in u_support
                               if (selected_mask >> root) & 1)
            selected_v = tuple(root for root in v_support
                               if (selected_mask >> root) & 1)
            for root in selected_u:
                require((state, root) in POINT_INDEX,
                        ("unrealized pointed source tail", left, state, root))
                for q_root in selected_v:
                    difference = (root - q_root) % P
                    require(difference != 0,
                            ("pointed same-root pair", left, root, q_root))
                    key = (cell_index, root, difference)
                    coefficients[tau_index][key] = (
                        coefficients[tau_index].get(key, 0)
                        + v_values[q_root] * jump
                    ) % PRIME
                    coefficient_updates += 1

    mask ^= events[positions[-1]]
    require(mask == 0, ("pointed two-digit endpoint mask", alpha, beta,
                        literal_tau, mask))
    banks = expand_coefficients(coefficients, cells)
    counts = (
        interval_count,
        mapped,
        active_segments,
        q_active_segments,
        weighted_segments,
        tuple(state_segments),
        coefficient_updates,
        tuple(len(row) for row in coefficients),
    )
    return banks, counts


def twist_joint_row(row, phase: int):
    return tuple(
        tuple(
            tuple(
                tuple(phase * value % PRIME for value in difference_row)
                for difference_row in r1_rows
            )
            for r1_rows in r0_rows
        )
        for r0_rows in row
    )


def gamma_two_digit_marginal(row):
    return tuple(
        tuple(
            tuple(
                sum(row[point][r0][r1][difference]
                    for point in STATE_POINTS[state] for difference in range(P))
                % PRIME
                for r1 in range(P)
            )
            for r0 in range(P)
        )
        for state in range(V)
    )


def gamma_one_digit_pointed(row):
    return tuple(
        tuple(
            tuple(sum(row[point][r0][r1][difference] for r1 in range(P)) % PRIME
                  for difference in range(P))
            for r0 in range(P)
        )
        for point in range(len(POINTS))
    )


def gamma_one_digit_state(row):
    pointed = gamma_one_digit_pointed(row)
    return tuple(
        tuple(
            tuple(sum(pointed[point][r0][difference]
                      for point in STATE_POINTS[state]) % PRIME
                  for difference in range(P))
            for r0 in range(P)
        )
        for state in range(V)
    )


def gamma_pointed_parent(row):
    return tuple(
        tuple(sum(row[point][r0][r1][difference]
                  for r0 in range(P) for r1 in range(P)) % PRIME
              for difference in range(P))
        for point in range(len(POINTS))
    )


def worker(alpha: int):
    zeta = C.context()["zeta"]
    tau_actual = [[[[[0 for _s in range(P)] for _r1 in range(P)]
                    for _r0 in range(P)] for _point in POINTS]
                  for _tau in range(P)]
    tau_support = [[[[[0 for _s in range(P)] for _r1 in range(P)]
                     for _r0 in range(P)] for _point in POINTS]
                   for _tau in range(P)]
    two_digit_alpha_rows = []
    one_digit_state_rows = []
    pointed_parent_rows = []
    actual_row_digests = []
    support_row_digests = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    update_count = 0
    key_counts = [0] * P
    control_rows = {}
    probe_rows = [[] for _probe in INVERSION_PROBES]

    for beta in range(P):
        banks, counts = integrate_joint(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        beta_two_digit = []
        for tau in range(P):
            actual = twist_joint_row(banks[0][tau], phase)
            support = twist_joint_row(banks[1][tau], phase)
            actual_row_digests.append(digest_json(actual))
            support_row_digests.append(digest_json(support))

            actual_parent = gamma_one_digit_pointed(actual)
            support_parent = gamma_one_digit_pointed(support)
            require(actual_parent == support_parent,
                    ("support gamma parent", alpha, beta, tau))
            beta_two_digit.append(gamma_two_digit_marginal(actual))
            one_digit_state_rows.append(tuple(
                tuple(
                    tuple(sum(actual_parent[point][r0][difference]
                              for point in STATE_POINTS[state]) % PRIME
                          for difference in range(P))
                    for r0 in range(P)
                )
                for state in range(V)
            ))
            pointed_parent_rows.append(gamma_pointed_parent(actual))

            if (alpha, beta, tau) in DIRECT_CONTROLS:
                control_rows[(alpha, beta, tau)] = (actual, support)
            for probe_index, (point, r0, r1, difference) in enumerate(INVERSION_PROBES):
                probe_rows[probe_index].append(actual[point][r0][r1][difference])

            for point in range(len(POINTS)):
                for r0 in range(P):
                    for r1 in range(P):
                        for difference in range(P):
                            tau_actual[tau][point][r0][r1][difference] = (
                                tau_actual[tau][point][r0][r1][difference]
                                + actual[point][r0][r1][difference]
                            ) % PRIME
                            tau_support[tau][point][r0][r1][difference] = (
                                tau_support[tau][point][r0][r1][difference]
                                + support[point][r0][r1][difference]
                            ) % PRIME
        two_digit_alpha_rows.append(tuple(beta_two_digit))
        scalar_counts = [left + right for left, right in zip(scalar_counts, counts[:5])]
        state_counts = [left + right for left, right in zip(state_counts, counts[5])]
        update_count += counts[6]
        key_counts = [left + right for left, right in zip(key_counts, counts[7])]

    require(digest_json(two_digit_alpha_rows) == T2.EXPECTED_ALPHA_GAMMA_SHA256[alpha],
            ("two-digit pre-integration marginal", alpha,
             digest_json(two_digit_alpha_rows)))

    alpha_phase = pow(zeta, -alpha % P, PRIME)
    for bank in (tau_actual, tau_support):
        for tau in range(P):
            for point in range(len(POINTS)):
                for r0 in range(P):
                    for r1 in range(P):
                        for difference in range(P):
                            bank[tau][point][r0][r1][difference] = (
                                bank[tau][point][r0][r1][difference] * alpha_phase
                            ) % PRIME

    return (
        alpha,
        freeze5(tau_actual),
        freeze5(tau_support),
        tuple(one_digit_state_rows),
        tuple(pointed_parent_rows),
        (
            digest_json(tuple(actual_row_digests)),
            digest_json(tuple(support_row_digests)),
        ),
        (tuple(scalar_counts), tuple(state_counts), update_count, tuple(key_counts)),
        tuple(sorted(control_rows.items())),
        tuple((INVERSION_PROBES[index], tuple(rows))
              for index, rows in enumerate(probe_rows)),
    )


def add_core(target, source):
    for tau in range(P):
        for point in range(len(POINTS)):
            for r0 in range(P):
                for r1 in range(P):
                    for difference in range(P):
                        target[tau][point][r0][r1][difference] = (
                            target[tau][point][r0][r1][difference]
                            + source[tau][point][r0][r1][difference]
                        ) % PRIME


def inverse_core(tau_core, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    roots = tuple(
        tuple(pow(zeta, -tau * relation % P, PRIME) for tau in range(P))
        for relation in range(P)
    )
    tensor = [[[[[0 for _t in range(P)] for _s in range(P)]
                for _r1 in range(P)] for _r0 in range(P)]
              for _point in POINTS]
    for point in range(len(POINTS)):
        for r0 in range(P):
            for r1 in range(P):
                for difference in range(P):
                    values = tuple(tau_core[tau][point][r0][r1][difference]
                                   for tau in range(P))
                    for relation in range(P):
                        tensor[point][r0][r1][difference][relation] = (
                            sum(values[tau] * roots[relation][tau] for tau in range(P))
                            * normalizer
                        ) % PRIME
    return freeze5(tensor)


def tensor_two_digit_marginal(tensor):
    return tuple(
        tuple(
            tuple(
                tuple(
                    sum(tensor[point][r0][r1][difference][relation]
                        for point in STATE_POINTS[state] for difference in range(P))
                    % PRIME
                    for relation in range(P)
                )
                for r1 in range(P)
            )
            for r0 in range(P)
        )
        for state in range(V)
    )


def tensor_one_digit_pointed(tensor):
    return tuple(
        tuple(
            tuple(
                tuple(sum(tensor[point][r0][r1][difference][relation]
                          for r1 in range(P)) % PRIME
                      for relation in range(P))
                for difference in range(P)
            )
            for r0 in range(P)
        )
        for point in range(len(POINTS))
    )


def tensor_one_digit_state(pointed):
    return tuple(
        tuple(
            tuple(
                tuple(sum(pointed[point][r0][difference][relation]
                          for point in STATE_POINTS[state]) % PRIME
                      for relation in range(P))
                for difference in range(P)
            )
            for r0 in range(P)
        )
        for state in range(V)
    )


def tensor_pointed_parent(pointed):
    return tuple(
        tuple(
            tuple(sum(pointed[point][r0][difference][relation]
                      for r0 in range(P)) % PRIME
                  for relation in range(P))
            for difference in range(P)
        )
        for point in range(len(POINTS))
    )


def r1_flat_lift(parent):
    inverse = pow(P, -1, PRIME)
    return tuple(
        tuple(
            tuple(
                tuple(
                    tuple(parent[point][r0][difference][relation] * inverse % PRIME
                          for relation in range(P))
                    for difference in range(P)
                )
                for _r1 in range(P)
            )
            for r0 in range(P)
        )
        for point in range(len(POINTS))
    )


def point_difference_flat_lift(two_digit):
    answer = [[[[[0 for _t in range(P)] for _s in range(P)]
                for _r1 in range(P)] for _r0 in range(P)]
              for _point in POINTS]
    for state in range(V):
        scale = pow(len(STATE_POINTS[state]) * P, -1, PRIME)
        for point in STATE_POINTS[state]:
            for r0 in range(P):
                for r1 in range(P):
                    for difference in range(P):
                        for relation in range(P):
                            answer[point][r0][r1][difference][relation] = (
                                two_digit[state][r0][r1][relation] * scale
                            ) % PRIME
    return freeze5(answer)


def rank_rows(matrix) -> int:
    if not matrix:
        return 0
    columns = len(matrix[0])
    basis = {}
    for source in matrix:
        row = [value % PRIME for value in source]
        require(len(row) == columns, "ragged rank matrix")
        for pivot, pivot_row in basis.items():
            factor = row[pivot]
            if factor:
                row[pivot:] = [
                    (value - factor * base) % PRIME
                    for value, base in zip(row[pivot:], pivot_row[pivot:])
                ]
        pivot = next((column for column, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], -1, PRIME)
        row[pivot:] = [value * inverse % PRIME for value in row[pivot:]]
        basis[pivot] = row
        if len(basis) == columns:
            break
    return len(basis)


def response_rows_base(base, points=None):
    points = tuple(range(len(POINTS))) if points is None else points
    return tuple(
        tuple(base[point][difference][relation]
              for difference in range(P) for relation in range(P))
        for point in points
    )


def response_rows_parent(parent, points=None, fixed_r0=None):
    points = tuple(range(len(POINTS))) if points is None else points
    r0_values = range(P) if fixed_r0 is None else (fixed_r0,)
    return tuple(
        tuple(parent[point][r0][difference][relation]
              for difference in range(P) for relation in range(P))
        for point in points for r0 in r0_values
    )


def response_rows_child(tensor, points=None, fixed_r0=None, fixed_r1=None):
    points = tuple(range(len(POINTS))) if points is None else points
    r0_values = range(P) if fixed_r0 is None else (fixed_r0,)
    r1_values = range(P) if fixed_r1 is None else (fixed_r1,)
    return tuple(
        tuple(tensor[point][r0][r1][difference][relation]
              for difference in range(P) for relation in range(P))
        for point in points for r0 in r0_values for r1 in r1_values
    )


def r1_conditional_rank(tensor):
    matrix = tuple(
        tuple(tensor[point][r0][r1][difference][relation]
              for point in range(len(POINTS)) for r0 in range(P)
              for difference in range(P) for relation in range(P))
        for r1 in range(P)
    )
    inverse = pow(P, -1, PRIME)
    mean = tuple(sum(matrix[r1][column] for r1 in range(P)) * inverse % PRIME
                 for column in range(len(matrix[0])))
    contrast = tuple(tuple((value - mean[column]) % PRIME
                           for column, value in enumerate(row)) for row in matrix)
    return rank_rows(matrix), rank_rows(contrast)


def response_rank_record(tensor):
    parent = tensor_one_digit_pointed(tensor)
    base = tensor_pointed_parent(parent)
    base_rows = response_rows_base(base)
    parent_rows = response_rows_parent(parent)
    child_rows = response_rows_child(tensor)
    global_record = (
        rank_rows(base_rows),
        rank_rows(parent_rows),
        rank_rows(base_rows + parent_rows),
        rank_rows(child_rows),
        rank_rows(parent_rows + child_rows),
        r1_conditional_rank(tensor),
    )
    state_record = []
    for state in range(V):
        points = STATE_POINTS[state]
        base_state = response_rows_base(base, points)
        parent_state = response_rows_parent(parent, points)
        child_state = response_rows_child(tensor, points)
        state_record.append((
            len(points),
            rank_rows(base_state),
            rank_rows(parent_state),
            rank_rows(child_state),
            rank_rows(base_state + parent_state + child_state),
        ))
    fixed_r0 = []
    for r0 in range(P):
        parent_block = response_rows_parent(parent, fixed_r0=r0)
        child_block = response_rows_child(tensor, fixed_r0=r0)
        fixed_r0.append((
            rank_rows(parent_block),
            rank_rows(child_block),
            rank_rows(parent_block + child_block),
        ))
    return global_record, tuple(state_record), tuple(fixed_r0)


def matrix_inverse(matrix):
    size = len(matrix)
    augmented = [
        [value % PRIME for value in row]
        + [int(column == row_index) for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    require(all(len(row) == 2 * size for row in augmented), "square inverse")
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if augmented[row][column]), None)
        require(pivot is not None, ("singular coordinate minor", column))
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], -1, PRIME)
        augmented[column] = [value * inverse % PRIME for value in augmented[column]]
        for row in range(size):
            if row == column or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                (value - factor * pivot_value) % PRIME
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return tuple(tuple(row[size:]) for row in augmented)


def pivot_columns(rows):
    work = [[value % PRIME for value in row] for row in rows]
    rank = 0
    pivots = []
    columns = len(work[0]) if work else 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, PRIME)
        work[rank] = [value * inverse % PRIME for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (value - factor * pivot_value) % PRIME
                for value, pivot_value in zip(work[row], work[rank])
            ]
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    return tuple(pivots)


def matrix_multiply(left, right):
    return tuple(
        tuple(sum(left[row][middle] * right[middle][column]
                  for middle in range(len(right))) % PRIME
              for column in range(len(right[0])))
        for row in range(len(left))
    )


def coordinate_matrix(parent_rows, child_rows):
    size = len(parent_rows)
    require(len(child_rows) == size, "carrier row count")
    rank = rank_rows(parent_rows)
    if rank != size:
        return None, ("parent-rank", rank)
    pivots = pivot_columns(parent_rows)
    require(len(pivots) == size, ("carrier pivots", pivots))
    minor = tuple(tuple(parent_rows[row][column] for column in pivots)
                  for row in range(size))
    child_minor = tuple(tuple(child_rows[row][column] for column in pivots)
                        for row in range(size))
    matrix = matrix_multiply(child_minor, matrix_inverse(minor))
    for row in range(size):
        for column in range(len(parent_rows[0])):
            expected = sum(matrix[row][middle] * parent_rows[middle][column]
                           for middle in range(size)) % PRIME
            if expected != child_rows[row][column]:
                return matrix, (row, column, expected, child_rows[row][column])
    return matrix, None


def trajectory_rows_parent(parent):
    return tuple(
        tuple(parent[point][r0][difference][relation]
              for r0 in range(P) for difference in range(P) for relation in range(P))
        for point in range(len(POINTS))
    )


def trajectory_rows_child(tensor, r1):
    return tuple(
        tuple(tensor[point][r0][r1][difference][relation]
              for r0 in range(P) for difference in range(P) for relation in range(P))
        for point in range(len(POINTS))
    )


def state_preserving(matrix) -> bool:
    return all(
        matrix[left][right] == 0
        for left in range(len(POINTS)) for right in range(len(POINTS))
        if POINTS[left][0] != POINTS[right][0]
    )


def matrix_rank_and_support(matrix):
    return rank_rows(matrix), sum(value != 0 for row in matrix for value in row)


def transition_record(tensor):
    parent = tensor_one_digit_pointed(tensor)
    parent_rows = trajectory_rows_parent(parent)
    parent_rank = rank_rows(parent_rows)
    matrices = []
    records = []
    for r1 in range(P):
        child_rows = trajectory_rows_child(tensor, r1)
        matrix, witness = coordinate_matrix(parent_rows, child_rows)
        exists = witness is None
        if matrix is None:
            rank_support = (0, 0)
            matrix_digest = "NONE"
            block = False
        else:
            rank_support = matrix_rank_and_support(matrix)
            matrix_digest = digest_json(matrix)
            block = exists and state_preserving(matrix)
        records.append((
            rank_rows(child_rows),
            rank_rows(parent_rows + child_rows),
            exists,
            block,
            rank_support,
            matrix_digest,
            witness,
        ))
        matrices.append(matrix if exists else None)

    all_exist = all(matrix is not None for matrix in matrices)
    all_block = all(record[3] for record in records)
    projective = False
    summed_digest = "NONE"
    if all_exist:
        summed = tuple(tuple(sum(matrices[r1][row][column] for r1 in range(P)) % PRIME
                             for column in range(len(POINTS)))
                       for row in range(len(POINTS)))
        identity = tuple(tuple(int(row == column) for column in range(len(POINTS)))
                         for row in range(len(POINTS)))
        projective = summed == identity
        summed_digest = digest_json(summed)
    matrix_bank_digest = digest_json(tuple(
        digest_json(matrix) if matrix is not None else "NONE" for matrix in matrices
    ))
    return (
        parent_rank,
        tuple(records),
        all_exist,
        all_block,
        projective,
        matrix_bank_digest,
        summed_digest,
    )


def address_dependent_record(tensor, source_kernel=None):
    parent = tensor_one_digit_pointed(tensor)
    parent_ranks = []
    contained = unique = state_block = diagonal = scalar = 0
    union_histogram = {}
    first_failure = None
    first_nondiagonal = None
    source_comparable = source_matches = 0
    first_source_mismatch = None
    matrix_digests = []
    rank_histogram = {}
    support_histogram = {}
    matrices_by_r0 = []
    for r0 in range(P):
        parent_rows = tuple(
            tuple(parent[point][r0][difference][relation]
                  for difference in range(P) for relation in range(P))
            for point in range(len(POINTS))
        )
        parent_rank = rank_rows(parent_rows)
        parent_ranks.append(parent_rank)
        r0_matrices = []
        for r1 in range(P):
            child_rows = tuple(
                tuple(tensor[point][r0][r1][difference][relation]
                      for difference in range(P) for relation in range(P))
                for point in range(len(POINTS))
            )
            union = rank_rows(parent_rows + child_rows)
            union_histogram[union] = union_histogram.get(union, 0) + 1
            if union == parent_rank:
                contained += 1
            matrix, witness = coordinate_matrix(parent_rows, child_rows)
            if witness is None:
                unique += 1
                matrix_rank, matrix_support = matrix_rank_and_support(matrix)
                rank_histogram[matrix_rank] = rank_histogram.get(matrix_rank, 0) + 1
                support_histogram[matrix_support] = (
                    support_histogram.get(matrix_support, 0) + 1
                )
                matrix_digests.append(digest_json(matrix))
                if state_preserving(matrix):
                    state_block += 1
                is_diagonal = all(
                    matrix[left][right] == 0
                    for left in range(len(POINTS)) for right in range(len(POINTS))
                    if left != right
                )
                if is_diagonal:
                    diagonal += 1
                elif first_nondiagonal is None:
                    first_nondiagonal = (r0, r1, digest_json(matrix))
                diagonal_values = tuple(matrix[index][index]
                                        for index in range(len(POINTS)))
                if is_diagonal and len(set(diagonal_values)) == 1:
                    scalar += 1
                if source_kernel is not None and is_diagonal:
                    for point in range(len(POINTS)):
                        expected = source_kernel[r0][r1][point]
                        if expected is None:
                            continue
                        source_comparable += 1
                        if matrix[point][point] == expected:
                            source_matches += 1
                        elif first_source_mismatch is None:
                            first_source_mismatch = (
                                r0, r1, point, expected, matrix[point][point]
                            )
                r0_matrices.append(matrix)
            else:
                r0_matrices.append(None)
            if union > parent_rank and first_failure is None:
                first_failure = (r0, r1, parent_rank, rank_rows(child_rows), union,
                                 witness)
        matrices_by_r0.append(tuple(r0_matrices))

    identity = tuple(tuple(int(row == column) for column in range(len(POINTS)))
                     for row in range(len(POINTS)))
    projective_r0 = 0
    for matrices in matrices_by_r0:
        if len(matrices) != P or any(matrix is None for matrix in matrices):
            continue
        summed = tuple(tuple(sum(matrices[r1][row][column] for r1 in range(P))
                             % PRIME for column in range(len(POINTS)))
                       for row in range(len(POINTS)))
        if summed == identity:
            projective_r0 += 1

    unique_global = len(set(matrix_digests))
    unique_by_r1 = tuple(
        len({digest_json(matrices_by_r0[r0][r1])
             for r0 in range(P) if matrices_by_r0[r0][r1] is not None})
        for r1 in range(P)
    )
    diagonal_profile_digest = digest_json(tuple(
        tuple(
            None if matrices_by_r0[r0][r1] is None else tuple(
                matrices_by_r0[r0][r1][point][point]
                for point in range(len(POINTS))
            )
            for r1 in range(P)
        )
        for r0 in range(P)
    ))
    return (
        tuple(parent_ranks),
        contained,
        unique,
        state_block,
        diagonal,
        scalar,
        projective_r0,
        tuple(sorted(rank_histogram.items())),
        tuple(sorted(support_histogram.items())),
        tuple(sorted(union_histogram.items())),
        unique_global,
        unique_by_r1,
        first_failure,
        first_nondiagonal,
        source_comparable,
        source_matches,
        first_source_mismatch,
        digest_json(tuple(matrix_digests)),
        diagonal_profile_digest,
    )


def direct_inversion_record(probe_chunks, tensor, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    records = []
    for probe_index, probe in enumerate(INVERSION_PROBES):
        values = []
        for relation in PROBE_RELATIONS:
            total = 0
            for alpha, probes in probe_chunks:
                recorded, rows = probes[probe_index]
                require(recorded == probe, ("probe order", recorded, probe))
                index = 0
                for _beta in range(P):
                    for tau in range(P):
                        phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                        total = (total + rows[index] * phase) % PRIME
                        index += 1
                require(index == P**2, ("probe row size", probe, index))
            value = total * normalizer % PRIME
            point, r0, r1, difference = probe
            require(value == tensor[point][r0][r1][difference][relation],
                    ("direct/compressed inversion", probe, relation, value,
                     tensor[point][r0][r1][difference][relation]))
            values.append(value)
        records.append((probe, tuple(values)))
    return tuple(records)


def main() -> None:
    _profiles, boundaries, cells, source_record = T2.two_digit_context()
    source_digest = digest_json((source_record, boundaries))
    require(source_digest == T2.EXPECTED_SOURCE_SHA256,
            ("two-digit source digest", source_digest))
    support_record = source_support_record(cells)
    source_kernels, source_ratio_record = source_ratio_kernels(cells)

    tau_actual = [[[[[0 for _s in range(P)] for _r1 in range(P)]
                    for _r0 in range(P)] for _point in POINTS]
                  for _tau in range(P)]
    tau_support = [[[[[0 for _s in range(P)] for _r1 in range(P)]
                     for _r0 in range(P)] for _point in POINTS]
                   for _tau in range(P)]
    one_digit_gamma = []
    pointed_gamma = []
    alpha_merkle = []
    control_rows = {}
    probe_chunks = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    update_count = 0
    key_counts = [0] * P

    with ProcessPoolExecutor(max_workers=4) as pool:
        for expected_alpha, chunk in enumerate(pool.map(worker, range(P))):
            alpha = chunk[0]
            require(alpha == expected_alpha, ("worker order", alpha, expected_alpha))
            add_core(tau_actual, chunk[1])
            add_core(tau_support, chunk[2])
            one_digit_gamma.extend(chunk[3])
            pointed_gamma.extend(chunk[4])
            alpha_merkle.append(chunk[5])
            control_rows.update(dict(chunk[7]))
            probe_chunks.append((alpha, chunk[8]))
            scalar_counts = [left + right
                             for left, right in zip(scalar_counts, chunk[6][0])]
            state_counts = [left + right
                            for left, right in zip(state_counts, chunk[6][1])]
            update_count += chunk[6][2]
            key_counts = [left + right
                          for left, right in zip(key_counts, chunk[6][3])]

    require(len(one_digit_gamma) == len(pointed_gamma) == P**3,
            ("parent gamma sizes", len(one_digit_gamma), len(pointed_gamma)))
    one_digit_gamma_digest = digest_json(tuple(one_digit_gamma))
    pointed_gamma_digest = digest_json(tuple(pointed_gamma))
    require(one_digit_gamma_digest == J1.EXPECTED_DIGESTS[0][0],
            ("one-digit/root gamma marginal", one_digit_gamma_digest))
    require(pointed_gamma_digest == J1.POINTED_GAMMA_SHA256,
            ("pointed gamma marginal", pointed_gamma_digest))

    zeta = C.context()["zeta"]
    actual = inverse_core(tau_actual, zeta)
    support = inverse_core(tau_support, zeta)
    two_digit = tensor_two_digit_marginal(actual)
    require(digest_json(two_digit) == T2.EXPECTED_DIGESTS[0],
            ("two-digit tensor marginal", digest_json(two_digit)))

    one_digit_pointed = tensor_one_digit_pointed(actual)
    one_digit_state = tensor_one_digit_state(one_digit_pointed)
    require(digest_json(one_digit_state) == J1.EXPECTED_DIGESTS[1][0],
            ("one-digit/root tensor marginal", digest_json(one_digit_state)))
    pointed_parent = tensor_pointed_parent(one_digit_pointed)
    require(digest_json(pointed_parent) == J1.POINTED_TENSOR_SHA256,
            ("pointed tensor marginal", digest_json(pointed_parent)))

    support_pointed = tensor_one_digit_pointed(support)
    require(support_pointed == one_digit_pointed,
            "support-normalized tensor misses one-digit pointed parent")
    require(all(actual[point][r0][r1][0][relation] == 0
                and support[point][r0][r1][0][relation] == 0
                for point in range(len(POINTS)) for r0 in range(P)
                for r1 in range(P) for relation in range(P)),
            "same-root slice")

    flat = r1_flat_lift(one_digit_pointed)
    require(tensor_one_digit_pointed(flat) == one_digit_pointed,
            "r1-flat one-digit pointed parent")
    point_difference_flat = point_difference_flat_lift(two_digit)
    require(tensor_two_digit_marginal(point_difference_flat) == two_digit,
            "point/difference-flat two-digit parent")
    tensors = (actual, support, flat, point_difference_flat)

    direct_records = []
    for alpha, beta, tau in DIRECT_CONTROLS:
        direct, counts = integrate_joint(alpha, beta, tau)
        phase = pow(zeta, beta, PRIME)
        rows = tuple(twist_joint_row(bank[0], phase) for bank in direct)
        require(rows == control_rows[(alpha, beta, tau)],
                ("literal joint guard", alpha, beta, tau))

        two_parent, _two_counts = T2.integrate_two_digits(alpha, beta, tau)
        two_row = T2.phase_bank(two_parent[0], phase)
        require(gamma_two_digit_marginal(rows[0]) == two_row,
                ("literal two-digit marginal", alpha, beta, tau))

        one_parent = J1.integrate_joint(alpha, beta, tau)[0]
        one_row = tuple(
            tuple(
                tuple(phase * value % PRIME for value in difference_row)
                for difference_row in branch_rows
            )
            for branch_rows in one_parent[0]
        )
        require(gamma_one_digit_state(rows[0]) == one_row,
                ("literal one-digit/root marginal", alpha, beta, tau))
        direct_records.append(((alpha, beta, tau), counts))

    inversion_record = direct_inversion_record(probe_chunks, actual, zeta)
    rank_records = tuple(response_rank_record(tensor) for tensor in tensors)
    transition_records = tuple(transition_record(tensor) for tensor in tensors)
    address_records = (
        address_dependent_record(actual, source_kernels[0]),
        address_dependent_record(support, source_kernels[1]),
        address_dependent_record(flat),
        address_dependent_record(point_difference_flat),
    )
    require(address_records[0][14:17] == (954, 954, None),
            ("actual source/tensor diagonal agreement", address_records[0][14:17]))
    require(address_records[1][14:17] == (1014, 1014, None),
            ("support source/tensor diagonal agreement", address_records[1][14:17]))

    # Incoming 2d52215a3 is the interpretation control: its rank-17 source
    # character bank has rank-four amplitude image but exactly rank-six
    # pointed carrier.  Here r1 has full 13/12 raw/contrast amplitude rank,
    # while every base/one/two response and union rank remains six.
    amplitude_carrier_record = (
        NESTED_AMPLITUDE_CARRIER_CONTROL,
        rank_records[0][0],
        rank_records[1][0],
        rank_records[2][0],
    )
    require(rank_records[0][0] == (6, 6, 6, 6, 6, (13, 12)),
            ("actual amplitude/carrier distinction", rank_records[0][0]))
    require(rank_records[1][0] == (6, 6, 6, 6, 6, (4, 3)),
            ("support amplitude/carrier distinction", rank_records[1][0]))
    require(rank_records[2][0] == (6, 6, 6, 6, 6, (1, 0)),
            ("flat amplitude/carrier distinction", rank_records[2][0]))

    inverse = pow(P, -1, PRIME)
    flat_transition = transition_records[2]
    require(flat_transition[2] and flat_transition[3] and flat_transition[4],
            ("flat transition closure", flat_transition))
    expected_flat_matrix = tuple(
        tuple(inverse if row == column else 0 for column in range(len(POINTS)))
        for row in range(len(POINTS))
    )
    expected_flat_digest = digest_json(expected_flat_matrix)
    require(all(record[5] == expected_flat_digest for record in flat_transition[1]),
            ("flat transition matrices", flat_transition[1]))

    tensor_digests = tuple(digest_json(tensor) for tensor in tensors)
    core_digests = (digest_json(freeze5(tau_actual)), digest_json(freeze5(tau_support)))
    work_record = (
        tuple(scalar_counts), tuple(state_counts), update_count, tuple(key_counts),
        support_record, source_ratio_record, tuple(direct_records), inversion_record,
    )

    verification_digests = (
        digest_json(tuple(alpha_merkle)),
        digest_json(rank_records),
        digest_json((transition_records, address_records)),
        digest_json(work_record),
    )
    if EXPECTED_JOINT_ALPHA_MERKLE_SHA256 != "TO_BE_PINNED":
        require(verification_digests[0] == EXPECTED_JOINT_ALPHA_MERKLE_SHA256,
                ("joint alpha Merkle drift", verification_digests[0]))
    if EXPECTED_TENSOR_DIGESTS != "TO_BE_PINNED":
        require((tensor_digests, core_digests) == EXPECTED_TENSOR_DIGESTS,
                ("joint tensor digests", tensor_digests, core_digests))
    if EXPECTED_RANK_RECORD_SHA256 != "TO_BE_PINNED":
        require(verification_digests[1] == EXPECTED_RANK_RECORD_SHA256,
                ("joint response ranks", verification_digests[1]))
    if EXPECTED_TRANSITION_RECORD_SHA256 != "TO_BE_PINNED":
        require(verification_digests[2] == EXPECTED_TRANSITION_RECORD_SHA256,
                ("joint carrier transitions", verification_digests[2]))
    if EXPECTED_WORK_RECORD_SHA256 != "TO_BE_PINNED":
        require(verification_digests[3] == EXPECTED_WORK_RECORD_SHA256,
                ("joint work record", verification_digests[3]))

    record = (
        TWO_DIGIT_SHA256,
        T2.EXPECTED_SEMANTIC_SHA256,
        TWO_DIGIT_AUDIT_SHA256,
        TWO_DIGIT_AUDIT_SEMANTIC,
        ONE_DIGIT_ROOT_SHA256,
        J1.EXPECTED_SEMANTIC_SHA256,
        POINTED_AUDIT_SHA256,
        POINTED_AUDIT_SEMANTIC,
        NESTED_SOURCE_CURRENT_SHA256,
        NESTED_SOURCE_CURRENT_SEMANTIC,
        amplitude_carrier_record,
        POINTS,
        STATE_POINTS,
        PRIME,
        C.JOINT_ROOT,
        zeta,
        source_digest,
        tuple(alpha_merkle),
        one_digit_gamma_digest,
        pointed_gamma_digest,
        work_record,
        tensor_digests,
        core_digests,
        rank_records,
        transition_records,
        address_records,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("joint semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    names = ("actual", "r1_support_normalized", "r1_flat", "point_difference_flat")
    print("== r=5 two inverse-cylinder digits x pointed-six root difference ==")
    print(f"parents=(two_digit_sha={TWO_DIGIT_SHA256},two_digit_independent_audit_sha={TWO_DIGIT_AUDIT_SHA256},one_digit_root_sha={ONE_DIGIT_ROOT_SHA256},pointed_independent_audit_sha={POINTED_AUDIT_SHA256},nested_amplitude_carrier_control_sha={NESTED_SOURCE_CURRENT_SHA256})")
    print(f"coordinates=(point=(state,u) in {POINTS},r0,r1,s=u-q!=0,relation=(1,0,t))")
    print("construction=ordered_root_pair reopened before endpoint multiplication;mod169 half-open profiles expanded only after sparse cell/root/difference aggregation")
    print(f"field=(prime={PRIME},root={C.JOINT_ROOT},zeta13={zeta});source_digest={source_digest}")
    print(f"r1_support_per_(cell,u,r0)=(values,patterns)={support_record}")
    print(f"source_cell_child_over_parent_proportionality_(actual,support)={source_ratio_record}")
    print(f"sparse_work=(counts={tuple(scalar_counts)},state_counts={tuple(state_counts)},updates={update_count},keys={tuple(key_counts)})")
    print(f"literal_controls={DIRECT_CONTROLS}: PASS;direct_2197_term_inversion={inversion_record}: PASS")
    print("preintegration_marginals=(sum_point,s->two_digit;sum_point,r1->one_digit_x_root;sum_r0,r1->independently_audited_pointed_root): PASS")
    print("inverse_marginals=(same three quotient maps): PASS;same_root_s0=0")
    print("hostiles=(support_normalized_same_one_digit_pointed_parent,r1_flat_same_parent,point_difference_flat_same_two_digit_parent): PASS")
    print(f"response_rank_order=(base,one_digit,base_union_one,two_digit,one_union_two,(r1_raw,r1_contrast));state_order=(point_count,base,one,two,union)")
    print(f"response_ranks={tuple(zip(names, rank_records))}")
    print(f"amplitude_vs_carrier_control=(incoming_source_current_(profile,amplitude,carrier,union),this_actual,this_support,this_flat)={amplitude_carrier_record}")
    print("interpretation=r1 enlarges amplitude coordinates to raw/contrast rank 13/12 but does not enlarge the pointed carrier: base=one=two=both unions=6")
    print("stationary_transition_record=(parent_trajectory_rank,per_r1(child_rank,union_rank,exists,state_block,(matrix_rank,nnz),matrix_digest,witness),all_exists,all_state_block,sum_K=I,matrix_bank_digest,sum_digest)")
    print(f"stationary_transitions={tuple(zip(names, transition_records))}")
    print("address_dependent_record=(parent_ranks,contained,unique,state_block,diagonal,scalar,projective_r0,rank_histogram,nnz_histogram,union_histogram,unique_global,unique_by_r1,first_failure,first_nondiagonal,source_comparable,source_matches,first_source_mismatch,matrix_digest,diagonal_profile_digest)")
    print(f"address_dependent_transitions={tuple(zip(names, address_records))}")
    print(f"alpha_gamma_Merkle_(actual,support)={tuple(alpha_merkle)}")
    print(f"digests=(tensors_{names}={tensor_digests},tau_cores={core_digests})")
    print(f"verification_digests=(alpha_Merkle,rank,transition,work)={verification_digests}")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT pointed two-cylinder response/transition probe on one owner base")
    print("typing=static finite response coordinates only;no chronology,no complete address,no arrival ancestry,no physical current,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
