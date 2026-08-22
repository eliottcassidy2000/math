#!/usr/bin/env python3
"""Exact bridge-only audit of the U_full endpoint atom-coupling debt.

The frozen refined endpoint bank computes, for every character ``ell``,

    gamma_ell = phase_ell * AX_ell * BY_ell

after ``AX_ell`` and ``BY_ell`` have each been summed.  This companion asks
whether either of two geometric identifications still visible in the interval
geometry can recover the frozen H-q5 bridge.  The coarser one restricts the
Cartesian product to pairs lying in the same maximal cyclic ``E`` component,
with counting weight on components.  The finer one identifies the two circle
points.  On the scaled circle,
the left field is supported on ``R E_ell intersect Q`` at frequency ``X``
and the right field is supported on ``R E_ell`` at frequency ``-Y``.  Their
pointwise product is therefore supported on ``R E_ell intersect Q`` at the
combined frequency ``X-Y=-w_c``.  We call its character bank the diagonal
bank.

These are hostile/cheap candidates only.  The component restriction keeps the
inherited two-factor unit but chooses only the unweighted equality relation on
geometric components; the point restriction has one endpoint denominator
instead of two.  THM-2471's genuine ancestry fibre instead uses linked nodes,
sheet labels, horizons, and one character-independent address map.  Equality
with the frozen bridge would not by itself prove a chronological arrival.
Inequality rules out only the named restrictions, not component-dependent
weights, an arbitrary rescaling of the point restriction, or a differently
typed ancestry lift.  The API/kernel certificate proves the stronger
current-data statement: the returned marginals do not select any such lift.
"""

from __future__ import annotations

import ast
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_endpoint_ufull_atom_bridge_kernel_probe_20260816.py"
OUTPUT = "05-knowledge/results/lrc_endpoint_ufull_atom_bridge_kernel_probe_20260816.out"

BRIDGE_PATH = ROOT / "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py"
BRIDGE_SHA256 = "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b"
REFINED_PATH = ROOT / "04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py"
REFINED_SHA256 = "ee2105742abee578a9c41ff7ec954a07ada324fccc2c643429e7ac6e6e6f8fc2"
REFINED_OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.out"
REFINED_OUTPUT_SHA256 = "10a98351cc59615a5b6d2b8f555e0936d1a39566d9906127edc2b0fbc3918e73"

THEOREM_PINS = (
    (
        "THM-2471",
        ROOT / "01-canon/theorems/THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary.md",
        "423882a4e3b06e2020fd7f3d0748fc9e904ced4d3214ef8cd77bfd5c6b9597eb",
    ),
    (
        "THM-2538",
        ROOT / "01-canon/theorems/THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary.md",
        "b3179e0849f581c6f5f5a44dab9b7d57d7c87f96d71259863315b23a241ae614",
    ),
    (
        "THM-2594",
        ROOT / "01-canon/theorems/THM-2594-realized-theta-slaved-contraction-at-the-r5-window.md",
        "b46e955c4dfc7fb1019791db82295643042c319eed52c02ed8d89573d4725cdc",
    ),
)

P = 13
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)
EXPECTED_MARGINAL_VALUES = (
    320618948602619577408,
    503604956476841920373,
)
EXPECTED_MARGINAL_BRIDGE = 389266878372286537904
EXPECTED_CARTESIAN_GAMMA_SHA256 = (
    "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
)
EXPECTED_DIAGONAL_GAMMA_SHA256 = (
    "771545a5cb1f0f03459b8d351de668ad950ece5fcb985fa61d599d643de3303f"
)
EXPECTED_SEGMENT_GAMMA_SHA256 = (
    "1270078f3c019a3fa5ab507128e5010149ba93bbb6f74d4c6967cc2753df21ea"
)
EXPECTED_SEGMENT_VALUES = (
    99143203042879994518,
    130742602587392835137,
    540653486701996040250,
)
EXPECTED_CYCLIC_GAMMA_SHA256 = (
    "1270078f3c019a3fa5ab507128e5010149ba93bbb6f74d4c6967cc2753df21ea"
)
EXPECTED_CYCLIC_VALUES = (
    99143203042879994518,
    130742602587392835137,
    540653486701996040250,
)
EXPECTED_CROSS_GAMMA_SHA256 = (
    "529cf77d326ada80ab99bfbf4b7b6444d41794820f1fc4552555846a5f49c978"
)
EXPECTED_CROSS_VALUES = (
    221475745559739582890,
    372862353889449085236,
    420866277916799378523,
)
EXPECTED_DIAGONAL_VALUES = (
    633668780131603861,
    405160484437854840264,
    167726070588785644466,
)
EXPECTED_DISCRETE_SCALE = 123588610788991450223
EXPECTED_DISCRETE_SCALED_BRIDGE = 24876622649422736677
EXPECTED_FREQUENCY_SCALE = 572122651179088191560
EXPECTED_FREQUENCY_SCALED_BRIDGE = 284514977757516176864
EXPECTED_SEMANTIC_SHA256 = (
    "b28f64e9beac1a9677d6562b606983b35824f54d9d2dc812d8678f8dbea47948"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_integers(values: tuple[int, ...]) -> str:
    return hashlib.sha256(
        ",".join(str(value) for value in values).encode("ascii")
    ).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_bridge_module():
    require(lf_sha256(BRIDGE_PATH) == BRIDGE_SHA256, "bridge source hash drift")
    spec = importlib.util.spec_from_file_location("thm3479_atom_bridge", BRIDGE_PATH)
    require(spec is not None and spec.loader is not None, "bridge module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_bridge_module()
CTX = None


def context():
    global CTX
    if CTX is None:
        word = M.to_current(M.U_FULL_REL)
        t_den = 182 * M.lcm_tuple(word)
        nn = M.R_DILATION * t_den
        prime, root, prime_factors, lucas_witnesses = M.FULL_EMBEDDINGS[0]
        M.verify_lucas_certificate(
            prime, prime_factors, lucas_witnesses, "U_full atom bridge prime"
        )
        M.verify_embedding(
            prime, root, nn, M.FULL_NN_FACTORS, "U_full atom bridge embedding"
        )
        zero = (0,) * 9
        q_intervals = M.fast_build_set(word, t_den, M.PATTERN_QA, zero)
        q_starts = [left for left, _right in q_intervals]
        embeddings = ((prime, root),)
        # AX has exponent -X*y and BY has exponent +Y*y, with
        # Y=X+w_target.  The pointwise product has exponent +w_target*y,
        # represented by fast_x_sweep at frequency -w_target.
        diagonal_frequency = -word[M.TARGET_B]
        diagonal_tabs = M.fast_make_tabs(
            q_intervals, diagonal_frequency, nn, embeddings
        )
        x_tabs = M.fast_make_tabs(
            q_intervals, M.X_FREQUENCY, nn, embeddings
        )
        CTX = (
            word,
            t_den,
            nn,
            prime,
            root,
            q_intervals,
            q_starts,
            embeddings,
            diagonal_frequency,
            diagonal_tabs,
            x_tabs,
        )
    return CTX


def maximal_linear_segments(
    intervals: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    output = []
    for left, right in intervals:
        require(left < right, ("empty interval", left, right))
        if output and output[-1][1] == left:
            output[-1] = (output[-1][0], right)
        else:
            if output:
                require(output[-1][1] < left,
                        ("overlapping intervals", output[-1], (left, right)))
            output.append((left, right))
    require(
        all(right < next_left
            for (_left, right), (next_left, _next_right)
            in zip(output, output[1:])),
        "nonmaximal internal segments",
    )
    return output


def joint_sweep(
    e_intervals: list[tuple[int, int]],
    q_intervals: list[tuple[int, int]],
    q_starts: list[int],
    t_den: int,
    nn: int,
    prime: int,
    root: int,
    diagonal_frequency: int,
    diagonal_tabs: tuple[tuple[list[int], list[int]], ...],
    x_tabs: tuple[tuple[list[int], list[int]], ...],
    y_frequency: int,
) -> tuple[int, int, int, int, int, int]:
    """AX, BY, point, cut-segment, cyclic-component, and overlap totals."""
    frequencies = (M.X_FREQUENCY, diagonal_frequency)
    tab_rows = (x_tabs[0], diagonal_tabs[0])
    wrap_forward = tuple(
        pow(root, (-frequency * t_den) % nn, prime)
        for frequency in frequencies
    )
    wrap_backward = tuple(
        pow(root, (frequency * t_den) % nn, prime)
        for frequency in frequencies
    )
    ax_total = 0
    by_total = 0
    point_total = 0
    component_total = 0
    overlap_total = 0
    first_component = None
    last_component = None
    for component_index, (e_left, e_right) in enumerate(e_intervals):
        scaled_left = M.R_DILATION * e_left
        start = scaled_left % t_den
        span = M.R_DILATION * (e_right - e_left)
        require(span < t_den, ("joint sweep span", span, t_den))
        stop = start + span
        index = bisect_right(q_starts, start) - 1
        offset = 0
        if index < 0:
            index = len(q_intervals) - 1
            offset = -t_den
        base_values = tuple(
            pow(root, (-frequency * (scaled_left - start)) % nn, prime)
            for frequency in frequencies
        )
        accumulators = [0, 0]
        wrap_values = [1, 1]
        if offset:
            wrap_values = list(wrap_backward)
        while True:
            q_left_0, q_right_0 = q_intervals[index]
            q_left = q_left_0 + offset
            q_right = q_right_0 + offset
            if q_left >= stop:
                break
            if q_right > start:
                left = max(start, q_left)
                right = min(stop, q_right)
                if left < right:
                    overlap_total += right - left
                    for position, frequency in enumerate(frequencies):
                        tabs_left, tabs_right = tab_rows[position]
                        if left == q_left:
                            value_left = (
                                tabs_left[index] * wrap_values[position]
                            ) % prime
                        else:
                            value_left = pow(
                                root, (-frequency * left) % nn, prime
                            )
                        if right == q_right:
                            value_right = (
                                tabs_right[index] * wrap_values[position]
                            ) % prime
                        else:
                            value_right = pow(
                                root, (-frequency * right) % nn, prime
                            )
                        accumulators[position] = (
                            accumulators[position] + value_left - value_right
                        ) % prime
            index += 1
            if index == len(q_intervals):
                index = 0
                offset += t_den
                for position in range(2):
                    wrap_values[position] = (
                        wrap_values[position] * wrap_forward[position]
                    ) % prime
        ax_component = base_values[0] * accumulators[0] % prime
        point_component = base_values[1] * accumulators[1] % prime
        by_component = (
            pow(root, (y_frequency * M.R_DILATION * e_left) % nn, prime)
            - pow(root, (y_frequency * M.R_DILATION * e_right) % nn, prime)
        ) % prime
        ax_total = (ax_total + ax_component) % prime
        by_total = (by_total + by_component) % prime
        point_total = (point_total + point_component) % prime
        component_total = (
            component_total + ax_component * by_component
        ) % prime
        if component_index == 0:
            first_component = (ax_component, by_component)
        last_component = (ax_component, by_component)
    cyclic_component_total = component_total
    if (
        len(e_intervals) > 1
        and e_intervals[0][0] == 0
        and e_intervals[-1][1] == t_den
    ):
        require(first_component is not None and last_component is not None,
                "cyclic endpoint components absent")
        first_ax, first_by = first_component
        last_ax, last_by = last_component
        cyclic_component_total = (
            cyclic_component_total
            + first_ax * last_by
            + last_ax * first_by
        ) % prime
    return (
        ax_total,
        by_total,
        point_total,
        component_total,
        cyclic_component_total,
        overlap_total,
    )


def diagonal_worker(
    alpha: int,
) -> tuple[
    int,
    int,
    int,
    int,
    int,
    tuple[int, ...],
    tuple[int, ...],
    tuple[int, ...],
    tuple[int, ...],
]:
    (
        word,
        t_den,
        nn,
        prime,
        root,
        q_intervals,
        q_starts,
        embeddings,
        diagonal_frequency,
        diagonal_tabs,
        x_tabs,
    ) = context()
    cartesian_rows = []
    point_rows = []
    component_rows = []
    cyclic_component_rows = []
    raw_interval_count = 0
    maximal_interval_count = 0
    cyclic_merge_count = 0
    overlap_length = 0
    for beta in range(P):
        for tau in range(P):
            ell = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
            e_raw = M.fast_build_set(word, t_den, M.PATTERN_E, ell)
            e_intervals = maximal_linear_segments(e_raw)
            raw_interval_count += len(e_raw)
            maximal_interval_count += len(e_intervals)
            cyclic_merge_count += int(
                len(e_intervals) > 1
                and e_intervals[0][0] == 0
                and e_intervals[-1][1] == t_den
            )
            ax, by, point, segment_component, cyclic_component, overlap = joint_sweep(
                e_intervals,
                q_intervals,
                q_starts,
                t_den,
                nn,
                prime,
                root,
                diagonal_frequency,
                diagonal_tabs,
                x_tabs,
                M.X_FREQUENCY + word[M.TARGET_B],
            )
            overlap_length += overlap
            phase = pow(root, beta * (nn // P) % nn, prime)
            cartesian_rows.append(phase * ax % prime * by % prime)
            point_rows.append(phase * point % prime)
            component_rows.append(phase * segment_component % prime)
            cyclic_component_rows.append(phase * cyclic_component % prime)
    return (
        alpha,
        raw_interval_count,
        maximal_interval_count,
        cyclic_merge_count,
        overlap_length,
        tuple(cartesian_rows),
        tuple(point_rows),
        tuple(component_rows),
        tuple(cyclic_component_rows),
    )


def inverse_value(
    gamma: tuple[int, ...], q: tuple[int, int, int], prime: int, zeta: int
) -> int:
    powers = tuple(pow(zeta, exponent, prime) for exponent in range(P))
    total = 0
    index = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                exponent = -(alpha * q[0] + beta * q[1] + tau * q[2]) % P
                total = (total + gamma[index] * powers[exponent]) % prime
                index += 1
    return total * pow(P**3, -1, prime) % prime


def reference_controls() -> tuple[object, ...]:
    (
        word,
        t_den,
        nn,
        prime,
        root,
        q_intervals_fast,
        q_starts,
        embeddings,
        diagonal_frequency,
        diagonal_tabs,
        x_tabs,
    ) = context()
    zero = (0,) * 9
    _v1, v2, _reps = M.target_representatives()
    q_intervals_reference = M.reference_build_set(
        word, t_den, M.PATTERN_QA, zero
    )
    rows = []
    for label, ell in (("zero", zero), ("v2", v2)):
        e_raw = M.fast_build_set(word, t_den, M.PATTERN_E, ell)
        e_fast = maximal_linear_segments(e_raw)
        e_reference = M.reference_build_set(word, t_den, M.PATTERN_E, ell)
        joint = joint_sweep(
            e_fast,
            q_intervals_fast,
            q_starts,
            t_den,
            nn,
            prime,
            root,
            diagonal_frequency,
            diagonal_tabs,
            x_tabs,
            M.X_FREQUENCY + word[M.TARGET_B],
        )
        fast_point, fast_overlap = M.fast_x_sweep(
            e_fast,
            q_intervals_fast,
            q_starts,
            diagonal_frequency,
            t_den,
            nn,
            embeddings,
            diagonal_tabs,
        )
        fast_ax, ax_overlap = M.fast_x_sweep(
            e_fast,
            q_intervals_fast,
            q_starts,
            M.X_FREQUENCY,
            t_den,
            nn,
            embeddings,
            x_tabs,
        )
        fast_by = M.fast_endpoint_sum(
            e_fast,
            -(M.X_FREQUENCY + word[M.TARGET_B]),
            nn,
            embeddings,
        )
        reference, reference_overlap, components = M.reference_marked_sum(
            e_reference,
            q_intervals_reference,
            diagonal_frequency,
            t_den,
            nn,
            embeddings,
        )
        public_component = 0
        public_component_pairs = []
        for interval in e_fast:
            ax_component, _component_overlap = M.fast_x_sweep(
                [interval],
                q_intervals_fast,
                q_starts,
                M.X_FREQUENCY,
                t_den,
                nn,
                embeddings,
                x_tabs,
            )
            by_component = M.fast_endpoint_sum(
                [interval],
                -(M.X_FREQUENCY + word[M.TARGET_B]),
                nn,
                embeddings,
            )
            public_component = (
                public_component + ax_component[0] * by_component[0]
            ) % prime
            public_component_pairs.append((ax_component[0], by_component[0]))
        public_cyclic_component = public_component
        if (
            len(e_fast) > 1
            and e_fast[0][0] == 0
            and e_fast[-1][1] == t_den
        ):
            first_ax, first_by = public_component_pairs[0]
            last_ax, last_by = public_component_pairs[-1]
            public_cyclic_component = (
                public_cyclic_component
                + first_ax * last_by
                + last_ax * first_by
            ) % prime
        require(
            joint[:4]
            == (fast_ax[0], fast_by[0], fast_point[0], public_component),
            (label, "joint/public mismatch", joint[:4]),
        )
        require(
            joint[4] == public_cyclic_component,
            (label, "joint/public cyclic mismatch", joint[4]),
        )
        require(fast_point == reference,
                (label, "diagonal fast/reference", fast_point, reference))
        require(fast_overlap == reference_overlap, (label, "overlap mismatch"))
        require(ax_overlap == reference_overlap, (label, "AX overlap mismatch"))
        require(joint[5] == reference_overlap, (label, "joint overlap mismatch"))
        rows.append(
            (
                label,
                len(e_raw),
                len(e_fast),
                len(e_reference),
                fast_overlap,
                components,
                fast_ax,
                fast_by,
                fast_point,
                public_component,
                public_cyclic_component,
            )
        )
    return tuple(rows)


def api_kernel_certificate() -> tuple[object, ...]:
    source = BRIDGE_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    functions = {
        node.name: node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    required = ("fast_x_sweep", "fast_endpoint_sum")
    require(all(name in functions for name in required), "endpoint API absent")
    rows = []
    for name in required:
        node = functions[name]
        args = tuple(argument.arg for argument in node.args.args)
        yields = tuple(
            type(item).__name__
            for item in ast.walk(node)
            if isinstance(item, (ast.Yield, ast.YieldFrom))
        )
        returns = tuple(
            ast.unparse(item.value)
            for item in ast.walk(node)
            if isinstance(item, ast.Return)
        )
        require(not yields, (name, "unexpected keyed/yield API", yields))
        require(
            not any(
                token in argument.lower()
                for argument in args
                for token in ("atom", "ancestry", "stalk", "joint", "key")
            ),
            (name, "unexpected ancestry argument", args),
        )
        rows.append((name, args, returns, yields))

    refined_source = REFINED_PATH.read_text(encoding="utf-8")
    require(
        "rows.append(phase * ax[0] % p * by[0] % p)" in refined_source,
        "refined marginal product call drift",
    )

    # The checkerboard is the minimal zero-margin direction.  It changes a
    # diagonal/common-atom functional but is invisible to every separately
    # summed endpoint bank.
    checkerboard = ((1, -1), (-1, 1))
    row_sums = tuple(sum(row) for row in checkerboard)
    column_sums = tuple(sum(checkerboard[i][j] for i in range(2)) for j in range(2))
    diagonal = checkerboard[0][0] + checkerboard[1][1]
    require(row_sums == column_sums == (0, 0), "checkerboard margins")
    require(diagonal == 2, "checkerboard diagonal")
    return (
        tuple(rows),
        (
            "composition",
            "phase*fast_x_sweep_total*fast_endpoint_sum_total",
        ),
        (
            "first_internal_common_label",
            "e_intervals list position; merged geometric component only",
        ),
        (
            "discarded_left_subkey",
            "q interval index and periodic wrap branch",
        ),
        (
            "absent_required_key",
            "THM-2471 linked (base,root,owner-sheet,word-sheet,source-sheet,horizon)",
        ),
        (checkerboard, row_sums, column_sums, diagonal),
        (
            "full_margin_kernel_model",
            "ker(e_U) tensor ker(e_V), dimension (m-1)(n-1)",
        ),
        (
            "actual_api_scope",
            "two character-weighted scalars; factors through a compatible full-margin observer only after an atomization is defined; kernel may be larger",
        ),
    )


def security_certificate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    forbidden_calls = {"eval", "exec", "compile", "__import__"}
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            bad.append("Assert")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in forbidden_calls:
                bad.append(node.func.id)
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), tuple(bad)


def main() -> None:
    require(lf_sha256(REFINED_PATH) == REFINED_SHA256, "refined source hash drift")
    require(
        lf_sha256(REFINED_OUTPUT_PATH) == REFINED_OUTPUT_SHA256,
        "refined output hash drift",
    )
    theorem_hashes = []
    for label, path, expected in THEOREM_PINS:
        actual = lf_sha256(path)
        require(actual == expected, (label, "hash drift", actual, expected))
        require("PROVED" in path.read_text(encoding="utf-8"), (label, "status drift"))
        theorem_hashes.append((label, actual))

    refined_output = REFINED_OUTPUT_PATH.read_text(encoding="utf-8")
    require(
        "factor_zero_counts=(bridge,left_K4,right_K4,product)=(0, 0, 0, 0)"
        in refined_output,
        "frozen bridge output drift",
    )
    require(
        str(EXPECTED_MARGINAL_BRIDGE) in refined_output,
        "frozen bridge residue absent",
    )

    controls = reference_controls()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(diagonal_worker, range(P)))
    require(tuple(row[0] for row in chunks) == tuple(range(P)), "worker order")
    raw_interval_counts = tuple(row[1] for row in chunks)
    maximal_interval_counts = tuple(row[2] for row in chunks)
    cyclic_merge_counts = tuple(row[3] for row in chunks)
    overlap_lengths = tuple(row[4] for row in chunks)
    cartesian_gamma = tuple(value for row in chunks for value in row[5])
    point_gamma = tuple(value for row in chunks for value in row[6])
    segment_gamma = tuple(value for row in chunks for value in row[7])
    cyclic_gamma = tuple(value for row in chunks for value in row[8])
    require(
        tuple(map(len, (cartesian_gamma, point_gamma, segment_gamma, cyclic_gamma)))
        == (P**3, P**3, P**3, P**3),
        "character bank sizes",
    )

    (
        word,
        t_den,
        nn,
        prime,
        root,
        q_intervals,
        _q_starts,
        _embeddings,
        diagonal_frequency,
        _diagonal_tabs,
        _x_tabs,
    ) = context()
    zeta = pow(root, nn // P, prime)
    require(pow(zeta, P, prime) == 1 and zeta != 1, "zeta order")
    cartesian_h = inverse_value(cartesian_gamma, Q_H, prime, zeta)
    cartesian_q5 = inverse_value(cartesian_gamma, Q_Q5, prime, zeta)
    cartesian_bridge = (cartesian_h - cartesian_q5) % prime
    diagonal_h = inverse_value(point_gamma, Q_H, prime, zeta)
    diagonal_q5 = inverse_value(point_gamma, Q_Q5, prime, zeta)
    diagonal_bridge = (diagonal_h - diagonal_q5) % prime
    segment_h = inverse_value(segment_gamma, Q_H, prime, zeta)
    segment_q5 = inverse_value(segment_gamma, Q_Q5, prime, zeta)
    segment_bridge = (segment_h - segment_q5) % prime
    cyclic_h = inverse_value(cyclic_gamma, Q_H, prime, zeta)
    cyclic_q5 = inverse_value(cyclic_gamma, Q_Q5, prime, zeta)
    cyclic_bridge = (cyclic_h - cyclic_q5) % prime
    cross_gamma = tuple(
        (cartesian - component) % prime
        for cartesian, component in zip(cartesian_gamma, cyclic_gamma)
    )
    cross_h = (cartesian_h - cyclic_h) % prime
    cross_q5 = (cartesian_q5 - cyclic_q5) % prime
    cross_bridge = (cartesian_bridge - cyclic_bridge) % prime
    cartesian_hash = digest_integers(cartesian_gamma)
    point_hash = digest_integers(point_gamma)
    segment_hash = digest_integers(segment_gamma)
    cyclic_hash = digest_integers(cyclic_gamma)
    cross_hash = digest_integers(cross_gamma)
    require(
        cartesian_hash == EXPECTED_CARTESIAN_GAMMA_SHA256,
        ("cartesian gamma hash", cartesian_hash),
    )
    require(
        (cartesian_h, cartesian_q5) == EXPECTED_MARGINAL_VALUES,
        ("cartesian inverse values", cartesian_h, cartesian_q5),
    )
    require(
        cartesian_bridge == EXPECTED_MARGINAL_BRIDGE,
        ("cartesian bridge", cartesian_bridge),
    )

    # Endpoint differences omit the geometric-series/Fourier denominators.
    # One exact discrete denominator restoration and one unsigned rational
    # frequency-magnitude probe are tested explicitly.  Neither is a semantic
    # ancestry selector; the latter is not an Archimedean Fourier normalization.
    x_frequency = M.X_FREQUENCY
    target_speed = word[M.TARGET_B]
    y_frequency = x_frequency + target_speed
    discrete_scale = (
        (1 - pow(root, (-x_frequency) % nn, prime))
        * (1 - pow(root, y_frequency % nn, prime))
        * pow(1 - pow(root, target_speed % nn, prime), -1, prime)
    ) % prime
    discrete_scaled_bridge = diagonal_bridge * discrete_scale % prime
    frequency_scale = (
        x_frequency * y_frequency * pow(target_speed, -1, prime)
    ) % prime
    frequency_scaled_bridge = diagonal_bridge * frequency_scale % prime
    require(discrete_scale == EXPECTED_DISCRETE_SCALE, "discrete scale drift")
    require(
        discrete_scaled_bridge == EXPECTED_DISCRETE_SCALED_BRIDGE,
        "discrete scaled bridge drift",
    )
    require(frequency_scale == EXPECTED_FREQUENCY_SCALE, "frequency scale drift")
    require(
        frequency_scaled_bridge == EXPECTED_FREQUENCY_SCALED_BRIDGE,
        "frequency scaled bridge drift",
    )

    require(
        point_hash == EXPECTED_DIAGONAL_GAMMA_SHA256,
        ("diagonal gamma hash", point_hash),
    )
    require(
        segment_hash == EXPECTED_SEGMENT_GAMMA_SHA256,
        ("segment gamma hash", segment_hash),
    )
    require(
        (segment_h, segment_q5, segment_bridge) == EXPECTED_SEGMENT_VALUES,
        ("segment values", segment_h, segment_q5, segment_bridge),
    )
    require(
        cyclic_hash == EXPECTED_CYCLIC_GAMMA_SHA256,
        ("cyclic gamma hash", cyclic_hash),
    )
    require(
        (cyclic_h, cyclic_q5, cyclic_bridge) == EXPECTED_CYCLIC_VALUES,
        ("cyclic values", cyclic_h, cyclic_q5, cyclic_bridge),
    )
    require(
        cross_hash == EXPECTED_CROSS_GAMMA_SHA256,
        ("cross gamma hash", cross_hash),
    )
    require(
        (cross_h, cross_q5, cross_bridge) == EXPECTED_CROSS_VALUES,
        ("cross values", cross_h, cross_q5, cross_bridge),
    )
    require(
        (diagonal_h, diagonal_q5, diagonal_bridge) == EXPECTED_DIAGONAL_VALUES,
        (
            "diagonal values",
            diagonal_h,
            diagonal_q5,
            diagonal_bridge,
        ),
    )

    api = api_kernel_certificate()
    security = security_certificate(ROOT / SCRIPT)
    consequence = (
        tuple(theorem_hashes),
        BRIDGE_SHA256,
        REFINED_SHA256,
        REFINED_OUTPUT_SHA256,
        (prime, root, t_den, nn),
        tuple(word),
        diagonal_frequency,
        len(q_intervals),
        raw_interval_counts,
        maximal_interval_counts,
        cyclic_merge_counts,
        overlap_lengths,
        controls,
        cartesian_hash,
        point_hash,
        segment_hash,
        cyclic_hash,
        cross_hash,
        (Q_H, cartesian_h, diagonal_h, segment_h, cyclic_h),
        (Q_Q5, cartesian_q5, diagonal_q5, segment_q5, cyclic_q5),
        cartesian_bridge,
        (Q_H, diagonal_h),
        (Q_Q5, diagonal_q5),
        diagonal_bridge,
        segment_bridge,
        cyclic_bridge,
        (Q_H, cross_h),
        (Q_Q5, cross_q5),
        cross_bridge,
        (discrete_scale, discrete_scaled_bridge),
        (frequency_scale, frequency_scaled_bridge),
        EXPECTED_MARGINAL_VALUES,
        EXPECTED_MARGINAL_BRIDGE,
        (
            diagonal_bridge == EXPECTED_MARGINAL_BRIDGE,
            segment_bridge == EXPECTED_MARGINAL_BRIDGE,
            cyclic_bridge == EXPECTED_MARGINAL_BRIDGE,
        ),
        api,
    )
    semantic_hash = hashlib.sha256(repr(consequence).encode("utf-8")).hexdigest()
    require(
        semantic_hash == EXPECTED_SEMANTIC_SHA256,
        ("semantic hash", semantic_hash),
    )

    print("LRC U_FULL ATOM-BRIDGE KERNEL PROBE")
    print("status=FINITE-EXACT bridge-only hostile/API result; LRC(14) OPEN")
    print(f"dependency_hashes={(BRIDGE_SHA256, REFINED_SHA256, REFINED_OUTPUT_SHA256)}")
    print(f"theorem_hashes={tuple(theorem_hashes)}")
    print(f"embedding=(prime={prime},root={root},order={nn})")
    print(f"universe=four primary F13^3 character banks plus derived cross bank, size={len(point_gamma)} each; no K4 evaluation")
    print(f"geometric_keys=(cut_segment:(ell,i),cyclic_component:(ell,[i]),point:(ell,y)); frequency_X_minus_Y={diagonal_frequency}")
    print(f"raw_interval_counts_by_alpha={raw_interval_counts} total={sum(raw_interval_counts)}")
    print(f"maximal_interval_counts_by_alpha={maximal_interval_counts} total={sum(maximal_interval_counts)}")
    print(f"cyclic_origin_merges_by_alpha={cyclic_merge_counts} total={sum(cyclic_merge_counts)}")
    print(f"overlap_lengths_by_alpha={overlap_lengths}")
    print(f"fast_reference_controls={controls}")
    print(f"gamma_sha256=(cartesian={cartesian_hash},cut_segment={segment_hash},cyclic_component={cyclic_hash},cross={cross_hash},point={point_hash})")
    print(f"cartesian_inverse_values=((q_H,{cartesian_h}),(q_q5,{cartesian_q5})); bridge={cartesian_bridge}")
    print(f"cut_segment_inverse_values=((q_H,{segment_h}),(q_q5,{segment_q5})); bridge={segment_bridge}")
    print(f"cyclic_component_inverse_values=((q_H,{cyclic_h}),(q_q5,{cyclic_q5})); bridge={cyclic_bridge}")
    print(f"cross_component_remainder=((q_H,{cross_h}),(q_q5,{cross_q5})); bridge={cross_bridge}")
    print(f"point_inverse_values=((q_H,{diagonal_h}),(q_q5,{diagonal_q5})); bridge={diagonal_bridge}")
    print(f"discrete_endpoint_scale_and_bridge={(discrete_scale, discrete_scaled_bridge)}")
    print(f"rational_frequency_magnitude_probe_and_bridge={(frequency_scale, frequency_scaled_bridge)}")
    print(f"frozen_marginal_values={EXPECTED_MARGINAL_VALUES}")
    print(f"frozen_marginal_bridge={EXPECTED_MARGINAL_BRIDGE}")
    print(f"geometric_recovery=(cut_segment={segment_bridge == EXPECTED_MARGINAL_BRIDGE},cyclic_component={cyclic_bridge == EXPECTED_MARGINAL_BRIDGE},point={diagonal_bridge == EXPECTED_MARGINAL_BRIDGE})")
    print(f"api_kernel_certificate={api}")
    print("cartesian_positive_control=recomputed phase*(sum_i AX_i)*(sum_j BY_j) bank exactly recovers the frozen gamma digest, role values, and bridge")
    print("first_failure=the E-component and point diagonals are different transforms; the returned API exposes neither component values nor a typed ancestry relation")
    print("normalization_scope=does not rule out arbitrary rescaling or a differently typed ancestry lift; those require the missing semantic key")
    print("kernel_scope=pinned proved THM-2538 gives the full-marginal mixed kernel after atomization; the actual two-scalar API is a further observer and may lose more")
    print("cheapest_engine_change=retain factor-labelled pre-merge cells, a shared THM-2471-style ancestry/horizon key, and one ell-independent address map r:Omega->F13^3; twist covariantly, multiply keyed factors, then sum and inverse DFT")
    print("replay_scope=one full normal exact sweep generated the frozen banks; no second full optimized sweep was run; frozen artifact receives compilation, static-security, cheap reference, and independent scope audits")
    print("nonconsequence=no chronological arrival, physical current, grouped exact-address coefficient, all-unit projector, row exclusion, or LRC(14)")
    print(f"security_ast={security}")
    print(f"semantic_sha256={semantic_hash}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
