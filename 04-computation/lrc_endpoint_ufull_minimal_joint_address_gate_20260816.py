#!/usr/bin/env python3
"""Exact minimum-joint-coordinate gate for the frozen U_full bridge.

This companion has two deliberately separated layers.

1.  It refines the existing U_full interval engine at the cheap ``ell=0``
    and ``ell=v2`` controls.  Surviving E and Q intervals retain the factor
    and branch that created each boundary; E--Q intersection cells retain
    both lineages and the periodic wrap turn.  Summing the labelled cells
    reproduces the existing AX and BY endpoint factors exactly.

2.  It proves the sharp 2x2 API statement suggested by the mixed-Haar
    kernel.  Row and column marginals lose one checkerboard coordinate.  One
    ell-independent joint address bit, mapped to the actual q_H and q_q5
    addresses, is necessary and sufficient to recover that coordinate.  A
    four-atom positive integer model recovers exactly the inherited frozen
    bridge, while a flat hostile has identical marginals and bridge zero.

The second layer is a minimum complete-record *schema witness*, not a U_full
ancestry construction.  The first layer shows why: the real endpoint cells
still expose no THM-2471 base/root/owner-sheet/word-sheet/source-sheet or
horizon projections.  Circle position is never called an ancestry stalk.
No K4, current, row exclusion, or LRC(14) conclusion is evaluated.
"""

from __future__ import annotations

import ast
from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass, fields
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_endpoint_ufull_minimal_joint_address_gate_20260816.py"
OUTPUT = "05-knowledge/results/lrc_endpoint_ufull_minimal_joint_address_gate_20260816.out"

BRIDGE_PATH = ROOT / "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py"
BRIDGE_SHA256 = "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b"
PREVIOUS_PATH = ROOT / "04-computation/lrc_endpoint_ufull_atom_bridge_kernel_probe_20260816.py"
PREVIOUS_SHA256 = "d1182de4d777bab20a8d423cf942151ac3149014b67d9c34883cbce37a7b0a9f"

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
        "5d3e14c9416ca52f6a523453ce3511ffb886fafe2b215b87ccf64a72e97ffd8b",
    ),
)

P = 13
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)
EXPECTED_BRIDGE = 389266878372286537904


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def load_bridge_module():
    require(lf_sha256(BRIDGE_PATH) == BRIDGE_SHA256, "bridge source hash drift")
    spec = importlib.util.spec_from_file_location("thm3479_minimal_joint", BRIDGE_PATH)
    require(spec is not None and spec.loader is not None, "bridge module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_bridge_module()


@dataclass(frozen=True)
class BoundaryLabel:
    family: str
    factor_index: int
    mode: str
    branch: int
    side: str
    turn: int = 0


@dataclass(frozen=True)
class LabelledInterval:
    left: int
    right: int
    left_boundary: BoundaryLabel
    right_boundary: BoundaryLabel
    parent_positive: tuple[int, int]
    cuts: tuple[int, ...]


@dataclass(frozen=True)
class VisibleFactorCell:
    e_index: int
    q_index: int
    wrap_turn: int
    left: int
    right: int
    e_parent_positive: tuple[int, int]
    q_parent_positive: tuple[int, int]
    e_cuts: tuple[int, ...]
    q_cuts: tuple[int, ...]
    left_sources: tuple[str, ...]
    right_sources: tuple[str, ...]
    left_factor: int


@dataclass(frozen=True)
class CompleteAtomRecord:
    """Minimum complete schema; values below are a formal finite model only."""

    base: int
    root: int
    owner_sheet: int
    word_sheet: int
    source_sheet: int
    left_horizon: int
    right_horizon: int
    address: tuple[int, int, int]
    left_factor_label: str
    right_factor_label: str
    measure: int


def boundary_text(label: BoundaryLabel) -> str:
    return (
        f"{label.family}:{label.factor_index}:{label.mode}:"
        f"{label.branch}:{label.side}:{label.turn}"
    )


def shifted_boundary(label: BoundaryLabel, turn: int) -> BoundaryLabel:
    return BoundaryLabel(
        label.family,
        label.factor_index,
        label.mode,
        label.branch,
        label.side,
        label.turn + turn,
    )


def labelled_in_comb(
    family: str,
    word: tuple[int, ...],
    t_den: int,
    index: int,
    shift: int,
    mode: str,
) -> list[LabelledInterval]:
    speed = word[index]
    unit = t_den // (182 * speed)
    low = (-13 - 14 * shift) % 182
    output = []
    for branch in range(speed):
        left = (low + 182 * branch) * unit
        right = left + 26 * unit
        enter = BoundaryLabel(family, index, mode, branch, "enter")
        leave = BoundaryLabel(family, index, mode, branch, "leave")
        parent = (index, branch)
        if right <= t_den:
            output.append(LabelledInterval(left, right, enter, leave, parent, ()))
        else:
            cut_high = BoundaryLabel(family, index, mode, branch, "circle-cut-high")
            cut_low = BoundaryLabel(family, index, mode, branch, "circle-cut-low", 1)
            output.append(
                LabelledInterval(left, t_den, enter, cut_high, parent, ())
            )
            output.append(
                LabelledInterval(0, right - t_den, cut_low, leave, parent, ())
            )
    output.sort(key=lambda interval: (interval.left, interval.right))
    return output


def labelled_subtract_comb(
    intervals: list[LabelledInterval],
    family: str,
    factor_index: int,
    mode: str,
    t_den: int,
    speed: int,
    period_denominator: int,
    low: int,
    high: int,
) -> list[LabelledInterval]:
    unit = t_den // (period_denominator * speed)
    step = period_denominator * unit
    window_length = (high - low) * unit
    base = (low % period_denominator) * unit
    output = []
    for interval in intervals:
        first_possible = interval.left - window_length + 1
        branch = -((base - first_possible) // step)
        window_left = base + branch * step
        cursor = interval.left
        cursor_boundary = interval.left_boundary
        cuts = interval.cuts
        while window_left < interval.right:
            window_right = window_left + window_length
            if window_right > cursor:
                if window_left > cursor:
                    cut_enter = BoundaryLabel(
                        family, factor_index, mode, branch, "excluded-enter"
                    )
                    output.append(
                        LabelledInterval(
                            cursor,
                            window_left,
                            cursor_boundary,
                            cut_enter,
                            interval.parent_positive,
                            cuts + (factor_index,),
                        )
                    )
                cursor = window_right
                cursor_boundary = BoundaryLabel(
                    family, factor_index, mode, branch, "excluded-leave"
                )
                cuts = cuts + (factor_index,)
                if cursor >= interval.right:
                    break
            branch += 1
            window_left += step
        if cursor < interval.right:
            output.append(
                LabelledInterval(
                    cursor,
                    interval.right,
                    cursor_boundary,
                    interval.right_boundary,
                    interval.parent_positive,
                    cuts,
                )
            )
    return output


def labelled_build_set(
    family: str,
    word: tuple[int, ...],
    t_den: int,
    pattern: dict[int, str],
    shift: tuple[int, ...],
) -> list[LabelledInterval]:
    positives = [index for index, mode in pattern.items() if mode == "in"]
    require(bool(positives), (family, "no positive comb"))
    start = min(positives, key=lambda index: word[index])
    intervals = labelled_in_comb(
        family, word, t_den, start, shift[start], pattern[start]
    )
    for index, mode in pattern.items():
        if mode == "guard_safe":
            intervals = labelled_subtract_comb(
                intervals,
                family,
                index,
                mode,
                t_den,
                word[index],
                91,
                -13 - 7 * shift[index],
                13 - 7 * shift[index],
            )
    rest = sorted(
        (word[index], index)
        for index, mode in pattern.items()
        if index != start and mode in ("in", "out")
    )
    for _speed, index in rest:
        mode = pattern[index]
        if mode == "out":
            low, high = -13 - 14 * shift[index], 13 - 14 * shift[index]
        else:
            low, high = 13 - 14 * shift[index], 169 - 14 * shift[index]
        intervals = labelled_subtract_comb(
            intervals,
            family,
            index,
            mode,
            t_den,
            word[index],
            182,
            low,
            high,
        )
    return intervals


def interval_projection(
    intervals: list[LabelledInterval],
) -> list[tuple[int, int]]:
    return [(interval.left, interval.right) for interval in intervals]


def labelled_left_cells(
    e_intervals: list[LabelledInterval],
    q_intervals: list[LabelledInterval],
    q_starts: list[int],
    t_den: int,
    nn: int,
    prime: int,
    root: int,
) -> tuple[tuple[VisibleFactorCell, ...], int, int, str]:
    cells = []
    overlap = 0
    total = 0
    digest = sha256()
    frequency = M.X_FREQUENCY
    for e_index, e_interval in enumerate(e_intervals):
        scaled_left = M.R_DILATION * e_interval.left
        start = scaled_left % t_den
        span = M.R_DILATION * (e_interval.right - e_interval.left)
        require(span < t_den, ("labelled sweep span", span, t_den))
        stop = start + span
        q_index = bisect_right(q_starts, start) - 1
        offset = 0
        if q_index < 0:
            q_index = len(q_intervals) - 1
            offset = -t_den
        base_value = pow(
            root, (-frequency * (scaled_left - start)) % nn, prime
        )
        while True:
            q_interval = q_intervals[q_index]
            q_left = q_interval.left + offset
            q_right = q_interval.right + offset
            if q_left >= stop:
                break
            if q_right > start:
                left = max(start, q_left)
                right = min(stop, q_right)
                if left < right:
                    overlap += right - left
                    left_value = pow(root, (-frequency * left) % nn, prime)
                    right_value = pow(root, (-frequency * right) % nn, prime)
                    contribution = base_value * (left_value - right_value) % prime
                    turn = offset // t_den
                    shifted_left = shifted_boundary(q_interval.left_boundary, turn)
                    shifted_right = shifted_boundary(q_interval.right_boundary, turn)
                    left_sources = []
                    right_sources = []
                    if left == start:
                        left_sources.append(
                            "E-scaled:" + boundary_text(e_interval.left_boundary)
                        )
                    if left == q_left:
                        left_sources.append("Q:" + boundary_text(shifted_left))
                    if right == stop:
                        right_sources.append(
                            "E-scaled:" + boundary_text(e_interval.right_boundary)
                        )
                    if right == q_right:
                        right_sources.append("Q:" + boundary_text(shifted_right))
                    cell = VisibleFactorCell(
                        e_index,
                        q_index,
                        turn,
                        left,
                        right,
                        e_interval.parent_positive,
                        q_interval.parent_positive,
                        e_interval.cuts,
                        q_interval.cuts,
                        tuple(left_sources),
                        tuple(right_sources),
                        contribution,
                    )
                    cells.append(cell)
                    total = (total + contribution) % prime
                    digest.update(
                        repr(
                            (
                                cell.e_index,
                                cell.q_index,
                                cell.wrap_turn,
                                cell.left,
                                cell.right,
                                cell.e_parent_positive,
                                cell.q_parent_positive,
                                cell.e_cuts,
                                cell.q_cuts,
                                cell.left_sources,
                                cell.right_sources,
                                cell.left_factor,
                            )
                        ).encode("ascii")
                    )
            q_index += 1
            if q_index == len(q_intervals):
                q_index = 0
                offset += t_den
    return tuple(cells), overlap, total, digest.hexdigest()


def right_endpoint_total(
    e_intervals: list[LabelledInterval],
    word: tuple[int, ...],
    nn: int,
    prime: int,
    root: int,
) -> tuple[int, str]:
    y_frequency = M.X_FREQUENCY + word[M.TARGET_B]
    total = 0
    digest = sha256()
    for e_index, interval in enumerate(e_intervals):
        left_value = pow(
            root, (y_frequency * M.R_DILATION * interval.left) % nn, prime
        )
        right_value = pow(
            root, (y_frequency * M.R_DILATION * interval.right) % nn, prime
        )
        contribution = (left_value - right_value) % prime
        total = (total + contribution) % prime
        digest.update(
            repr(
                (
                    e_index,
                    interval.left,
                    interval.right,
                    boundary_text(interval.left_boundary),
                    boundary_text(interval.right_boundary),
                    interval.parent_positive,
                    interval.cuts,
                    contribution,
                )
            ).encode("ascii")
        )
    return total, digest.hexdigest()


def boundary_census(intervals: list[LabelledInterval]) -> tuple[object, ...]:
    counter = Counter()
    for interval in intervals:
        for label in (interval.left_boundary, interval.right_boundary):
            counter[(label.factor_index, label.mode, label.side)] += 1
    return tuple(sorted((key, value) for key, value in counter.items()))


def actual_factor_controls() -> tuple[object, ...]:
    word = M.to_current(M.U_FULL_REL)
    t_den = 182 * M.lcm_tuple(word)
    nn = M.R_DILATION * t_den
    prime, root, prime_factors, lucas_witnesses = M.FULL_EMBEDDINGS[0]
    M.verify_lucas_certificate(
        prime, prime_factors, lucas_witnesses, "U_full minimal-address prime"
    )
    M.verify_embedding(
        prime, root, nn, M.FULL_NN_FACTORS, "U_full minimal-address embedding"
    )
    embeddings = ((prime, root),)
    zero = (0,) * 9
    _v1, v2, _reps = M.target_representatives()
    labelled_q = labelled_build_set(
        "Q", word, t_den, M.PATTERN_QA, zero
    )
    projected_q = interval_projection(labelled_q)
    fast_q = M.fast_build_set(word, t_den, M.PATTERN_QA, zero)
    require(projected_q == fast_q, "labelled Q projection drift")
    q_starts = [interval.left for interval in labelled_q]
    tabs = M.fast_make_tabs(fast_q, M.X_FREQUENCY, nn, embeddings)
    rows = []
    for label, ell in (("zero", zero), ("v2", v2)):
        labelled_e = labelled_build_set(
            "E", word, t_den, M.PATTERN_E, ell
        )
        projected_e = interval_projection(labelled_e)
        fast_e = M.fast_build_set(word, t_den, M.PATTERN_E, ell)
        require(projected_e == fast_e, (label, "labelled E projection drift"))
        cells, overlap, ax, cell_digest = labelled_left_cells(
            labelled_e,
            labelled_q,
            q_starts,
            t_den,
            nn,
            prime,
            root,
        )
        by, right_digest = right_endpoint_total(
            labelled_e, word, nn, prime, root
        )
        fast_ax, fast_overlap = M.fast_x_sweep(
            fast_e,
            fast_q,
            [left for left, _right in fast_q],
            M.X_FREQUENCY,
            t_den,
            nn,
            embeddings,
            tabs,
        )
        fast_by = M.fast_endpoint_sum(
            fast_e,
            -(M.X_FREQUENCY + word[M.TARGET_B]),
            nn,
            embeddings,
        )
        require((ax,) == fast_ax, (label, "AX cell reconstruction"))
        require((by,) == fast_by, (label, "BY endpoint reconstruction"))
        require(overlap == fast_overlap, (label, "overlap reconstruction"))
        require(len(cells) > len(labelled_e), (label, "no factor refinement"))
        rows.append(
            (
                label,
                len(labelled_e),
                len(labelled_q),
                len(cells),
                overlap,
                ax,
                by,
                cell_digest,
                right_digest,
                boundary_census(labelled_e),
            )
        )
    return (
        prime,
        root,
        nn,
        boundary_census(labelled_q),
        tuple(rows),
    )


def address_of(record: CompleteAtomRecord) -> tuple[int, int, int]:
    """The one address map; crucially it has no ell argument."""

    return record.address


def make_records(matrix: tuple[tuple[int, int], tuple[int, int]]) -> tuple[CompleteAtomRecord, ...]:
    records = []
    for owner_sheet in range(2):
        for word_sheet in range(2):
            source_sheet = owner_sheet ^ word_sheet
            address = Q_H if source_sheet == 0 else Q_Q5
            records.append(
                CompleteAtomRecord(
                    base=2 * owner_sheet + word_sheet,
                    root=0,
                    owner_sheet=owner_sheet,
                    word_sheet=word_sheet,
                    source_sheet=source_sheet,
                    left_horizon=0,
                    right_horizon=1,
                    address=address,
                    left_factor_label=f"L{owner_sheet}",
                    right_factor_label=f"R{word_sheet}",
                    measure=matrix[owner_sheet][word_sheet],
                )
            )
    return tuple(records)


def matrix_margins(
    matrix: tuple[tuple[int, int], tuple[int, int]]
) -> tuple[tuple[int, int], tuple[int, int], int, int]:
    rows = tuple(sum(row) for row in matrix)
    columns = tuple(matrix[0][j] + matrix[1][j] for j in range(2))
    total = sum(rows)
    mixed_haar = matrix[0][0] - matrix[0][1] - matrix[1][0] + matrix[1][1]
    return rows, columns, total, mixed_haar


def character_bank(
    records: tuple[CompleteAtomRecord, ...],
    prime: int,
    zeta: int,
) -> tuple[int, ...]:
    rows = []
    fixed_addresses = tuple(address_of(record) for record in records)
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                ell = (alpha, beta, tau)
                total = 0
                for record, address in zip(records, fixed_addresses):
                    exponent = sum(a * b for a, b in zip(ell, address)) % P
                    # Both labelled factors are evaluated on the record before
                    # this sum.  They are Boolean 1 on their own joint atom.
                    total = (
                        total + record.measure * pow(zeta, exponent, prime)
                    ) % prime
                rows.append(total)
    return tuple(rows)


def inverse_value(
    bank: tuple[int, ...],
    address: tuple[int, int, int],
    prime: int,
    zeta: int,
) -> int:
    total = 0
    index = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                exponent = -(
                    alpha * address[0]
                    + beta * address[1]
                    + tau * address[2]
                ) % P
                total = (total + bank[index] * pow(zeta, exponent, prime)) % prime
                index += 1
    return total * pow(P**3, -1, prime) % prime


def direct_address_mass(
    records: tuple[CompleteAtomRecord, ...],
    address: tuple[int, int, int],
    prime: int,
) -> int:
    return sum(
        record.measure
        for record in records
        if address_of(record) == address
    ) % prime


def reconstruct_2x2(
    rows: tuple[int, int],
    columns: tuple[int, int],
    mixed_haar: int,
    prime: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    total = sum(rows) % prime
    entry_00 = (
        (mixed_haar - total + 2 * rows[0] + 2 * columns[0])
        * pow(4, -1, prime)
    ) % prime
    return (
        (entry_00, (rows[0] - entry_00) % prime),
        ((columns[0] - entry_00) % prime,
         (rows[1] - columns[0] + entry_00) % prime),
    )


def minimal_joint_address_gate(
    prime: int, root: int, nn: int
) -> tuple[object, ...]:
    require(EXPECTED_BRIDGE % 4 == 0, "bridge is not an integral checkerboard")
    lam = EXPECTED_BRIDGE // 4
    positive = ((2 * lam, 0), (0, 2 * lam))
    hostile = ((lam, lam), (lam, lam))
    positive_margins = matrix_margins(positive)
    hostile_margins = matrix_margins(hostile)
    require(
        positive_margins[:3] == hostile_margins[:3],
        "positive/hostile margins differ",
    )
    require(positive_margins[3] == EXPECTED_BRIDGE, "positive Haar bridge")
    require(hostile_margins[3] == 0, "hostile Haar bridge")

    positive_records = make_records(positive)
    hostile_records = make_records(hostile)
    schema = tuple(field.name for field in fields(CompleteAtomRecord))
    required_schema = (
        "base",
        "root",
        "owner_sheet",
        "word_sheet",
        "source_sheet",
        "left_horizon",
        "right_horizon",
        "address",
    )
    require(all(name in schema for name in required_schema), "incomplete atom schema")
    source_tree = ast.parse((ROOT / SCRIPT).read_text(encoding="utf-8"))
    address_nodes = tuple(
        node
        for node in source_tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "address_of"
    )
    require(len(address_nodes) == 1, "address function ambiguity")
    address_node = address_nodes[0]
    address_signature = tuple(argument.arg for argument in address_node.args.args)
    address_names = tuple(
        sorted(
            node.id
            for node in ast.walk(address_node)
            if isinstance(node, ast.Name)
        )
    )
    require(address_signature == ("record",), "address map depends on ell")
    require("ell" not in address_names, "address body depends on ell")

    zeta = pow(root, nn // P, prime)
    require(pow(zeta, P, prime) == 1 and zeta != 1, "bad order-13 root")
    positive_bank = character_bank(positive_records, prime, zeta)
    hostile_bank = character_bank(hostile_records, prime, zeta)
    positive_values = tuple(
        inverse_value(positive_bank, address, prime, zeta)
        for address in (Q_H, Q_Q5)
    )
    hostile_values = tuple(
        inverse_value(hostile_bank, address, prime, zeta)
        for address in (Q_H, Q_Q5)
    )
    require(
        positive_values
        == tuple(direct_address_mass(positive_records, q, prime) for q in (Q_H, Q_Q5)),
        "positive DFT/direct mismatch",
    )
    require(
        hostile_values
        == tuple(direct_address_mass(hostile_records, q, prime) for q in (Q_H, Q_Q5)),
        "hostile DFT/direct mismatch",
    )
    positive_bridge = (positive_values[0] - positive_values[1]) % prime
    hostile_bridge = (hostile_values[0] - hostile_values[1]) % prime
    require(positive_bridge == EXPECTED_BRIDGE, "frozen bridge not recovered")
    require(hostile_bridge == 0, "hostile address bridge")

    reconstructed_positive = reconstruct_2x2(
        positive_margins[0], positive_margins[1], positive_margins[3], prime
    )
    reconstructed_hostile = reconstruct_2x2(
        hostile_margins[0], hostile_margins[1], hostile_margins[3], prime
    )
    require(
        reconstructed_positive
        == tuple(tuple(value % prime for value in row) for row in positive),
        "positive reconstruction",
    )
    require(
        reconstructed_hostile
        == tuple(tuple(value % prime for value in row) for row in hostile),
        "hostile reconstruction",
    )

    checkerboard = ((1, -1), (-1, 1))
    checker_rows = tuple(sum(row) for row in checkerboard)
    checker_columns = tuple(
        checkerboard[0][j] + checkerboard[1][j] for j in range(2)
    )
    address_pairing = (
        checkerboard[0][0]
        - checkerboard[0][1]
        - checkerboard[1][0]
        + checkerboard[1][1]
    )
    require(checker_rows == checker_columns == (0, 0), "checker margins")
    require(address_pairing == 4, "joint address misses checkerboard")
    row_only = tuple(
        sum(checkerboard[i][j] * (1 if i == 0 else -1)
            for i in range(2) for j in range(2))
        for _choice in range(1)
    )[0]
    column_only = sum(
        checkerboard[i][j] * (1 if j == 0 else -1)
        for i in range(2) for j in range(2)
    )
    require(row_only == column_only == 0, "one-factor address sees kernel")

    return (
        lam,
        positive,
        hostile,
        positive_margins,
        hostile_margins,
        schema,
        address_signature,
        address_names,
        positive_values,
        hostile_values,
        positive_bridge,
        hostile_bridge,
        digest_json(positive_bank),
        digest_json(hostile_bank),
        checkerboard,
        checker_rows,
        checker_columns,
        address_pairing,
        row_only,
        column_only,
        reconstructed_positive,
        reconstructed_hostile,
    )


def security_certificate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    forbidden_calls = {"eval", "exec", "compile", "__import__"}
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            bad.append("Assert")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in forbidden_calls
        ):
            bad.append(node.func.id)
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), tuple(bad)


def main() -> None:
    require(lf_sha256(PREVIOUS_PATH) == PREVIOUS_SHA256, "previous probe drift")
    theorem_hashes = []
    for theorem_id, path, expected_hash in THEOREM_PINS:
        actual_hash = lf_sha256(path)
        require(actual_hash == expected_hash, (theorem_id, "hash drift", actual_hash))
        theorem_hashes.append((theorem_id, actual_hash))

    geometry = actual_factor_controls()
    prime, root, nn = geometry[:3]
    gate = minimal_joint_address_gate(prime, root, nn)
    security = security_certificate(ROOT / SCRIPT)

    required_actual = (
        "base",
        "root",
        "owner_sheet",
        "word_sheet",
        "source_sheet",
        "left_horizon",
        "right_horizon",
        "address",
    )
    visible_actual = (
        "circle_interval",
        "E_factor_lineage",
        "Q_factor_lineage",
        "E_component_index",
        "Q_component_index",
        "periodic_wrap_turn",
        "boundary_provenance",
    )
    actual_missing = tuple(field for field in required_actual if field not in visible_actual)
    require(actual_missing == required_actual, "manufactured ancestry projection")

    semantic = (
        "FINITE-EXACT minimum 2x2 bridge-coordinate/API theorem",
        geometry,
        gate,
        visible_actual,
        actual_missing,
        "one joint parity/address scalar resolves the 1D checkerboard kernel",
        "the actual U_full endpoint geometry supplies no map to that scalar",
        "schema witness only; no K4/current/row/LRC consequence",
    )

    print("LRC U_FULL MINIMUM JOINT-ADDRESS BRIDGE GATE")
    print("status=FINITE-EXACT API/minimal-kernel result; not a U_full ancestry realization; LRC(14) OPEN")
    print(f"dependency_hashes={(BRIDGE_SHA256, PREVIOUS_SHA256, tuple(theorem_hashes))}")
    print(f"embedding=(prime={prime},root={root},order={nn})")
    print(f"actual_factor_controls=(q_boundary_census={geometry[3]},rows={geometry[4]})")
    print(f"visible_actual_cell_fields={visible_actual}")
    print(f"required_ancestry_projection_fields={required_actual}")
    print(f"actual_projection_gate=FAIL_MISSING_{actual_missing}")
    print(f"complete_schema_fields={gate[5]}")
    print(f"ell_independent_address=(signature={gate[6]},body_names={gate[7]})")
    print(f"checkerboard_minimality=(matrix={gate[14]},row_sums={gate[15]},column_sums={gate[16]},joint_address_value={gate[17]},row_only_value={gate[18]},column_only_value={gate[19]},kernel_dimension=1)")
    print(f"frozen_bridge_quarter={gate[0]}")
    print(f"positive_matrix={gate[1]}; hostile_matrix={gate[2]}")
    print(f"shared_margins=(rows={gate[3][0]},columns={gate[3][1]},total={gate[3][2]})")
    print(f"mixed_haar=(positive={gate[3][3]},hostile={gate[4][3]})")
    print(f"inverse_address_values=(positive_H_q5={gate[8]},hostile_H_q5={gate[9]})")
    print(f"bridge_only_test=(positive={gate[10]},hostile={gate[11]},expected={EXPECTED_BRIDGE})")
    print(f"character_bank_sha256=(positive={gate[12]},hostile={gate[13]})")
    print(f"reconstruction_from_margins_plus_one_joint_scalar=(positive={gate[20]},hostile={gate[21]})")
    print("positive_control=linked factors are multiplied on four complete records before summation; one fixed parity address map recovers the exact frozen q_H-q_q5 bridge")
    print("hostile_control=the flat coupling has identical row/column marginals and every separate-factor scalar, but its q_H-q_q5 bridge is zero")
    print("minimum_coordinate=one joint pairing bit, equivalently one mixed-Haar coefficient nonzero on [[1,-1],[-1,1]]; no left-only or right-only address can see it")
    print("calibration_scope=the inherited bridge fixes the positive model's checkerboard coefficient; this proves API sufficiency/minimality, not independent recovery from U_full geometry")
    print("first_unresolved_map=visible factor-labelled circle cells -> THM-2471 (base,root,owner-sheet,word-sheet,source-sheet,horizons,address)")
    print("terminology=the actual cells share circle geometry only; the formal four-atom model is a common-base schema witness, never an asserted U_full ancestry stalk")
    print("no_K4_evaluated=True")
    print("nonconsequence=no physical current, chronology, grouped C(a;X,m), all-unit B(q), scalar-row exclusion, or LRC(14)")
    print(f"security_ast={security}")
    print(f"semantic_sha256={digest_json(semantic)}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
