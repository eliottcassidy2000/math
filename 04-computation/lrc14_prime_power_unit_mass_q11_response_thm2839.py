#!/usr/bin/env python3
"""Exact companion for the THM-2839 prime-power unit-mass q11 theorem.

The physical atom, semantic-cell witness, and endpoint constructors are the
pinned THM-2806/2782 objects.  Only the depth-five outer Boolean word is
varied through the three canonical patterns QA, QB, QAB.  Besides the three
same-word profiles, the script records the complete directed cross-word
matrix and every union coarsening of these three disjoint atoms.

It reconstructs the inherited THM-2835 fixture and verifies the complete
q11 ancestry response/provenance census.  It makes no row-decrement or
LRC(14) claim.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter
from hashlib import sha256
import importlib.util
from itertools import combinations
from pathlib import Path
from struct import pack
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    COMP / "lrc14_169_twist_variance_opus_20260726.py":
        "c0aa9c55c3e7d38dc28b4705e58c776a979f17d2874d1b762f96b97d2364e5e9",
    COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py":
        "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation
import lrc14_169_twist_variance_opus_20260726 as atlas


# One inherited loader pins raw bytes and therefore rejects a CRLF checkout.
# The LF hashes above are immutable; load the verified source directly.
def crlf_safe_present_loader():
    path = COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    spec = importlib.util.spec_from_file_location("q11_outer_present", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


allocation.physical.target.load_present_module = crlf_safe_present_loader

P = 13
Q = 11
DEPTH = P**5
WORD_ORDER = ("QA", "QB", "QAB")


def circular_shift_interval(interval, shift, period):
    left = (interval[0] + shift) % period
    right = (interval[1] + shift) % period
    require(left < right, "selected interval crossed the circle seam")
    return left, right


def containing_window(interval, offset, ancestry, period, pieces, starts):
    """Return the x-window induced by one whole-cylinder word hit."""
    low = interval[0] + offset + ancestry * period
    high = interval[1] + offset + ancestry * period
    midpoint = (low + high) // 2
    index = bisect_right(starts, midpoint // DEPTH) - 1
    require(index >= 0, "declared whole-cylinder hit has no interval")
    left, right = pieces[index]
    require(
        left * DEPTH <= low and high <= right * DEPTH,
        "declared whole-cylinder hit crosses a word seam",
    )
    return (
        left * DEPTH - ancestry * period - offset,
        right * DEPTH - ancestry * period - offset,
    )


def whole_cylinder_hit(interval, offset, ancestry, period, pieces, starts):
    """Test one label and certify that midpoint membership is not a seam hit."""
    low = interval[0] + offset + ancestry * period
    high = interval[1] + offset + ancestry * period
    midpoint = (low + high) // 2
    index = bisect_right(starts, midpoint // DEPTH) - 1
    midpoint_hit = (
        index >= 0
        and pieces[index][0] * DEPTH <= midpoint
        < pieces[index][1] * DEPTH
    )
    whole_hit = (
        midpoint_hit
        and pieces[index][0] * DEPTH <= low
        and high <= pieces[index][1] * DEPTH
    )
    require(
        midpoint_hit == whole_hit,
        "midpoint/whole-cylinder word census differs",
    )
    return whole_hit


def word_flags(interval, offset, period, pieces, starts):
    """All depth-five labels; assert midpoint and whole-cylinder flags agree."""
    result = bytearray(DEPTH)
    for ancestry in range(DEPTH):
        result[ancestry] = whole_cylinder_hit(
            interval, offset, ancestry, period, pieces, starts
        )
    return result


def overlap_count(left, right):
    return sum(a & b for a, b in zip(left, right))


def interval_overlap_measure(left, right):
    """Exact overlap length of two sorted half-open interval unions."""
    i = j = total = 0
    while i < len(left) and j < len(right):
        a0, a1 = left[i]
        b0, b1 = right[j]
        total += max(0, min(a1, b1) - max(a0, b0))
        if a1 <= b1:
            i += 1
        else:
            j += 1
    return total


def safe_comb_contains(interval, module, speed, denominator, lo, hi):
    """Exact containment in one periodic safe-comb complement."""
    require(
        module.T % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = module.T // (denominator * speed)
    step = denominator * unit
    length = (hi - lo) * unit
    base = (lo % denominator) * unit
    left, right = interval
    k0 = (left - base) // step
    for k in range(k0 - 1, k0 + 3):
        forbidden_left = base + k * step
        forbidden_right = forbidden_left + length
        if max(left, forbidden_left) < min(right, forbidden_right):
            return False
    return True


def comb_intervals(period, speed, denominator, lo, hi):
    """Materialize one strict periodic comb on an exactly resolving grid."""
    require(
        period % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = period // (denominator * speed)
    step = denominator * unit
    length = (hi - lo) * unit
    base = (lo % denominator) * unit
    pieces = []
    for index in range(speed):
        left = base + index * step
        right = left + length
        if right <= period:
            pieces.append((left, right))
        else:
            pieces.extend(((left, period), (0, right - period)))
    return tuple(sorted(pieces))


def union_flags(banks):
    return bytearray(any(values) for values in zip(*banks))


def strongly_connected_sizes(nodes, edges):
    """Exact Kosaraju census for one small finite directed relation."""
    reverse = {node: set() for node in nodes}
    for source, targets in edges.items():
        for target in targets:
            reverse[target].add(source)

    seen = set()
    order = []
    for root in nodes:
        if root in seen:
            continue
        stack = [(root, False)]
        seen.add(root)
        while stack:
            node, closing = stack.pop()
            if closing:
                order.append(node)
                continue
            stack.append((node, True))
            for target in edges[node]:
                if target not in seen:
                    seen.add(target)
                    stack.append((target, False))

    seen.clear()
    sizes = []
    for root in reversed(order):
        if root in seen:
            continue
        size = 0
        stack = [root]
        seen.add(root)
        while stack:
            node = stack.pop()
            size += 1
            for target in reverse[node]:
                if target not in seen:
                    seen.add(target)
                    stack.append(target)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def main():
    _module, full_module, details, e3, clocks, q_pairs = (
        allocation.build_geometry()
    )
    period = full_module.T
    allocation_unit = period // P
    atom = allocation.ATOM_INTERVAL
    target0 = tuple(x + allocation.physical.SHIFT for x in atom)
    target11 = circular_shift_interval(
        target0, Q * allocation_unit, period
    )

    old = allocation.physical.relative.lift.m.core.old
    require(
        period == atlas.T_DEN == old.T == 297836897838480
        and old.rail.DEPTH == DEPTH
        and old.base.W == atlas.W
        and old.base.PAT_QA == atlas.PAT_QA
        and old.host.PAT_QB == atlas.PAT_QB,
        "outer-word constructors no longer share the pinned typed row",
    )
    patterns = {
        "QA": atlas.PAT_QA,
        "QB": atlas.PAT_QB,
        "QAB": atlas.PAT_QAB,
    }
    require(
        tuple(
            (patterns[name][7], patterns[name][8])
            for name in WORD_ORDER
        ) == (("in", "out"), ("out", "in"), ("in", "in")),
        "canonical a/b word signatures changed",
    )

    # The physical q11 cell and all data downstream of the local six factors
    # are unchanged before the outer word is imposed.
    witness = (0, 5, 1)
    factors = allocation.section_factors(
        full_module, e3, clocks, *witness
    )
    require(
        all(
            allocation.contained(interval, factor)
            for interval in (atom, target0, target11)
            for factor in factors
        ),
        "q11 semantic-reselected physical cell changed",
    )

    target7 = circular_shift_interval(
        target0, 7 * allocation_unit, period
    )

    def physical_signature(interval, s, t, clock):
        return (
            allocation.contained(interval, e3),
            allocation.contained(interval, clocks[clock]),
            safe_comb_contains(
                interval,
                full_module,
                full_module.W[1],
                182,
                -14 * s - 13,
                -14 * s + 13,
            ),
            safe_comb_contains(
                interval,
                full_module,
                full_module.W[2],
                182,
                -14 * t - 13,
                -14 * t + 13,
            ),
            safe_comb_contains(
                interval,
                full_module,
                full_module.C2,
                182,
                14 * s - 13,
                14 * s + 13,
            ),
            safe_comb_contains(
                interval,
                full_module,
                full_module.C3,
                182,
                14 * t - 13,
                14 * t + 13,
            ),
        )

    q11_cells = tuple(
        (s, t, clock)
        for s in allocation.COMMON_S
        for t in allocation.COMMON_T
        for clock in range(7)
        if all(physical_signature(atom, s, t, clock))
        and all(physical_signature(target0, s, t, clock))
        and all(physical_signature(target11, s, t, clock))
    )
    q7_failure_census = Counter(
        tuple(
            index
            for index, bit in enumerate(
                physical_signature(target7, *cell)
            )
            if not bit
        )
        for cell in q11_cells
    )
    q7_fully_physical_cells = tuple(
        (s, t, clock)
        for s in allocation.COMMON_S
        for t in allocation.COMMON_T
        for clock in range(7)
        if all(physical_signature(atom, s, t, clock))
        and all(physical_signature(target0, s, t, clock))
        and all(physical_signature(target7, s, t, clock))
    )
    require(
        q11_cells
        == tuple(
            (s, t, 1)
            for s in allocation.COMMON_S
            for t in range(5, 12)
        )
        and q7_failure_census
        == {(0,): 35, (0, 3): 14, (0, 2): 10, (0, 2, 3): 4}
        and q7_fully_physical_cells == ()
        and tuple(
            allocation.contained(target7, clock)
            for clock in clocks
        )
        == (False, True, False, False, False, False, False),
        "q7 physical failure stratification changed",
    )

    source_semantic = (
        allocation.physical.relative.private.delayed_carry_pair(
            ((*atom, 1),), q_pairs[1], {}
        )
    )
    target_semantic = (
        allocation.physical.relative.private.delayed_carry_pair(
            ((*target11, 1),), q_pairs[1], {}
        )
    )
    require(
        source_semantic[12] == (0, allocation.SEMANTIC_VALUE)
        and target_semantic[6] == (0, allocation.SEMANTIC_VALUE),
        "delayed semantic values changed",
    )

    endpoint_base = allocation.endpoint_base
    endpoint_value = (
        231164267889491750,
        630230755085920022,
    )
    endpoint_rows = {}
    gauge_endpoint_rows = {}
    target_atom = ((*target0, 1),)
    for address in (
        allocation.RIGHT_ORIGIN,
        allocation.add(
            allocation.RIGHT_ORIGIN, allocation.TARGET_STEP
        ),
    ):
        ell = endpoint_base.REPS[address]
        present = tuple(
            endpoint_base.endpoint.build_set(endpoint_base.PAT_E3, ell)
        )
        starts = tuple(left for left, _right in present)
        rows = []
        for q in (0, Q):
            shifted = allocation.physical.overlap.shift_weighted(
                target_atom, q * allocation_unit
            )
            restricted = allocation.indexed_weighted_intersection(
                shifted, present, starts
            )
            values = tuple(
                allocation.endpoint_sum(
                    restricted, -endpoint_base.Y0, embedding
                )
                for embedding in endpoint_base.endpoint.MODS
            )
            mass = sum(
                (right - left) * weight
                for left, right, weight in restricted
            )
            rows.append((q, len(restricted), mass, values))
        endpoint_rows[address] = tuple(rows)
        ell4 = tuple(
            (
                ell[index]
                + 4 * endpoint_base.WMOD[index]
            ) % P
            for index in range(9)
        )
        gauge_present = tuple(
            endpoint_base.endpoint.build_set(
                endpoint_base.PAT_E3, ell4
            )
        )
        require(
            gauge_present
            == allocation.physical.overlap.shift_union(
                present, -4 * allocation_unit
            ),
            "four-step representative gauge covariance changed",
        )
        gauge_starts = tuple(left for left, _right in gauge_present)
        gauge_carrier = allocation.physical.overlap.shift_weighted(
            target_atom, 7 * allocation_unit
        )
        gauge_restricted = allocation.indexed_weighted_intersection(
            gauge_carrier, gauge_present, gauge_starts
        )
        gauge_endpoint_rows[address] = (
            ell4,
            len(gauge_restricted),
            sum(
                (right - left) * weight
                for left, right, weight in gauge_restricted
            ),
            tuple(
                allocation.endpoint_sum(
                    gauge_restricted, -endpoint_base.Y0, embedding
                )
                for embedding in endpoint_base.endpoint.MODS
            ),
        )
    expected_endpoint = (
        (0, 1, 26444880, endpoint_value),
        (11, 1, 26444880, endpoint_value),
    )
    require(
        all(row == expected_endpoint for row in endpoint_rows.values()),
        "two-right-origin q11 endpoint data changed",
    )
    require(
        all(
            row[1:] == (1, 26444880, endpoint_value)
            for row in gauge_endpoint_rows.values()
        ),
        "four-step marked-gauge q7 endpoint copy changed",
    )

    pieces = {
        name: tuple(old.base.build_set(pattern, old.base.ZELL))
        for name, pattern in patterns.items()
    }
    require(
        all(
            interval_overlap_measure(pieces[left], pieces[right]) == 0
            for left, right in combinations(WORD_ORDER, 2)
        ),
        "canonical outer-word atoms ceased to be disjoint",
    )
    starts = {
        name: tuple(left for left, _right in pieces[name])
        for name in WORD_ORDER
    }
    source = {
        name: word_flags(atom, 0, period, pieces[name], starts[name])
        for name in WORD_ORDER
    }
    target = {
        (name, q): word_flags(
            atom,
            allocation.physical.SHIFT + q * allocation_unit,
            period,
            pieces[name],
            starts[name],
        )
        for name in WORD_ORDER
        for q in range(P)
    }

    source_counts = {
        name: sum(source[name]) for name in WORD_ORDER
    }
    target_counts = {
        name: tuple(sum(target[name, q]) for q in range(P))
        for name in WORD_ORDER
    }
    common_counts = {
        name: tuple(
            overlap_count(source[name], target[name, q])
            for q in range(P)
        )
        for name in WORD_ORDER
    }
    require(
        source_counts == {"QA": 0, "QB": 66099, "QAB": 10786}
        and target_counts["QA"]
        == (0, 10787, 10787, 10787, 10786, 10783, 11232,
            0, 10782, 10785, 10785, 10787, 10788)
        and common_counts["QA"] == (0,) * P
        and target_counts["QB"]
        == (66099, 0, 0, 0, 0, 0, 0, 65652, 0, 0, 0, 0, 0)
        and common_counts["QB"]
        == (66099, 0, 0, 0, 0, 0, 0, 65612, 0, 0, 0, 0, 0)
        and target_counts["QAB"]
        == (10786, 0, 0, 0, 0, 0, 0, 11232, 0, 0, 0, 0, 0)
        and common_counts["QAB"]
        == (10786, 0, 0, 0, 0, 0, 0, 10779, 0, 0, 0, 0, 0),
        "same-word ancestry census changed",
    )

    directed = {}
    for source_name in WORD_ORDER:
        for target_name in WORD_ORDER:
            directed[source_name, target_name] = tuple(
                overlap_count(
                    source[source_name], target[target_name, q]
                )
                for q in range(P)
            )
    require(
        directed["QB", "QA"]
        == (0, 0, 0, 0, 0, 0, 449, 0, 449, 449, 449, 449, 449)
        and directed["QB", "QAB"]
        == (0, 0, 0, 0, 0, 0, 0, 449, 0, 0, 0, 0, 0)
        and directed["QAB", "QA"]
        == (0, 10786, 10785, 10785, 10783, 10780, 10779,
            0, 10329, 10329, 10328, 10328, 10328)
        and all(
            directed["QA", target_name] == (0,) * P
            for target_name in WORD_ORDER
        )
        and directed["QAB", "QB"] == (0,) * P,
        "directed cross-word ancestry census changed",
    )

    # Myhill--Nerode-style finite response state: record the three source
    # word bits, followed by all 13 x 3 target word bits in residue-major
    # order.  The code uses only whole-cylinder predicates already checked
    # above, so it cannot merge a midpoint seam with an open word cell.
    response_codes = []
    for ancestry in range(DEPTH):
        code = sum(
            source[name][ancestry] << index
            for index, name in enumerate(WORD_ORDER)
        )
        bit = len(WORD_ORDER)
        for q in range(P):
            for name in WORD_ORDER:
                code |= target[name, q][ancestry] << bit
                bit += 1
        response_codes.append(code)
    response_state_counts = Counter(response_codes)
    expected_response_state_counts = {
        0: 294357,
        18: 38,
        100: 1,
        4708: 2,
        37476: 3,
        299620: 1,
        33554432: 40,
        33554450: 65612,
        69505636: 450,
        1277465188: 1,
        549755813888: 1,
        618475290624: 2,
        627065225216: 1,
        628138967040: 3,
        628342390784: 1,
        628342390802: 449,
        628342685696: 1,
        628342690304: 1,
        628342690368: 1,
        628342690404: 10328,
    }
    require(
        response_state_counts == expected_response_state_counts,
        "complete source/target response-state census changed",
    )
    adjacent_response_edges = Counter(
        (
            response_codes[ancestry],
            response_codes[(ancestry + 1) % DEPTH],
        )
        for ancestry in range(DEPTH)
    )
    response_successors = {
        code: {
            target_code
            for (source_code, target_code) in adjacent_response_edges
            if source_code == code
        }
        for code in response_state_counts
    }
    response_scc_sizes = strongly_connected_sizes(
        tuple(response_state_counts), response_successors
    )
    require(
        sum(len(targets) for targets in response_successors.values()) == 44
        and sum(
            len(targets) == 1
            for targets in response_successors.values()
        ) == 15
        and max(map(len, response_successors.values())) == 12
        and response_scc_sizes == (20,),
        "response-code successor graph changed",
    )

    # The normal h=9 map is addition by nine on the combined depth-six
    # digit n=13*a+q.  On quotient states it keeps the response code for
    # q=0,1,2,3 and uses the adjacent-code relation for q=4,...,12.
    natural_nodes = tuple(
        (q, code) for q in range(P) for code in response_state_counts
    )
    natural_edges = {}
    for q, code in natural_nodes:
        target_q = (q + 9) % P
        carry = (q + 9) // P
        natural_edges[q, code] = {
            (
                target_q,
                target_code if carry else code,
            )
            for target_code in (
                response_successors[code] if carry else (code,)
            )
        }
    natural_scc_sizes = strongly_connected_sizes(
        natural_nodes, natural_edges
    )
    require(
        len(natural_nodes) == 260
        and sum(len(targets) for targets in natural_edges.values()) == 476
        and natural_scc_sizes == (260,),
        "natural-extension response automaton changed",
    )

    qb_to_qa_q11 = tuple(
        ancestry
        for ancestry, (source_hit, target_hit) in enumerate(
            zip(source["QB"], target["QA", Q])
        )
        if source_hit and target_hit
    )
    source_qb_residue_counts = tuple(
        sum(
            source["QB"][ancestry]
            for ancestry in range(residue, DEPTH, P**2)
        )
        for residue in range(P**2)
    )
    target_qa_q11_residue_counts = tuple(
        sum(
            target["QA", Q][ancestry]
            for ancestry in range(residue, DEPTH, P**2)
        )
        for residue in range(P**2)
    )
    source_qb_residue_support = tuple(
        residue
        for residue, count in enumerate(source_qb_residue_counts)
        if count
    )
    target_qa_q11_residue_support = tuple(
        residue
        for residue, count in enumerate(target_qa_q11_residue_counts)
        if count
    )
    label_digest = sha256()
    for ancestry in qb_to_qa_q11:
        label_digest.update(pack(">I", ancestry))
    gap_census = tuple(sorted(Counter(
        right - left
        for left, right in zip(
            qb_to_qa_q11, qb_to_qa_q11[1:]
        )
    ).items()))
    require(
        len(qb_to_qa_q11) == 449
        and qb_to_qa_q11[0] == 59306
        and qb_to_qa_q11[-1] == 311961
        and source_qb_residue_support == tuple(range(12, 157))
        and target_qa_q11_residue_support
        == tuple(range(0, 11)) + tuple(range(156, 169))
        and (
            set(source_qb_residue_support)
            & set(target_qa_q11_residue_support)
        ) == {156}
        and source_qb_residue_counts[156] == 455
        and target_qa_q11_residue_counts[156] == 449
        and {ancestry % 169 for ancestry in qb_to_qa_q11} == {156}
        and label_digest.hexdigest()
        == "2bba5b07079eec3c359e155efa930c204c92cc599260e4aacf0f652e48629212"
        and gap_census
        == ((169, 378), (845, 2), (1014, 7), (1183, 10),
            (1352, 1), (1521, 7), (1690, 2), (1859, 2),
            (2028, 3), (2197, 11), (2366, 2), (2704, 1),
            (2873, 2), (3211, 2), (3549, 2), (3887, 6),
            (4394, 1), (4732, 2), (6422, 2), (7267, 1),
            (7436, 2), (8619, 1), (11323, 1)),
        "q11 QB-to-QA whole-cylinder label progression changed",
    )

    # Follow exactly these 449 cross-word labels along THM-2829's natural
    # q11,h9 lift.  Since 11+9=20=7+13, the target residue is q=7 and the
    # ancestry label receives the compulsory +1 borrow.  This does not
    # manufacture the missing physical action; it only types its exact
    # label-level semantic landing if that action is supplied.
    borrow_word_profiles = {
        delta: {
            name: sum(
                target[name, 7][(ancestry + delta) % DEPTH]
                for ancestry in qb_to_qa_q11
            )
            for name in WORD_ORDER
        }
        for delta in (-2, -1, 0, 1, 2)
    }
    reverse_exceptions = tuple(
        ancestry
        for ancestry in qb_to_qa_q11
        if not target["QB", 7][(ancestry - 1) % DEPTH]
    )
    require(
        borrow_word_profiles
        == {
            -2: {"QA": 0, "QB": 446, "QAB": 0},
            -1: {"QA": 0, "QB": 447, "QAB": 0},
            0: {"QA": 0, "QB": 0, "QAB": 449},
            1: {"QA": 0, "QB": 0, "QAB": 449},
            2: {"QA": 0, "QB": 0, "QAB": 449},
        }
        and reverse_exceptions == (107978, 154622),
        "q11 cross-word borrow semantic profile changed",
    )
    horn_state = 628342390802
    horn_successor_state = 628342690404
    horn_state_counts = Counter(
        response_codes[ancestry] for ancestry in qb_to_qa_q11
    )
    horn_successor_counts = Counter(
        response_codes[(ancestry + 1) % DEPTH]
        for ancestry in qb_to_qa_q11
    )
    horn_successor_predecessors = {
        source_code: multiplicity
        for (source_code, target_code), multiplicity
        in adjacent_response_edges.items()
        if target_code == horn_successor_state
    }
    require(
        horn_state_counts == {horn_state: 449}
        and response_state_counts[horn_state] == 449
        and horn_successor_counts == {horn_successor_state: 449}
        and response_state_counts[horn_successor_state] == 10328
        and adjacent_response_edges[
            horn_state, horn_successor_state
        ] == 449
        and response_successors[horn_state] == {horn_successor_state}
        and horn_successor_predecessors
        == {
            549755813888: 1,
            618475290624: 2,
            627065225216: 1,
            628138967040: 3,
            628342390784: 1,
            628342390802: 449,
            628342685696: 1,
            628342690304: 1,
            628342690368: 1,
            628342690404: 9868,
        }
        and all(
            response_codes[ancestry]
            != response_codes[(ancestry + 1) % DEPTH]
            for ancestry in qb_to_qa_q11
        ),
        "449-horn response-state isolation changed",
    )

    # A target-state predecessor bit is the minimal set-theoretic memory for
    # the horn.  The full two-digit ancestry residue makes it intrinsic on
    # this finite bank: C at residue 157 is exactly H+1.  The low residue
    # alone does not suffice, and a globally right-resolving response lift
    # does not stay small.
    horn_successor_labels = {
        (ancestry + 1) % DEPTH for ancestry in qb_to_qa_q11
    }
    successor_state_labels = tuple(
        ancestry
        for ancestry, code in enumerate(response_codes)
        if code == horn_successor_state
    )
    successor_residue_census = Counter(
        ancestry % (P**2) for ancestry in successor_state_labels
    )
    expected_successor_residue_census = {
        0: 449,
        1: 449,
        2: 449,
        3: 448,
        4: 448,
        5: 449,
        6: 450,
        7: 450,
        8: 450,
        9: 449,
        10: 449,
        157: 449,
        158: 449,
        159: 450,
        160: 450,
        161: 450,
        162: 449,
        163: 448,
        164: 448,
        165: 448,
        166: 449,
        167: 449,
        168: 449,
    }
    successor_low_residue_census = Counter(
        ancestry % P for ancestry in successor_state_labels
    )
    require(
        successor_residue_census
        == expected_successor_residue_census
        and horn_successor_labels
        == {
            ancestry
            for ancestry in successor_state_labels
            if ancestry % (P**2) == 157
        }
        and successor_low_residue_census[1] == 898,
        "two-digit horn-provenance selector changed",
    )

    pair_successors = {}
    pair_frequencies = Counter()
    for ancestry in range(DEPTH):
        pair = (
            response_codes[(ancestry - 1) % DEPTH],
            response_codes[ancestry],
        )
        pair_frequencies[pair] += 1
        pair_successors.setdefault(pair, set()).add(
            response_codes[(ancestry + 1) % DEPTH]
        )
    require(
        len(pair_successors) == 44
        and sum(map(len, pair_successors.values())) == 68
        and sum(
            len(targets) == 1
            for targets in pair_successors.values()
        ) == 39
        and max(map(len, pair_successors.values())) == 12
        and pair_frequencies[horn_state, horn_successor_state] == 449
        and pair_successors[horn_state, horn_successor_state]
        == {horn_successor_state},
        "one-step-memory response lift changed",
    )

    maximal_proper_period = P**4
    proper_period_matches = sum(
        response_codes[ancestry]
        == response_codes[(ancestry + maximal_proper_period) % DEPTH]
        for ancestry in range(DEPTH)
    )
    first_period_failure = next(
        ancestry
        for ancestry in range(DEPTH)
        if response_codes[ancestry]
        != response_codes[
            (ancestry + maximal_proper_period) % DEPTH
        ]
    )
    require(
        proper_period_matches == 281510
        and first_period_failure == 30601
        and (
            response_codes[first_period_failure],
            response_codes[
                (
                    first_period_failure
                    + maximal_proper_period
                ) % DEPTH
            ],
        )
        == (0, 33554450),
        "response word unexpectedly descended below depth five",
    )

    # Exact cyclic-period audit.  Since DEPTH=13^5, every period divides one
    # of 1,13,...,13^5.  The horn indicator is exactly the H-state indicator,
    # and H has only the H->C successor, so the marked-edge indicator agrees
    # with it pointwise.  Full least period for any one of these observables
    # rules out every smaller deterministic invertible phase realization:
    # if two phase states on the distinguished orbit agreed, their future
    # response words would agree and their separation would be a period.
    period_candidates = tuple(P**exponent for exponent in range(6))
    horn_indicator = tuple(
        code == horn_state for code in response_codes
    )
    marked_edge_indicator = tuple(
        response_codes[ancestry] == horn_state
        and response_codes[(ancestry + 1) % DEPTH]
        == horn_successor_state
        for ancestry in range(DEPTH)
    )

    def cyclic_mismatch_table(sequence):
        return tuple(
            (
                shift,
                sum(
                    sequence[index]
                    != sequence[(index + shift) % DEPTH]
                    for index in range(DEPTH)
                ),
                next(
                    (
                        index
                        for index in range(DEPTH)
                        if sequence[index]
                        != sequence[(index + shift) % DEPTH]
                    ),
                    None,
                ),
            )
            for shift in period_candidates
        )

    response_period_table = cyclic_mismatch_table(response_codes)
    horn_period_table = cyclic_mismatch_table(horn_indicator)
    edge_period_table = cyclic_mismatch_table(marked_edge_indicator)
    require(
        response_period_table
        == (
            (1, 2052, 59161),
            (13, 14600, 59149),
            (169, 24857, 58993),
            (2197, 116235, 56965),
            (28561, 89783, 30601),
            (371293, 0, None),
        )
        and horn_period_table
        == (
            (1, 898, 59305),
            (13, 898, 59293),
            (169, 142, 59137),
            (2197, 686, 57109),
            (28561, 532, 30745),
            (371293, 0, None),
        )
        and horn_indicator == marked_edge_indicator
        and edge_period_table == horn_period_table,
        "response/horn/edge least cyclic period changed",
    )

    # The horn is spectrally stronger than merely having full period.
    # Every a in L is 156+13^2*k.  Thus its depth-five mask polynomial is
    #
    #     F(z)=z^156 G(z^(13^2)),        G(1)=449.
    #
    # At a depth-five character m, G is evaluated at a root whose order is
    # one of 1,13,169,2197.  If a nontrivial p-power root were a zero, the
    # corresponding cyclotomic polynomial would divide G over Z, forcing
    # p=Phi_(p^j)(1) to divide G(1)=449, contradiction.  The order-one
    # value is 449.  Hence every one of the 13^5 Fourier characters,
    # including all phi(13^5)=342732 unit characters, survives.
    horn_reduced_indices = tuple(
        (ancestry - 156) // (P**2) for ancestry in qb_to_qa_q11
    )
    horn_spectral_strata = (
        ("13-unit", DEPTH - DEPTH // P, P**3),
        ("valuation-1", DEPTH // P - DEPTH // (P**2), P**2),
        ("valuation-2", DEPTH // (P**2) - DEPTH // (P**3), P),
        ("valuation-at-least-3", DEPTH // (P**3), 1),
    )
    require(
        min(horn_reduced_indices) == 350
        and max(horn_reduced_indices) == 1845
        and len(set(horn_reduced_indices)) == 449
        and 449 % P == 7
        and horn_spectral_strata
        == (
            ("13-unit", 342732, 2197),
            ("valuation-1", 26364, 169),
            ("valuation-2", 2028, 13),
            ("valuation-at-least-3", 169, 1),
        )
        and sum(row[1] for row in horn_spectral_strata) == DEPTH,
        "horn cyclotomic full-spectrum certificate changed",
    )

    # Prime-power orbit-mark unit lemma, specialized to this horn.
    # In F_p[C_(p^5)] put e=z-1.  Then e^(p^5)=0, so the augmentation
    # ideal is nilpotent and an element is a unit iff its augmentation is.
    # Here h=7+n with n in (e), hence
    #
    #   h^(-1)=7^(-1) sum_{j=0}^{p^5-1}(-7^(-1)n)^j.
    #
    # The circulant determinant is therefore a p-unit.  Over Q, translates
    # of h form a basis of the whole regular module, so every equivariant
    # linear realization of this depth-five horn observable has dimension
    # at least p^5.  This is unrelated to the 13/26 ranks of the q-response
    # quotient in THM-2835.  The inverse cannot be nonnegative: a positive
    # convolution inverse of a positive element with 449 support points
    # would create at least two positive support sums.
    horn_augmentation_mod_p = len(qb_to_qa_q11) % P
    horn_circulant_determinant_mod_p = pow(
        horn_augmentation_mod_p, DEPTH, P
    )
    require(
        horn_augmentation_mod_p == 7
        and pow(horn_augmentation_mod_p, -1, P) == 2
        and horn_circulant_determinant_mod_p == 7
        and len(qb_to_qa_q11) > 1,
        "horn group-algebra unit certificate changed",
    )

    def rotate_target_response(code, shift):
        source_bits = code & 7
        rows = tuple(
            (code >> (3 + 3 * q)) & 7 for q in range(P)
        )
        rotated = source_bits
        for q in range(P):
            rotated |= rows[(q - shift) % P] << (3 + 3 * q)
        return rotated

    horn_q_rotations = tuple(
        rotate_target_response(horn_state, shift)
        for shift in range(P)
    )
    horn_q_rotation_frequencies = tuple(
        response_state_counts.get(code, 0)
        for code in horn_q_rotations
    )
    require(
        horn_q_rotations
        == (
            628342390802,
            628692615306,
            631494411338,
            653908779594,
            833223725642,
            2267743294026,
            549760307810,
            35951370,
            287610946,
            2300887554,
            18407100418,
            147256803330,
            1178054426626,
        )
        and horn_q_rotation_frequencies == (449,) + (0,) * 12,
        "horn response acquired a physical cyclic q-orbit",
    )

    # Compare the word gauge with the endpoint/harmonic gauge.  Four forward
    # representative steps are the same residue displacement as inverse
    # normal label +9: target harmonic 11 becomes 7.  The word gauge alone
    # reindexes ancestry by 4*13^4 and preserves the horn.  The endpoint
    # gauge also moves the source carrier by -4*T/13; on that moved source
    # every transported horn label lies in the all-safe no-word chamber.
    gauge4_ell = tuple(
        4 * (speed % P) % P for speed in old.base.W
    )
    gauge4_shift = 4 * P**4
    gauge4_horn_labels = tuple(
        (ancestry - gauge4_shift) % DEPTH
        for ancestry in qb_to_qa_q11
    )
    gauge4_pieces = {
        name: tuple(old.base.build_set(pattern, gauge4_ell))
        for name, pattern in patterns.items()
    }
    gauge4_starts = {
        name: tuple(left for left, _right in gauge4_pieces[name])
        for name in WORD_ORDER
    }

    def gauge4_word_count(interval, offset, name, delta):
        flags = word_flags(
            interval,
            offset,
            period,
            gauge4_pieces[name],
            gauge4_starts[name],
        )
        return sum(
            flags[(ancestry + delta) % DEPTH]
            for ancestry in gauge4_horn_labels
        )

    gauge4_source_interval = (
        atom[0] - 4 * allocation_unit,
        atom[1] - 4 * allocation_unit,
    )
    gauge4_word_census = {
        "fixed_source": tuple(
            gauge4_word_count(atom, 0, name, 0)
            for name in WORD_ORDER
        ),
        "q11_target": tuple(
            gauge4_word_count(
                atom,
                allocation.physical.SHIFT + Q * allocation_unit,
                name,
                0,
            )
            for name in WORD_ORDER
        ),
        "q7_borrow_target": tuple(
            gauge4_word_count(
                atom,
                allocation.physical.SHIFT + 7 * allocation_unit,
                name,
                1,
            )
            for name in WORD_ORDER
        ),
        "moved_source": tuple(
            gauge4_word_count(
                gauge4_source_interval, 0, name, 0
            )
            for name in WORD_ORDER
        ),
    }
    gauge4_factor_pieces = {}
    for index, speed in enumerate(old.base.W):
        if index == 0:
            factor = comb_intervals(
                period,
                speed,
                91,
                -13 - 7 * gauge4_ell[index],
                13 - 7 * gauge4_ell[index],
            )
        else:
            factor = tuple(
                old.base.in_comb(index, gauge4_ell[index])
            )
        gauge4_factor_pieces[index] = (
            factor,
            tuple(left for left, _right in factor),
        )
    gauge4_moved_danger_vectors = Counter(
        tuple(
            int(
                whole_cylinder_hit(
                    gauge4_source_interval,
                    0,
                    ancestry,
                    period,
                    *gauge4_factor_pieces[index],
                )
            )
            for index in range(9)
        )
        for ancestry in gauge4_horn_labels
    )
    require(
        gauge4_ell == (4, 4, 4, 4, 4, 4, 0, 0, 0)
        and gauge4_word_census
        == {
            "fixed_source": (0, 449, 0),
            "q11_target": (449, 0, 0),
            "q7_borrow_target": (0, 0, 449),
            "moved_source": (0, 0, 0),
        }
        and gauge4_moved_danger_vectors
        == {(0, 0, 0, 0, 0, 0, 0, 0, 0): 449},
        "word/endpoint representative-gauge clutch boundary changed",
    )

    # The two reverse failures are not failures of the a/b switch.  In each
    # case slot 7 is safe and slot 8 is dangerous, as QB requests, but one
    # of the common low ordinary safeties turns dangerous: speed 14 for the
    # first label and speed 53 for the second.
    reverse_bad_patterns = {}
    for slot in (1, 4):
        pattern = dict(patterns["QB"])
        pattern[slot] = "in"
        reverse_bad_patterns[slot] = pattern
    reverse_bad_pieces = {
        slot: tuple(old.base.build_set(pattern, old.base.ZELL))
        for slot, pattern in reverse_bad_patterns.items()
    }
    reverse_bad_starts = {
        slot: tuple(left for left, _right in reverse_bad_pieces[slot])
        for slot in reverse_bad_pieces
    }
    reverse_bad_factor = {}
    q7_offset = allocation.physical.SHIFT + 7 * allocation_unit
    for ancestry in reverse_exceptions:
        hits = tuple(
            slot
            for slot in (1, 4)
            if whole_cylinder_hit(
                atom,
                q7_offset,
                ancestry - 1,
                period,
                reverse_bad_pieces[slot],
                reverse_bad_starts[slot],
            )
        )
        require(len(hits) == 1, "reverse exception lacks a unique low bit")
        reverse_bad_factor[ancestry] = (hits[0], old.base.W[hits[0]])
    require(
        reverse_bad_factor
        == {107978: (1, 14), 154622: (4, 53)},
        "reverse exception low-factor mechanism changed",
    )

    first = qb_to_qa_q11[0]
    source_window = containing_window(
        atom, 0, first, period, pieces["QB"], starts["QB"]
    )
    target_window = containing_window(
        atom,
        allocation.physical.SHIFT + Q * allocation_unit,
        first,
        period,
        pieces["QA"],
        starts["QA"],
    )
    common_window = (
        max(source_window[0], target_window[0]),
        min(source_window[1], target_window[1]),
    )
    margins = (
        atom[0] - common_window[0],
        common_window[1] - atom[1],
    )
    require(
        common_window == (138281416853580, 159555480984900)
        and margins == (3723575735880, 17550461950560),
        "first q11 QB-to-QA whole-cylinder margins changed",
    )

    # Boolean unions are the closest lawful coarsenings.  They forget the
    # exact semantic word; the census quantifies precisely what that costs.
    union_rows = {}
    for size in range(1, len(WORD_ORDER) + 1):
        for subset in combinations(WORD_ORDER, size):
            source_union = union_flags(
                tuple(source[name] for name in subset)
            )
            target_row = []
            common_row = []
            for q in range(P):
                target_union = union_flags(
                    tuple(target[name, q] for name in subset)
                )
                target_row.append(sum(target_union))
                common_row.append(
                    overlap_count(source_union, target_union)
                )
            union_rows[subset] = (
                sum(source_union),
                tuple(target_row),
                tuple(common_row),
            )
    require(
        union_rows[("QA", "QB")][2][Q] == 449
        and union_rows[("QA", "QAB")][2][Q] == 10328
        and union_rows[("QB", "QAB")][2][Q] == 0
        and union_rows[("QA", "QB", "QAB")][2][Q] == 10777,
        "q11 union-coarsening census changed",
    )

    print("Q11 OUTER-WORD ANCESTRY SCAN")
    print("status=PROVED VERIFIED-EXACT; no row decrement or LRC14")
    print(
        f"physical_cell=(s,t,clock)={witness}; "
        f"I={atom}; J11={target11}; all_six_factors=1"
    )
    print(
        "q7_physical_stratification="
        "(factor_order=(E3,clock,q1src,q2src,q1target,q2target),"
        f"q11_valid_cells={len(q11_cells)},"
        f"q7_failure_census={tuple(sorted(q7_failure_census.items()))},"
        "q7_fully_physical_cells=0,"
        "clock_q7_support=(0,1,0,0,0,0,0),"
        "best_horn_cell_products=35*449=15715_all_killed_by_E3)"
    )
    print(
        "delayed_semantic="
        f"(source12={source_semantic[12]},target6={target_semantic[6]}); "
        f"right_endpoint_rows={endpoint_rows}"
    )
    print(
        "word_signatures_(a,b)="
        f"{tuple((name, patterns[name][7], patterns[name][8]) for name in WORD_ORDER)}; "
        "canonical_THM2584_current=QB"
    )
    for name in WORD_ORDER:
        print(
            f"{name}: source={source_counts[name]}; "
            f"target_by_q={target_counts[name]}; "
            f"same_word_common_by_q={common_counts[name]}"
        )
    for source_name in WORD_ORDER:
        for target_name in WORD_ORDER:
            print(
                f"directed_{source_name}_to_{target_name}="
                f"{directed[source_name, target_name]}"
            )
    print(
        "QB_to_QA_q11="
        f"(count={len(qb_to_qa_q11)},first={first},"
        f"last={qb_to_qa_q11[-1]},residue_mod169=156,"
        "source_QB_residue_arc=12..156,"
        "target_QA_residue_arc=156..168,0..10,"
        "residue156_counts=455/449,"
        f"digest={label_digest.hexdigest()},gap_census={gap_census},"
        f"whole_atom_window={common_window},margins={margins})"
    )
    print(
        "QB_to_QA_natural_q7_borrow="
        f"(word_profiles_delta_-2..2={borrow_word_profiles},"
        "plus1_lands_QAB=449,plus1_lands_QB=0,"
        "reverse_minus1_lands_QB=447,"
        f"reverse_exceptions={reverse_exceptions},"
        f"exception_low_factors={reverse_bad_factor})"
    )
    print(
        "complete_response_automaton="
        f"(state_count={len(response_state_counts)},"
        f"state_frequencies={tuple(sorted(response_state_counts.items()))},"
        "adjacent_edges=44,deterministic_states=15,max_outdegree=12,"
        f"code_SCC_sizes={response_scc_sizes},"
        "natural_h9_nodes=260,natural_h9_edges=476,"
        f"natural_h9_SCC_sizes={natural_scc_sizes})"
    )
    print(
        "horn_response_state="
        f"(source_state={horn_state},frequency=449,"
        f"plus1_state={horn_successor_state},"
        "plus1_state_global_frequency=10328,"
        "horn_edge_multiplicity=449,"
        f"plus1_state_predecessor_multiplicities="
        f"{tuple(sorted(horn_successor_predecessors.items()))},"
        "single_q7_word_is_carry_blind=1,"
        "full_response_code_detects_carry=1)"
    )
    print(
        "horn_provenance_lift="
        "(one_predecessor_bit_suffices_set_theoretically=1,"
        "two_digit_target_residue=157_selects_exactly_449,"
        "low_residue1_selects=898,"
        f"C_residue_mod169={tuple(sorted(successor_residue_census.items()))},"
        "pair_states=44,pair_edges=68,pair_deterministic=39,"
        "pair_max_outdegree=12,HC_pair_frequency=449,"
        f"maximal_proper_period_shift={maximal_proper_period},"
        f"period_matches={proper_period_matches},"
        f"first_period_failure={first_period_failure}:0->33554450,"
        "fundamental_response_period=371293)"
    )
    print(
        "cyclic_observable_periods="
        f"(candidates={period_candidates},"
        f"response_shift_mismatch_first={response_period_table},"
        f"horn_shift_mismatch_first={horn_period_table},"
        "marked_HC_edge_equals_horn="
        f"{int(marked_edge_indicator == horn_indicator)},"
        "least_periods=(371293,371293,371293),"
        "deterministic_invertible_observer_lower_bound=371293)"
    )
    print(
        "horn_cyclotomic_full_spectrum="
        "(mask=Z^156*G(Z^169),"
        "reduced_index_range=350..1845,G(1)=449=7_mod13,"
        "effective_orders=(2197,169,13,1),"
        f"spectral_strata={horn_spectral_strata},"
        "nonzero_characters=371293,nonzero_13_units=342732,"
        "mechanism=Phi_(13^j)(1)=13_cannot_divide_G(1))"
    )
    print(
        "horn_group_algebra_unit="
        "(mod13_ring=F13[e]/(e^371293),augmentation=449=7_mod13,"
        "augmentation_inverse=2,nilpotence_index=371293,"
        f"circulant_determinant_mod13={horn_circulant_determinant_mod_p},"
        "rational_translate_span=371293,"
        "linear_equivariant_observer_lower_bound=371293,"
        "unique_rational_signed_convolution_inverse=1,"
        "nonnegative_inverse=0_for_support449,"
        "THM2835_q_response_ranks_13_26_are_a_different_quotient=1)"
    )
    print(
        "horn_cyclic_q_orbit="
        f"(rotated_codes={horn_q_rotations},"
        f"actual_frequencies={horn_q_rotation_frequencies},"
        "verdict=only_base_profile_is_physical)"
    )
    print(
        "marked_gauge_vertical_clutch="
        f"(ell4={gauge4_ell},ancestry_reindex=-{gauge4_shift},"
        f"word_census_QA_QB_QAB={gauge4_word_census},"
        "q7_endpoint_gauge_rows="
        f"{gauge_endpoint_rows},"
        "moved_source_danger_vector="
        f"{tuple(gauge4_moved_danger_vectors.items())},"
        "verdict=word_gauge_retains_H_to_C_but_endpoint_gauge_moves_"
        "source_to_all_safe_no_word)"
    )
    for subset, (source_count, target_row, common_row) in union_rows.items():
        print(
            f"union_{'+'.join(subset)}="
            f"(source={source_count},target_by_q={target_row},"
            f"common_by_q={common_row})"
        )
    print(
        "CONCLUSION: no canonical fixed word has a nonzero q11 common "
        "ancestry label.  QA alone has 10,787 q11 target labels but zero "
        "source labels; QB and QAB have no q11 target labels.  The nearest "
        "positive repair is the directed QB->QA cospan: exactly 449 "
        "whole-cylinder labels, all in residue 156 modulo 169.  It "
        "preserves the physical q11 cell, "
        "owner row, delayed values, and endpoint scalars, but changes the "
        "outer semantic word from {b} to {a}; hence it is not the same "
        "THM-2584/2806 semantic current.  Union coarsening makes the q11 "
        "intersection positive (449 for QA+QB; 10,777 for all three) "
        "only by forgetting that word label.  A closure therefore needs "
        "a lawful word-transition arrow/natural transformation, not a "
        "different static outer mask.  Conditional on the missing q11,h9 "
        "action, its compulsory +1 ancestry borrow sends all 449 labels "
        "to the q7 QAB={a,b} block.  The reverse -1 shift returns 447 "
        "labels to QB; its two failures retain the requested a/b bits but "
        "lose respectively the speed-14 and speed-53 common safety.  This "
        "types an exact label-level semantic horn QB->QA->QAB.  The full "
        "response automaton isolates its source as one 449-label state "
        "and its +1 image as one deterministic edge, but the image state "
        "has 10 predecessor states and the global quotient is one SCC.  "
        "On the 35 best physical cells E3 still vanishes identically at "
        "q7.  The marked gauge can preserve either the word horn or the "
        "q7 endpoint copy: applying the endpoint transport moves the "
        "source into the all-safe no-word chamber.  Hence none of these "
        "finite selectors supplies the physical current-level horn filler."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
