#!/usr/bin/env python3
"""Exact all-packet three-twist census on THM-2334's typed control row.

This is deliberately a positive control, not a covering-row computation.  It
imports the independently refereed interval engine of THM-2334 and combines
it with THM-3666's physical pair-swap coordinates.  The 120 owner packets
collapse to 30 ordered graft pairs; every resulting support-minimal
three-twist defect is certified nonzero in two exact cyclotomic embeddings.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENTS = (
    ROOT / "04-computation/lrc14_169_twist_two_twist_referee_thm2334.py",
    ROOT / "05-knowledge/results/lrc14_169_twist_two_twist_referee_thm2334.out",
    ROOT / "04-computation/lrc_owner_pivot_dual_pair_swap_twists_thm3666.py",
    ROOT / "05-knowledge/results/lrc_owner_pivot_dual_pair_swap_twists_thm3666.out",
)
EXPECTED_PARENT_HASHES = (
    "0e4a9e181263647e13d2a6738b6996c45df901d9d2b37d4d589dfddfbdd91480",
    "fabae42d085cd9c338a142cbd0898116fde48f27e9c14aa1c3bcc3abfdb99130",
    "8871c9d8a6eff3d1be33d7944197cffbb8a11c60f6d61132be27a89ce22ff96a",
    "1015cff6d84f2571cabaaf3e66fad92bee2cba4ba9d616769643dc4809e5c9b9",
)

P = 13
UNITS = tuple(range(6))
EXPECTED_SEMANTIC_SHA256 = "a33d438898b3c714557bc652add55c82e7d648ded6e9416a5ddb0a674282e27f"
EXPECTED_COMPACT_PACKET_SHA256 = (
    "9521c4654c9fddb43094e1194b36ef75f6118b8b216918c3f80787270c6db55a"
)


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def digest_json(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(encoded).hexdigest()


def valuation_13(value: int) -> int:
    require(value != 0, "valuation requested at zero")
    value = abs(value)
    exponent = 0
    while value % P == 0:
        value //= P
        exponent += 1
    return exponent


def load_referee():
    path = PARENTS[0]
    spec = importlib.util.spec_from_file_location("thm2334_referee", path)
    require(spec is not None and spec.loader is not None, "cannot load referee")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dipole(target: int, graft: int, coefficient: int) -> tuple[int, ...]:
    vector = [0] * 9
    vector[target] = coefficient % P
    vector[graft] = -coefficient % P
    return tuple(vector)


def marked_intervals(ref, present, word):
    """Materialize E intersect {t:R*t mod 1 in Q} on the NN grid."""
    scaled_present = [(ref.R * left, ref.R * right) for left, right in present]
    present_index = 0
    output = []
    for q_left, q_right in ref.word_preimages(word):
        while (present_index < len(scaled_present)
               and scaled_present[present_index][1] <= q_left):
            present_index += 1
        scan = present_index
        while scan < len(scaled_present) and scaled_present[scan][0] < q_right:
            e_left, e_right = scaled_present[scan]
            left = max(e_left, q_left)
            right = min(e_right, q_right)
            if left < right:
                if output and output[-1][1] == left:
                    output[-1] = (output[-1][0], right)
                else:
                    output.append((left, right))
            scan += 1
    return output


def jump_map(intervals):
    jumps = {}
    for left, right in intervals:
        jumps[left] = jumps.get(left, 0) + 1
        jumps[right] = jumps.get(right, 0) - 1
    return {endpoint: jump for endpoint, jump in jumps.items() if jump}


def main() -> None:
    parent_hashes = tuple(file_sha256(path) for path in PARENTS)
    require(parent_hashes == EXPECTED_PARENT_HASHES,
            ("parent hashes", parent_hashes, EXPECTED_PARENT_HASHES))
    ref = load_referee()
    require(ref.W == (1, 14, 27, 40, 53, 66, 13, 2197, 742586), ref.W)
    require(ref.R == 169 and ref.X == 13 and ref.M == 1 and ref.Y == 742599,
            (ref.R, ref.X, ref.M, ref.Y))
    require(tuple((prime, root) for prime, root in ref.EMBEDDINGS) == (
        (352341050142921841, 435817657216),
        (956354278959359281, 153943385426666320),
    ), ref.EMBEDDINGS)

    zero = (0,) * 9
    q_a = ref.build_boolean_set(ref.PATTERN_QA, zero)

    # Only thirteen present sets are needed: the base and the six negative
    # alpha / six negative beta dipoles.  Omitted-unit choices do not alter
    # either dual character.
    shifts = {"zero": zero}
    for graft in UNITS:
        shifts[f"a{graft}"] = dipole(ref.TARGET_A, graft, -1)
        shifts[f"b{graft}"] = dipole(ref.TARGET_B, graft, -1)

    values = {}
    selected_present = {}
    for label, shift in shifts.items():
        present = ref.build_boolean_set(ref.PATTERN_E, shift)
        if label in ("zero", "a1", "b2"):
            selected_present[label] = present
        ax_x, overlap, components = ref.marked_endpoint_sum(present, q_a, ref.X)
        by_y = ref.endpoint_sum_on_t_den(present, -ref.Y)
        ax_one, _, _ = ref.marked_endpoint_sum(present, q_a, 1)
        gammas = []
        for index, (prime, root) in enumerate(ref.EMBEDDINGS):
            phase = pow(root,
                        (ref.M * shift[ref.TARGET_B] * (ref.NN // P)) % ref.NN,
                        prime)
            gammas.append(phase * ax_x[index] % prime * by_y[index] % prime)
        require(overlap > 0 and components > 0, ("empty marked product", label))
        require(all(value != 0 for value in ax_x + by_y),
                ("zero endpoint factor", label, ax_x, by_y))
        values[label] = {
            "shift": shift,
            "gamma": tuple(gammas),
            "overlap": overlap,
            "components": components,
            "marked_jet_1": ax_one,
        }

    chart_rows = []
    packet_rows = []
    base = values["zero"]
    for graft_a in UNITS:
        for graft_b in UNITS:
            if graft_a == graft_b:
                continue
            va = values[f"a{graft_a}"]
            vb = values[f"b{graft_b}"]
            residual = tuple(
                (base["gamma"][j] + va["gamma"][j] - 2 * vb["gamma"][j]) % prime
                for j, (prime, _root) in enumerate(ref.EMBEDDINGS)
            )
            require(all(residual), ("vanishing three-twist defect", graft_a, graft_b))
            mass_residual = base["overlap"] + va["overlap"] - 2 * vb["overlap"]
            jet_residual = tuple(
                (base["marked_jet_1"][j] + va["marked_jet_1"][j]
                 - 2 * vb["marked_jet_1"][j]) % prime
                for j, (prime, _root) in enumerate(ref.EMBEDDINGS)
            )
            chart = (
                graft_a,
                graft_b,
                residual,
                mass_residual,
                jet_residual,
                (base["overlap"], va["overlap"], vb["overlap"]),
                (base["components"], va["components"], vb["components"]),
            )
            chart_rows.append(chart)
            for omitted in UNITS:
                if omitted not in (graft_a, graft_b):
                    packet_rows.append((omitted,) + chart)

    require(len(chart_rows) == 30, len(chart_rows))
    require(len(packet_rows) == 120, len(packet_rows))
    require(len({(row[1], row[2]) for row in packet_rows}) == 30,
            "packet-to-chart collapse failed")

    first = chart_rows[0]
    require(first[:2] == (0, 1), first[:2])
    # Also report the historical q1/q2 chart used by THM-2334 and the audit.
    canonical = next(row for row in chart_rows if row[:2] == (1, 2))
    c_a = values["a1"]
    c_b = values["b2"]
    require((base["overlap"], c_a["overlap"], c_b["overlap"]) ==
            (60084076348296, 54135630512964, 61887542465528),
            "canonical overlap controls changed")
    require((base["components"], c_a["components"], c_b["components"]) ==
            (188056, 169431, 186674), "canonical component controls changed")
    require(canonical[2] ==
            (200189390529282589, 368827291397979410),
            ("canonical residual", canonical[2]))
    require(canonical[3] == -9555378069796,
            ("canonical mass residual", canonical[3]))
    require(canonical[4] ==
            (212565399985231344, 197131116389696829),
            ("canonical first jet", canonical[4]))

    mass_valuations = tuple(
        (valuation_13(base["overlap"]), valuation_13(values[f"a{a}"]["overlap"]),
         valuation_13(values[f"b{b}"]["overlap"]), valuation_13(row[3]))
        for row in chart_rows
        for a, b in (row[:2],)
    )

    # Reproduce the independent audit's compact byte contract exactly.
    compact_charts = tuple(
        (row[0], row[1], row[2][0], row[2][1], row[3],
         valuation_13(base["overlap"]),
         valuation_13(values[f"a{row[0]}"]["overlap"]),
         valuation_13(values[f"b{row[1]}"]["overlap"]),
         valuation_13(row[3]))
        for row in chart_rows
    )
    compact_packets = tuple(
        (omitted,) + next(
            row for row in compact_charts if row[0] == graft_a and row[1] == graft_b
        )
        for omitted in UNITS
        for graft_a in UNITS if graft_a != omitted
        for graft_b in UNITS if graft_b not in (omitted, graft_a)
    )
    compact_bytes = json.dumps(compact_packets, separators=(",", ":")).encode("ascii")
    compact_digest = sha256(compact_bytes).hexdigest()
    require(len(compact_packets) == 120 and len(compact_bytes) == 8269,
            (len(compact_packets), len(compact_bytes)))
    require(compact_digest == EXPECTED_COMPACT_PACKET_SHA256,
            (compact_digest, EXPECTED_COMPACT_PACKET_SHA256))
    require(compact_packets[0] ==
            (0, 1, 2, 200189390529282589, 368827291397979410,
             -9555378069796, 2, 1, 1, 1), compact_packets[0])
    require(compact_packets[-1] ==
            (5, 4, 3, 87332794336807381, 384335766739000838,
             -10444758045280, 2, 1, 1, 1), compact_packets[-1])

    # A local discontinuity witness: choose the first endpoint, in term order,
    # owned by exactly one of the three canonical marked products.
    interval_bank = {
        label: marked_intervals(ref, selected_present[label], q_a)
        for label in ("zero", "a1", "b2")
    }
    require(tuple(len(interval_bank[label]) for label in ("zero", "a1", "b2")) ==
            (188056, 169431, 186674), "marked component census")
    jump_bank = {label: jump_map(interval_bank[label]) for label in interval_bank}
    exclusive = None
    for label in ("zero", "a1", "b2"):
        for endpoint in sorted(jump_bank[label]):
            owners = tuple(name for name in ("zero", "a1", "b2")
                           if endpoint in jump_bank[name])
            if len(owners) == 1:
                weights = {"zero": 1, "a1": 1, "b2": -2}
                exclusive = (endpoint, label, jump_bank[label][endpoint],
                             weights[label] * jump_bank[label][endpoint])
                break
        if exclusive is not None:
            break
    require(exclusive == (7791347595131100, "zero", 1, 1), exclusive)
    exclusive_time = Fraction(exclusive[0], ref.NN)
    require(exclusive_time == Fraction(1609245, 10396204), exclusive_time)

    semantic_payload = (
        parent_hashes,
        tuple((label, values[label]) for label in sorted(values)),
        tuple(chart_rows),
        tuple(packet_rows),
        compact_packets,
        mass_valuations,
        exclusive,
    )
    semantic = digest_json(semantic_payload)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3669 typed-control all-packet three-twist census ==")
    print(f"parents_sha256={parent_hashes}")
    print(f"row={ref.W};word={{a}};R={ref.R};triangle=({ref.X},{ref.M},{ref.Y})")
    print(f"cyclotomic_order={ref.NN};embeddings={ref.EMBEDDINGS}")
    print("unique_present_sets=13;ordered_graft_charts=30;owner_packets=120")
    print("all_30_chart_defects_nonzero_in_both_embeddings=True")
    print("all_120_packet_defects_nonzero_in_both_embeddings=True")
    print(f"chart_ledger_sha256={digest_json(tuple(chart_rows))}")
    print(f"packet_ledger_sha256={digest_json(tuple(packet_rows))}")
    print(f"independent_compact_packet_sha256={compact_digest};bytes={len(compact_bytes)}")
    print("canonical_packet=(omitted=0,graft_a=1,graft_b=2)")
    print(f"canonical_gamma=(H0={base['gamma']},Ha={c_a['gamma']},Hb={c_b['gamma']})")
    print(f"canonical_overlap_numerators={(base['overlap'], c_a['overlap'], c_b['overlap'])}")
    print(f"canonical_components={(base['components'], c_a['components'], c_b['components'])}")
    print(f"canonical_three_twist_residual={canonical[2]}")
    print(f"canonical_mass_residual={canonical[3]}/{ref.NN}={Fraction(canonical[3], ref.NN)}")
    print(f"canonical_overlap_v13={(valuation_13(base['overlap']), valuation_13(c_a['overlap']), valuation_13(c_b['overlap']))};residual_v13={valuation_13(canonical[3])}")
    print(f"canonical_marked_product_first_jet={canonical[4]}")
    print(f"canonical_first_exclusive_boundary={exclusive};time={exclusive_time}")
    print(f"mass_valuation_patterns={tuple(sorted(set(mass_valuations)))}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256={file_sha256(source)}")
    print("scope=typed non-cover positive control;no covering-row transfer or LRC(14) closure")
    print("PASS")


if __name__ == "__main__":
    main()
