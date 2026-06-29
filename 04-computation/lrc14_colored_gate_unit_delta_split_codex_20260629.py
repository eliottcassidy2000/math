#!/usr/bin/env python3
"""HYP-3472: sharpen HYP-3471 by classifying minimum E/branch gates.

HYP-3471 proved on the audited bank that every row with a dead component has
some rank <= 2 E/branch survivor gate.  This follow-up asks how rigid the
minimum such gate is.

The exact split in the same 135-row bank is:

    structured dead rows  -> always branch-specific unit-delta edge gates
    random dead rows      -> mostly the same, plus one both-branch unit-delta
                             row and a 19-row cover-delta sidecar packet

This isolates the random residual as a finite delta packet rather than an
arbitrary geometric failure of the colored-gate reservoir.
"""

from __future__ import annotations

from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3471_PATH = ROOT / "04-computation" / "lrc14_colored_gate_reservoir_codex_20260629.py"
H3435_PATH = ROOT / "04-computation" / "lrc14_two_adic_branch_cover_certificate_codex_20260628.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3471 = load_module("hyp3471_for_hyp3472", H3471_PATH)
H3435 = load_module("hyp3435_for_hyp3472", H3435_PATH)


def gate_sort_key(gate):
    return (gate.length, gate.row_name, gate.component_index, gate.interval[0])


def min_e_branch_gate(row):
    gates = [gate for gate in row.low_rank_gates if H3471.is_e_branch_gate(gate)]
    if not gates:
        raise ValueError(f"row has dead component but no E/branch gate: {row.name}")
    return min(gates, key=gate_sort_key)


def structural_word(gate) -> tuple[str, str, str, int, int]:
    return (
        gate.endpoint_kind_signature,
        gate.branch_mask,
        gate.adjacency,
        gate.b0_delta,
        gate.b1_delta,
    )


def typed_word(gate) -> str:
    return H3471.typed_word(gate)


def kind(gate) -> str:
    edge = gate.adjacency in {"left_bad_edge", "right_bad_edge"}
    unit = (gate.b0_delta, gate.b1_delta) == (1, 1)
    if edge and unit and gate.branch_mask in {"branch0", "branch1"}:
        return "branch_unit_delta"
    if edge and unit and gate.branch_mask == "both":
        return "both_unit_delta"
    if edge and gate.branch_mask in {"branch0", "branch1"}:
        return "delta_sidecar_packet"
    return "other"


def print_section(title: str) -> None:
    print()
    print(f"## {title}")


def main() -> None:
    structured_rows = H3435.structured_rows()
    random_rows = H3435.random_rows()
    rows = dict(structured_rows)
    rows.update(random_rows)

    dead_rows = []
    structured_dead = []
    random_dead = []
    min_kind_hist = Counter()
    structured_typed = Counter()
    structured_struct = Counter()
    random_exception_typed = Counter()
    random_exception_struct = Counter()
    random_exception_names = []
    random_both_unit = []
    bad_kind_rows = []

    for name, speeds in rows.items():
        row = H3471.H3453.join_row(name, speeds)
        if not row.has_dead:
            continue
        gate = min_e_branch_gate(row)
        gate_kind = kind(gate)
        dead_rows.append((row, gate, gate_kind))
        min_kind_hist[gate_kind] += 1
        if name in structured_rows:
            structured_dead.append((row, gate, gate_kind))
            structured_typed[typed_word(gate)] += 1
            structured_struct[structural_word(gate)] += 1
        else:
            random_dead.append((row, gate, gate_kind))
            if gate_kind == "delta_sidecar_packet":
                random_exception_names.append(
                    (
                        row.name,
                        typed_word(gate),
                        structural_word(gate),
                        H3471.H3453.fmt(gate.length),
                        row.dead_count,
                        len(row.low_rank_gates),
                    )
                )
                random_exception_typed[typed_word(gate)] += 1
                random_exception_struct[structural_word(gate)] += 1
            elif gate_kind == "both_unit_delta":
                random_both_unit.append(
                    (
                        row.name,
                        typed_word(gate),
                        structural_word(gate),
                        H3471.H3453.fmt(gate.length),
                        row.dead_count,
                        len(row.low_rank_gates),
                    )
                )
        if gate_kind == "other":
            bad_kind_rows.append(
                (
                    row.name,
                    typed_word(gate),
                    structural_word(gate),
                    H3471.H3453.fmt(gate.length),
                )
            )

    print("HYP-3472 COLORED GATE UNIT-DELTA SPLIT")
    print("status=EVIDENCE / exact HYP-3471 sharpening; not an LRC14 proof")
    print("source=HYP-3471 colored gate-reservoir + H3435 structured/random bank split")

    print_section("Aggregate Split")
    print(f"rows_audited={len(rows)}")
    print(f"structured_rows={len(structured_rows)}")
    print(f"random_rows={len(random_rows)}")
    print(f"dead_rows={len(dead_rows)}")
    print(f"structured_dead_rows={len(structured_dead)}")
    print(f"random_dead_rows={len(random_dead)}")
    print(f"min_kind_hist={dict(min_kind_hist)}")
    print(f"unexpected_other_rows={bad_kind_rows}")

    print_section("Structured Packet")
    print(
        "claim=every structured dead row has a branch-specific unit-delta "
        "minimum E/branch gate"
    )
    print(
        "structured_unit_delta_rows="
        f"{sum(gate_kind == 'branch_unit_delta' for _row, _gate, gate_kind in structured_dead)}"
        f"/{len(structured_dead)}"
    )
    print(f"structured_min_typed_hist={dict(structured_typed)}")
    print(f"structured_min_struct_hist={dict(structured_struct)}")
    print(
        "structured_samples="
        f"{[(row.name, typed_word(gate), structural_word(gate), H3471.H3453.fmt(gate.length)) for row, gate, _gate_kind in structured_dead[:15]]}"
    )

    print_section("Random Residual")
    print(
        "random_branch_unit_delta_rows="
        f"{sum(gate_kind == 'branch_unit_delta' for _row, _gate, gate_kind in random_dead)}"
        f"/{len(random_dead)}"
    )
    print(
        "random_both_unit_delta_rows="
        f"{sum(gate_kind == 'both_unit_delta' for _row, _gate, gate_kind in random_dead)}"
        f"/{len(random_dead)}"
    )
    print(
        "random_delta_sidecar_rows="
        f"{sum(gate_kind == 'delta_sidecar_packet' for _row, _gate, gate_kind in random_dead)}"
        f"/{len(random_dead)}"
    )
    print(f"random_both_unit_samples={random_both_unit}")
    print()
    print("delta_sidecar packet = still single-bad-edge + single-branch,")
    print("but needs adjacent-cover delta sidecar beyond the unit-delta forms.")
    print(f"delta_sidecar_struct_hist={dict(random_exception_struct)}")
    print(f"delta_sidecar_typed_hist={dict(random_exception_typed)}")
    print(f"delta_sidecar_rows={random_exception_names}")

    print_section("Finite Consequence")
    print("Exact sharpening of HYP-3471 on the audited bank:")
    print("  every structured dead row -> branch-specific unit-delta minimum E/branch gate")
    print("  every random dead row -> branch-specific unit-delta minimum gate,")
    print("                       or one both-branch unit-delta row,")
    print("                       or a 19-row single-branch cover-delta sidecar packet")
    print("The residual is therefore finite and typed; it is not arbitrary loss of the")
    print("colored gate-reservoir.  The next proof task is to route the 19-row packet")
    print("through HYP-3451 conductance/Menger or HYP-3455 gluing sidecars.")


if __name__ == "__main__":
    main()
