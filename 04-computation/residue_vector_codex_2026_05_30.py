#!/usr/bin/env python3
"""
codex-2026-05-30

Compact residue-vector table for the projection/transport synthesis.

This is intentionally small: it reuses the exact helpers from the existing
projection-defect bridge and emits one row per representative example.  The
goal is not a general feature extractor yet, but a reproducible seed table for
HYP-1779/HYP-1781.
"""

from __future__ import annotations

from projection_defect_bridge_s12 import (
    circulant_rows,
    cycle_projection_even_graph,
    deletion_projection_profile,
    gentourng_reps,
    hamiltonian_paths,
    odd_cycle_support_stats,
    rows_from_nauty_bits,
    score_sequence,
    thm025_rows,
    transitive_rows,
)


def residue_row(label: str, rows: tuple[int, ...]) -> dict[str, object]:
    stats = odd_cycle_support_stats(rows)
    deletion = deletion_projection_profile(rows)
    even = cycle_projection_even_graph(rows)
    max_loss = max((float(r["loss_frac"]) for r in deletion), default=0.0)
    kill_vertices = tuple(int(r["v"]) for r in deletion if int(r["kept"]) == 0)
    top = sorted(deletion, key=lambda r: (-float(r["loss_frac"]), int(r["v"])))[:3]
    top_signature = tuple((int(r["v"]), int(r["lost"]), int(r["kept"])) for r in top)

    return {
        "label": label,
        "n": len(rows),
        "H": hamiltonian_paths(rows),
        "scores": score_sequence(rows),
        "cycles": len(stats["masks"]),
        "supports": int(stats["support_count"]),
        "support_excess": int(stats["support_excess"]),
        "max_support_mult": int(stats["max_support_mult"]),
        "alpha": tuple(int(x) for x in stats["alpha"]),
        "complete_Omega": bool(stats["complete"]),
        "max_loss_frac": round(max_loss, 6),
        "kill_vertices": kill_vertices,
        "top_deletion_v_lost_kept": top_signature,
        "even_edges": None if even is None else int(even["edges"]),
        "even_degree_seq": None if even is None else tuple(int(x) for x in even["degrees"]),
    }


def main() -> None:
    n8_reps = gentourng_reps(8)
    examples = [
        ("transitive n=7", transitive_rows(7)),
        ("Paley T7", circulant_rows(7, {1, 2, 4})),
        ("interval T7", circulant_rows(7, {1, 2, 3})),
        ("H=63 cid=2519", rows_from_nauty_bits(8, n8_reps[2519])),
        ("H=63 cid=3285", rows_from_nauty_bits(8, n8_reps[3285])),
        ("THM-025 n=9", thm025_rows()),
    ]

    rows = [residue_row(label, tournament) for label, tournament in examples]
    keys = [
        "label",
        "n",
        "H",
        "scores",
        "cycles",
        "supports",
        "support_excess",
        "max_support_mult",
        "alpha",
        "complete_Omega",
        "max_loss_frac",
        "kill_vertices",
        "top_deletion_v_lost_kept",
        "even_edges",
        "even_degree_seq",
    ]

    print("Residue-vector seed table (codex-2026-05-30)")
    print("=" * 72)
    print(" | ".join(keys))
    print("-" * 72)
    for row in rows:
        print(" | ".join(str(row[k]) for k in keys))

    print()
    print("Reading guide:")
    print("  max_loss_frac=1 means deletion of one vertex kills all odd cycles.")
    print("  complete_Omega=True with cycles=r gives H=1+2r by OCF.")
    print("  even_edges=None means the cycle-space projection is ambiguous at even n.")


if __name__ == "__main__":
    main()
