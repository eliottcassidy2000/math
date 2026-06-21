#!/usr/bin/env python3
"""
codex-2026-06-21-S75

Classify sorted-W-cell leaks for the HYP-2780/HYP-2781 joint-coupling route.

This is a focused follow-up to lrc14_cell_width_majorization_codex_s75.py.
It reuses that exact Fraction width engine and asks whether sorted-prefix
violations in modest extended bounded banks are finite one-hole AP-extension
events, and at which sorted prefix the leak is repaid.

Tournament Analysis:
  Vertices are leak explanations, not runners or cells.
  Pair observable is

      (explains_k9, explains_k10, survives_wider_span, finite_family,
       proof_ready)

  The switch is lexicographic comparison.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


sys.stdout.reconfigure(line_buffering=True)


ROOT = Path(__file__).resolve().parents[1]
ENGINE_PATH = ROOT / "04-computation" / "lrc14_cell_width_majorization_codex_s75.py"


def load_engine():
    spec = spec_from_file_location("cell_width_engine_s75", ENGINE_PATH)
    mod = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


engine = load_engine()


@lru_cache(None)
def sorted_prefixes(E):
    return engine.sorted_top_sums(engine.Wvec(tuple(E)))


def total_width(E):
    return sum(engine.Wvec(tuple(E)))


def ap_hole_extra_signature(E):
    Eset = set(E)
    k = len(E)
    ap = set(range(k))
    holes = tuple(sorted(ap - Eset))
    extras = tuple(sorted(x for x in E if x >= k))
    below = tuple(sorted(x for x in E if x < k))
    gaps = tuple(b - a for a, b in zip(E, E[1:]))
    return holes, extras, below, gaps


def leak_records(k, span):
    consec = tuple(range(k))
    base_prefix = sorted_prefixes(consec)
    base_total = total_width(consec)
    records = []
    for E in engine.full_residue_bank(k, span):
        pref = sorted_prefixes(E)
        gains = tuple(pref[i] - base_prefix[i] for i in range(5))
        leak_levels = tuple(i + 1 for i, g in enumerate(gains) if g > 0)
        if not leak_levels:
            continue
        repay = None
        for i, g in enumerate(gains, start=1):
            if i >= min(leak_levels) and g <= 0:
                repay = i
                break
        if repay is None:
            dt = total_width(E) - base_total
            repay = 6 if dt <= 0 else None
        holes, extras, below, gaps = ap_hole_extra_signature(E)
        records.append(
            {
                "E": E,
                "levels": leak_levels,
                "gains": gains,
                "dtotal": total_width(E) - base_total,
                "repay": repay,
                "holes": holes,
                "extras": extras,
                "below": below,
                "gaps": gaps,
            }
        )
    return records


def family_label(rec):
    holes = rec["holes"]
    extras = rec["extras"]
    k = len(rec["E"])
    if len(holes) == 1 and len(extras) == 1 and extras[0] == k:
        return f"one_hole_extend_h{holes[0]}"
    if len(holes) == 0 and len(extras) == 0:
        return "consec"
    if len(holes) == 1 and len(extras) == 1:
        return "one_hole_one_extra"
    if len(extras) == 0:
        return "inside_ap_window"
    return f"holes{len(holes)}_extras{len(extras)}"


LENSES = {
    "one_hole_extension": (3, 3, 2, 3, 2),
    "single_largest_cell_sink": (3, 3, 3, 2, 2),
    "sorted_prefix_lemma": (3, 2, 3, 2, 3),
    "fixed_window_lemma": (1, 1, 1, 1, 2),
    "conductance_order": (1, 1, 1, 1, 1),
    "raw_per_cell_max": (0, 0, 0, 2, 0),
}


def tournament():
    names = list(LENSES)
    adj = {u: set() for u in names}
    for i, u in enumerate(names):
        for v in names[i + 1 :]:
            if (LENSES[u], u) >= (LENSES[v], v):
                adj[u].add(v)
            else:
                adj[v].add(u)
    ranking = sorted(names, key=lambda n: (LENSES[n], n), reverse=True)
    return ranking, {n: len(adj[n]) for n in names}


def main():
    spans = {8: 15, 9: 14, 10: 14, 11: 14, 12: 14}
    print("LRC14 sorted-width leak classifier (codex S75)")
    print("method: exact W_a Fractions, sorted-prefix leaks vs consec/AP")
    print()
    all_records = {}
    for k, span in spans.items():
        records = leak_records(k, span)
        all_records[k] = records
        by_family = Counter(family_label(r) for r in records)
        by_repay = Counter(r["repay"] for r in records)
        print(f"== k={k}, span<={span}: leaks={len(records)} ==")
        print(f"family histogram: {dict(sorted(by_family.items()))}")
        print(f"repay-prefix histogram: {dict(sorted(by_repay.items(), key=lambda kv: str(kv[0])))}")
        for rec in records[:12]:
            level_gain = ", ".join(f"top{L}:{rec['gains'][L-1]}" for L in rec["levels"])
            print(
                "  "
                f"E={rec['E']} family={family_label(rec)} "
                f"levels={rec['levels']} gains=({level_gain}) "
                f"repay={rec['repay']} dtotal={rec['dtotal']} "
                f"holes={rec['holes']} extras={rec['extras']} gaps={rec['gaps']}"
            )
        if len(records) > 12:
            print(f"  ... {len(records) - 12} more leak records omitted")
        print()

    print("cross-k synthesis")
    union_families = Counter()
    positive_total = []
    for k, records in all_records.items():
        union_families.update(family_label(r) for r in records)
        positive_total.extend((k, r) for r in records if r["dtotal"] > 0)
    print(f"all leak families: {dict(sorted(union_families.items()))}")
    print(f"positive-total leak records: {len(positive_total)}")
    for k, rec in positive_total[:10]:
        print(f"  k={k} E={rec['E']} family={family_label(rec)} dtotal={rec['dtotal']}")

    print()
    print("Tournament Analysis over leak explanations")
    print("pair observable = (explains_k9, explains_k10, survives_wider_span, finite_family, proof_ready)")
    ranking, scores = tournament()
    for name in ranking:
        print(f"  {name}: score={scores[name]}, observable={LENSES[name]}")


if __name__ == "__main__":
    main()
