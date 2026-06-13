#!/usr/bin/env python3
"""
goodcut_transport_excess_s15.py

opus-2026-05-29-S15

Cross the ordered half-line balance view of THM-346 with the good-cut
coordinate g(t).  For selected Hamming mask families, every source tiling
contributes one ordered half-line per mask.  We stratify these half-lines by
the source good-cut bucket and ask:

  * how much transport stays inside the merged tournament bucket;
  * how cross-bucket transport decomposes into spine/ribs/sea;
  * where tournament-only and even-only projection defects sit;
  * whether g-neutral transport is the hiding place for quotient defects.

Exact for n=3..6.
"""

from collections import Counter, defaultdict

from goodcut_projection_defect_s14 import LABELS, defect_label, good_cut_count
from projection_defect_waggly_layers_s1 import tile_pairs
from tiling_quotient_bucket_balance_s5 import (
    build_context,
    masks_at_distance,
    selected_mask_families,
)


GEOMETRY = ("self", "spine", "ribs", "sea")


def pct(x, total):
    return 100 * x / total if total else 0.0


def geometry_label(ctx, bits, other):
    c1 = ctx["t_class"][bits]
    c2 = ctx["t_class"][other]
    if c1 == c2:
        return "self"
    sc1 = ctx["t_info"][c1]["SC"]
    sc2 = ctx["t_info"][c2]["SC"]
    if sc1 and sc2:
        return "spine"
    if sc1 != sc2:
        return "ribs"
    return "sea"


def add_halfline(counter, ctx, g, bits, other):
    geom = geometry_label(ctx, bits, other)
    t_changed = ctx["t_class"][bits] != ctx["t_class"][other]
    e_changed = ctx["e_class"][bits] != ctx["e_class"][other]
    label = defect_label(t_changed, e_changed)
    dg = g[other] - g[bits]

    counter["total"] += 1
    counter[f"geom:{geom}"] += 1
    counter[f"defect:{label}"] += 1
    counter["dg_abs_sum"] += abs(dg)
    counter["dg_max_abs"] = max(counter["dg_max_abs"], abs(dg))
    if dg == 0:
        counter["dg_zero"] += 1
    elif dg > 0:
        counter["dg_up"] += 1
    else:
        counter["dg_down"] += 1

    if dg != 0 and label in ("silent_both", "even_only"):
        counter["bad_nonzero_dg_without_t_change"] += 1


def mask_families_for_report(m):
    families = []
    for label, masks in selected_mask_families(m):
        if label == "all waggly":
            continue
        families.append((label, masks))
    return families


def print_bucket_rows(by_g):
    print(
        "    g  half-lines  self%  spine%   ribs%    sea%  "
        "T-only% E-only% joint%  dg0%  dg_up% dg_dn% mean|dg| bad"
    )
    for g_bucket in sorted(by_g):
        c = by_g[g_bucket]
        total = c["total"]
        mean_dg = c["dg_abs_sum"] / total if total else 0.0
        print(
            f"    {g_bucket:1d} {total:10d} "
            f"{pct(c['geom:self'], total):6.2f} "
            f"{pct(c['geom:spine'], total):7.2f} "
            f"{pct(c['geom:ribs'], total):7.2f} "
            f"{pct(c['geom:sea'], total):7.2f} "
            f"{pct(c['defect:tournament_only'], total):8.2f} "
            f"{pct(c['defect:even_only'], total):7.2f} "
            f"{pct(c['defect:joint'], total):6.2f} "
            f"{pct(c['dg_zero'], total):5.1f} "
            f"{pct(c['dg_up'], total):6.1f} "
            f"{pct(c['dg_down'], total):6.1f} "
            f"{mean_dg:8.3f} "
            f"{c['bad_nonzero_dg_without_t_change']:3d}"
        )


def summarize_family(ctx, g, label, masks):
    by_g = defaultdict(Counter)
    by_dg_geom = defaultdict(Counter)
    for bits in range(ctx["N"]):
        source_g = g[bits]
        for mask in masks:
            other = bits ^ mask
            add_halfline(by_g[source_g], ctx, g, bits, other)
            dg_abs = abs(g[other] - source_g)
            by_dg_geom[dg_abs][geometry_label(ctx, bits, other)] += 1

    print(f"\n  {label}: moves={len(masks)}")
    print_bucket_rows(by_g)
    print("    geometry by |Delta g|:")
    for dg in sorted(by_dg_geom):
        c = by_dg_geom[dg]
        total = sum(c.values())
        bits = ", ".join(f"{name}={pct(c[name], total):.1f}%" for name in GEOMETRY)
        print(f"      |dg|={dg}: {total} half-lines; {bits}")


def analyze_n(n):
    ctx = build_context(n)
    tiles = tile_pairs(n)
    g = {bits: good_cut_count(bits, tiles) for bits in range(ctx["N"])}
    g_sizes = Counter(g.values())

    print("\n" + "=" * 100)
    print(
        f"n={n}, m={ctx['m']}, tilings={ctx['N']}, "
        f"merged buckets={len(set(ctx['t_class'].values()))}, "
        f"even buckets={len(set(ctx['e_class'].values()))}, "
        f"g sizes={dict(sorted(g_sizes.items()))}"
    )

    for label, masks in mask_families_for_report(ctx["m"]):
        summarize_family(ctx, g, label, masks)

    # The full all-waggly family is useful at n<=5: it is the complete graph
    # on tilings and reveals the background sea/SC mass without swamping output.
    if ctx["m"] <= 6:
        summarize_family(ctx, g, "all waggly", list(range(1, 1 << ctx["m"])))


def main():
    print("GOOD-CUT TRANSPORT EXCESS")
    print("opus-2026-05-29-S15")
    print("Rows are ordered half-lines from source good-cut bucket g.")
    print("'bad' counts nonzero Delta g half-lines that do not change merged tournament class.")
    for n in range(3, 7):
        analyze_n(n)


if __name__ == "__main__":
    main()
