#!/usr/bin/env python3
"""
endpoint_sc_collision_s95.py

Probe the collision layer behind HYP-1778.

For complement-merged endpoint transfer, odd column boundary is exactly the
self-complementary child set.  HYP-1778 says non-private SC columns occur only
in the strongly connected child bucket.  This script measures the SC boundary
block, its private/decomposable part, and the residual top-SC collision block.
"""

from collections import Counter, defaultdict

import endpoint_private_goodcut_s95 as epg
import endpoint_transfer_bucket_recursion_s95 as et


def odd_column_owners(rows, ncols):
    owners = [[] for _ in range(ncols)]
    for i, row in enumerate(rows):
        for j, weight in row.items():
            if weight % 2:
                owners[j].append(i)
    return owners


def row_masks_for_columns(rows, columns):
    col_index = {col: i for i, col in enumerate(columns)}
    masks = []
    for row in rows:
        mask = 0
        for col, weight in row.items():
            if weight % 2 and col in col_index:
                mask |= 1 << col_index[col]
        masks.append(mask)
    return masks


def profile_bucket(profile):
    if len(profile) != 1:
        return "mixed"
    return next(iter(profile))


def node_label(level, mid):
    node = level["merged"][mid]
    pieces = [f"mid={mid}", node["type"], f"M={node['M']}"]
    if len(node["members"]) == 1:
        cid = node["members"][0]
        info = level["classes"][level["cans"][cid]]
        pieces.extend([f"H={info['H']}", f"F={info['F']}", f"aut={info['aut']}"])
    else:
        hs = []
        fs = []
        for cid in node["members"]:
            info = level["classes"][level["cans"][cid]]
            hs.append(info["H"])
            fs.append(info["F"])
        pieces.extend([f"Hs={hs}", f"Fs={fs}"])
    return ":".join(pieces)


def summarize_transition(n, parent, child):
    transfer = et.analyze_transfer(parent, child)
    rows = transfer["merged_rows"]
    child_profiles = epg.class_goodcut_profiles(child, merged=True)
    parent_profiles = epg.class_goodcut_profiles(parent, merged=True)
    owners = odd_column_owners(rows, len(child["merged"]))

    child_top_bucket = child["n"] - 1
    groups = {
        "all": list(range(len(child["merged"]))),
        "sc": [],
        "decomp_sc": [],
        "top_sc": [],
        "nsc": [],
        "private_sc": [],
        "nonprivate_sc": [],
    }

    support_spectrum_by_group = defaultdict(Counter)
    nonprivate_examples = []
    owner_private_decomp = Counter()
    owner_private_top = Counter()
    owner_nonprivate_top_incidence = Counter()
    nonprivate_owner_bucket_tuples = Counter()
    nonprivate_child_h_mod8 = Counter()
    nonprivate_child_h = []

    for mid, node in enumerate(child["merged"]):
        bucket = profile_bucket(child_profiles[mid])
        support = len(owners[mid])
        if node["type"] == "SC":
            groups["sc"].append(mid)
            if bucket == child_top_bucket:
                groups["top_sc"].append(mid)
            else:
                groups["decomp_sc"].append(mid)
            if support == 1:
                groups["private_sc"].append(mid)
                support_spectrum_by_group["private_sc"][support] += 1
                owner = owners[mid][0]
                if bucket == child_top_bucket:
                    owner_private_top[owner] += 1
                else:
                    owner_private_decomp[owner] += 1
            else:
                groups["nonprivate_sc"].append(mid)
                support_spectrum_by_group["nonprivate_sc"][support] += 1
                owner_buckets = tuple(
                    sorted(profile_bucket(parent_profiles[owner]) for owner in owners[mid])
                )
                nonprivate_owner_bucket_tuples[owner_buckets] += 1
                cid = node["members"][0]
                h_value = child["classes"][child["cans"][cid]]["H"]
                nonprivate_child_h_mod8[h_value % 8] += 1
                nonprivate_child_h.append(h_value)
                if bucket == child_top_bucket:
                    for owner in owners[mid]:
                        owner_nonprivate_top_incidence[owner] += 1
                if len(nonprivate_examples) < 20:
                    nonprivate_examples.append((mid, bucket, support, owners[mid]))
        else:
            groups["nsc"].append(mid)

        if node["type"] == "SC":
            name = "top_sc" if bucket == child_top_bucket else "decomp_sc"
            support_spectrum_by_group[name][support] += 1
            support_spectrum_by_group["sc"][support] += 1
        else:
            support_spectrum_by_group["nsc"][support] += 1

    print("\n" + "-" * 78)
    print(f"n={n}->{n+1}: parent merged={len(parent['merged'])}, child merged={len(child['merged'])}")
    print(f"  merged GF2 rank={transfer['merged_parity']['rank']}")
    print(f"  child top bucket={child_top_bucket}")

    for name in ["sc", "decomp_sc", "top_sc", "nsc", "private_sc", "nonprivate_sc"]:
        columns = groups[name]
        rank = et.gf2_rank(row_masks_for_columns(rows, columns))
        print(
            f"  block {name}: cols={len(columns)}, rank={rank}, "
            f"support_spectrum={dict(sorted(support_spectrum_by_group[name].items()))}"
        )

    sc_by_bucket = Counter()
    private_sc_by_bucket = Counter()
    nonprivate_sc_by_bucket = Counter()
    for mid in groups["sc"]:
        bucket = profile_bucket(child_profiles[mid])
        sc_by_bucket[bucket] += 1
        if len(owners[mid]) == 1:
            private_sc_by_bucket[bucket] += 1
        else:
            nonprivate_sc_by_bucket[bucket] += 1

    print(f"  SC by bucket={dict(sorted(sc_by_bucket.items()))}")
    print(f"  private SC by bucket={dict(sorted(private_sc_by_bucket.items()))}")
    print(f"  non-private SC by bucket={dict(sorted(nonprivate_sc_by_bucket.items()))}")
    print(f"  parent owners with decomp-private SC={len(owner_private_decomp)}")
    print(f"  parent owners with top-private SC={len(owner_private_top)}")
    print(f"  parent owners incident to top non-private SC={len(owner_nonprivate_top_incidence)}")
    if nonprivate_owner_bucket_tuples:
        print(f"  non-private SC owner-bucket tuples={dict(sorted(nonprivate_owner_bucket_tuples.items()))}")
        print(f"  non-private SC child H mod 8={dict(sorted(nonprivate_child_h_mod8.items()))}")
        print(f"  non-private SC child H values={sorted(nonprivate_child_h)}")

    if nonprivate_examples:
        print("  non-private SC examples:")
        for mid, bucket, support, owner_list in nonprivate_examples:
            owner_labels = [
                f"{owner}:g={dict(parent_profiles[owner])}:M={parent['merged'][owner]['M']}"
                for owner in owner_list
            ]
            print(f"    {node_label(child, mid)}:g={bucket}:support={support}")
            print(f"      owners={owner_labels}")

    return {
        "rank": transfer["merged_parity"]["rank"],
        "sc": len(groups["sc"]),
        "decomp_sc": len(groups["decomp_sc"]),
        "top_sc": len(groups["top_sc"]),
        "private_sc": len(groups["private_sc"]),
        "nonprivate_sc": len(groups["nonprivate_sc"]),
        "nonprivate_sc_by_bucket": nonprivate_sc_by_bucket,
        "decomp_rank": et.gf2_rank(row_masks_for_columns(rows, groups["decomp_sc"])),
        "top_rank": et.gf2_rank(row_masks_for_columns(rows, groups["top_sc"])),
        "sc_rank": et.gf2_rank(row_masks_for_columns(rows, groups["sc"])),
        "nonprivate_rank": et.gf2_rank(row_masks_for_columns(rows, groups["nonprivate_sc"])),
    }


def main():
    print("=" * 78)
    print("ENDPOINT SC COLLISION S95")
    print("=" * 78)
    print("Merged endpoint-transfer SC boundary: private pivots vs collision block.")

    levels = {n: et.build_level(n) for n in range(2, 8)}
    summaries = []
    for n in range(2, 7):
        summaries.append(summarize_transition(n, levels[n], levels[n + 1]))

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    for key in [
        "rank",
        "sc",
        "decomp_sc",
        "top_sc",
        "private_sc",
        "nonprivate_sc",
        "decomp_rank",
        "top_rank",
        "sc_rank",
        "nonprivate_rank",
    ]:
        print(f"  {key}: {[summary[key] for summary in summaries]}")
    print(
        "  nonprivate_sc_by_bucket: "
        f"{[dict(sorted(summary['nonprivate_sc_by_bucket'].items())) for summary in summaries]}"
    )

    print("\n" + "=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print("""
The merged SC boundary splits into a decomposable part and a strongly connected
top bucket.  Through 6->7, every non-private SC column lies in the top bucket.
Thus endpoint-transfer rank has a triangular decomposable boundary plus a
smaller strongly connected collision block.
""".strip())


if __name__ == "__main__":
    main()
