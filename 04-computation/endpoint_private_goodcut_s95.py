#!/usr/bin/env python3
"""
endpoint_private_goodcut_s95.py

Probe whether endpoint-transfer private odd child witnesses are visible in the
good-cut/Morse coordinate.

The endpoint-transfer scripts use the fixed path 0 -> 1 -> ... -> n-1.  The
good-cut scripts use the reflected path n-1 -> ... -> 0.  Reflection preserves
the staircase, so a bit set in endpoint coordinates contributes the reflected
interval of cuts.
"""

from collections import Counter, defaultdict

import endpoint_transfer_bucket_recursion_s95 as et


def good_cut_count_endpoint_bits(bits, n):
    cuts = set()
    for idx, (i, j) in enumerate(et.tile_pairs(n)):
        if (bits >> idx) & 1:
            hi = n - 1 - i
            lo = n - 1 - j
            cuts.update(range(lo + 1, hi + 1))
    return len(cuts)


def class_goodcut_profiles(level, merged=False):
    profiles = defaultdict(Counter)
    for bits, cid in enumerate(level["tiling_cid"]):
        key = level["cid_to_mid"][cid] if merged else cid
        profiles[key][good_cut_count_endpoint_bits(bits, level["n"])] += 1
    return profiles


def odd_support(rows):
    col_support = Counter()
    owners = defaultdict(list)
    for i, row in enumerate(rows):
        for j, weight in row.items():
            if weight % 2:
                col_support[j] += 1
                owners[j].append(i)
    private = {j: owners[j][0] for j, s in col_support.items() if s == 1}
    return col_support, private


def summarize_private(
    name,
    parent,
    child,
    rows,
    child_profiles,
    parent_profiles,
    label_child,
    child_kind=None,
):
    _, private = odd_support(rows)
    by_owner = defaultdict(list)
    for child_id, owner in private.items():
        by_owner[owner].append(child_id)

    private_bucket_spectrum = Counter()
    private_bucket_columns = Counter()
    private_profile_width = Counter()
    private_kind_spectrum = Counter()
    owner_private_bucket_span = Counter()
    owner_no_private = 0
    examples = []

    for owner in range(len(rows)):
        child_ids = by_owner.get(owner, [])
        if not child_ids:
            owner_no_private += 1
            continue
        buckets = set()
        for child_id in child_ids:
            profile = child_profiles[child_id]
            if child_kind is not None:
                private_kind_spectrum[child_kind(child_id)] += 1
            private_profile_width[max(profile) - min(profile)] += 1
            if len(profile) == 1:
                private_bucket_columns[next(iter(profile))] += 1
            else:
                private_bucket_columns["mixed"] += 1
            for g, count in profile.items():
                private_bucket_spectrum[g] += count
                buckets.add(g)
        owner_private_bucket_span[max(buckets) - min(buckets)] += 1
        if len(examples) < 8:
            examples.append((owner, dict(parent_profiles[owner]), [
                (label_child(child_id), dict(child_profiles[child_id]))
                for child_id in child_ids[:4]
            ]))

    print(f"  {name}: private child columns={len(private)}, owners covered={len(by_owner)}/{len(rows)}")
    print(f"    owners without private child={owner_no_private}")
    print(f"    private child good-cut mass={dict(sorted(private_bucket_spectrum.items()))}")
    print(f"    private child columns by good-cut={dict(sorted(private_bucket_columns.items(), key=lambda kv: kv[0] if isinstance(kv[0], int) else 10**9))}")
    print(f"    private child profile width spectrum={dict(sorted(private_profile_width.items()))}")
    if child_kind is not None:
        print(f"    private child kind spectrum={dict(sorted(private_kind_spectrum.items()))}")
        all_sc_buckets = Counter()
        private_sc_buckets = Counter()
        for child_id, profile in child_profiles.items():
            bucket = next(iter(profile)) if len(profile) == 1 else "mixed"
            if child_kind(child_id) == "SC":
                all_sc_buckets[bucket] += 1
                if child_id in private:
                    private_sc_buckets[bucket] += 1
        missing_sc_buckets = all_sc_buckets.copy()
        missing_sc_buckets.subtract(private_sc_buckets)
        missing_sc_buckets = +missing_sc_buckets
        print(f"    all SC child columns by good-cut={dict(sorted(all_sc_buckets.items()))}")
        print(f"    private SC child columns by good-cut={dict(sorted(private_sc_buckets.items()))}")
        print(f"    non-private SC child columns by good-cut={dict(sorted(missing_sc_buckets.items()))}")
    print(f"    owner private-bucket span spectrum={dict(sorted(owner_private_bucket_span.items()))}")
    print("    examples:")
    for owner, parent_prof, child_rows in examples:
        print(f"      owner={owner}, parent_g={parent_prof}")
        for child_label, child_prof in child_rows:
            print(f"        child={child_label}, child_g={child_prof}")

    return {
        "private_cols": len(private),
        "owners": len(by_owner),
        "owner_no_private": owner_no_private,
        "private_bucket_spectrum": private_bucket_spectrum,
        "private_bucket_columns": private_bucket_columns,
        "private_profile_width": private_profile_width,
        "private_kind_spectrum": private_kind_spectrum,
        "owner_private_bucket_span": owner_private_bucket_span,
    }


def main():
    print("=" * 78)
    print("ENDPOINT PRIVATE GOOD-CUT S95")
    print("=" * 78)
    print("Do private endpoint-transfer child witnesses have good-cut structure?")

    levels = {n: et.build_level(n) for n in range(2, 8)}
    sequences = {
        "t_private_cols": [],
        "t_owners": [],
        "t_width0": [],
        "m_private_cols": [],
        "m_owners": [],
        "m_width0": [],
    }

    for n in range(2, 7):
        print("\n" + "-" * 78)
        print(f"n={n}->{n+1}")
        parent = levels[n]
        child = levels[n + 1]
        transfer = et.analyze_transfer(parent, child)
        parent_class_g = class_goodcut_profiles(parent, merged=False)
        child_class_g = class_goodcut_profiles(child, merged=False)
        parent_merged_g = class_goodcut_profiles(parent, merged=True)
        child_merged_g = class_goodcut_profiles(child, merged=True)

        class_summary = summarize_private(
            "unmerged tournament",
            parent,
            child,
            transfer["class_rows"],
            child_class_g,
            parent_class_g,
            lambda cid: f"cid={cid}:H={child['classes'][child['cans'][cid]]['H']}:F={child['classes'][child['cans'][cid]]['F']}",
        )
        merged_summary = summarize_private(
            "merged tournament",
            parent,
            child,
            transfer["merged_rows"],
            child_merged_g,
            parent_merged_g,
            lambda mid: f"mid={mid}:type={child['merged'][mid]['type']}:M={child['merged'][mid]['M']}",
            child_kind=lambda mid: child["merged"][mid]["type"],
        )

        sequences["t_private_cols"].append(class_summary["private_cols"])
        sequences["t_owners"].append(class_summary["owners"])
        sequences["t_width0"].append(class_summary["private_profile_width"].get(0, 0))
        sequences["m_private_cols"].append(merged_summary["private_cols"])
        sequences["m_owners"].append(merged_summary["owners"])
        sequences["m_width0"].append(merged_summary["private_profile_width"].get(0, 0))

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    for key, seq in sequences.items():
        print(f"  {key}: {seq}")

    print("\n" + "=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print("""
1. If private child columns are good-cut pure, then the triangular rank witness
   respects the good-cut Morse coordinate.

2. If private witnesses concentrate in high good-cut buckets, endpoint
   parity-injectivity is tied to crossing many cut barriers; if they spread
   across all buckets, the rank mechanism is independent of the Morse height.

3. Merged private witnesses are expected to be less pure because complement
   merging hides chirality.  Any persistent good-cut purity there strengthens
   the conjecture that good-cut count descends to G_n/Z_2.
""")


if __name__ == "__main__":
    main()
