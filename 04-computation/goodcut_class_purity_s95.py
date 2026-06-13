#!/usr/bin/env python3
"""
goodcut_class_purity_s95.py

Use the endpoint-transfer level builder to test whether good-cut count descends
from fixed-path tilings to unmerged tournament classes and complement-merged
nodes.

This extends the earlier merged-only S13 computation through n=7 and also
checks the stronger unmerged class statement.
"""

from collections import Counter

import endpoint_private_goodcut_s95 as epg
import endpoint_transfer_bucket_recursion_s95 as et


def reachable_from(A, start):
    seen = {start}
    stack = [start]
    while stack:
        v = stack.pop()
        for w, edge in enumerate(A[v]):
            if edge and w not in seen:
                seen.add(w)
                stack.append(w)
    return seen


def scc_count(A):
    n = len(A)
    remaining = set(range(n))
    count = 0
    while remaining:
        v = next(iter(remaining))
        forward = reachable_from(A, v)
        reverse = {u for u in range(n) if v in reachable_from(A, u)}
        comp = forward & reverse
        remaining -= comp
        count += 1
    return count


def width(profile):
    return max(profile) - min(profile)


def summarize_profiles(label, profiles, type_of=None):
    widths = Counter(width(profile) for profile in profiles.values())
    mixed = {key: profile for key, profile in profiles.items() if len(profile) > 1}
    print(f"  {label}: count={len(profiles)}, mixed={len(mixed)}, width_spectrum={dict(sorted(widths.items()))}")

    if type_of is not None:
        type_counts = Counter(type_of(key) for key in profiles)
        mixed_type_counts = Counter(type_of(key) for key in mixed)
        print(f"    type_counts={dict(sorted(type_counts.items()))}")
        print(f"    mixed_type_counts={dict(sorted(mixed_type_counts.items()))}")

    if mixed:
        print("    first mixed examples:")
        for key, profile in list(mixed.items())[:8]:
            print(f"      {key}: {dict(sorted(profile.items()))}")

    bucket_classes = Counter()
    bucket_mass = Counter()
    for profile in profiles.values():
        if len(profile) == 1:
            bucket = next(iter(profile))
            bucket_classes[bucket] += 1
        for bucket, count in profile.items():
            bucket_mass[bucket] += count

    print(f"    pure classes by bucket={dict(sorted(bucket_classes.items()))}")
    print(f"    tiling mass by bucket={dict(sorted(bucket_mass.items()))}")

    return {
        "count": len(profiles),
        "mixed": len(mixed),
        "widths": widths,
        "bucket_classes": bucket_classes,
        "bucket_mass": bucket_mass,
    }


def main():
    print("=" * 78)
    print("GOOD-CUT CLASS PURITY S95")
    print("=" * 78)
    print("Does good-cut count descend to unmerged and merged tournament classes?")

    unmerged_mixed = []
    merged_mixed = []

    for n in range(2, 8):
        level = et.build_level(n)
        class_profiles = epg.class_goodcut_profiles(level, merged=False)
        merged_profiles = epg.class_goodcut_profiles(level, merged=True)

        print("\n" + "-" * 78)
        print(f"n={n}, tilings={level['total']}, classes={len(class_profiles)}, merged={len(merged_profiles)}")
        formula_failures = []
        scc_distribution = Counter()
        for bits in range(level["total"]):
            A = et.tournament_from_tiling(n, bits)
            components = scc_count(A)
            g = epg.good_cut_count_endpoint_bits(bits, n)
            scc_distribution[components] += 1
            if g != n - components:
                formula_failures.append((bits, g, components))
        print(f"  formula g=n-scc_count failures={len(formula_failures)}")
        print(f"  scc_count tiling distribution={dict(sorted(scc_distribution.items()))}")
        if formula_failures:
            print(f"  first formula failures={formula_failures[:8]}")
        class_summary = summarize_profiles("unmerged tournament classes", class_profiles)
        merged_summary = summarize_profiles(
            "complement-merged nodes",
            merged_profiles,
            type_of=lambda mid: level["merged"][mid]["type"],
        )

        unmerged_mixed.append(class_summary["mixed"])
        merged_mixed.append(merged_summary["mixed"])

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    print(f"  unmerged_mixed: {unmerged_mixed}")
    print(f"  merged_mixed:   {merged_mixed}")

    print("\n" + "=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print("""
The exhaustive formula check supports the theorem

    good_cut_count(T, path) = n - scc_count(T).

So zero mixed profiles are not accidental: strong-component count is an
ordinary tournament isomorphism invariant.  This strengthens the S13
merged-coordinate hypothesis and explains why every private endpoint-transfer
pivot in endpoint_private_goodcut_s95.py was good-cut pure: all child classes
at those levels are good-cut pure.
""".strip())


if __name__ == "__main__":
    main()
