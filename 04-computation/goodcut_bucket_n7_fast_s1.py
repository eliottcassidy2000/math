#!/usr/bin/env python3
"""
goodcut_bucket_n7_fast_s1.py

kind-pasteur-2026-05-30-S1

Fast n=7 probe for the good-cut descent hypothesis:

  Does g(tau) = #good cuts descend from fixed-base tilings to merged
  tournament isomorphism classes G_7/Z_2?

The earlier exact script `goodcut_bucket_merged_s13.py` verifies n=3..6
using full permutation canonization.  This script keeps the same fixed
base-path tiling model but accelerates n=7 by hashing tilings with several
isomorphism invariants before canonicalizing within the remaining buckets.
"""

from collections import Counter, defaultdict
from itertools import permutations
import time


def tile_pairs(n):
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            tiles.append((x, y))
    return tiles


def mask_to_adj(mask, n, tiles):
    verts = list(range(n, 0, -1))
    vert_idx = {v: i for i, v in enumerate(verts)}
    adj = [[0] * n for _ in range(n)]

    # Fixed base path: n -> n-1 -> ... -> 1.
    for i in range(n - 1):
        adj[i][i + 1] = 1

    for k, (hi, lo) in enumerate(tiles):
        i = vert_idx[hi]
        j = vert_idx[lo]
        if (mask >> k) & 1:
            adj[j][i] = 1
        else:
            adj[i][j] = 1
    return adj


def complement_adj(adj):
    n = len(adj)
    return [[adj[j][i] for j in range(n)] for i in range(n)]


def hamiltonian_paths(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            val = dp[mask][v]
            if not val:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += val
    return sum(dp[(1 << n) - 1])


def c3_count(adj):
    n = len(adj)
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    total += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    total += 1
    return total


def deletion_fingerprint(adj):
    n = len(adj)
    fp = []
    for skip in range(n):
        sub = []
        for i in range(n):
            if i == skip:
                continue
            row = []
            for j in range(n):
                if j != skip:
                    row.append(adj[i][j])
            sub.append(row)
        fp.append(hamiltonian_paths(sub))
    return tuple(sorted(fp))


def score_sequence(adj):
    return tuple(sorted((sum(row) for row in adj), reverse=True))


def neighborhood_signature(adj):
    n = len(adj)
    out = [sum(row) for row in adj]
    sigs = []
    for v in range(n):
        succ = tuple(sorted((out[u] for u in range(n) if adj[v][u]), reverse=True))
        pred = tuple(sorted((out[u] for u in range(n) if adj[u][v]), reverse=True))
        sigs.append((out[v], succ, pred))
    return tuple(sorted(sigs, reverse=True))


def invariant_key(adj):
    return (
        score_sequence(adj),
        c3_count(adj),
        hamiltonian_paths(adj),
        deletion_fingerprint(adj),
        neighborhood_signature(adj),
    )


def grouped_permutations(groups):
    if not groups:
        yield []
        return
    first, *rest_groups = groups
    for p in permutations(first):
        for rest in grouped_permutations(rest_groups):
            yield list(p) + rest


def canonical_form(adj):
    n = len(adj)
    scores = [sum(row) for row in adj]
    by_score = defaultdict(list)
    for v, score in enumerate(scores):
        by_score[score].append(v)
    groups = [by_score[s] for s in sorted(by_score)]

    best = None
    for perm in grouped_permutations(groups):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best


def good_cut_set(mask, tiles):
    cuts = set()
    for k, (hi, lo) in enumerate(tiles):
        if (mask >> k) & 1:
            cuts.update(range(lo + 1, hi + 1))
    return frozenset(cuts)


def build_classes(n):
    tiles = tile_pairs(n)
    total = 1 << len(tiles)
    t0 = time.time()

    print(f"n={n}, tiles={len(tiles)}, tilings={total}")
    print("Phase 1: invariant hashing")
    hash_groups = defaultdict(list)
    for mask in range(total):
        adj = mask_to_adj(mask, n, tiles)
        hash_groups[invariant_key(adj)].append((mask, adj))
        if (mask + 1) % 4096 == 0:
            elapsed = time.time() - t0
            print(f"  {mask + 1}/{total} tilings, groups={len(hash_groups)}, elapsed={elapsed:.1f}s")

    print(
        f"  groups={len(hash_groups)}, largest={max(len(v) for v in hash_groups.values())}, "
        f"elapsed={time.time() - t0:.1f}s"
    )

    print("Phase 2: canonical split inside hash groups")
    canon_to_cid = {}
    class_info = []
    mask_to_cid = {}
    ambiguous_groups = 0

    for group in hash_groups.values():
        local = defaultdict(list)
        for mask, adj in group:
            local[canonical_form(adj)].append((mask, adj))
        if len(local) > 1:
            ambiguous_groups += 1

        for canon, members in local.items():
            if canon not in canon_to_cid:
                cid = len(canon_to_cid)
                canon_to_cid[canon] = cid
                rep_adj = members[0][1]
                class_info.append(
                    {
                        "cid": cid,
                        "canon": canon,
                        "rep_mask": members[0][0],
                        "rep_adj": rep_adj,
                        "H": hamiltonian_paths(rep_adj),
                        "score": score_sequence(rep_adj),
                        "size": 0,
                    }
                )
            cid = canon_to_cid[canon]
            class_info[cid]["size"] += len(members)
            for mask, _ in members:
                mask_to_cid[mask] = cid

    print(
        f"  unmerged classes={len(class_info)}, ambiguous_hash_groups={ambiguous_groups}, "
        f"elapsed={time.time() - t0:.1f}s"
    )

    print("Phase 3: complement merge and good-cut profiles")
    for info in class_info:
        comp_canon = canonical_form(complement_adj(info["rep_adj"]))
        info["comp_cid"] = canon_to_cid[comp_canon]
        info["SC"] = info["comp_cid"] == info["cid"]

    pair_to_mid = {}
    cid_to_mid = {}
    for info in class_info:
        pair = tuple(sorted((info["cid"], info["comp_cid"])))
        if pair not in pair_to_mid:
            pair_to_mid[pair] = len(pair_to_mid)
        cid_to_mid[info["cid"]] = pair_to_mid[pair]

    merged_profiles = defaultdict(Counter)
    merged_H = {}
    bucket_counts = Counter()
    cutset_counts = Counter()
    for mask in range(total):
        g = len(good_cut_set(mask, tiles))
        mid = cid_to_mid[mask_to_cid[mask]]
        merged_profiles[mid][g] += 1
        bucket_counts[g] += 1
        cutset_counts[good_cut_set(mask, tiles)] += 1
        merged_H.setdefault(mid, class_info[mask_to_cid[mask]]["H"])

    mixed = {
        mid: profile for mid, profile in merged_profiles.items()
        if len(profile) > 1
    }
    max_span = max(max(p) - min(p) for p in merged_profiles.values())
    pure_class_buckets = Counter(next(iter(profile)) for profile in merged_profiles.values() if len(profile) == 1)

    print("\nRESULT")
    print(f"  unmerged classes={len(class_info)}")
    print(f"  SC unmerged classes={sum(1 for info in class_info if info['SC'])}")
    print(f"  merged classes={len(merged_profiles)}")
    print(f"  good-cut tiling counts={dict(sorted(bucket_counts.items()))}")
    print(f"  distinct good-cut sets={len(cutset_counts)} of {2 ** (n - 1)} possible")
    print(
        f"  merged purity: pure={len(merged_profiles) - len(mixed)}, "
        f"mixed={len(mixed)}, max_span={max_span}"
    )
    print(f"  pure merged classes by g={dict(sorted(pure_class_buckets.items()))}")

    if mixed:
        print("  HYP-1764 status at n=7: REFUTED")
        print("  first mixed merged profiles:")
        for mid, profile in sorted(mixed.items())[:20]:
            print(f"    mid={mid:3d} H={merged_H[mid]:4d} profile={dict(sorted(profile.items()))}")
    else:
        print("  HYP-1764 status at n=7: CONFIRMED")

    print(f"  total elapsed={time.time() - t0:.1f}s")


def main():
    print("GOOD-CUT DESCENT TO MERGED CLASSES: FAST n=7 PROBE")
    print("kind-pasteur-2026-05-30-S1")
    build_classes(7)


if __name__ == "__main__":
    main()
