"""Independent exact small-universe check for the historical five-core bridge.

Run from any directory: python check_five_cores.py
No third-party packages or inherited counting code. Writes receipt.json here.
"""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json


CORES = [(0,), (2, 4, 1), (2, 4, 9, 3),
         (2, 4, 9, 19, 7), (2, 4, 25, 3, 11)]


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def strong(rows):
    n = len(rows)
    full = (1 << n) - 1
    for root in range(n):
        seen = frontier = 1 << root
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            fresh = rows[bit.bit_length() - 1] & ~seen
            seen |= fresh
            frontier |= fresh
        if seen != full:
            return False
    return True


def paths_dp(rows):
    n = len(rows)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            available = rows[v] & ~mask
            while available:
                bit = available & -available
                available ^= bit
                w = bit.bit_length() - 1
                dp[mask | bit][w] += count
    return sum(dp[full])


def paths_permutation(rows):
    return sum(all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def relabel(rows, p):
    return tuple(sum(((rows[p[i]] >> p[j]) & 1) << j
                     for j in range(len(rows))) for i in range(len(rows)))


def canonical(rows):
    return min(relabel(rows, p) for p in permutations(range(len(rows))))


def converse(rows):
    full = (1 << len(rows)) - 1
    return tuple(full ^ (1 << v) ^ row for v, row in enumerate(rows))


def from_bits(n, bits):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if bits >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def rooted_profiles(rows):
    """Independent permutation/bin enumeration, with one empty path on empty set."""
    n = len(rows)
    full = (1 << n) - 1
    paths = {mask: [p for p in permutations([i for i in range(n) if mask >> i & 1])
                    if all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))]
             for mask in range(full + 1)}
    h = len(paths[full])
    check(h == paths_dp(rows), "rooted permutation count disagrees with subset DP")
    result = []
    for r in range(n):
        s = sum(p[0] == r for p in paths[full])
        q = len(paths[full ^ (1 << r)])
        f = sum(sum(p[-1] == r for p in paths[mask]) * len(paths[full ^ mask])
                for mask in range(full + 1) if mask >> r & 1)
        result.append({"root": r, "s_start": s, "q_delete": q,
                       "F_ordered_two_bins": f, "ready_gap": f-h,
                       "joint_delta": f-h-q+s})
    return result


def birth_beaten_only_by(rows, r):
    n = len(rows)
    return tuple(row | ((1 << n) if v == r else 0) for v, row in enumerate(rows)) + (
        ((1 << n) - 1) ^ (1 << r),)


def ready_word(word):
    """Literal chart r=0,a=1,b=2, chain=3,..., with a->b."""
    n = len(word) + 3
    rows = [0] * n
    rows[1] |= 1 | (1 << 2)
    rows[2] |= 1
    for j, letter in enumerate(word):
        u = j + 3
        rows[0] |= 1 << u
        for v in range(u + 1, n):
            rows[u] |= 1 << v
        for port, bit in enumerate(letter, start=1):
            if bit == "1":
                rows[u] |= 1 << port
            else:
                rows[port] |= 1 << u
    return tuple(rows)


def main():
    census = []
    small_classes = {}
    independent_cases = 0
    for n in range(1, 7):
        histogram = Counter()
        strong_count = 0
        universe = 1 << (n * (n - 1) // 2)
        for bits in range(universe):
            rows = from_bits(n, bits)
            h = paths_dp(rows)
            if n <= 5:
                check(h == paths_permutation(rows), f"count mismatch {rows}")
                independent_cases += 1
            if strong(rows):
                strong_count += 1
                histogram[h] += 1
                if h <= 9:
                    key = canonical(rows)
                    small_classes[key] = h
        census.append({"n": n, "labelled_universe": universe,
                       "labelled_strong": strong_count,
                       "strong_minimum": min(histogram, default=None),
                       "strong_H_histogram": dict(sorted(histogram.items()))})
    expected = {canonical(rows): paths_dp(rows) for rows in CORES}
    check(small_classes == expected, "five cores fail complete bounded-H census")
    details = []
    for i, rows in enumerate(CORES):
        aut = sum(relabel(rows, p) == rows for p in permutations(range(len(rows))))
        anti = [list(p) for p in permutations(range(len(rows)))
                if relabel(rows, p) == converse(rows)]
        dual_images = [j for j, target in enumerate(CORES)
                       if len(target) == len(rows)
                       and canonical(target) == canonical(converse(rows))]
        check(dual_images == [i], "current bank is not pointwise converse-fixed")
        details.append({"index": i, "rows": list(rows), "order": len(rows),
                        "H": paths_dp(rows), "scores": sorted(r.bit_count() for r in rows),
                        "automorphism_order": aut, "converse_isomorphisms": anti,
                        "converse_bank_images": dual_images,
                        "canonical_rows": list(canonical(rows)),
                        "rooted_profiles": rooted_profiles(rows)})
    check(census[-1]["strong_minimum"] == 15, "six-vertex strong floor mismatch")
    check([d["H"] for d in details] == [1, 3, 5, 9, 9], "wrong count list")
    check([d["automorphism_order"] for d in details] == [1, 3, 1, 1, 1], "aut mismatch")
    check(details[3]["scores"] == details[4]["scores"], "H9 scores differ")
    check(details[3]["canonical_rows"] != details[4]["canonical_rows"], "H9 types collapse")
    ready_checks = []
    for i, r, gap, child_h in [(2, 1, 2, 19), (2, 3, 1, 17),
                              (3, 2, 2, 31), (4, 4, 1, 29)]:
        rows = CORES[i]
        profile = details[i]["rooted_profiles"][r]
        child = birth_beaten_only_by(birth_beaten_only_by(rows, r), len(rows))
        check(profile["ready_gap"] == gap, "wrong ready gap")
        check(paths_dp(child) == child_h, "reset identity/control mismatch")
        check(child_h == 3 * paths_dp(rows) + 2 * gap, "reset formula mismatch")
        ready_checks.append({"core_index": i, "root": r, "gap": gap,
                             "child_rows": list(child), "child_H": child_h,
                             "child_rooted_profiles": rooted_profiles(child)})
    check(1 not in [p["ready_gap"] for p in details[3]["rooted_profiles"]],
          "first H9 core unexpectedly 1-ready")
    check(2 not in [p["ready_gap"] for p in details[4]["rooted_profiles"]],
          "second H9 core unexpectedly 2-ready")
    for reset in ready_checks:
        check(all(p["ready_gap"] not in (1, 2) for p in reset["child_rooted_profiles"]),
              "standard reset unexpectedly retains a ready root")
    for record, blocks in [(ready_checks[2], [[5], [6], [3], [4], [0], [1]]),
                           (ready_checks[3], [[4, 5, 6], [3], [0], [1]])]:
        child = record["child_rows"]
        check(sorted(v for block in blocks for v in block) == [0, 1, 3, 4, 5, 6],
              "wrong apex deletion vertex set")
        for i, left in enumerate(blocks):
            for right in blocks[i+1:]:
                check(all(child[a] >> b & 1 for a in left for b in right),
                      "ordinal decomposition cross edge mismatch")
        for block in blocks:
            induced = tuple(sum(((child[a] >> b) & 1) << j for j, b in enumerate(block))
                            for a in block)
            check(strong(induced), "deletion component is not strong")
        record["one_apex_bank_return"] = {
            "apex_old_vertex": 2,
            "remaining_sccs_in_forward_order": blocks,
            "scope": "Literal one-apex bank membership; no ready-root return"}
    ready_embeddings = []
    for word, core_index, root in [(["11"], 2, 1), (["11", "01"], 3, 2)]:
        rows = ready_word(word)
        target = CORES[core_index]
        maps = [p for p in permutations(range(len(rows)))
                if relabel(rows, p) == target and p[root] == 0]
        check(bool(maps), "literal old ready-word/core embedding failed")
        ready_embeddings.append({"word": word, "word_rows": list(rows),
                                 "core_index": core_index, "core_root": root,
                                 "old_vertices_in_core_order": list(maps[0])})
    receipt = {
        "status": "FINITE-EXACT; all-order core classification additionally uses inherited strong subtournament theorem",
        "baseline_commit": "0f96e904a5",
        "universe": "All labelled tournaments on vertex sets {0,...,n-1}, 1<=n<=6, singleton counted as strong",
        "inherited_filters": "None; exhaustive orientations, then exact reachability and H<=9",
        "positive_controls": "Five supplied row-mask cores are recovered; H=1,3,5,9,9",
        "hostile_controls": "H9 cores same score multiset but nonisomorphic; converse fixes both individually; strong order6 floor15",
        "independent_path": "All n<=5 Hamilton counts checked by all vertex permutations versus subset dynamic programming",
        "independent_cases": independent_cases,
        "labelled_total": sum(item["labelled_universe"] for item in census),
        "reproduction": "python 05-knowledge/results/chanpak_five_core_bridges_20261110/history/check_five_cores.py",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "core_details": details,
        "ready_reset_checks": ready_checks,
        "literal_ready_word_embeddings": ready_embeddings,
        "census": census,
    }
    target = Path(__file__).with_name("receipt.json")
    target.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"labelled_total": receipt["labelled_total"],
                      "independent_cases": independent_cases,
                      "core_H": [d["H"] for d in details],
                      "aut_orders": [d["automorphism_order"] for d in details],
                      "strong_minima": [d["strong_minimum"] for d in census],
                      "all_cores_individually_self_dual": True,
                      "receipt": str(target)}, indent=2))


if __name__ == "__main__":
    main()
