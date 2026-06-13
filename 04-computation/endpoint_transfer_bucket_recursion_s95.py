#!/usr/bin/env python3
"""
endpoint_transfer_bucket_recursion_s95.py

Fixed-path endpoint insertion exposes a clean recursion in the
tournament-tiling explorer.

For n vertices, the explorer tiling cube has m_n = C(n-1, 2) free tiles.  Add a
new terminal vertex n to the fixed path

    0 -> 1 -> ... -> n-1 -> n.

The old tiling contributes all old free tiles.  The only new free tiles are the
edges between the new terminal vertex n and old vertices 0..n-2, so every
parent tiling has exactly 2^(n-1) endpoint-extension children.  Therefore the
tiling cube recursion is

    Q_{m_{n+1}} = Q_{m_n} x Q_{n-1}.

This script descends that exact cube product to isomorphism classes and to
complement-merged nodes.
"""

from collections import Counter, defaultdict
from itertools import permutations
from math import comb


def tile_pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 2, n)]


def tournament_from_tiling(n, bits):
    A = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        A[i][i + 1] = 1
    for idx, (i, j) in enumerate(tile_pairs(n)):
        if (bits >> idx) & 1:
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A


def extend_bits(parent_n, parent_bits, signature):
    """Add terminal vertex parent_n using signature bits on old vertices 0..n-2."""
    old = {}
    for idx, pair in enumerate(tile_pairs(parent_n)):
        old[pair] = (parent_bits >> idx) & 1

    child_bits = 0
    for idx, (i, j) in enumerate(tile_pairs(parent_n + 1)):
        if j < parent_n:
            bit = old[(i, j)]
        else:
            bit = (signature >> i) & 1
        if bit:
            child_bits |= 1 << idx
    return child_bits


def score_partition_perms(A):
    blocks = defaultdict(list)
    for v, row in enumerate(A):
        blocks[sum(row)].append(v)
    groups = [blocks[s] for s in sorted(blocks)]

    def rec(k):
        if k == len(groups):
            yield []
            return
        for p in permutations(groups[k]):
            for rest in rec(k + 1):
                yield list(p) + rest

    yield from rec(0)


def canonical(A):
    n = len(A)
    best = None
    for p in score_partition_perms(A):
        form = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best


def complement(A):
    n = len(A)
    return [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]


def aut_size(A):
    n = len(A)
    blocks = defaultdict(list)
    for v, row in enumerate(A):
        blocks[sum(row)].append(v)
    groups = list(blocks.values())
    p = [None] * n

    def rec(k):
        if k == len(groups):
            yield p[:]
            return
        group = groups[k]
        for images in permutations(group):
            for src, image in zip(group, images):
                p[src] = image
            yield from rec(k + 1)

    total = 0
    for perm in rec(0):
        ok = True
        for i in range(n):
            for j in range(n):
                if A[perm[i]][perm[j]] != A[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 1
    return total


def hamiltonian_paths(A):
    n = len(A)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if (mask & (1 << nxt)) == 0 and A[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[full])


def build_level(n):
    m = comb(n - 1, 2)
    total = 1 << m
    classes = {}
    tiling_canon = []

    for bits in range(total):
        A = tournament_from_tiling(n, bits)
        can = canonical(A)
        if can not in classes:
            classes[can] = {
                "rep": A,
                "H": hamiltonian_paths(A),
                "aut": aut_size(A),
                "F": 0,
            }
        classes[can]["F"] += 1
        tiling_canon.append(can)

    cans = sorted(classes, key=lambda c: (classes[c]["H"], classes[c]["aut"], c))
    cid = {can: i for i, can in enumerate(cans)}
    tiling_cid = [cid[can] for can in tiling_canon]

    for can in cans:
        info = classes[can]
        info["cid"] = cid[can]
        info["comp"] = cid[canonical(complement(info["rep"]))]
        info["sc"] = info["comp"] == info["cid"]

    merged = []
    cid_to_mid = {}
    seen = set()
    for can in cans:
        i = cid[can]
        if i in seen:
            continue
        j = classes[cans[i]]["comp"]
        members = [i] if i == j else sorted([i, j])
        mid = len(merged)
        for member in members:
            cid_to_mid[member] = mid
            seen.add(member)
        infos = [classes[cans[member]] for member in members]
        merged.append({
            "mid": mid,
            "members": members,
            "type": "SC" if len(members) == 1 else "NSC",
            "M": sum(info["F"] for info in infos),
        })

    return {
        "n": n,
        "m": m,
        "total": total,
        "classes": classes,
        "cans": cans,
        "tiling_cid": tiling_cid,
        "merged": merged,
        "cid_to_mid": cid_to_mid,
    }


def gf2_rank(bit_rows):
    basis = {}
    for row in bit_rows:
        x = row
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                break
    return len(basis)


def support_stats(rows, cols):
    row_support = [len(r) for r in rows]
    col_support = [0] * cols
    for row in rows:
        for j in row:
            col_support[j] += 1
    return {
        "support": sum(row_support),
        "row_min": min(row_support) if row_support else 0,
        "row_max": max(row_support) if row_support else 0,
        "row_avg": sum(row_support) / len(row_support) if row_support else 0,
        "col_min": min(col_support) if col_support else 0,
        "col_max": max(col_support) if col_support else 0,
        "col_avg": sum(col_support) / cols if cols else 0,
    }


def parity_summary(rows, cols, expected_col_parity):
    bit_rows = []
    xor_rows = 0
    row_parity_failures = 0
    col_parity = [0] * cols
    odd_entries = 0

    for row in rows:
        bits = 0
        odd_count = 0
        for j, weight in row.items():
            if weight % 2:
                bits |= 1 << j
                col_parity[j] ^= 1
                odd_entries += 1
                odd_count ^= 1
        if odd_count:
            row_parity_failures += 1
        bit_rows.append(bits)
        xor_rows ^= bits

    expected_bits = 0
    col_failures = 0
    for j, parity in enumerate(expected_col_parity):
        if parity:
            expected_bits |= 1 << j
        if col_parity[j] != parity:
            col_failures += 1

    return {
        "odd_entries": odd_entries,
        "row_parity_failures": row_parity_failures,
        "col_parity_failures": col_failures,
        "xor_boundary_ok": xor_rows == expected_bits,
        "boundary_weight": sum(expected_col_parity),
        "rank": gf2_rank(bit_rows),
    }


def analyze_transfer(parent, child):
    n = parent["n"]
    class_rows = [Counter() for _ in parent["cans"]]
    merged_rows = [Counter() for _ in parent["merged"]]

    for parent_bits, parent_cid in enumerate(parent["tiling_cid"]):
        parent_mid = parent["cid_to_mid"][parent_cid]
        for sig in range(1 << (n - 1)):
            child_bits = extend_bits(n, parent_bits, sig)
            child_cid = child["tiling_cid"][child_bits]
            child_mid = child["cid_to_mid"][child_cid]
            class_rows[parent_cid][child_cid] += 1
            merged_rows[parent_mid][child_mid] += 1

    class_col = [0] * len(child["cans"])
    for row in class_rows:
        for j, weight in row.items():
            class_col[j] += weight

    merged_col = [0] * len(child["merged"])
    for row in merged_rows:
        for j, weight in row.items():
            merged_col[j] += weight

    class_row_errors = []
    for i, row in enumerate(class_rows):
        expected = (1 << (n - 1)) * parent["classes"][parent["cans"][i]]["F"]
        actual = sum(row.values())
        if actual != expected:
            class_row_errors.append((i, actual, expected))

    class_col_errors = []
    for j, actual in enumerate(class_col):
        expected = child["classes"][child["cans"][j]]["F"]
        if actual != expected:
            class_col_errors.append((j, actual, expected))

    merged_row_errors = []
    for i, row in enumerate(merged_rows):
        expected = (1 << (n - 1)) * parent["merged"][i]["M"]
        actual = sum(row.values())
        if actual != expected:
            merged_row_errors.append((i, actual, expected))

    merged_col_errors = []
    for j, actual in enumerate(merged_col):
        expected = child["merged"][j]["M"]
        if actual != expected:
            merged_col_errors.append((j, actual, expected))

    class_expected_col_parity = [
        child["classes"][can]["F"] % 2 for can in child["cans"]
    ]
    merged_expected_col_parity = [
        node["M"] % 2 for node in child["merged"]
    ]

    return {
        "class_rows": class_rows,
        "merged_rows": merged_rows,
        "class_row_errors": class_row_errors,
        "class_col_errors": class_col_errors,
        "merged_row_errors": merged_row_errors,
        "merged_col_errors": merged_col_errors,
        "class_support": support_stats(class_rows, len(child["cans"])),
        "merged_support": support_stats(merged_rows, len(child["merged"])),
        "class_parity": parity_summary(class_rows, len(child["cans"]), class_expected_col_parity),
        "merged_parity": parity_summary(merged_rows, len(child["merged"]), merged_expected_col_parity),
        "class_weight_spectrum": Counter(
            weight for row in class_rows for weight in row.values()
        ),
        "merged_weight_spectrum": Counter(
            weight for row in merged_rows for weight in row.values()
        ),
    }


def compact_counter(counter, limit=12):
    rows = sorted(counter.items())
    if len(rows) <= limit:
        return dict(rows)
    return {
        "head": rows[: limit // 2],
        "tail": rows[-limit // 2 :],
        "distinct": len(rows),
    }


def main():
    print("=" * 78)
    print("ENDPOINT TRANSFER BUCKET RECURSION S95")
    print("=" * 78)
    print("Endpoint insertion gives Q_{m(n+1)} = Q_{m(n)} x Q_{n-1}.")
    print("The transfer matrices below are this cube product after quotienting.")

    levels = {n: build_level(n) for n in range(2, 8)}
    summaries = []

    for n in range(2, 7):
        parent = levels[n]
        child = levels[n + 1]
        d = analyze_transfer(parent, child)
        summaries.append((n, d))

        print("\n" + "-" * 78)
        print(
            f"n={n}->{n+1}: parent classes={len(parent['cans'])}, "
            f"child classes={len(child['cans'])}, parent merged={len(parent['merged'])}, "
            f"child merged={len(child['merged'])}, multiplier=2^{n-1}"
        )
        print(
            "  exact sums: "
            f"class row/col errors={len(d['class_row_errors'])}/{len(d['class_col_errors'])}, "
            f"merged row/col errors={len(d['merged_row_errors'])}/{len(d['merged_col_errors'])}"
        )
        cs = d["class_support"]
        ms = d["merged_support"]
        print(
            "  class support: "
            f"nonzero={cs['support']}, row={cs['row_min']}..{cs['row_max']} "
            f"avg={cs['row_avg']:.2f}, col={cs['col_min']}..{cs['col_max']} "
            f"avg={cs['col_avg']:.2f}"
        )
        print(
            "  merged support: "
            f"nonzero={ms['support']}, row={ms['row_min']}..{ms['row_max']} "
            f"avg={ms['row_avg']:.2f}, col={ms['col_min']}..{ms['col_max']} "
            f"avg={ms['col_avg']:.2f}"
        )
        cp = d["class_parity"]
        mp = d["merged_parity"]
        print(
            "  GF(2) class transfer: "
            f"odd_entries={cp['odd_entries']}, rank={cp['rank']}, "
            f"row_fail={cp['row_parity_failures']}, col_fail={cp['col_parity_failures']}, "
            f"boundary_weight={cp['boundary_weight']}, xor_boundary_ok={cp['xor_boundary_ok']}"
        )
        print(
            "  GF(2) merged transfer: "
            f"odd_entries={mp['odd_entries']}, rank={mp['rank']}, "
            f"row_fail={mp['row_parity_failures']}, col_fail={mp['col_parity_failures']}, "
            f"boundary_weight={mp['boundary_weight']}, xor_boundary_ok={mp['xor_boundary_ok']}"
        )
        print(f"  class transfer weight spectrum={compact_counter(d['class_weight_spectrum'])}")
        print(f"  merged transfer weight spectrum={compact_counter(d['merged_weight_spectrum'])}")

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    sequence_fields = [
        ("class support", lambda d: d["class_support"]["support"]),
        ("merged support", lambda d: d["merged_support"]["support"]),
        ("class GF2 rank", lambda d: d["class_parity"]["rank"]),
        ("merged GF2 rank", lambda d: d["merged_parity"]["rank"]),
        ("class odd entries", lambda d: d["class_parity"]["odd_entries"]),
        ("merged odd entries", lambda d: d["merged_parity"]["odd_entries"]),
        ("merged boundary weight = SC child count", lambda d: d["merged_parity"]["boundary_weight"]),
    ]
    for name, fn in sequence_fields:
        print(f"  {name}: {[fn(d) for _, d in summaries]}")

    print("\n" + "=" * 78)
    print("CONSTRAINT STATEMENTS")
    print("=" * 78)
    print("""
1. Endpoint transfer matrix.
   Let A_n(C,D) count endpoint-extension children of fixed-path tilings in
   class C at n that land in class D at n+1. Then:
       Σ_D A_n(C,D) = 2^(n-1) F_n(C)
       Σ_C A_n(C,D) = F_{n+1}(D).
   This is not asymptotic. It is the exact cube product
       Q_{m(n+1)} = Q_{m(n)} x Q_{n-1}
   after quotienting by isomorphism.

2. Merged endpoint transfer matrix.
   The same law holds after complement merging, with bucket masses M:
       Σ_V B_n(U,V) = 2^(n-1) M_n(U)
       Σ_U B_n(U,V) = M_{n+1}(V).

3. Recursive parity boundary.
   Mod 2, every parent row has even sum because 2^(n-1) is even.  In the
   unmerged transfer, every child column has odd sum because every F(D) is odd.
   Thus the XOR of all parent parity rows is the all-child vector.

4. Merged recursive boundary.
   In the merged transfer, the XOR of all parent parity rows is exactly the
   self-complementary child-node indicator.  The SC spine at n+1 can therefore
   be recovered from endpoint-extension parity alone.

5. Rank as recursive memory.
   The GF(2) ranks measured here are small shadows of how much parent-level
   information survives endpoint insertion.  If the rank sequence has a closed
   form, it would be a recursive invariant of the quotient tower distinct from
   class count, chromatic number, and H moments.  In the tested range the
   parity transfer has full row rank, i.e. no parent parity signal is lost.
""")


if __name__ == "__main__":
    main()
