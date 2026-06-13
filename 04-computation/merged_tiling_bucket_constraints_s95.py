#!/usr/bin/env python3
"""
merged_tiling_bucket_constraints_s95.py

Exact small-n constraints for the tournament-tiling-explorer after merging
complement/isomorphism-class pairs.

The explorer fixes one Hamiltonian path P0, so its tiling cube has
    m = C(n-1, 2)
    2^m tilings.

For an isomorphism class C:
    F(C) = number of fixed-P0 tilings in C = H(C)/|Aut(C)|.

This script studies how the power of two is split after complement merging.
The key parity mechanism is:
    H(C) is odd (Redei)
    |Aut(C)| is odd (a tournament has no involutory automorphism)
so every unmerged class bucket F(C) is odd.

After complement merging:
    SC node:   M = F(C), odd
    NSC node:  M = F(C) + F(C^op) = 2F(C), twice an odd

Thus the merged tiling cube is a partition of 2^m into exactly SC_n odd
buckets and (V_n - SC_n)/2 even buckets, each even bucket having v2 exactly 1.
"""

from collections import Counter, defaultdict
from itertools import permutations
from math import comb


def tournament_from_tiling(n, bits):
    A = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        A[i][i + 1] = 1
    idx = 0
    for i in range(n):
        for j in range(i + 2, n):
            if (bits >> idx) & 1:
                A[j][i] = 1
            else:
                A[i][j] = 1
            idx += 1
    return A


def score_partition_perms(A):
    n = len(A)
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


def v2(x):
    if x == 0:
        return 10**9
    c = 0
    while x % 2 == 0:
        c += 1
        x //= 2
    return c


def compact_counter(counter, limit=14):
    rows = sorted(counter.items())
    if len(rows) <= limit:
        return dict(rows)
    head = rows[: limit // 2]
    tail = rows[-limit // 2 :]
    return {"head": head, "tail": tail, "distinct": len(rows)}


def analyze(n):
    m = comb(n - 1, 2)
    tiling_count = 1 << m
    classes = {}
    tiling_class = []

    for bits in range(tiling_count):
        A = tournament_from_tiling(n, bits)
        can = canonical(A)
        if can not in classes:
            H = hamiltonian_paths(A)
            aut = aut_size(A)
            classes[can] = {
                "rep": A,
                "H": H,
                "aut": aut,
                "F": 0,
            }
        classes[can]["F"] += 1
        tiling_class.append(can)

    cans = sorted(classes, key=lambda c: (classes[c]["H"], classes[c]["aut"], c))
    cid = {can: i for i, can in enumerate(cans)}
    for can in cans:
        info = classes[can]
        info["cid"] = cid[can]
        info["comp"] = cid[canonical(complement(info["rep"]))]
        info["sc"] = info["comp"] == info["cid"]

    fiber_mismatches = [
        (info["cid"], info["F"], info["H"], info["aut"])
        for info in classes.values()
        if info["F"] * info["aut"] != info["H"]
    ]

    merged = []
    mid_of_cid = {}
    seen = set()
    for can in cans:
        i = cid[can]
        if i in seen:
            continue
        j = classes[cans[i]]["comp"]
        members = [i] if i == j else sorted([i, j])
        mid = len(merged)
        for member in members:
            mid_of_cid[member] = mid
            seen.add(member)
        infos = [classes[cans[member]] for member in members]
        merged.append(
            {
                "mid": mid,
                "members": members,
                "type": "SC" if len(members) == 1 else "NSC",
                "M": sum(info["F"] for info in infos),
                "H": infos[0]["H"],
                "aut": infos[0]["aut"],
                "F_half": infos[0]["F"],
            }
        )

    class_to_mid = {can: mid_of_cid[cid[can]] for can in cans}
    tiling_mid = [class_to_mid[can] for can in tiling_class]

    loops = Counter()
    cross = Counter()
    for bits, mid in enumerate(tiling_mid):
        for k in range(m):
            nb = bits ^ (1 << k)
            if bits < nb:
                other = tiling_mid[nb]
                if mid == other:
                    loops[mid] += 1
                else:
                    cross[tuple(sorted((mid, other)))] += 1

    cross_incidence = Counter()
    odd_cross_incidence = Counter()
    for (a, b), weight in cross.items():
        cross_incidence[a] += weight
        cross_incidence[b] += weight
        if weight % 2:
            odd_cross_incidence[a] += 1
            odd_cross_incidence[b] += 1

    incidence_errors = []
    for node in merged:
        mid = node["mid"]
        lhs = 2 * loops[mid] + cross_incidence[mid]
        rhs = m * node["M"]
        if lhs != rhs:
            incidence_errors.append((mid, lhs, rhs))

    odd_nodes = [node for node in merged if node["M"] % 2 == 1]
    even_nodes = [node for node in merged if node["M"] % 2 == 0]
    bad_even_v2 = [node for node in even_nodes if v2(node["M"]) != 1]
    bad_odd_type = [node for node in odd_nodes if node["type"] != "SC"]
    bad_even_type = [node for node in even_nodes if node["type"] != "NSC"]
    bad_unmerged_parity = [
        info["cid"] for info in classes.values() if info["F"] % 2 == 0
    ]
    bad_order_parity = [
        info["cid"] for info in classes.values() if info["aut"] % 2 == 0
    ]

    merged_bucket_spectrum = Counter(node["M"] for node in merged)
    odd_bucket_spectrum = Counter(node["M"] for node in odd_nodes)
    even_half_spectrum = Counter(node["M"] // 2 for node in even_nodes)
    loop_spectrum = Counter(loops[node["mid"]] for node in merged)
    cross_weight_spectrum = Counter(cross.values())
    odd_cross_edges = sum(1 for weight in cross.values() if weight % 2)

    sc_mass = sum(node["M"] for node in odd_nodes)
    nsc_mass = sum(node["M"] for node in even_nodes)
    all_f2 = sum(info["F"] ** 2 for info in classes.values())
    merged_m2 = sum(node["M"] ** 2 for node in merged)
    quotient_representative_sum = sc_mass + nsc_mass // 2

    odd_weighted_cross_incidence_nodes = 0
    parity_incidence_failures = []
    if m % 2 == 1:
        for node in merged:
            expected = node["M"] % 2
            actual = cross_incidence[node["mid"]] % 2
            odd_weighted_cross_incidence_nodes += actual
            if actual != expected:
                parity_incidence_failures.append((node["mid"], actual, expected))
    else:
        for node in merged:
            odd_weighted_cross_incidence_nodes += cross_incidence[node["mid"]] % 2
            if cross_incidence[node["mid"]] % 2:
                parity_incidence_failures.append((node["mid"], 1, 0))

    return {
        "n": n,
        "m": m,
        "tiling_count": tiling_count,
        "classes": classes,
        "cans": cans,
        "merged": merged,
        "fiber_mismatches": fiber_mismatches,
        "bad_unmerged_parity": bad_unmerged_parity,
        "bad_order_parity": bad_order_parity,
        "bad_even_v2": bad_even_v2,
        "bad_odd_type": bad_odd_type,
        "bad_even_type": bad_even_type,
        "incidence_errors": incidence_errors,
        "parity_incidence_failures": parity_incidence_failures,
        "sc_mass": sc_mass,
        "nsc_mass": nsc_mass,
        "quotient_representative_sum": quotient_representative_sum,
        "all_f2": all_f2,
        "merged_m2": merged_m2,
        "loops": loops,
        "cross": cross,
        "loop_spectrum": loop_spectrum,
        "cross_weight_spectrum": cross_weight_spectrum,
        "odd_cross_edges": odd_cross_edges,
        "odd_weighted_cross_incidence_nodes": odd_weighted_cross_incidence_nodes,
        "merged_bucket_spectrum": merged_bucket_spectrum,
        "odd_bucket_spectrum": odd_bucket_spectrum,
        "even_half_spectrum": even_half_spectrum,
        "cross_incidence": cross_incidence,
        "odd_cross_incidence": odd_cross_incidence,
    }


def main():
    print("=" * 78)
    print("MERGED TILING BUCKET CONSTRAINTS S95")
    print("=" * 78)
    print("Convention: merged NSC nodes receive both complementary fixed-path fibers,")
    print("so merged bucket masses still sum to the full explorer cube 2^m.")
    print()

    summaries = []
    for n in range(3, 8):
        d = analyze(n)
        summaries.append(d)
        V = len(d["cans"])
        SC = sum(1 for can in d["cans"] if d["classes"][can]["sc"])
        Mnodes = len(d["merged"])
        NSC = Mnodes - SC
        total_loop = sum(d["loops"].values())
        total_cross = sum(d["cross"].values())
        total_edges = d["m"] * (1 << (d["m"] - 1)) if d["m"] else 0

        print("-" * 78)
        print(
            f"n={n}, m={d['m']}, 2^m={d['tiling_count']}, "
            f"classes V={V}, SC={SC}, merged=(V+SC)/2={Mnodes}"
        )
        print(
            "  theorem checks: "
            f"fiber={len(d['fiber_mismatches'])}, "
            f"F-even={len(d['bad_unmerged_parity'])}, "
            f"Aut-even={len(d['bad_order_parity'])}, "
            f"merged parity/type={len(d['bad_odd_type']) + len(d['bad_even_type'])}, "
            f"even-v2-not-1={len(d['bad_even_v2'])}"
        )
        print(
            "  bucket split: "
            f"odd SC buckets={SC}, even NSC buckets={NSC}, "
            f"SC odd mass={d['sc_mass']}, NSC even mass={d['nsc_mass']}, "
            f"sum={d['sc_mass'] + d['nsc_mass']}"
        )
        print(
            "  one-representative quotient sum: "
            f"SC_mass + NSC_mass/2 = {d['quotient_representative_sum']} "
            f"(loss={d['tiling_count'] - d['quotient_representative_sum']})"
        )
        print(
            "  second moments: "
            f"Σ_unmerged F²={d['all_f2']}, "
            f"Σ_merged M²={d['merged_m2']}, "
            f"merge bonus={d['merged_m2'] - d['all_f2']}"
        )
        print(
            "  cube-edge buckets: "
            f"loops={total_loop}, cross={total_cross}, total={total_edges}, "
            f"incidence_errors={len(d['incidence_errors'])}, "
            f"parity_incidence_failures={len(d['parity_incidence_failures'])}"
        )
        print(
            "  edge parity: "
            f"odd cross buckets={d['odd_cross_edges']}/{len(d['cross'])}, "
            f"cross weight spectrum={compact_counter(d['cross_weight_spectrum'])}"
        )
        print(
            "  bucket spectra: "
            f"M={compact_counter(d['merged_bucket_spectrum'])}, "
            f"SC odd={compact_counter(d['odd_bucket_spectrum'])}, "
            f"NSC half={compact_counter(d['even_half_spectrum'])}"
        )
        print(f"  loop spectrum={compact_counter(d['loop_spectrum'])}")

    print("=" * 78)
    print("SEQUENCES n=3..7")
    print("=" * 78)
    fields = [
        ("V classes", lambda d: len(d["cans"])),
        ("SC classes = odd merged buckets", lambda d: sum(1 for node in d["merged"] if node["type"] == "SC")),
        ("merged nodes", lambda d: len(d["merged"])),
        ("SC odd mass", lambda d: d["sc_mass"]),
        ("NSC even mass", lambda d: d["nsc_mass"]),
        ("one-representative quotient sum", lambda d: d["quotient_representative_sum"]),
        ("Σ merged M²", lambda d: d["merged_m2"]),
        ("merge bonus Σ_M M² - Σ_C F²", lambda d: d["merged_m2"] - d["all_f2"]),
        ("merged silent loop edges", lambda d: sum(d["loops"].values())),
        ("merged cross edges", lambda d: sum(d["cross"].values())),
        ("odd cross-edge buckets", lambda d: d["odd_cross_edges"]),
        ("odd weighted cross-incidence nodes", lambda d: d["odd_weighted_cross_incidence_nodes"]),
    ]
    for name, fn in fields:
        print(f"  {name}: {[fn(d) for d in summaries]}")

    print("=" * 78)
    print("CONSTRAINT STATEMENTS")
    print("=" * 78)
    print("""
1. Odd-bucket theorem.
   Every unmerged explorer fiber F(C)=H(C)/|Aut(C)| is odd. Redei gives H(C)
   odd, and tournament automorphism groups have odd order because an involution
   would reverse the edge inside one of its 2-cycles.

2. Merged parity theorem.
   In the complement-merged explorer, a node has odd mass exactly when it is
   self-complementary. NSC nodes have mass 2F(C), hence v2(mass)=1. Therefore:
       2^m = sum_{SC nodes} odd + sum_{NSC nodes} 2*odd.
   Since 2^m and every NSC bucket are even, the SC mass is even; because it is
   a sum of odd buckets, the number SC_n of self-complementary classes is even
   for every n with a nonempty fixed-path cube.

3. Self-complementary mass identity.
   Let S_n be the sum of merged odd bucket masses. Then the sum obtained by
   keeping only one representative from every NSC complement pair is
       Q_n = S_n + (2^m - S_n)/2 = 2^(m-1) + S_n/2.
   The previous parity theorem forces S_n even, so Q_n is always integral.

4. Merged second-moment identity.
   If NSC pairs have half-fiber f, merging changes the second moment by
       (2f)^2 - f^2 - f^2 = 2f^2.
   Hence
       Σ_merged M^2 = Σ_unmerged F^2 + 2 Σ_NSC-pairs F^2.
   This is the pair-count cost of hiding complement chirality.

5. Weighted quotient-edge balance.
   Let lambda_u be the number of one-tile cube edges staying inside merged node
   u, and tau_uv the number of cube edges between merged nodes u and v. Then
       2 lambda_u + Σ_{v!=u} tau_uv = m M_u.
   Loops consume two flip incidences, cross buckets consume one incidence at
   each endpoint. Summing over u recovers m*2^m = 2E(cube).

6. Incidence parity constraint.
   Because 2 lambda_u is even,
       Σ_{v!=u} tau_uv == m M_u (mod 2).
   Thus when m is odd, precisely the SC merged nodes have odd weighted
   cross-incidence; when m is even, every merged node has even weighted
   cross-incidence. This turns self-complementarity into an edge-parity
   boundary condition on the merged graph.
""")


if __name__ == "__main__":
    main()
