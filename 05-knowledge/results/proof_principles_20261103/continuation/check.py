"""Exact bounded audit for PROOF.md; standard library only, -O safe."""
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

from scout import add, cyclic, endpoints, h, inflate, reverse, root_profiles, tournament


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def hp_permutations(rows):
    return sum(all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def kernel(rows):
    n = len(rows)
    e = endpoints(rows)
    b = endpoints(reverse(rows))
    full = (1 << n) - 1
    c = [[0] * n for _ in rows]
    for s in range(1, full):
        for u in range(n):
            if not e[s][u]:
                continue
            for v in range(n):
                c[u][v] += e[s][u] * b[full ^ s][v]
    w = [[0] * n for _ in rows]
    for u, v in combinations(range(n), 2):
        require((c[u][v] + c[v][u]) % 2 == 0, "nonintegral inherited weight")
        w[u][v] = w[v][u] = (c[u][v] + c[v][u]) // 2
        require(w[u][v] > 0, "nonpositive inherited weight")
    return sum(e[-1]), w


def signature(rows):
    base, w = kernel(rows)
    n = len(rows)
    d = tuple(sum(w[u][v] for v in range(n) if rows[u] >> v & 1)
              for u in range(n))
    return (base, tuple(w[u][v] for u, v in combinations(range(n), 2)), d), w


def reaches(rows):
    reach = [r | 1 << u for u, r in enumerate(rows)]
    for k in range(len(rows)):
        for u in range(len(rows)):
            if reach[u] >> k & 1:
                reach[u] |= reach[k]
    return reach


def is_strong(rows):
    return all(r == (1 << len(rows)) - 1 for r in reaches(rows))


def check_profiles():
    tested = []
    # Every cyclic Cayley tournament obtained by choosing one member of
    # each pair {d,-d}, at orders 3,5,7. Translation proves transitivity.
    for n in (3, 5, 7):
        half = (n - 1) // 2
        for choice in range(1 << half):
            steps = [d if choice >> (d - 1) & 1 else n - d
                     for d in range(1, half + 1)]
            rows = [sum(1 << ((u + d) % n) for d in steps) for u in range(n)]
            b, e = root_profiles(rows)
            require(b == e, ("vertex-transitive root balance", n, choice))
            require(all(b[u] == b[0] for u in range(n)), "root transitivity")
            tested.append([n, choice])
    return tested


def check_ambient():
    totals = {}
    permutation_checks = 0
    for module_order, largest_quotient in ((3, 4), (5, 3)):
        m = cyclic(module_order)
        total = 0
        quotient_count = 0
        for qn in range(1, largest_quotient + 1):
            for qm in range(1 << (qn * (qn - 1) // 2)):
                q = tournament(qn, qm)
                a = inflate(q, m)
                b = inflate(q, reverse(m))
                require(h(a) == h(b), "ambient base H")
                quotient_count += 1
                for mask in range(1 << len(a)):
                    aa, bb = add(a, mask), add(b, mask)
                    ha, hb = h(aa), h(bb)
                    require(ha == hb, ("module reversal", module_order, qn, qm, mask))
                    if len(aa) <= 6:
                        require(ha == hp_permutations(aa), "first independent count")
                        require(hb == hp_permutations(bb), "second independent count")
                        permutation_checks += 2
                    total += 1
        totals[str(module_order)] = {
            "quotient_orders": [1, largest_quotient],
            "quotients": quotient_count,
            "all_probe_masks_compared": total,
        }
    return totals, permutation_checks


def check_collisions():
    summary = []
    for n in range(1, 6):
        groups = {}
        for m in range(1 << (n * (n - 1) // 2)):
            rows = tournament(n, m)
            sig, w = signature(rows)
            groups.setdefault(sig, []).append((rows, w))
        pairs = 0
        smallest_difference = None
        for group in groups.values():
            for (a, wa), (b, wb) in combinations(group, 2):
                require(wa == wb, "signature weights disagree")
                ra, rb = reaches(a), reaches(b)
                flip_count = 0
                divergences = [0] * n
                for u, v in combinations(range(n), 2):
                    if (a[u] >> v & 1) == (b[u] >> v & 1):
                        continue
                    flip_count += 1
                    x, y = (u, v) if a[u] >> v & 1 else (v, u)
                    divergences[x] += wa[u][v]
                    divergences[y] -= wa[u][v]
                    require(ra[x] >> y & 1 and ra[y] >> x & 1, "flipped arc crosses SCC")
                require(all(x == 0 for x in divergences), "collision not balanced")
                require(flip_count >= 3, "collision changes fewer than 3 arcs")
                require(all(bool(ra[u] >> v & 1 and ra[v] >> u & 1) ==
                            bool(rb[u] >> v & 1 and rb[v] >> u & 1)
                            for u in range(n) for v in range(n)), "SCC partitions differ")
                for mask in range(1 << n):
                    require(h(add(a, mask)) == h(add(b, mask)), "signature collision not full")
                smallest_difference = flip_count if smallest_difference is None else min(smallest_difference, flip_count)
                pairs += 1
        summary.append({"n": n, "labelled_bases": 1 << (n * (n - 1) // 2),
                        "response_classes": len(groups), "collision_pairs": pairs,
                        "smallest_arc_difference": smallest_difference})
    return summary


def controls():
    # Strong arbitrary ambient quotient, followed by identical successive probes.
    a = [10, 12, 9, 16, 7]
    b = [12, 9, 10, 16, 7]
    require(is_strong(a) and is_strong(b), "ambient not strong")
    require([h(a), h(b)] == [15, 15], "ambient H")
    require(all(h(add(a, x)) == h(add(b, x)) for x in range(32)), "ambient one-step")
    aa, bb = add(a, 1), add(b, 1)
    require([h(aa), h(bb)] == [31, 31], "ambient first attachment")
    aaa, bbb = add(aa, 2), add(bb, 2)
    require([h(aaa), h(bbb)] == [81, 87], "ambient two-step witness")
    require([hp_permutations(aaa), hp_permutations(bbb)] == [81, 87], "independent strong witness")

    # Regular degree is a strictly weaker sidecar.
    regular = [22, 88, 50, 69, 104, 11, 37]
    require(all(r.bit_count() == 3 for r in regular), "hostile not regular")
    starts, ends = root_profiles(regular)
    require((starts[0][1], ends[0][1]) == (24, 25), "regular endpoint imbalance")
    require([h(regular), h(reverse(regular))] == [171, 171], "regular base")
    require([h(add(regular, 1)), h(add(reverse(regular), 1))] == [403, 407], "regular separation")

    # Original triangle ceases to be a module after a nonuniform attachment.
    triangle = [2, 4, 1]
    nonmod_a, nonmod_b = add(triangle, 4), add(reverse(triangle), 4)
    require([h(nonmod_a), h(nonmod_b)] == [5, 5], "nonmodule first stage")
    require([h(add(nonmod_a, 2)), h(add(nonmod_b, 2))] == [13, 11], "nonmodule hostile")

    # A circulation in the old kernel alone is insufficient.
    c = [24, 1, 3, 6, 14]
    d = [10, 16, 3, 6, 13]
    hc, wc = kernel(c)
    hd, wd = kernel(d)
    tri = (0, 1, 4)
    require(all(sum(c[u] >> v & 1 for v in tri) == 1 for u in tri), "not directed triangle")
    require({wc[u][v] for u, v in combinations(tri, 2)} == {4}, "old weights not balanced")
    require((hc, hd) == (11, 15) and wc != wd, "old-flow hostile")
    require([hp_permutations(c), hp_permutations(d)] == [11, 15], "independent flow hostile")

    return {"strong_ambient": {"base": [15, 15], "first": [31, 31], "second": [81, 87]},
            "regular_seven": {"base": [171, 171], "endpoint_at_0": [24, 25], "probe": [403, 407]},
            "nonmodule_triangle": {"first": [5, 5], "second": [13, 11]},
            "old_weight_circulation": {"old_triangle_weights": [4, 4, 4], "counts": [11, 15]}}


def main():
    here = Path(__file__).resolve().parent
    frozen = {name: sha256((here / name).read_bytes()).hexdigest()
              for name in ("check.py", "scout.py")}
    profiles = check_profiles()
    ambient, permutation_checks = check_ambient()
    collisions = check_collisions()
    witnesses = controls()
    require(frozen == {name: sha256((here / name).read_bytes()).hexdigest() for name in frozen},
            "source changed during computation")
    print(json.dumps({"status": "PASS", "source_sha256": frozen,
                      "vertex_transitive_profile_controls": profiles,
                      "ambient_universe": ambient,
                      "independent_ambient_permutation_counts": permutation_checks,
                      "collision_census": collisions, "witnesses": witnesses}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
