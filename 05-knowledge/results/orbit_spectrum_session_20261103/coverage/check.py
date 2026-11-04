#!/usr/bin/env python3
"""Exact certificates for the marked-clone cone and moving-mark construction."""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
from math import gcd, isqrt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / "05-knowledge/results/session134_constructive_spectrum_frontier_20261101/two_apex/word_recovery.py"
spec = importlib.util.spec_from_file_location("inherited_two_apex", SOURCE)
word = importlib.util.module_from_spec(spec)
spec.loader.exec_module(word)


def need(test, message):
    if not test:
        raise RuntimeError(message)


def graph(w):
    n = len(w) + 2
    a = [0] * n
    a[0] |= 2
    for i, c in enumerate(w, 2):
        for p in range(2):
            if int(c) >> p & 1:
                a[p] |= 1 << i
            else:
                a[i] |= 1 << p
        for j in range(2, i):
            a[j] |= 1 << i
    return tuple(a)


def induced(a, vertices):
    return tuple(sum(1 << j for j, y in enumerate(vertices)
                     if a[x] >> y & 1) for x in vertices)


def tables(a):
    n = len(a)
    end = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for v in range(n):
        end[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not mask >> v & 1:
                continue
            before = mask ^ (1 << v)
            if before:
                end[mask][v] = sum(end[before][u] for u in range(n)
                                   if before >> u & 1 and a[u] >> v & 1)
        h[mask] = sum(end[mask])
    return h, end


def profile(a, r):
    h, end = tables(a)
    full = (1 << len(a)) - 1
    q = h[full ^ (1 << r)]
    # Starts are counted independently by reversal of the adjacency relation.
    reverse = tuple(sum(1 << j for j in range(len(a)) if a[j] >> i & 1)
                    for i in range(len(a)))
    _, rev_end = tables(reverse)
    s = rev_end[full][r]
    f = sum(end[mask][r] * h[full ^ mask]
            for mask in range(1, full + 1) if mask >> r & 1)
    return h[full], s, q, f


def perm_counts(a, r):
    n = len(a)
    h = s = f = 0
    for order in itertools.permutations(range(n)):
        ok = all(a[u] >> v & 1 for u, v in zip(order, order[1:]))
        h += ok
        s += ok and order[0] == r
        k = order.index(r) + 1
        f += all(a[order[i]] >> order[i+1] & 1
                 for i in range(n - 1) if i != k - 1)
    rest = induced(a, [v for v in range(n) if v != r])
    q = sum(all(rest[u] >> v & 1 for u, v in zip(p, p[1:]))
            for p in itertools.permutations(range(n - 1)))
    return h, s, q, f


def strong(a):
    n = len(a)
    reach = [a[i] | (1 << i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(x == (1 << n) - 1 for x in reach)


def clone(a, r, k):
    n = len(a)
    b = list(a) + [0] * k
    for j in range(n, n + k):
        b[r] |= 1 << j
        b[j] |= ((1 << n) - 1) ^ (1 << r)
        for i in range(n, j):
            b[i] |= 1 << j
    return tuple(b)


def birth(a, r):
    return clone(a, r, 1)


def marked_words(a):
    n = len(a)
    answer = {}
    for x, y in itertools.combinations(range(n), 2):
        if a[y] >> x & 1:
            x, y = y, x
        rest = [i for i in range(n) if i not in (x, y)]
        order = sorted(rest, key=lambda v: -sum(a[v] >> u & 1 for u in rest))
        if not all(a[u] >> v & 1 for i, u in enumerate(order) for v in order[i+1:]):
            continue
        w = "".join(str((a[x] >> v & 1) + 2 * (a[y] >> v & 1)) for v in order)
        answer[(x, y)] = (w, tuple(order))
    return answer


def floor_strong(n):
    if n == 3:
        return 3
    if n % 3 == 0:
        return 3 * 5 ** ((n-3)//3)
    if n % 3 == 1:
        return 5 ** ((n-1)//3)
    return 9 * 5 ** ((n-5)//3)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def main():
    controls = {"suffix_profiles": 0, "clone_formula_checks": 0,
                "literal_permutation_profiles": 0, "marked_pairs": 0,
                "complete_suffix_queries": 0, "direct_named_repairs": 0}
    ext = graph("232")
    exterior_profile = profile(ext, 2)
    need(exterior_profile == (9, 3, 5, 11), "five-vertex exterior profile")
    need(perm_counts(ext, 2) == exterior_profile, "literal exterior profile")
    controls["literal_permutation_profiles"] += 1
    # Complete short suffix universe, without a strongness or favorable-cut filter.
    cones = {}
    for length in range(5):
        for v in map("".join, itertools.product("0123", repeat=length)):
            a = graph(v)
            h, s, q, f = profile(a, 1)
            delta = f-h-q+s
            A, B = f-q, q-delta
            need(delta >= 0 and delta % 2 == 0, "rooted defect")
            need(h-s >= q >= 1, "positive insertion inequality")
            need(A >= 1 and 2-A <= B <= A and B % 2 == 1, "clone cone")
            need((A, B) != (4, -1), "missing numerical pair")
            controls["suffix_profiles"] += 1
            cones.setdefault((A, B), v)
            for m in range(1, 5):
                value = word.hamilton("2" * m + v)
                need(value == A * 2**m + B, "suffix clone formula")
                controls["clone_formula_checks"] += 1
                if m >= 3:
                    need(value not in (3 * 2**m + 5, 3 * 2**m + 7),
                         "fixed-mark ternary obstruction")
                if length <= 2 and m <= 3:
                    need(tables(graph("2" * m + v))[0][-1] == value,
                         "independent suffix subset DP")
            if length <= 2:
                need(perm_counts(a, 1) == (h, s, q, f), "permutation short profile")
                controls["literal_permutation_profiles"] += 1
    # Every labelled exterior of order 2..4, at every root with an in-neighbor.
    arbitrary = 0
    for n in range(2, 5):
        pairs = list(itertools.combinations(range(n), 2))
        for bits in range(1 << len(pairs)):
            a = [0] * n
            for j, (u, v) in enumerate(pairs):
                if bits >> j & 1:
                    a[u] |= 1 << v
                else:
                    a[v] |= 1 << u
            for r in range(n):
                if not any(a[u] >> r & 1 for u in range(n)):
                    continue
                h, s, q, f = profile(a, r)
                A, B = f-q, h-s-f+2*q
                need(A >= 1 and 2-A <= B <= A and (A, B) != (4, -1),
                     "arbitrary exterior cone")
                for k in (1, 2):
                    need(tables(clone(a, r, k))[0][-1] == A * 2**k + B,
                         "arbitrary exterior clone law")
                arbitrary += 1
    rows = []
    for m in range(1, 11):
        source = graph("2" * m)
        source_h = 2**m + 1
        need(tables(source)[0][-1] == source_h and strong(source), "strong source")
        need(word.recover_full_profile(word.full_profile("2" * m)) == "2" * m,
             "source word decode and replay")
        pairs = marked_words(source)
        controls["marked_pairs"] += len(pairs)
        actual_words = {x[0] for x in pairs.values()}
        expected = {"2" * m}
        expected |= {"2" + "0"*i + "1"*(m-1-i) for i in range(m)}
        expected |= {"1"*i + "3"*(m-1-i) + "2" for i in range(m)}
        if m >= 3:
            need(actual_words == expected, "complete feedback-pair classification")
        target = 3 * source_h + 4
        nmax = 3
        while floor_strong(nmax + 1) <= target:
            nmax += 1
        cap = nmax - len(source)
        hits = {2: [], 4: []}
        for w in sorted(actual_words):
            for length in range(cap + 1):
                for suffix in map("".join, itertools.product("0123", repeat=length)):
                    if suffix.endswith("3"):
                        continue
                    h = word.hamilton(w + suffix)
                    need(strong(graph(w + suffix)), "strong completed suffix")
                    controls["complete_suffix_queries"] += 1
                    for b in (2, 4):
                        if h == 3 * source_h + b:
                            hits[b].append([w, suffix])
        if m >= 3:
            need(not hits[4], "all-marked plus-four hostile")
            need(hits[2] == sorted(hits[2]), "deterministic hit order")
        if m >= 2:
            w = "2" + "0"*(m-1) + "32"
            repaired = graph(w)
            need(word.recover_full_profile(word.full_profile(w)) == w,
                 "moving-mark child word decode and replay")
            source_profile = profile(source, 2)
            need(source_profile[3] - source_profile[0] == 1,
                 "inherited one-ready source boundary")
            # Explicit map from old source labels to re-marked word labels.
            # old x=0,y=1,c_i=i+1; new pair is (c_m,x).
            old_to_new = [1, 2] + list(range(3, m+2)) + [0]
            need(induced(repaired, old_to_new) == source, "all old arcs retained")
            old_names = ["x", "y"] + ["c_"+str(i) for i in range(1, m+1)]
            repaired_names = [None]*len(repaired)
            for name, vertex in zip(old_names, old_to_new):
                repaired_names[vertex] = name
            repaired_names[-2:] = ["birth_a", "birth_b"]
            need(len(set(repaired_names)) == len(repaired), "unique physical names")
            need(strong(repaired), "strong repaired child")
            hp = profile(repaired, 1)
            need(hp[0] == 3*2**m+5 and hp[3] == 6*2**m-3,
                 "repaired H and retained-root F")
            if m <= 4:
                need(perm_counts(repaired, 1) == hp, "literal clone skeleton values")
                controls["literal_permutation_profiles"] += 1
            if m in (2, 3, 4, 6, 9):
                first = birth(repaired, 1)
                final = birth(first, len(repaired))
                bp = profile(final, len(repaired))
                b, t = 15*2**m-1, 3*2**m+5
                need(bp == (b, t, t, b), "actual balanced reset")
                need(strong(first) and strong(final), "strong reset stages")
                controls["direct_named_repairs"] += 1
                balanced_graph = {"rows": final,
                                  "names": repaired_names+["reset_root", "reset_child"],
                                  "root": len(repaired), "child": len(repaired)+1}
            else:
                bp = None
                balanced_graph = None
            need(gcd(15*2**m-1, 3*2**m+5) == (13 if m % 12 == 11 else 1),
                 "balanced primitive selector")
            repaired_data = {"word": w, "profile_at_old_x": list(hp),
                             "old_to_new": old_to_new, "balanced_profile": bp,
                             "rows": repaired, "names": repaired_names,
                             "word_continuation_profile": word.full_profile(w),
                             "balanced_graph": balanced_graph}
        else:
            repaired_data = {"word": "220", "H": word.hamilton("220"),
                             "excluded_naive_word": "232",
                             "excluded_naive_H": word.hamilton("232")}
            need(repaired_data["H"] == 11 and repaired_data["excluded_naive_H"] == 9,
                 "moving-mark m=1 boundary")
        rows.append({"m": m, "source_H": source_h,
                     "marked_pair_count": len(pairs), "distinct_words": len(actual_words),
                     "maximum_suffix_length": cap,
                     "hits_plus_two": hits[2], "hits_plus_four": hits[4],
                     "moving_mark_child": repaired_data})
    for m, suffix, expected in ((1, "20", 11), (1, "12", 13),
                                (2, "22", 17), (2, "02", 19)):
        need(word.hamilton("2"*m + suffix) == expected, "sharp base cases")
    seven_rows = []
    for j in range(1, 8):
        b = 15 * 8**j - 1
        q = b // 7
        t = 3 * 8**j + 5
        need(b % 7 == 0 and gcd(b, t) == 1, "primitive seven lane")
        need((b % 49 == 0) == (j % 7 == 5), "mod 49 selector")
        if j <= 3:
            need(prime(q), "named seven-prime examples")
        seven_rows.append({"j": j, "H": b, "quotient": q,
                           "quotient_prime": prime(q), "t": t})
    need(all(gcd(15*2**m-1, 3*2**m+5) ==
             (13 if m % 12 == 11 else 1) for m in range(2, 122)),
         "whole modular period repeated")
    source_paths = [
        SOURCE,
        ROOT / "05-knowledge/results/return6325_completion_20261103/CLONE_LAW.md",
        ROOT / "05-knowledge/results/adaptive_metagraph_reformulation_20261103/ADAPTIVE.md",
        ROOT / "05-knowledge/results/return18979_palindrome_20261103/STABLE.md",
        ROOT / "05-knowledge/results/prime_lane_ports_five_groups_20261102/realization.md",
    ]
    print(json.dumps({
        "status": "PASS",
        "scope": "all-size proofs in NOTE.md; finite exact implementation controls",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in source_paths},
        "five_vertex_exterior_profile": exterior_profile,
        "controls": controls,
        "arbitrary_labelled_rooted_exteriors": arbitrary,
        "small_coefficient_pairs": [{"A": A, "B": B, "first_suffix": v}
                                   for (A, B), v in sorted(cones.items()) if A <= 8],
        "clone_rows": rows,
        "seven_lane": seven_rows,
        "positive_bin_skeleton_coefficients": [21, 12, 0],
        "word_continuation_suffixes": word.FULL_SUFFIXES,
        "joint_data_scope": "Seventeen two-apex suffix queries at feedback pair (c_m,x), plus full named graph rows; balanced child has its own named root and child. These are not identified with an unrestricted birth observer.",
        "open": "all-target coverage; all-m plus-four with moving feedback marks",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

