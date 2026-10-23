"""Exact controls for marked-overlap gluing and zero-slot information loss.

Run with Python normally and with -O.  No assertion is a validity gate.
The finite universe and hostiles are described in REPORT.md.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import importlib.util
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise ValueError(message)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    return tuple(p.index(i) for i in range(len(p)))


def affine(a, b, p):
    return tuple((a * x + b) % p for x in range(p))


def fixed(p, marks):
    return all(p[x] == x for x in marks)


def automorphism(p, adjacency):
    return all(adjacency[i][j] == adjacency[p[i]][p[j]]
               for i in range(len(p)) for j in range(len(p)))


def chart_controls():
    """Independent 7! automorphism census; literal three-cell hostiles."""
    p = 7
    paley = [[int((v - u) % p in (1, 2, 4)) for v in range(p)]
             for u in range(p)]
    group = tuple(affine(a, b, p) for a in (1, 2, 4) for b in range(p))
    actual = tuple(g for g in permutations(range(p)) if automorphism(g, paley))
    need(set(actual) == set(group), "affine H21 is not the full P7 automorphism group")
    ident = tuple(range(p))
    dilation = affine(2, 0, p)
    translation = affine(1, 1, p)
    for g, h in product(group, repeat=2):
        need(compose(g, h) in group, "H21 is not closed")
    base_profile = {}
    for k in (0, 1, 2):
        hist = Counter()
        for marks in combinations(range(p), k):
            hist[len([g for g in group if fixed(g, marks)])] += 1
        base_profile[str(k)] = dict(sorted(hist.items()))
    need(base_profile == {"0": {21: 1}, "1": {3: 7}, "2": {1: 21}},
         "H21 base profile")
    need(fixed(dilation, [0]) and not fixed(dilation, [0, 1]),
         "one-owner dilation hostile")

    def host(twist):
        out = [[0] * 21 for _ in range(21)]
        for i in range(3):
            for u in range(7):
                for v in range(7):
                    out[7*i+u][7*i+v] = paley[u][v]
        for i, j in combinations(range(3), 2):
            transport = twist if (i, j) == (0, 2) else ident
            for u in range(7):
                for v in range(7):
                    out[7*i+u][7*j+v] = int(v == transport[u])
                    out[7*j+v][7*i+u] = 1 - out[7*i+u][7*j+v]
        return out

    source, target = host(ident), host(dilation)
    for matrix in (source, target):
        need(all(matrix[i][j]+matrix[j][i] == 1
                 for i in range(21) for j in range(i+1, 21)), "not a tournament")
        need([sum(matrix[7*i]) for i in range(3)] == [5, 10, 15],
             "three actual cells lost degree distinction")

    def compatible(f, g, i, j, left, right):
        return all(left[7*i+u][7*j+v] == right[7*i+f[u]][7*j+g[v]]
                   for u in range(7) for v in range(7))

    counts = {}
    for marks in ((), (0,), (0, 1)):
        allowed = [g for g in group if fixed(g, marks)]
        pairs = {}
        for i, j in combinations(range(3), 2):
            pairs[(i, j)] = {(f, g) for f, g in product(allowed, repeat=2)
                             if compatible(f, g, i, j, source, target)}
        solutions = [charts for charts in product(allowed, repeat=3)
                     if all((charts[i], charts[j]) in pairs[(i, j)]
                            for i, j in pairs)]
        need(not solutions, "twisted target admitted a global chart")
        pair_counts = [len(pairs[e]) for e in sorted(pairs)]
        expected = [21, 21, 21] if not marks else (
            [3, 3, 3] if len(marks) == 1 else [1, 0, 1])
        need(pair_counts == expected, "pairwise chart count differs")
        counts[str(len(marks))] = {"pair_counts_01_02_12": pair_counts,
                                   "global_chart_count": len(solutions)}
    owner_sections = [v for v in product(range(7), repeat=3)
                      if v[0] == v[1] == v[2] and v[2] == dilation[v[0]]]
    need(owner_sections == [(0, 0, 0)], "owner section differs from chart section")
    # Existence of chart comparisons and uniqueness of charts have different
    # groups.  Tree compatibility makes every automorphism diagonal here.
    centralizers = {}
    for name, matrix in (("flat", source), ("twisted", target)):
        actual_diagonal = []
        for g in group:
            full = tuple(7*i+g[u] for i in range(3) for u in range(7))
            if automorphism(full, matrix):
                actual_diagonal.append(g)
        hist = [sum(fixed(g, marks) for g in actual_diagonal)
                for marks in ((), (0,), (1,), (0, 1))]
        expected = [21, 3, 3, 1] if name == "flat" else [3, 3, 1, 1]
        need(hist == expected, "centralizer/marking profile differs")
        centralizers[name] = hist
    commutator = compose(dilation, compose(translation,
                  compose(inverse(dilation), inverse(translation))))
    need(commutator == translation, "H21 commutator is not translation one")
    need(not any(translation[x] == x for x in range(7)), "translation acquired a fixed owner")
    # The universal base criterion is checked for every subgroup of H21.
    subgroups = {frozenset([ident])}
    todo = list(subgroups)
    while todo:
        old = todo.pop()
        for g in group:
            enlarged = set(old) | {g}
            change = True
            while change:
                nxt = enlarged | {compose(a, b) for a in enlarged for b in enlarged}
                change = len(nxt) != len(enlarged)
                enlarged = nxt
            enlarged = frozenset(enlarged)
            if enlarged not in subgroups:
                subgroups.add(enlarged)
                todo.append(enlarged)
    base_tests = 0
    for subgroup in subgroups:
        for mask in range(1 << 7):
            marks = [i for i in range(7) if mask >> i & 1]
            stabilizer = [g for g in subgroup if fixed(g, marks)]
            # Identity is the only undetected loop exactly at bases.
            need((len(stabilizer) == 1) ==
                 all(g == ident for g in subgroup if fixed(g, marks)),
                 "base criterion")
            base_tests += 1
    return {"P7_automorphisms": len(actual), "base_profile": base_profile,
            "literal_order21_chart_hostile": counts,
            "owner_sections": owner_sections,
            "automorphisms_none_0_1_01": centralizers,
            "commutator": commutator,
            "H21_subgroups": len(subgroups), "all_subset_base_tests": base_tests}


def fibonacci_controls():
    # Import only the inherited finite grammar, not its result file.
    path = HERE.parents[1] / "odd_cell_gluing_20261023" / "check_fibonacci.py"
    spec = importlib.util.spec_from_file_location("inherited_zero_slot", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    word = module.U
    need(word.startswith("CABCACABCABCACAB"), "user's sixteen-letter prefix")
    fib = [0, 1]
    while len(fib) < len(module.H):
        fib.append(fib[-1] + fib[-2])
    for k in range(4, len(module.W)):
        need(module.W[k].count("C") == fib[k-3], "full-block C count")
    counts = []
    for n in range(2001):
        exact = sum(1 if i == 3 else fib[i-3] for i in module.zeck(n))
        need(exact == word[:n].count("C"), "Zeckendorf-prefix C count")
        counts.append(exact)

    # Distinct C/B replacements are distinct coherent histories.  Check every
    # one of the 256 choices through21, with all shared shorter prefixes.
    n = 21
    positions = [i for i, c in enumerate(word[:n], 1) if c == "C"]
    need(len(positions) == 8, "21-prefix should have eight independent positions")
    permutations3 = tuple(permutations(range(3)))
    identity3 = tuple(range(3))
    # Assign all three letters different permutations; the test is not a
    # constant-product collision manufactured by identical letter maps.
    letter_map = {"A": (1, 2, 0), "B": (0, 2, 1), "C": (1, 0, 2)}
    need(len(set(letter_map.values())) == 3, "letter transports must differ")
    buckets = defaultdict(list)
    outputs = set()
    for mask in range(1 << len(positions)):
        changed = {pos for j, pos in enumerate(positions) if mask >> j & 1}
        text = "".join("B" if i in changed else c for i, c in enumerate(word[:n], 1))
        for size in range(n+1):
            need(module.row(size, changed) == text[:size], "recolouring lost prefix compatibility")
        result = identity3
        for c in text:
            result = compose(letter_map[c], result)
        need(result in permutations3, "not a permutation product")
        buckets[(len(changed), result)].append((mask, text))
        outputs.add(text)
    need(len(outputs) == 2**len(positions), "independent family has merged words")
    collision = next((key, members[:2]) for key, members in sorted(buckets.items()) if len(members) > 1)
    key, members = collision
    first_difference = next(i+1 for i, (a, b) in enumerate(zip(members[0][1], members[1][1])) if a != b)
    need(key[0] > 0 and members[0][1] != members[1][1], "nontrivial charge/product collision")

    # Given fixed boundary locations, all interval products are necessary
    # and sufficient for preserving every boundary pin query on bijective
    # chains.  Exhaust all S3 words of length3 and a shared port after step1.
    pinned = defaultdict(set)
    for a, b, c in product(permutations3, repeat=3):
        products = (a, compose(c, b))
        relation = frozenset((x, a[x], c[b[a[x]]]) for x in range(3))
        pinned[products].add(relation)
    need(len(pinned) == 36 and all(len(v) == 1 for v in pinned.values()),
         "interval products do not determine pinned boundary relation")
    need(len({next(iter(v)) for v in pinned.values()}) == 36,
         "distinct interval products have merged boundary relations")
    return {"prefix_limit": 2000, "user_prefix": word[:16],
            "independent_choices_at_21": len(positions),
            "coherent_prefixes_at_21": len(outputs),
            "same_charge_product_bucket_count": len(buckets),
            "largest_same_charge_product_fibre": max(map(len, buckets.values())),
            "collision": {"replacement_count": key[0], "product": key[1],
                          "histories": members, "first_different_owner": first_difference},
            "interval_product_relations": len(pinned),
            "large_block": {"index": len(module.W)-1,
                            "length": len(module.W[-1]),
                            "independent_bits": module.W[-1].count("C")}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = {"charts": chart_controls(), "fibonacci": fibonacci_controls(), "checks": CHECKS}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    destination = HERE / "receipt.json"
    if args.write:
        destination.write_text(text, encoding="utf-8")
    else:
        need(json.loads(destination.read_text(encoding="utf-8")) == json.loads(text),
             "receipt mismatch; inspect before replacing")
    print(json.dumps({"checks": CHECKS, "receipt": str(destination), "status": "PASS"}))


if __name__ == "__main__":
    main()
