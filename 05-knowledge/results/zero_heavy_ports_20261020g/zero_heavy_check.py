"""Exact anchor-word quotient, one-feedback count, and bounded re-anchoring controls."""
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path


def load(name, filename):
    path = Path(__file__).resolve().parents[1] / "metagraph_growth_20261020" / filename
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E, G = load("anchor_enum", "enumeration.py"), load("anchor_growth", "growth_check.py")


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def transitive_order(rows, vertices):
    order = tuple(sorted(vertices, key=lambda v: -sum(rows[v] >> w & 1 for w in vertices)))
    return order if all(rows[u] >> v & 1 for u, v in combinations(order, 2)) else None


def from_word(a, header, word):
    n = a + len(word)
    rows = list(E.matrix(a, header)) + [0] * len(word)
    for p, letter in enumerate(word, a):
        for anchor in range(a):
            u, v = (anchor, p) if letter >> anchor & 1 else (p, anchor)
            rows[u] |= 1 << v
        for q in range(p + 1, n):
            rows[p] |= 1 << q
    return tuple(rows)


def record(rows, anchors):
    base = transitive_order(rows, tuple(v for v in range(len(rows)) if v not in anchors))
    if base is None:
        return None
    header = sum((rows[anchors[i]] >> anchors[j] & 1) << p
                 for p, (i, j) in enumerate(combinations(range(len(anchors)), 2)))
    word = tuple(sum((rows[v] >> z & 1) << i for i, v in enumerate(anchors)) for z in base)
    return (header, word), base


def representations(rows, a):
    for anchors in permutations(range(len(rows)), a):
        value = record(rows, anchors)
        if value is not None:
            yield anchors, value[0], value[1]


def pair_free_word(a, header, word):
    """Finite-language homogeneous-pair predicate, separate from clone_chains."""
    if any(x == y for x, y in zip(word, word[1:])):
        return False
    anchor_rows = E.matrix(a, header)
    for i, j in combinations(range(a), 2):
        if all((anchor_rows[i] >> z & 1) == (anchor_rows[j] >> z & 1)
               for z in range(a) if z not in (i, j)):
            if all((letter >> i & 1) == (letter >> j & 1) for letter in word):
                return False
    for i in range(a):
        for p, letter in enumerate(word):
            if (all(not (x >> i & 1) for x in word[:p]) and
                    all(x >> i & 1 for x in word[p + 1:]) and
                    all((letter >> j & 1) == 1 - (anchor_rows[i] >> j & 1)
                        for j in range(a) if j != i)):
                return False
    return True


def transport(rows, old, new):
    """Re-encode using only the old header/word, selected positions, and insertion order."""
    (header, word), base = record(rows, old)
    expected, newbase = record(rows, new)
    a = len(old)
    old_index, base_index = {v: i for i, v in enumerate(old)}, {v: i for i, v in enumerate(base)}
    anchor_rows = E.matrix(a, header)

    def arc(u, v):
        if u in old_index and v in old_index:
            return anchor_rows[old_index[u]] >> old_index[v] & 1
        if u in old_index:
            return word[base_index[v]] >> old_index[u] & 1
        if v in old_index:
            return 1 - (word[base_index[u]] >> old_index[v] & 1)
        return int(base_index[u] < base_index[v])

    newheader = sum(arc(new[i], new[j]) << p
                    for p, (i, j) in enumerate(combinations(range(a), 2)))
    output, index, consumed, maxlag = [], 0, 0, 0
    for v in base:
        if v in new:
            consumed += 1
            maxlag = max(maxlag, abs(len(output) - consumed))
            continue
        while index < len(newbase) and newbase[index] in old_index:
            x = newbase[index]
            output.append(sum(arc(z, x) << i for i, z in enumerate(new)))
            index += 1
            maxlag = max(maxlag, abs(len(output) - consumed))
        need(index < len(newbase) and newbase[index] == v, "old base order retained")
        output.append(sum(arc(z, v) << i for i, z in enumerate(new)))
        consumed += 1
        index += 1
        maxlag = max(maxlag, abs(len(output) - consumed))
    while index < len(newbase):
        x = newbase[index]
        need(x in old_index, "only displaced old anchors remain")
        output.append(sum(arc(z, x) << i for i, z in enumerate(new)))
        index += 1
        maxlag = max(maxlag, abs(len(output) - consumed))
    need((newheader, tuple(output)) == expected, "metadata transport reconstructs full new record")
    need(maxlag <= a and len(output) == len(word), "bounded-delay equal-length re-anchoring")
    return maxlag


def classes_through(maximum):
    layers = {1: {(0,)}}
    for n in range(2, maximum + 1):
        codes = {E.canonical(E.extension(rows, mask))
                 for rows in layers[n - 1] for mask in range(1 << (n - 1))}
        layers[n] = {E.decode(n, code) for code in codes}
    need(tuple(len(layers[n]) for n in range(1, maximum + 1)) ==
         (1, 1, 2, 4, 12, 56, 456)[:maximum], "independent class universe")
    return layers


def joins_of_triangles(m):
    return tuple(sum((int(i // 3 < j // 3) if i // 3 != j // 3 else
                      int((j - i) % 3 == 1)) << j for j in range(3 * m) if i != j)
                 for i in range(3 * m))


def separating(rows, anchors):
    signatures = [tuple(rows[a] >> v & 1 for a in anchors)
                  for v in range(len(rows)) if v not in anchors]
    return len(signatures) == len(set(signatures))


def main():
    layers = classes_through(7)
    old_counts, two_counts, comparisons, maxlag = [], [], 0, 0
    root_hist = Counter()
    for n, rows_set in sorted(layers.items()):
        one, two = 0, 0
        for rows in rows_set:
            one_reps = list(representations(rows, 1))
            if one_reps:
                one += 1
                roots, rooted_classes = len(one_reps), len({rep for _, rep, _ in one_reps})
                if transitive_order(rows, tuple(range(n))) is not None:
                    need(roots == rooted_classes == n, "transitive root correction")
                else:
                    need((roots, rooted_classes) in ((1, 1), (2, 2), (3, 1)), "nontransitive root classification")
                    root_hist[n, roots, rooted_classes] += 1
            for a in (1, 2):
                if a > n:
                    continue
                reps = one_reps if a == 1 else list(representations(rows, a))
                if a == 2 and reps:
                    two += 1
                for anchors, encoded, _ in reps:
                    need(pair_free_word(a, *encoded) == (len(G.clone_chains(rows)) == n),
                         "full port-word clone predicate")
                if reps:
                    old = reps[0][0]
                    for new, _, _ in reps:
                        maxlag = max(maxlag, transport(rows, old, new))
                        comparisons += 1
        old_counts.append(one)
        if n >= 2:
            need(one == 2 ** (n - 1) - (n - 1) - comb(n - 2, 2), "one-feedback class formula")
            two_counts.append(two)
    print("all-class universe n=1..7=", tuple(len(layers[n]) for n in layers), flush=True)
    print("feedback<=1 class counts=", tuple(old_counts), flush=True)
    print("feedback<=2 class counts n=2..7=", tuple(two_counts), flush=True)
    print("eligible-root histogram (n,vertices,rooted_classes)=", dict(sorted(root_hist.items())), flush=True)

    binary_counts, core_counts, word_tests = [], [], 0
    for n in range(1, 11):
        ordinary, cores = set(), set()
        for word in product(range(2), repeat=n - 1):
            rows = from_word(1, 0, word)
            code = E.canonical(rows)
            ordinary.add(code)
            pf = pair_free_word(1, 0, word)
            need(pf == (len(G.clone_chains(rows)) == n), "binary clone classifier")
            if pf:
                cores.add(code)
            word_tests += 1
        expected = 1 if n == 1 else 2 ** (n - 1) - (n - 1) - comb(n - 2, 2)
        expected_core = (1, 0, 1)[n - 1] if n <= 3 else 2
        need(len(ordinary) == expected and len(cores) == expected_core, "binary quotient/formal counts")
        binary_counts.append(len(ordinary))
        core_counts.append(len(cores))
    print("independent binary-word quotient n=1..10=", tuple(binary_counts), flush=True)
    print("pair-free feedback<=1 counts=", tuple(core_counts), "words_checked=", word_tests, flush=True)

    owners = []
    for n in range(2, 8):
        accepted = 0
        for header in range(2):
            for word in product(range(4), repeat=n - 2):
                rows = from_word(2, header, word)
                smallest = min(rep for _, rep, _ in representations(rows, 2))
                accepted += (header, word) == smallest
        need(accepted == two_counts[n - 2], "two-anchor canonical-owner language equals class census")
        owners.append(accepted)
    print("two-anchor least-word owners=", tuple(owners), flush=True)
    print("exact re-anchoring controls=", comparisons, "maximum observed lag=", maxlag, flush=True)

    separating_minima = []
    for m in range(1, 5):
        rows = joins_of_triangles(m)
        need(len(G.clone_chains(rows)) == len(rows), "anchor hostile is pair free")
        minimum = next(a for a in range(len(rows) + 1)
                       if any(separating(rows, anchors) for anchors in combinations(range(len(rows)), a)))
        need(minimum == m + m // 2, "sharp linear signature-anchor requirement")
        separating_minima.append(minimum)
    print("ordered-C3-join separating-anchor minima m=1..4=", tuple(separating_minima), flush=True)
    print("PASS exact one-feedback quotient, full clone tests, bounded port transport, and anchor hostile")


if __name__ == "__main__":
    main()
