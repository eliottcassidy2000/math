"""Independent direct-word / complete-curl audit of the information packet.

No imports from either producer or parent auditor. Always-active checks.
"""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import json

CHECKS = 0


def need(p, message):
    global CHECKS
    CHECKS += 1
    if not p:
        raise RuntimeError(message)


def tournament(n, mask):
    a = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = (-1 if mask >> k & 1 else 1)
        a[j][i] = -a[i][j]
    return a


def encode(a):
    return sum((a[i][j] < 0) << k for k, (i, j) in enumerate(combinations(range(len(a)), 2)))


def all_curl(a):
    return tuple(a[i][j] + a[j][k] + a[k][i] for i, j, k in combinations(range(len(a)), 3))


def permute(a, order):
    return [[a[i][j] for j in order] for i in order]


def direct_coefficients():
    records = []
    for n in range(3, 9):
        edges = list(combinations(range(n), 2))
        index = {e:i for i, e in enumerate(edges)}
        coeff = defaultdict(int)
        for word in permutations(range(n)):
            signed = [(index[tuple(sorted((u, v)))], 1 if u < v else -1)
                      for u, v in zip(word, word[1:])]
            for (e, s), (f, t) in combinations(signed, 2):
                coeff[tuple(sorted((e, f)))] += s*t
        nonzero = 0
        for e, f in combinations(range(len(edges)), 2):
            meet = set(edges[e]) & set(edges[f])
            expected = 0
            if meet:
                j = next(iter(meet))
                i = next(v for v in edges[e] if v != j)
                k = next(v for v in edges[f] if v != j)
                sign = (1 if i < j else -1) * (1 if j < k else -1)
                expected = sign*2*factorial(n-2)
                nonzero += 1
            need(coeff[e, f] == expected, "direct word coefficient including orientation sign")
        records.append({"n": n, "words": factorial(n), "adjacent_pair_coefficients": nonzero,
                        "normalized_magnitude": str(Fraction(2*factorial(n-2), 2**(n-1)))})
    return records


def direct_response():
    records = []
    one_bit_hostile = None
    for n in range(1, 6):
        m = n*(n-1)//2
        size, full = 2**m, 2**m-1
        orders = list(permutations(range(n)))
        matrices = [tournament(n, z) for z in range(size)]
        valid_words = [[w for w in orders if all(a[u][v] == 1 for u, v in zip(w, w[1:]))]
                       for a in matrices]
        hs = [len(words) for words in valid_words]
        rs = [tuple(hs[z ^ f] for f in range(size)) for z in range(size)]
        curls = [all_curl(a) for a in matrices]
        response_classes, curl_classes, start_classes = defaultdict(list), defaultdict(list), defaultdict(set)
        joint = set()
        for z in range(size):
            response_classes[rs[z]].append(z)
            curl_classes[curls[z]].append(z)
            joint.add((curls[z], rs[z]))
            for w in valid_words[z]:
                start_classes[curls[z], w[0]].add(z)
        need(all(set(group) == {group[0], group[0] ^ full} for group in response_classes.values()),
             "literal full-response converse fibres")
        if n >= 3:
            need(len(joint) == size, "full-curl/full-response aligned injectivity")
        need(all(len(group) == 1 for group in start_classes.values()), "valid path-start selector")
        if n >= 3 and one_bit_hostile is None:
            for group in curl_classes.values():
                sides = defaultdict(list)
                for z in group:
                    sides[z & 1].append(z)
                for sub in sides.values():
                    if len(sub) > 1:
                        one_bit_hostile = {"n": n, "masks": sub[:2], "known_arc": [0, 1],
                                           "common_arc_bit": sub[0] & 1, "full_curl": curls[sub[0]]}
                        break
                if one_bit_hostile is not None:
                    break
        # Independently check the entire SCC-rotation description using
        # transitive closure, not the curl decoder or a score reconstruction.
        for z, a in enumerate(matrices):
            reach = [[i == j or a[i][j] == 1 for j in range(n)] for i in range(n)]
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
            components = []
            unseen = set(range(n))
            while unseen:
                i = min(unseen)
                block = {j for j in unseen if reach[i][j] and reach[j][i]}
                components.append(block)
                unseen -= block
            component_snapshot = tuple(components)
            components.sort(key=lambda block: -sum(a[next(iter(block))][next(iter(other))] == 1
                                                    for other in component_snapshot if other is not block))
            rotations = set()
            for shift in range(len(components)):
                blocks = components[shift:] + components[:shift]
                b = [row[:] for row in a]
                for l, r in combinations(range(len(blocks)), 2):
                    for u in blocks[l]:
                        for v in blocks[r]:
                            b[u][v], b[v][u] = 1, -1
                rotations.add(encode(b))
            need(rotations == set(curl_classes[curls[z]]), "all full-curl fibres are exactly SCC rotations")
        record = {"n": n, "tournaments": size, "literal_response_fibres": len(response_classes),
                  "full_curl_fibres": len(curl_classes), "aligned_fibres": len(joint),
                  "valid_start_keys": len(start_classes)}
        if n <= 4:
            ordinary, separate, aligned = set(), set(), set()
            for z, a in enumerate(matrices):
                moved = [encode(permute(a, p)) for p in orders]
                ordinary.add(min(moved))
                # Canonicalize ACTUAL full response arrays, not a proxy
                # ordinary/converse canonical code based on the theorem.
                separate.add((min(curls[t] for t in moved), min(rs[t] for t in moved)))
                aligned.add(min((curls[t], rs[t]) for t in moved))
            need(len(ordinary) == len(aligned), "aligned isomorphism classification")
            if n <= 3:
                need(len(ordinary) == len(separate), "minimal alignment hostile boundary")
            else:
                need(len(ordinary) == 4 and len(separate) == 3, "order-four separate quotient loss")
            record.update({"tournament_isotypes": len(ordinary), "separate_isotypes": len(separate),
                           "aligned_isotypes": len(aligned)})
        records.append(record)
    need(one_bit_hostile is not None and one_bit_hostile["n"] == 3, "curl plus one bit is insufficient already at3")
    return records, one_bit_hostile


def main():
    coefficients = direct_coefficients()
    responses, hostile = direct_response()
    output = {"status": "PASS", "coefficient_controls": coefficients,
              "literal_response_controls": responses, "curl_plus_one_arc_hostile": hostile,
              "checks": CHECKS, "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
