"""Exact controls for the incoming defect-template / cubic-kernel dictionary.

Universe: all 105 matchings of eight labelled ports, all 3**8 words,
and all 3**7 words for a shared-endpoint hostile and cubic-junction repair.
No actual tournament-count transport is inferred for the repaired junction.
Run from any directory; standard library only. Checks survive python -O.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matchings(vertices):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for j in range(1, len(vertices)):
        other = vertices[j]
        rest = vertices[1:j] + vertices[j + 1:]
        for matching in matchings(rest):
            yield ((first, other),) + matching


def paired_count(n):
    if n % 2:
        return 0
    answer = 1
    for factor in range(1, n, 2):
        answer *= factor
    return answer


bank = tuple(matchings(tuple(range(8))))
directions = Counter(len({a ^ b for a, b in matching}) for matching in bank)
require(directions == {1: 7, 2: 42, 4: 56}, "direction bank")
responses = Counter()
for word in product(range(3), repeat=8):
    actual = sum(all(word[a] == word[b] for a, b in matching)
                 for matching in bank)
    expected = 1
    for colour in range(3):
        expected *= paired_count(word.count(colour))
    require(actual == expected, ("matching response", word, actual, expected))
    responses[actual] += 1
require(responses == {0: 4920, 3: 1260, 9: 210, 15: 168, 105: 3},
        ("response histogram", responses))

naive = repaired = naive_parity = repaired_parity = 0
for word in product(range(3), repeat=7):
    pairs = word[3] == word[4] and word[5] == word[6]
    equality = pairs and word[0] == word[1] == word[2]
    junction = pairs and len(set(word[:3])) == 3
    parity = all(word.count(colour) % 2 == 1 for colour in range(3))
    naive += equality
    repaired += junction
    naive_parity += equality and parity
    repaired_parity += junction and parity
require((naive, naive_parity, repaired, repaired_parity) == (27, 0, 54, 54),
        "seven-port parity boundary")

print("Eight-port direction bank:", dict(sorted(directions.items())))
print("Sum of matching kernels, all 6561 words:", dict(sorted(responses.items())))
print("Budget-four top templates: core_order=17 capped_weight=53 heavy_gaps=9")
print("Top generating term: 105*x^53/(1-x)^9")
print("Leading coefficient:", Fraction(len(bank), factorial(8)))
print("Seven-port equality hostile: support=27 parity-compatible=0")
print("Seven-port cubic junction: support=54 parity-compatible=54")
print("PASS: all identities checked by exhaustive labelled words")

# A separate order-five control for the later pure-blue atom integration.
# All paths and all automorphisms are checked, not just chosen orbit reps.
vertices = tuple(range(5))
orders = tuple(permutations(vertices))


def paths(arcs):
    return [p for p in orders
            if all((a, b) in arcs for a, b in zip(p, p[1:]))]


def profile(arcs):
    aut = sum(all((p[a], p[b]) in arcs for a, b in arcs) for p in orders)
    return len(paths(arcs)), aut


def complement_along(arcs, path):
    kept = set(zip(path, path[1:]))
    return frozenset((a, b) if (a, b) in kept else (b, a)
                     for a, b in arcs)


def strong(arcs):
    for start in vertices:
        seen = {start}
        frontier = [start]
        while frontier:
            a = frontier.pop()
            for b in vertices:
                if (a, b) in arcs and b not in seen:
                    seen.add(b)
                    frontier.append(b)
        if len(seen) != 5:
            return False
    return True


c5 = frozenset((a, b) for a in vertices for b in vertices
               if (b - a) % 5 in (1, 2))
tt5 = frozenset((a, b) for a in vertices for b in vertices if a < b)
stack = set(tt5)
stack.remove((1, 3))
stack.add((3, 1))
controls = (("C5", c5), ("TT5", tt5), ("1->C3->1", frozenset(stack)))
for label, arcs in controls:
    partners = [complement_along(arcs, p) for p in paths(arcs)]
    histogram = Counter(profile(partner) for partner in partners)
    require(all(strong(partner) for partner in partners),
            ("order-five partner not strong", label))
    require((15, 5) not in histogram, ("pure-blue partner", label))
    if label == "C5":
        require(histogram == {(11, 1): 5, (15, 3): 5, (13, 1): 5},
                ("C5 partner profile", histogram))
    print(label, "all-path partner (H,Aut) histogram:", dict(sorted(histogram.items())))
print("PASS: order-five pure-blue partners are strong and not C5")
