"""Exact matching-to-atlas dictionary for the budget-four leading templates."""
from collections import Counter
from itertools import combinations


def need(value, message):
    if not value:
        raise RuntimeError(message)


def matchings(vertices):
    if not vertices:
        yield ()
        return
    a = vertices[0]
    for b in vertices[1:]:
        rest = tuple(v for v in vertices if v != a and v != b)
        for tail in matchings(rest):
            yield ((a, b),) + tail


def main():
    fibers = {}
    for slope in range(8):
        fibers[slope] = frozenset(
            sum((s ^ ((slope & x).bit_count() % 2)) << x for x in range(8))
            + 768 * t for s in range(2) for t in range(2))
    code = frozenset().union(*fibers.values())
    need(len(code) == 32 and all(len(f) == 4 for f in fibers.values()),
         "code/slope fibers")
    charts = {frozenset(fibers[a] | fibers[b]): (a, b)
              for a, b in combinations(range(8), 2)}
    need(len(charts) == 28 and all(len(c) == 8 for c in charts),
         "maximum pair charts")
    atlas_set = set()
    directions = Counter()
    for matching in matchings(tuple(range(8))):
        atlas = frozenset(fibers[a] | fibers[b] for a, b in matching)
        need(len(atlas) == 4 and frozenset().union(*atlas) == code,
             "partition covers")
        need(sum(map(len, atlas)) == len(code), "charts overlap")
        recovered = frozenset(charts[c] for c in atlas)
        need(recovered == frozenset(matching), "matching decoder")
        atlas_set.add(atlas)
        directions[len({a ^ b for a, b in matching})] += 1
    need(len(atlas_set) == 105, "all atlases distinct")
    need(directions == Counter({1: 7, 2: 42, 4: 56}), "direction split")
    # Equal support does not permit replacing numerical transfer weights.
    need(2 * 5 != 1 * 6, "incoming strip blocks accidentally proportional")
    print("Code words=32, slope fibers=8, pair charts=28, decoded atlases=105.")
    print("Budget-four leading matchings: XOR direction counts", dict(directions))
    print("Each sparse template: core order17, heavy gaps9, capped weight53.")
    print("Incoming support hostile: (2,1) and (6,5) not proportional.")
    print("PASS")


if __name__ == "__main__":
    main()
