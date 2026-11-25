"""Exact rational and integral lattice certificates after deleting forbidden triples at (1,0)."""
from collections import Counter
from math import gcd
from pathlib import Path

from joint_characters import determinant, square_minor
from three_stars import A, B, IDENTITY, PIVOT, SLOTS, make_chart, routes


def legal_triples():
    chart = make_chart()
    centers = (PIVOT, A, B)
    branches = [h for h in SLOTS if h != 1]
    allowed = [set(IDENTITY) - {r[0] for r in routes(chart, 1, j)} for j in centers]
    triples = {}
    for h in branches:
        for t in IDENTITY:
            targets = tuple(chart[h, j][t] for j in centers)
            if all(y in accepted for y, accepted in zip(targets, allowed)):
                triples[h, t] = targets
    return branches, allowed, triples


def lattice_minors(branches, allowed, triples):
    labels = [{y: z for z, y in enumerate(sorted(accepted))} for accepted in allowed]
    differences = []
    for h in branches:
        fiber = [targets for (k, t), targets in triples.items() if k == h]
        for point in fiber[1:]:
            row = [0] * 156
            for j, (v, ref) in enumerate(zip(point, fiber[0])):
                a, b = labels[j][v], labels[j][ref]
                if a:
                    row[52*j + a - 1] += 1
                if b:
                    row[52*j + b - 1] -= 1
            differences.append(row)
    assert len(differences) == 2452
    minors = [abs(determinant(square_minor(differences, seed))) for seed in (0, 1)]
    assert minors == [47304122656582883, 16329197407783105] and gcd(*minors) == 1
    return minors


def interior(branches, allowed, triples):
    baseline = 354158094286698
    weights = dict.fromkeys(triples, baseline)
    changed = set()
    for line in Path(__file__).with_name("joint_slot_weights.txt").read_text().splitlines():
        if line.startswith("#"):
            continue
        h, t, extra = map(int, line.split())
        assert extra > 0 and (h, t) not in changed
        weights[h, t] += extra
        changed.add((h, t))
    assert len(changed) == 209
    sums = [Counter() for _ in range(4)]
    for (h, t), weight in weights.items():
        for histogram, label in zip(sums, (h, *triples[h, t])):
            histogram[label] += weight
    for histogram, labels in zip(sums, (branches, *allowed)):
        assert set(histogram) == set(labels) and set(histogram.values()) == {100 * baseline}
    return 100 * baseline


if __name__ == "__main__":
    branches, allowed, triples = legal_triples()
    assert len(triples) == 2505
    print("Positive integral cover degree", interior(branches, allowed, triples))
    print("Legal-support maximal minors", lattice_minors(branches, allowed, triples))
