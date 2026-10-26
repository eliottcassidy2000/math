"""Independent audit of the parent's planar whole-frame transport theorem.

No producer import. Barycentric Gaussian elimination and multi-loop networks
provide controls beyond the producer's cofactor and triangle calculations.
"""
from fractions import Fraction as F
from itertools import permutations, product, combinations
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
PROOF = HERE.parent / "metric_transport/RESULT.md"
PRODUCER = HERE.parent / "metric_transport/check.py"


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


PERMS = tuple(permutations(range(3)))
MAPS = tuple(tuple(a ^ sum(((x >> i) & 1) << p[i] for i in range(3))
                   for x in range(8)) for p in PERMS for a in range(8))


def sign(frame, x):
    return tuple(tuple((-1 if x & (1 << i) else 1)*v for v in u)
                 for i, u in enumerate(frame))


def barycentric_origin(vertices):
    rows = [[F(v[j]) for v in vertices] + [F(0)] for j in range(2)]
    rows.append([F(1)] * 4)
    for c in range(3):
        pivots = [r for r in range(c, 3) if rows[r][c]]
        if not pivots:
            return None
        r = pivots[0]
        rows[r], rows[c] = rows[c], rows[r]
        scale = rows[c][c]
        rows[c] = [z/scale for z in rows[c]]
        for r in range(3):
            if r != c:
                scale = rows[r][c]
                rows[r] = [z-scale*w for z, w in zip(rows[r], rows[c])]
    return tuple(rows[i][-1] for i in range(3))


def balanced(frame):
    answer = set()
    for x in range(8):
        coeff = barycentric_origin(sign(frame, x))
        if coeff is not None and all(z > 0 for z in coeff):
            answer.add(x)
    need(len(answer) == 2 and {x ^ 7 for x in answer} == answer, "unique balanced antipodal pair")
    return answer


def points(frame):
    return tuple(tuple(s*z for z in u) for u in frame for s in (1, -1))


def sqdistance(a, b, metric):
    return sum(g*(x-y)**2 for x, y, g in zip(a, b, metric))


def isometry(source, target, map8, metric):
    # Recover the six-endpoint permutation from the affine cube map by
    # differentiating its three coordinate flips, not by calling a Gram test.
    a = map8[0]
    image = []
    for i in range(3):
        direction = map8[1 << i] ^ a
        need(direction in (1, 2, 4), "signed permutation derivative")
        j = direction.bit_length()-1
        image.extend((2*j + ((a >> j) & 1), 2*j + (((a >> j) & 1) ^ 1)))
    s, t = points(source), points(target)
    return all(sqdistance(s[i], s[j], metric) == sqdistance(t[image[i]], t[image[j]], metric)
               for i, j in combinations(range(6), 2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "metric-independent-normal.json")
    args = parser.parse_args()
    flower = ((F(1), F(0)), (F(-1, 2), F(1, 2)), (F(-1, 2), F(-1, 2)))
    frames = [(flower, (1, 3)),
              (((1, 0), (F(3, 5), F(4, 5)), (F(-5, 13), F(12, 13))), (1, 1)),
              (((1, 0), (2, 3), (-1, 4)), (1, 1))]  # unequal lengths
    transport_tests = 0
    for frame, metric in frames:
        marked = [sign(frame, x) for x in range(8)]
        balanced_sets = [balanced(f) for f in marked]
        for a, b in product(range(8), repeat=2):
            for g in MAPS:
                if isometry(marked[a], marked[b], g, metric):
                    need({g[x] for x in balanced_sets[a]} == balanced_sets[b], "barycentric metric covariance")
                transport_tests += 1
    # Network has vertices0,1, two parallel0->1 edges and one loop at0.
    # The two independent cycles test a common section, not separate loops.
    markings = (0, 3, 5, 6)
    f = {a: sign(flower, a) for a in markings}
    seams = {(a, b): [g for g in MAPS if isometry(f[a], f[b], g, (1, 3))]
             for a, b in product(markings, repeat=2)}
    histogram = Counter()
    for a, b in product(markings, repeat=2):
        for edge1, edge2, loop in product(seams[a, b], seams[a, b], seams[a, a]):
            sections = [(x, y) for x, y in product(range(8), repeat=2)
                        if edge1[x] == y == edge2[x] and loop[x] == x]
            parity_sections = [(i, j) for i, j in product(range(2), repeat=2)
                               if j == (i ^ (edge1[0].bit_count() & 1))
                               and j == (i ^ (edge2[0].bit_count() & 1))
                               and (loop[0].bit_count() & 1) == 0]
            need(bool(sections) == bool(parity_sections), "multi-cycle network iff")
            if sections:
                need({x for x, y in sections} >= balanced(f[a]), "shared balanced section on both cycles")
            histogram[len(sections)] += 1
    cycle = (1, 2, 0)
    g0 = tuple(sum(((x >> i) & 1) << cycle[i] for i in range(3)) for x in range(8))
    g5 = tuple(x ^ 5 for x in g0)
    need(not any(g0[x] == x == g5[x] for x in range(8)), "abstract two-loop hostile")
    axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    need(isometry(axes, axes, g0, (1, 1, 1)) and isometry(axes, axes, g5, (1, 1, 1)),
         "three-dimensional whole-frame hostile")
    need(not (g0[0].bit_count() & 1) and not (g5[0].bit_count() & 1), "3D parity passes")
    report = {"status": "PASS", "metric_transport_tests": transport_tests,
              "two_vertex_two_cycle_networks": sum(histogram.values()),
              "full_section_histogram": dict(sorted(histogram.items())),
              "three_dimensional_counterexample": "both incompatible parity-zero loops are unit-frame isometries of the coordinate octahedron",
              "producer_sha256": hashlib.sha256(PRODUCER.read_bytes()).hexdigest(),
              "proof_sha256": hashlib.sha256(PROOF.read_bytes()).hexdigest(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
