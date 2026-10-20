#!/usr/bin/env python3
"""Independent bounded exact controls for the all-order equivariance obstruction.

No repository modules, external packages, random choices, or nauty are used.
Universe: all labelled simple graphs/tournaments on {0,...,n-1}; switching
classes use their unique representative with vertex 0 a source. A tournament
mask bit is 1 precisely when its increasing-label reference arc is reversed.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def valuation(a, p):
    s = 0
    while a and a % p == 0:
        a //= p
        s += 1
    return s


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    r = [0] * len(a)
    for i, j in enumerate(a):
        r[j] = i
    return tuple(r)


def cycle(n, entries):
    p = list(range(n))
    for i, j in zip(entries, entries[1:] + entries[:1]):
        p[i] = j
    return tuple(p)


def double_transposition(p):
    return sum(i != p[i] for i in range(len(p))) == 4 and all(p[p[i]] == i for i in range(len(p)))


def overlapping_controls(n):
    cycles = []
    for support in itertools.combinations(range(n), 3):
        cycles += [(set(support), cycle(n, support)), (set(support), cycle(n, support[::-1]))]
    counts = {"overlap_one": 0, "overlap_two": 0}
    for (sa, a), (sb, b) in itertools.combinations(cycles, 2):
        overlap = len(sa & sb)
        if overlap not in (1, 2):
            continue
        counts["overlap_one" if overlap == 1 else "overlap_two"] += 1
        c = b if overlap == 2 else compose(compose(b, a), inverse(b))
        assert len(sa & {i for i in range(n) if c[i] != i}) == 2
        assert any(double_transposition(compose(a, d)) for d in (c, inverse(c)))
    return counts


class Universe:
    def __init__(self, n):
        self.n = n
        self.edges = list(itertools.combinations(range(n), 2))
        self.m = len(self.edges)
        self.size = 1 << self.m
        self.perms = list(itertools.permutations(range(n)))
        edge_index = {e: i for i, e in enumerate(self.edges)}
        self.maps = []
        self.domain_flips = []
        self.image_flips = []
        for p in self.perms:
            bits, dom, img = [], 0, 0
            for i, (a, b) in enumerate(self.edges):
                j = edge_index[tuple(sorted((p[a], p[b])))]
                bits.append(1 << j)
                if p[a] > p[b]:
                    dom |= 1 << i
                    img |= 1 << j
            self.maps.append(bits)
            self.domain_flips.append(dom)
            self.image_flips.append(img)
        self.normalized = [self.normalize(mask) for mask in range(self.size)]

    def graph_action(self, mask, pi):
        ans = 0
        bits = self.maps[pi]
        while mask:
            low = mask & -mask
            ans |= bits[low.bit_length() - 1]
            mask ^= low
        return ans

    def tournament_action(self, mask, pi):
        return self.graph_action(mask, pi) ^ self.image_flips[pi]

    def normalize(self, mask):
        if self.n == 0:
            return mask
        switched = {b for i, (a, b) in enumerate(self.edges) if a == 0 and mask >> i & 1}
        cut = sum(1 << i for i, (a, b) in enumerate(self.edges) if (a in switched) != (b in switched))
        return mask ^ cut

    def class_action(self, mask, pi):
        return self.normalized[self.tournament_action(mask, pi)]

    def catalogue(self, kind):
        action = {"graph": self.graph_action, "tournament": self.tournament_action, "class": self.class_action}[kind]
        remaining = set(self.normalized if kind == "class" else range(self.size))
        out = []
        while remaining:
            rep = min(remaining)
            images = [action(rep, pi) for pi in range(len(self.perms))]
            orbit = set(images)
            stabilizer = [pi for pi, image in enumerate(images) if image == rep]
            assert len(orbit) * len(stabilizer) == math.factorial(self.n)
            remaining.difference_update(orbit)
            out.append({"representative": rep, "orbit": orbit, "stabilizer": stabilizer, "order": len(stabilizer)})
        return out

    def even(self, graph):
        mask = graph["representative"]
        return all((mask & self.domain_flips[pi]).bit_count() % 2 == 0 for pi in graph["stabilizer"])

    def euler(self, graph):
        mask = graph["representative"]
        degrees = [0] * self.n
        for i, (a, b) in enumerate(self.edges):
            if mask >> i & 1:
                degrees[a] += 1
                degrees[b] += 1
        return all(d % 2 == 0 for d in degrees)

    def matching(self, source, target):
        adjacency = []
        for a in source:
            row = []
            for j, b in enumerate(target):
                if any(all(self.graph_action(copy, pi) == copy for pi in a["stabilizer"]) for copy in b["orbit"]):
                    row.append(j)
            adjacency.append(row)
        assigned = {}
        def augment(i, visited):
            for j in adjacency[i]:
                if j in visited:
                    continue
                visited.add(j)
                if j not in assigned or augment(assigned[j], visited):
                    assigned[j] = i
                    return True
            return False
        for i in range(len(source)):
            augment(i, set())
        return len(assigned), adjacency

    def controls(self):
        allgraphs = self.catalogue("graph")
        even = [x for x in allgraphs if self.even(x)]
        eveneuler = [x for x in even if self.euler(x)]
        tournaments = self.catalogue("tournament")
        classes = self.catalogue("class")
        result = {"n": self.n, "graph_types": len(allgraphs)}
        for name, source, target in [("tournaments", tournaments, even), ("switching_classes", classes, eveneuler)]:
            assert len(source) == len(target)
            maxmatch, adjacency = self.matching(source, target)
            totals = {}
            sylow_counts = {}
            for p in (3, 5, 7):
                a = sum(valuation(x["order"], p) for x in source)
                b = sum(valuation(x["order"], p) for x in target)
                assert a == b, (self.n, name, p, a, b)
                totals[str(p)] = a
                ambient = valuation(math.factorial(self.n), p)
                source_count = sum(valuation(x["order"], p) == ambient for x in source)
                target_count = sum(valuation(x["order"], p) == ambient for x in target)
                assert source_count == target_count
                sylow_counts[str(p)] = source_count
            if self.n >= 5:
                for x in source:
                    assert not any(double_transposition(self.perms[pi]) for pi in x["stabilizer"])
            for x in source:
                for pi in x["stabilizer"]:
                    perm = self.perms[pi]
                    involution = perm != tuple(range(self.n)) and compose(perm, perm) == tuple(range(self.n))
                    if name == "tournaments":
                        assert not involution
                    elif involution:
                        assert all(perm[v] != v for v in range(self.n))
            joint_source = sum(all(valuation(x["order"], p) == valuation(math.factorial(self.n), p) for p in (3, 5)) for x in source)
            joint_target = sum(all(valuation(x["order"], p) == valuation(math.factorial(self.n), p) for p in (3, 5)) for x in target)
            if self.n >= 5:
                assert joint_source == 0 and joint_target >= 1
                assert len(source) - maxmatch >= joint_target
            result[name] = {
                "full_ambient_sylow_counts_both_sides": sylow_counts,
                "full_ambient_sylow_3_5_joint_source_count": joint_source,
                "full_ambient_sylow_3_5_joint_target_count": joint_target,
                "source_types": len(source), "target_types": len(target), "maximum_matching": maxmatch,
                "deficit": len(source) - maxmatch, "odd_valuation_totals_both_sides": totals,
                "source_orders": [x["order"] for x in source], "target_orders": [x["order"] for x in target],
                "source_canonical_masks": [x["representative"] for x in source],
                "target_canonical_masks": [x["representative"] for x in target], "adjacency": adjacency,
            }
            assert maxmatch == len(source) if self.n <= 4 else maxmatch < len(source)
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=6, choices=range(0, 7))
    parser.add_argument("--out", default="exact_controls.json")
    args = parser.parse_args()
    report = {"scope": "Complete labelled universes through max_n; all 3-cycle pairs for n=5,6,7; arithmetic gap through 10000", "max_n": args.max_n}
    report["orders"] = []
    for n in range(args.max_n + 1):
        row = Universe(n).controls()
        report["orders"].append(row)
        print(json.dumps({"n": n, "tournament_types": row["tournaments"]["source_types"], "tournament_deficit": row["tournaments"]["deficit"], "class_types": row["switching_classes"]["source_types"], "class_deficit": row["switching_classes"]["deficit"]}, sort_keys=True))
    report["overlap_controls"] = {str(n): overlapping_controls(n) for n in (5, 6, 7)}
    def factorial_valuation(n, p):
        total = 0
        while n:
            n //= p
            total += n
        return total
    for n in range(5, 10001):
        k = n // 3
        assert factorial_valuation(k, 5) < factorial_valuation(n, 5)
    report["arithmetic_gap"] = {"first": 5, "last": 10000, "checks": 9996, "passed": True}
    report["all_checks_passed"] = True
    Path(args.out).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("ALL CHECKS PASSED")

if __name__ == "__main__":
    main()
