#!/usr/bin/env python3
"""Independent brute controls and exhaustive declared finite universes."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import importlib.util
import json
import sys

from enumerate import (need, rows, encode, relabel, canonical, enumerate_classes,
                       hamilton_count, chart, unchart, clone_core, inflate,
                       weighted_key, vertex_orbits, weighted_core_allocations,
                       port_coordinates, marked_key, hamilton_paths,
                       chart_complement, chart_reverse_converse, reframe_ports)


def brute(n, code):
    """Independent matrix/permutation path; no refinement or canonical call."""
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = (code >> k) & 1
        a[j][i] = 1 - a[i][j]
    best, autos = None, []
    for p in permutations(range(n)):
        value = 0
        shift = 0
        for i in range(n):
            for j in range(i + 1, n):
                value += a[p[i]][p[j]] << shift
                shift += 1
        best = value if best is None else min(best, value)
        if all(a[i][j] == a[p[i]][p[j]] for i in range(n) for j in range(n)):
            autos.append(p)
    return best, set(autos)


def load_reverse():
    path = Path(__file__).resolve().parents[2] / "lca_ports_reverse_20261025/counting/reverse.py"
    spec = importlib.util.spec_from_file_location("independent_reverse_ordinary", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run():
    levels, growth = enumerate_classes(8)
    counts = [len(level) for level in levels]
    need(counts == [1, 1, 1, 2, 4, 12, 56, 456, 6880], "known A000568 counts")
    brute_cases = relabel_cases = 0
    for n in range(6):
        can_to_brute, brute_to_can = {}, {}
        for code in range(1 << (n * (n - 1) // 2)):
            cf = canonical(n, code)
            bc, ba = brute(n, code)
            need(set(cf.automorphisms()) == ba, "independent full automorphism group")
            need(cf.code in levels[n], "labeled census missing class")
            need(can_to_brute.setdefault(cf.code, bc) == bc, "canonical false merge")
            need(brute_to_can.setdefault(bc, cf.code) == cf.code, "canonical false split")
            need(relabel(n, code, cf.orders[0]) == cf.code, "canonical transport")
            need(canonical(n, cf.code).code == cf.code, "idempotence")
            brute_cases += 1
        need(len(can_to_brute) == len(levels[n]), "complete small class census")
        for rec in levels[n].values():
            for p in permutations(range(n)):
                cp = relabel(n, rec.code, p)
                need(canonical(n, cp).code == rec.code, "all-label equivariance")
                relabel_cases += 1

    chart_cases = 0
    chart_receipts = []
    for n in range(8):
        width = max(0, (n - 1) * (n - 2) // 2) if n else 0
        fibers = Counter()
        for word in range(1 << width):
            code = unchart(n, word)
            need(chart(n, code, tuple(range(n)))[0] == word, "literal chart roundtrip")
            fibers[canonical(n, code).code] += 1
            chart_cases += 1
        need(set(fibers) == set(levels[n]), "chart family is not onto the class set")
        for code, size in fibers.items():
            rec = levels[n][code]
            h = hamilton_count(n, code)
            need(h == size * len(rec.automorphisms), "chart fiber H/Aut")
        chart_receipts.append({"n": n, "words": 1 << width, "classes": len(fibers)})

    class_cases = clone_cases = allocations = 0
    class_receipts = []
    clone_counts = []
    for n, level in enumerate(levels):
        labelled_mass = Fraction(0)
        chart_mass = 0
        root_mass = 0
        weighted_keys = set()
        for code, rec in level.items():
            group = rec.automorphisms
            need(len(group) % 2 == 1, "even tournament automorphism order")
            a = rows(n, code)
            for g in group:
                need(encode(a, g) == code, "transported automorphism is not literal")
            h = hamilton_count(n, code)
            need(h % len(group) == 0, "nonintegral chart fiber")
            labelled_mass += Fraction(factorial(n), len(group))
            chart_mass += h // len(group)
            root_mass += len(vertex_orbits(n, group))
            word, path = chart(n, code)
            need(unchart(n, word) == relabel(n, code, path), "class chart transport")
            blocks, q, weights = clone_core(n, code)
            flat = tuple(v for block in blocks for v in block)
            need(inflate(q, weights) == relabel(n, code, flat), "literal clone reconstruction")
            need(len(clone_core(len(weights), q)[0]) == len(weights), "core has surviving homogeneous pair")
            key = weighted_key(q, weights)
            need(key not in weighted_keys, "weighted core false merger")
            weighted_keys.add(key)
            # Unique transitive internal orders: the weighted core has the
            # same automorphism order, without symmetric internal factors.
            core_group = canonical(len(weights), q, weights).automorphisms()
            need(len(core_group) == len(group), "weighted-core symmetry changed")
            class_cases += 1
            clone_cases += 1
        need(labelled_mass == 1 << (n * (n - 1) // 2), "orbit-stabilizer labeled mass")
        width = max(0, (n - 1) * (n - 2) // 2) if n else 0
        need(chart_mass == 1 << width, "Hamiltonian-chart mass")
        if n:
            need(root_mass == growth[n - 1]["rooted_extensions"], "rooted extension/deletion fiber identity")
        reconstructed = set()
        core_layers = Counter()
        for qn, q, weights in weighted_core_allocations(levels, n):
            code = canonical(n, inflate(q, weights)).code
            need(code not in reconstructed, "weighted-core allocation duplicated a class")
            reconstructed.add(code)
            core_layers[qn] += 1
            allocations += 1
        need(reconstructed == set(level), "weighted-core route not complete")
        clone_counts.append({"n": n, "by_core_order": dict(sorted(core_layers.items()))})
        class_receipts.append({"n": n, "classes": len(level), "rooted_classes": root_mass,
                               "chart_words": chart_mass, "labelled_tournaments": int(labelled_mass)})

    reverse = load_reverse()
    exact_counts = []
    for n in range(9):
        value, _ = reverse.OrdinaryDAG(n).exact()
        need(value == counts[n], "independent reverse counting recurrence")
        exact_counts.append(value)

    path_cases = marked_inventory_cases = 0
    frame_dependence = None
    for n in range(7):
        for code, rec in levels[n].items():
            paths = tuple(hamilton_paths(n, code))
            need(len(paths) == hamilton_count(n, code), "literal path inventory vs subset DP")
            word_multiplicities = Counter()
            complement_owners = set()
            for path in paths:
                word, _ = chart(n, code, path)
                word_multiplicities[word] += 1
                complement = chart_complement(n, word)
                reflected = chart_reverse_converse(n, word)
                need(chart_complement(n, complement) == word, "complement involution")
                need(chart_reverse_converse(n, reflected) == word, "reflection involution")
                need(chart_reverse_converse(n, complement) == chart_complement(n, reflected), "chart involutions commute")
                # Independent literal arc flip in the original label frame.
                adjacent = {frozenset((u, v)) for u, v in zip(path, path[1:])}
                changed = code
                for k, (u, v) in enumerate(combinations(range(n), 2)):
                    if frozenset((u, v)) not in adjacent:
                        changed ^= 1 << k
                need(relabel(n, changed, path) == unchart(n, complement), "path complement arc transport")
                opposite = code ^ ((1 << (n * (n - 1) // 2)) - 1)
                need(relabel(n, opposite, path[::-1]) == unchart(n, reflected), "reverse-converse arc transport")
                ports = (0, n - 1, 0) if n else ()
                frame_ports = reframe_ports(path, ports)
                need(marked_key(n, code, ports) == marked_key(n, unchart(n, word), frame_ports), "literal path-port transport")
                need(marked_key(n, changed, ports) == marked_key(n, unchart(n, complement), frame_ports), "complement port transport")
                need(marked_key(n, opposite, ports) == marked_key(n, unchart(n, reflected),
                     tuple(n - 1 - v for v in frame_ports)), "reverse-converse port transport")
                complement_owners.add(canonical(n, changed).code)
                path_cases += 1
            need(set(word_multiplicities.values()) == {len(rec.automorphisms)}, "full path/word orbit multiplicities")
            if frame_dependence is None and len(complement_owners) > 1:
                frame_dependence = {"n": n, "source": code, "target_classes": sorted(complement_owners)}
            if n <= 5:
                lists = [(u, v, u) for u in range(n) for v in range(n)] if n else [()]
                for ports in lists:
                    marked_group = tuple(g for g in rec.automorphisms if all(g[v] == v for v in ports))
                    inventory = Counter((chart(n, code, p)[0], reframe_ports(p, ports)) for p in paths)
                    need(set(inventory.values()) == {len(marked_group)}, "marked path inventory stabilizer")
                    need(len(inventory) * len(marked_group) == len(paths), "marked chart fiber H/Aut(T,ports)")
                    marked_inventory_cases += 1
    need(frame_dependence is not None, "unmarked path-complement false descent not exposed")

    # Hostile 1: root orbits, not a uniform n-to-one multiplicity.
    need(growth[2]["rooted_extensions"] == 4 and Fraction(4, 3) != counts[3], "false divide by n not exposed")
    # Hostile 2: canonicalizing only graph shape or H loses a class at n4.
    same_h = {}
    collision = None
    for code, rec in levels[4].items():
        h = hamilton_count(4, code)
        if h in same_h:
            collision = (same_h[h], code, h)
        same_h[h] = code
    need(collision is not None and collision[0] != collision[1], "H collision missing")
    # Hostile 3: arbitrary transitive blocks merge; they are not unique cores.
    merged = inflate(1, (1, 2))  # 0 -> 1 at quotient order 2.
    need(canonical(3, merged).code == canonical(3, inflate(0, (3,))).code, "clone merger witness")
    need(len(clone_core(3, merged)[0]) == 1, "nonmaximal clone decomposition accepted")
    # Hostile 4: weights without incidence misidentify quotient vertices.
    c3 = canonical(3, 5).code  # directed 0->1->2->0 has bits 1,0,1.
    ga = canonical(3, c3).automorphisms()
    need(len(ga) == 3, "C3 control")
    need(weighted_key(c3, (1, 2, 3)) != weighted_key(c3, (1, 3, 2)), "weight multiset wrongly accepted")
    need(canonical(6, inflate(c3, (1, 2, 3))).code != canonical(6, inflate(c3, (1, 3, 2))).code,
         "opposite weighted C3 controls merged")
    # Hostile 5: a pin's position inside a transitive module must survive.
    tt3 = (1 << 3) - 1
    blocks, _, _ = clone_core(3, tt3)
    ports1, ports2 = (0, 2, 0), (0, 1, 0)
    need(port_coordinates(blocks, ports1) == ((0, 0), (0, 2), (0, 0)), "port internal ranks")
    need(marked_key(3, tt3, ports1) != marked_key(3, tt3, ports2), "port rank loss")
    port_cases = 0
    for p in permutations(range(3)):
        inverse = {v: i for i, v in enumerate(p)}
        need(marked_key(3, relabel(3, tt3, p), tuple(inverse[v] for v in ports1)) == marked_key(3, tt3, ports1),
             "marked port transport")
        port_cases += 1

    return {"status": "PASS", "scope": "exact class enumeration through n8; explicit finite gates",
            "counts": counts, "growth": growth, "brute_labeled_cases_n0_5": brute_cases,
            "all_labelings_of_representatives_n0_5": relabel_cases,
            "literal_chart_cases_n0_7": chart_cases, "charts": chart_receipts,
            "class_and_clone_cases_n0_8": class_cases, "weighted_core_allocations": allocations,
            "masses": class_receipts, "clone_layers": clone_counts,
            "independent_reverse_counts": exact_counts,
            "literal_path_and_operation_cases_n0_6": path_cases,
            "marked_port_inventories_n0_5": marked_inventory_cases,
            "hostiles": {"false_divide_by_n": "rooted=4, n=3, classes=2",
                         "first_H_collision": collision, "arbitrary_clone_partition_merges": True,
                         "C3_weight_multiset_loses_orientation": True,
                         "repeated_ports_and_internal_ranks": port_cases,
                         "path_complement_depends_on_retained_spine": frame_dependence}}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
