#!/usr/bin/env python3
"""Exact independent controls; Python standard library only, never bare assert."""
import argparse
import csv
import hashlib
import itertools as it
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def rational(value):
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pad(poly, length):
    require(len(poly) <= length, "polynomial padding too short")
    return poly + [0] * (length - len(poly))


class Universe:
    def __init__(self, n):
        self.n = n
        self.pairs = list(it.combinations(range(n), 2))
        self.index = {pair: i for i, pair in enumerate(self.pairs)}
        self.N = len(self.pairs)
        self.full = (1 << self.N) - 1
        self.permutations = list(it.permutations(range(n)))
        self.path_data = [(p, self.edges(p)) for p in self.permutations]

    def edge_index(self, u, v):
        return self.index[tuple(sorted((u, v)))]

    def edges(self, path):
        mask = bits = 0
        for u, v in zip(path, path[1:]):
            bit = 1 << self.edge_index(u, v)
            mask |= bit
            if u < v:
                bits |= bit
        return mask, bits

    def arc(self, source, u, v):
        return bool((source >> self.edge_index(u, v)) & 1) == (u < v)

    def cell_mask(self, cell):
        return sum(1 << self.index[pair] for pair in it.combinations(cell, 2))

    def moments(self, source):
        degrees = [0] * self.n
        for j, (u, v) in enumerate(self.pairs):
            degrees[u if (source >> j) & 1 else v] += 1
        score_square = sum((2 * d - self.n + 1) ** 2 for d in degrees)
        # Direct directed-triple test, independent of the score identity.
        cycles = sum(self.arc(source, u, v) == self.arc(source, v, w)
                     == self.arc(source, w, u)
                     for u, v, w in it.combinations(range(self.n), 3))
        return score_square, cycles

    def relabel(self, source, permutation):
        out = 0
        for j, (u, v) in enumerate(self.pairs):
            a, b = (u, v) if (source >> j) & 1 else (v, u)
            a, b = permutation[a], permutation[b]
            if a < b:
                out |= 1 << self.edge_index(a, b)
        return out

    def records(self):
        records, hosts = [], []
        for source in range(1 << self.N):
            square, cycles = self.moments(source)
            require(24 * cycles == self.n * (self.n ** 2 - 1) - 3 * square,
                    "direct triangle census disagrees with score identity")
            paths = [p for p, (mask, bits) in self.path_data
                     if source & mask == bits]
            hosts.append((source, len(paths), square, cycles))
            for path in paths:
                mask, _ = self.edges(path)
                partner = source ^ (self.full ^ mask)
                records.append((source, path, partner, square, cycles))
        return records, hosts


def orient_free(base, positions):
    for choices in range(1 << len(positions)):
        yield base | sum(1 << position for j, position in enumerate(positions)
                         if (choices >> j) & 1)


def census(n):
    u = Universe(n)
    records, hosts = u.records()
    f = u.N - n + 1
    count = math.factorial(n) * (1 << f)
    expected_square = n * n - 3 * n + 4 if n >= 2 else 0
    expected_cycles = Fraction(math.comb(n, 3) + n - 2, 4) if n >= 2 else 0
    require(len(records) == count, "marked record count")
    require(sum(r[3] for r in records) == count * expected_square,
            "marked score-square expectation")
    require(sum(r[4] for r in records) == count * expected_cycles,
            "marked cyclic-triple expectation")

    # Independent traversal: observer first, then its independent nonpath arcs.
    observer_first = set()
    for path, (mask, bits) in u.path_data:
        positions = [j for j in range(u.N) if not (mask >> j) & 1]
        for source in orient_free(bits, positions):
            observer_first.add((source, path))
    require(observer_first == {(r[0], r[1]) for r in records},
            "two independent enumerations of marked records")
    require(len({(r[0], r[2]) for r in records}) == len(records),
            "labelled source/partner pair recovers its observer uniquely")
    for source, path, partner, _, _ in records:
        mask, bits = u.edges(path)
        require((u.full ^ (source ^ partner)) == mask and partner & mask == bits,
                "source/partner agreement arcs recover the directed path")

    eh = Fraction(count, 1 << u.N)
    ec = Fraction(sum(h[3] for h in hosts), 1 << u.N)
    ehc = Fraction(sum(h[1] * h[3] for h in hosts), 1 << u.N)
    covariance = ehc - eh * ec
    require(eh == Fraction(math.factorial(n), 1 << (n - 1)), "uniform E H")
    require(ec == Fraction(math.comb(n, 3), 4), "uniform E cyclic triples")
    require(covariance == (eh * Fraction(n - 2, 4) if n >= 2 else 0),
            "uniform host covariance")

    # Small-order orbit grouping is only an audit of orbit-stabilizer weights.
    groups, canonical_of = defaultdict(list), {}
    for host in hosts:
        canonical = min(u.relabel(host[0], p) for p in u.permutations)
        groups[canonical].append(host)
        canonical_of[host[0]] = canonical
    fixed_path_frequencies = Counter((r[1], canonical_of[r[0]]) for r in records)
    weighted = [Fraction(0), Fraction(0), Fraction(0)]
    by_class = {}
    for representative, group in groups.items():
        aut = sum(u.relabel(representative, p) == representative
                  for p in u.permutations)
        require(len(group) * aut == math.factorial(n), "orbit-stabilizer")
        _, h, square, cycles = group[0]
        require(all(g[1:] == (h, square, cycles) for g in group),
                "invariant statistics under relabelling")
        require(h % aut == 0, "automorphisms act freely on ordered Hamiltonian paths")
        for path in u.permutations:
            require(fixed_path_frequencies[path, representative] == h // aut,
                    "fixed-observer realization count H/Aut")
        by_class[representative] = (h // aut, square, cycles)
        for j, statistic in enumerate((1, square, cycles)):
            weighted[j] += Fraction(h * statistic, aut)
    require(weighted == [Fraction(1 << f), Fraction((1 << f) * expected_square),
                         (1 << f) * expected_cycles], "automorphism-weighted sums")
    visited, merged_sums = set(), [0, 0, 0]
    for representative, (multiplicity, square, cycles) in by_class.items():
        if representative in visited:
            continue
        opposite = canonical_of[u.full ^ representative]
        require(by_class[opposite] == (multiplicity, square, cycles),
                "converse preserves invariant statistics and H/Aut")
        mass = multiplicity if opposite == representative else 2 * multiplicity
        visited.update((representative, opposite))
        for j, statistic in enumerate((1, square, cycles)):
            merged_sums[j] += mass * statistic
    require(merged_sums == weighted, "converse-merged fixed-observer checksums")
    return u, records, {
        "n": n, "labelled_hosts": 1 << u.N, "marked_records": count,
        "sum_H_score_square": sum(r[3] for r in records),
        "sum_H_cyclic_triangles": sum(r[4] for r in records),
        "E_H_score_square": rational(expected_square),
        "E_H_cyclic_triangles": rational(expected_cycles),
        "uniform_E_H": rational(eh), "uniform_E_cyclic_triangles": rational(ec),
        "uniform_Cov_H_cyclic_triangles": rational(covariance),
        "automorphism_weighted_sums_H_Hscore_Hcycles": list(map(rational, weighted)),
        "fixed_observer_class_counts_verified": True,
        "converse_merged_sums_mass_massscore_masscycles": merged_sums,
    }


def path_blocks(n, local_paths):
    edges = set(edge for path in local_paths for edge in zip(path, path[1:]))
    out, incoming = defaultdict(set), defaultdict(set)
    for a, b in edges:
        out[a].add(b)
        incoming[b].add(a)
    if any(len(out[v]) > 1 or len(incoming[v]) > 1 for v in range(n)):
        return None
    visited, blocks = set(), []
    for head in range(n):
        if incoming[head]:
            continue
        block, current = [], head
        while current not in visited:
            visited.add(current)
            block.append(current)
            if not out[current]:
                break
            current = next(iter(out[current]))
        blocks.append(tuple(block))
    return tuple(blocks) if len(visited) == n else None


def allowed(path, pins):
    return all(0 <= clock < len(path) and path[clock] == owner
               for owner, clock in pins)


def block_predictions(u, source, covered, local_paths):
    blocks = path_blocks(u.n, local_paths)
    if blocks is None:
        return None, []
    f = u.N - covered.bit_count() - len(blocks) + 1
    require(f >= 0, "negative free exponent")
    orders = []
    for order in it.permutations(blocks):
        path = tuple(v for block in order for v in block)
        forced = [(a, b) if (source >> j) & 1 else (b, a)
                  for j, (a, b) in enumerate(u.pairs) if (covered >> j) & 1]
        joins = [(a[-1], b[0]) for a, b in zip(order, order[1:])]
        require(all(not (covered >> u.edge_index(a, b)) & 1 for a, b in joins),
                "a joining pair must be uncovered")
        forced += joins
        b = [0] * u.n
        for tail, head in forced:
            b[tail] += 1
            b[head] -= 1
        orders.append((path, sum(a < b for a, b in joins), sum(x * x for x in b)))
    return f, orders


def local_frames(universes):
    stats = Counter()
    for n, (u, records) in universes.items():
        if n > 4:
            continue
        cells = [cell for k in range(1, n + 1)
                 for cell in it.combinations(range(n), k)]
        families = {(), *((cell,) for cell in cells), *it.combinations(cells, 2)}
        if n >= 3:
            families.add(tuple((0, v) for v in range(1, n)))
            families.add(((0, 1), (0, 2), (1, 2)))
        for family in sorted(families):
            masks = [u.cell_mask(cell) for cell in family]
            covered = 0
            for mask in masks:
                covered |= mask
            positions = [j for j in range(u.N) if (covered >> j) & 1]
            for base in orient_free(0, positions):
                menus = [[path for path in it.permutations(cell)
                          if all(u.arc(base, a, b) for a, b in zip(path, path[1:]))]
                         for cell in family]
                hosts = [record for record in records if record[0] & covered == base]
                for paths in it.product(*menus):
                    # Compare the entire local K-partner, not merely path subsequences.
                    targets = [(base ^ (mask ^ u.edges(path)[0])) & mask
                               for mask, path in zip(masks, paths)]
                    actual = [r for r in hosts
                              if all(r[2] & mask == target
                                     for mask, target in zip(masks, targets))]
                    f, orders = block_predictions(u, base, covered, paths)
                    stats["full_local_record_instances"] += 1
                    stats["forest_failures"] += f is None
                    stats["marked_completions_without_pins"] += len(actual)
                    pin_cases = [()] + [((v, k),) for v in range(n) for k in range(n)]
                    pin_cases += [((0, n),)]
                    if n >= 2:
                        pin_cases += [((0, 0), (1, 0)), ((0, 0), (0, 1))]
                    for pins in pin_cases:
                        accepted = [r for r in actual if allowed(r[1], pins)]
                        accepted_orders = [o for o in orders if allowed(o[0], pins)]
                        stats["pin_conditioned_cases"] += 1
                        stats["zero_completion_cases"] += not accepted
                        polynomial = Counter((r[0] & (u.full ^ covered)).bit_count()
                                             for r in accepted)
                        expected = Counter()
                        square_sum = 0
                        for _, a, bnorm in accepted_orders:
                            for k in range(f + 1):
                                expected[a + k] += math.comb(f, k)
                            square_sum += (1 << f) * (2 * f + bnorm)
                        require(polynomial == expected, "full-frame completion polynomial")
                        require(sum(r[3] for r in accepted) == square_sum,
                                "full-frame imbalance-square sum")
                        require(24 * sum(r[4] for r in accepted)
                                == n * (n * n - 1) * len(accepted) - 3 * square_sum,
                                "full-frame cyclic-triple sum")
                        if accepted:
                            average = Fraction(square_sum, len(accepted))
                            predicted = 2 * f + Fraction(sum(o[2] for o in accepted_orders),
                                                        len(accepted_orders))
                            require(average == predicted, "nonempty conditional mean")
                        else:
                            require(not accepted_orders, "empty case must not invent a measure")
    require(path_blocks(4, [(0, 1), (0, 2), (0, 3)]) is None, "branch hostile")
    require(path_blocks(3, [(0, 1), (1, 2), (2, 0)]) is None, "cycle hostile")
    return dict(stats)


def hostile():
    blocks = [(0, 3), (6, 1), (2, 5), (4, 7)]
    n, e, f, N = 8, 24, 21, 28
    histogram = Counter(sum(a[-1] < b[0] for a, b in zip(order, order[1:]))
                        for order in it.permutations(blocks))
    reverse_blocks = [block[::-1] for block in blocks]
    reverse_histogram = Counter(sum(a[-1] < b[0] for a, b in zip(order, order[1:]))
                                for order in it.permutations(reverse_blocks))
    B = [histogram[k] for k in range(4)]
    Br = [reverse_histogram[k] for k in range(4)]
    require(B == [4, 10, 10, 0], "hostile block-order polynomial")
    require(Br == [0, 10, 10, 4], "literal-converse block-order polynomial")
    binomial = [math.comb(f, k) for k in range(f + 1)]
    # Independent repeated binary convolution, representing every free-arc bit.
    binary_choices = [1]
    for _ in range(f):
        binary_choices = multiply(binary_choices, [1, 1])
    require(binary_choices == binomial, "independent free-arc multiplicities")
    G = pad(multiply(B[:-1], binary_choices), e + 1)
    Gr = pad(multiply(Br, binomial), e + 1)
    require(Gr == G[::-1], "uncovered converse reciprocity")
    covered_forward = sum(a < b for a, b in blocks)
    require(covered_forward == 3, "original covered forward arc count")
    F = pad([0] * covered_forward + G, N + 1)
    Fr = pad([0] * (4 - covered_forward) + Gr, N + 1)
    union = [a + b for a, b in zip(F, Fr)]
    factorized = [0, 0] + multiply([2 * math.comb(22, k) for k in range(23)], [5, 2, 5])
    require(union == pad(factorized, N + 1), "full-source factorization")
    require(union == union[::-1], "28-reciprocal converse-symmetric union")
    require(sum(G) == math.factorial(4) * (1 << f), "single input total")
    require(sum(union) == 2 * math.factorial(4) * (1 << f), "union total")
    require(10 * 10 - 4 * 10 * 4 == -60, "original nonreal-root discriminant")
    require(2 * 2 - 4 * 5 * 5 == -96, "symmetric nonreal-root discriminant")
    require(F[::-1] == Fr, "all-source converse reciprocity")
    for family in (blocks, reverse_blocks):
        for order in it.permutations(family):
            divergence = [0] * n
            arcs = list(family) + [(a[-1], b[0]) for a, b in zip(order, order[1:])]
            for tail, head in arcs:
                divergence[tail] += 1
                divergence[head] -= 1
            require(sum(x * x for x in divergence) == 2,
                    "hostile fixed arcs form a single Hamiltonian path")
            require(2 * f + sum(x * x for x in divergence) == 44,
                    "hostile conditional imbalance-square expectation")
    require(Fraction(n * (n * n - 1), 24) - Fraction(44, 8) == Fraction(31, 2),
            "hostile conditional cyclic-triangle expectation")
    return {
        "n": n, "blocks": blocks, "p": 4, "e": e, "free_arcs_per_order": f,
        "block_ascent_coefficients": B, "converse_block_ascent_coefficients": Br,
        "covered_forward_arcs_original": covered_forward,
        "covered_forward_arcs_converse": 4 - covered_forward,
        "uncovered_coefficients_original": G, "uncovered_coefficients_converse": Gr,
        "all_forward_coefficients_original": F, "all_forward_coefficients_converse": Fr,
        "all_forward_coefficients_disjoint_union": union,
        "marked_completions_original": sum(G), "marked_completions_union": sum(union),
        "factor_original": "(1+z)^21 (4+10z+10z^2)",
        "factor_union_all_forward": "2 z^2 (1+z)^22 (5+2z+5z^2)",
        "quadratic_discriminants": [-60, -96],
        "reciprocity_exponent_all_forward": N,
        "E_score_square_both_inputs": 44, "E_cyclic_triangles_both_inputs": "31/2",
    }


def eulerian_controls():
    rows, recurrence = [], [1]
    for n in range(1, 9):
        direct = Counter(sum(a < b for a, b in zip(p, p[1:]))
                         for p in it.permutations(range(n)))
        require([direct[k] for k in range(n)] == recurrence, "Eulerian insertion recurrence")
        rows.append({"n": n, "coefficients": recurrence})
        following = [0] * (n + 1)
        for k, coefficient in enumerate(recurrence):
            following[k] += (k + 1) * coefficient
            following[k + 1] += (n - k) * coefficient
        recurrence = following
    return rows


def metagraph_atlas():
    """Audit frozen original fibres, never infer them from Laplacian diagonals."""
    here = Path(__file__).resolve().parent
    source = here / "sources" / "merged_nodes.tsv.txt"
    pins = json.loads((here / "sources" / "SOURCE_PINS.json").read_text())
    require(hashlib.sha256(source.read_bytes()).hexdigest() == pins["atlas"]["sha256"],
            "frozen original atlas byte pin")
    rows = list(csv.DictReader(source.read_text().splitlines(), delimiter="\t"))
    require(len(rows) == 3528, "frozen order-eight atlas size")
    sums = [0, 0, 0]
    maximum_aut = 0

    def matrix(code):
        a = [[0] * 8 for _ in range(8)]
        for (v, w), bit in zip(it.combinations(range(8), 2), f"{code:028b}"):
            a[v][w] = int(bit)
            a[w][v] = 1 - int(bit)
        return a

    def isomorphism_count(a, b):
        da, db = list(map(sum, a)), list(map(sum, b))
        if sorted(da) != sorted(db):
            return 0
        order = sorted(range(8), key=lambda v: (da.count(da[v]), v))
        image = {}
        used = set()

        def search(depth):
            if depth == 8:
                return 1
            v = order[depth]
            total = 0
            for w in range(8):
                if w in used or da[v] != db[w]:
                    continue
                if any(a[v][old] != b[w][mapped] for old, mapped in image.items()):
                    continue
                image[v] = w
                used.add(w)
                total += search(depth + 1)
                used.remove(w)
                del image[v]
            return total
        return search(0)

    def paths(a):
        # Independent exact subset DP, not path-permutation brute force.
        dp = [[0] * 8 for _ in range(256)]
        for v in range(8):
            dp[1 << v][v] = 1
        for mask in range(1, 256):
            for last in range(8):
                if not dp[mask][last]:
                    continue
                for following in range(8):
                    if not (mask >> following) & 1 and a[last][following]:
                        dp[mask | (1 << following)][following] += dp[mask][last]
        return sum(dp[-1])

    matrix_packet = json.loads((here / "sources" / "prior_four_source_packet.json.txt").read_text())
    matrix_controls = []
    for slot in matrix_packet["slots"]:
        row = rows[slot["source"]]
        require(int(row["canon_a"]) == slot["canon"], "prior actual source code control")
        require(matrix(slot["canon"]) == slot["canonical_matrix"],
                "canonical encoding agrees with frozen explicit source matrix")
        matrix_controls.append(slot["source"])
    for index, row in enumerate(rows):
        require(int(row["merged_id"]) == index, "atlas row identity")
        a, b = matrix(int(row["canon_a"])), matrix(int(row["canon_b"]))
        aut = isomorphism_count(a, a)
        maximum_aut = max(maximum_aut, aut)
        h = paths(a)
        require(aut > 0 and h % aut == 0, "atlas path/aut integrality")
        require(isomorphism_count(a, list(map(list, zip(*b)))) == aut,
                "each atlas pair is literally converse-isomorphic")
        self_converse = row["class_a"] == row["class_b"]
        require((int(row["canon_a"]) == int(row["canon_b"])) == self_converse,
                "atlas self-converse class flag")
        multiplicity = h // aut * (1 if self_converse else 2)
        require(multiplicity == int(row["fibre"]), "every original fibre equals H/Aut or 2H/Aut")
        square = sum((2 * d - 7) ** 2 for d in map(sum, a))
        cycles = sum(a[u][v] == a[v][w] == a[w][u]
                     for u, v, w in it.combinations(range(8), 3))
        require(24 * cycles == 8 * 63 - 3 * square, "atlas independent triangle/score identity")
        for j, statistic in enumerate((1, square, cycles)):
            sums[j] += multiplicity * statistic
    require(sums == [1 << 21, 92274688, 32505856], "order-eight metagraph checksums")
    return {
        "rows_checked": len(rows), "every_fibre_verified_by_H_over_Aut": True,
        "H_method": "exact Hamiltonian-subset dynamic program",
        "Aut_method": "degree-filtered exact adjacency-preserving bijection search",
        "largest_automorphism_order_observed": maximum_aut,
        "explicit_matrix_encoding_controls": matrix_controls,
        "sum_fibre": sums[0], "sum_fibre_score_square": sums[1],
        "sum_fibre_cyclic_triangles": sums[2],
        "no_laplacian_degree_or_stationary_walk_identification": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    universes, censuses = {}, []
    for n in range(1, 6):
        u, records, row = census(n)
        universes[n] = (u, records)
        censuses.append(row)
    local = local_frames(universes)
    nonreal = hostile()
    eulerian = eulerian_controls()
    atlas = metagraph_atlas()
    report = {
        "status": "PASS", "checks": CHECKS, "weighted_host_census": censuses,
        "full_local_frame_controls": local, "nonreal_root_controls": nonreal,
        "Eulerian_recurrence_controls": eulerian,
        "original_order8_metagraph_atlas": atlas,
        "scope": "All labelled tournaments and marked Hamiltonian observers through n=5; "
                 "all coherent zero/one/two-distinct-cell families of arbitrary nonempty cell sizes "
                 "on n<=4, all local Hamiltonian records, all single global pins, impossible "
                 "double/out-of-range pins, branch/cycle controls; exact 24-order hostile "
                 "polynomials with all free-arc multiplicities. No large-host brute census, "
                 "unmarked isomorphism enumeration formula, or spectral assertion.",
    }
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
