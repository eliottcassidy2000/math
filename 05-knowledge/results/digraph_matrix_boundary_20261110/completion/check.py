#!/usr/bin/env python3
"""Exact small controls for named path contractions and resolution overlap.

All checks use explicit exceptions, including under python -O. No external
graph library, imported research implementation, or random seed is used.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def h_dp(rows):
    n = len(rows)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for u, weight in enumerate(dp[mask]):
            if not weight:
                continue
            available = rows[u] & ~mask
            while available:
                bit = available & -available
                dp[mask | bit][bit.bit_length() - 1] += weight
                available -= bit
    return sum(dp[-1])


def words(rows):
    return [p for p in permutations(range(len(rows)))
            if all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))]


def from_code(n, code, radix=4, fixed_pair=None):
    rows = [0] * n
    for u, v in combinations(range(n), 2):
        if (u, v) == fixed_pair:
            digit = 3
        else:
            digit = code % radix
            code //= radix
            if radix == 2:
                digit += 1
        if digit & 1:
            rows[u] |= 1 << v
        if digit & 2:
            rows[v] |= 1 << u
    need(code == 0, "Unused input digits")
    return tuple(rows)


def forest_blocks(n, selected):
    neighbors = [set() for _ in range(n)]
    for u, v in selected:
        need(u != v and 0 <= u < n and 0 <= v < n, "Invalid selected pair")
        neighbors[u].add(v)
        neighbors[v].add(u)
    if any(len(v) > 2 for v in neighbors):
        return None
    unseen = set(range(n))
    blocks = []
    while unseen:
        first = min(unseen)
        stack = [first]
        component = set()
        while stack:
            v = stack.pop()
            if v in component:
                continue
            component.add(v)
            stack.extend(neighbors[v] - component)
        unseen -= component
        if len(component) == 1:
            blocks.append((first,))
            continue
        ends = sorted(v for v in component if len(neighbors[v]) == 1)
        if len(ends) != 2:
            return None
        path, previous, current = [], None, ends[0]
        while current is not None:
            path.append(current)
            nxt = neighbors[current] - {previous}
            previous, current = current, (next(iter(nxt)) if nxt else None)
        need(set(path) == component, "Path component reconstruction")
        blocks.append(tuple(path))
    return blocks


def endpoint_quotient(rows, directed_blocks):
    result = [0] * len(directed_blocks)
    for i, block in enumerate(directed_blocks):
        for j, other in enumerate(directed_blocks):
            if i != j and rows[block[-1]] >> other[0] & 1:
                result[i] |= 1 << j
    return tuple(result)


def contraction_count(rows, selected):
    blocks = forest_blocks(len(rows), selected)
    if blocks is None:
        return 0
    choices = []
    for block in blocks:
        directions = [block] if len(block) == 1 else [block, block[::-1]]
        choices.append([b for b in directions
                        if all(rows[u] >> v & 1 for u, v in zip(b, b[1:]))])
    return sum(h_dp(endpoint_quotient(rows, orientation))
               for orientation in product(*choices))


def all_pair_responses_by_words(rows):
    n = len(rows)
    pairs = list(combinations(range(n), 2))
    indices = {p: i for i, p in enumerate(pairs)}
    exact = [0] * (1 << len(pairs))
    paths = words(rows)
    for path in paths:
        used = sum(1 << indices[tuple(sorted((u, v)))]
                   for u, v in zip(path, path[1:]))
        exact[used] += 1
    contains = exact.copy()
    for i in range(len(pairs)):
        for mask in range(1 << len(pairs)):
            if not mask >> i & 1:
                contains[mask] += contains[mask | (1 << i)]
    return pairs, contains, paths


def tournament_embedding(rows):
    m = len(rows)
    n = 2 * m
    out = [0] * n
    for a, b in combinations(range(n), 2):
        i, kind_a = divmod(a, 2)
        j, kind_b = divmod(b, 2)
        if i == j or kind_a == kind_b:
            u, v = a, b
        elif kind_a == 1:
            u, v = (a, b) if rows[i] >> j & 1 else (b, a)
        else:
            u, v = (b, a) if rows[j] >> i & 1 else (a, b)
        out[u] |= 1 << v
    result = tuple(out)
    need(all(((result[u] >> v) ^ (result[v] >> u)) & 1
             for u, v in combinations(range(n), 2)), "Embedding is a tournament")
    blocks = [(2 * i, 2 * i + 1) for i in range(m)]
    need(endpoint_quotient(result, blocks) == rows, "Embedding recovers original rows")
    return result, blocks


def forced_arcs_inclusion_exclusion(rows, arcs):
    result = 0
    for deleted in range(1 << len(arcs)):
        altered = list(rows)
        for i, (u, v) in enumerate(arcs):
            if deleted >> i & 1:
                altered[u] &= ~(1 << v)
        result += (-1 if deleted.bit_count() % 2 else 1) * h_dp(tuple(altered))
    return result


def matrix(rows):
    return [''.join(str(row >> v & 1) for v in range(len(rows))) for row in rows]


def strong(rows):
    n = len(rows)
    for start in range(n):
        reached, fringe = {start}, [start]
        while fringe:
            u = fringe.pop()
            for v in range(n):
                if rows[u] >> v & 1 and v not in reached:
                    reached.add(v)
                    fringe.append(v)
        if len(reached) != n:
            return False
    return True


def canonical(rows):
    images = []
    for p in permutations(range(len(rows))):
        transformed = tuple(sum(((rows[p[u]] >> p[v]) & 1) << v
                                for v in range(len(rows)))
                            for u in range(len(rows)))
        images.append((transformed, p))
    return min(images)


def one_pair_record(rows):
    need(rows[0] & 2 and rows[1] & 1, "Marked digon exists")
    resolutions = []
    for remove in ((1, 0), (0, 1)):
        r = list(rows)
        r[remove[0]] &= ~(1 << remove[1])
        r = tuple(r)
        c, p = canonical(r)
        resolutions.append(dict(rows=r, H=h_dp(r), canonical=c,
                                canonical_order=p, words=words(r)))
    missing = list(rows)
    missing[0] &= ~2
    missing[1] &= ~1
    path_sets = [set(tuple(p) for p in record['words']) for record in resolutions]
    need(len(path_sets[0] | path_sets[1]) == h_dp(rows), "Union count")
    need(len(path_sets[0] & path_sets[1]) == h_dp(tuple(missing)), "Intersection count")
    return dict(rows=rows, matrix=matrix(rows), H=h_dp(rows), strong=strong(rows),
                missing_H=h_dp(tuple(missing)), resolutions=resolutions,
                intersection=sorted(path_sets[0] & path_sets[1]), words=words(rows))


def padded(rows, before, after):
    n = len(rows)
    total = before + n + after
    result = []
    for u in range(total):
        if before <= u < before + n:
            result.append((rows[u - before] << before) |
                          sum(1 << v for v in range(before + n, total)))
        else:
            result.append(sum(1 << v for v in range(u + 1, total)))
    return tuple(result)


def main():
    graph_count = response_count = forbidden_support_count = 0
    by_order = []
    for n in range(5):
        count = 4 ** (n * (n - 1) // 2)
        local_responses = 0
        for code in range(count):
            rows = from_code(n, code)
            pairs, contains, paths = all_pair_responses_by_words(rows)
            need(h_dp(rows) == len(paths), "Subset DP versus permutation enumeration")
            for selected_mask, expected in enumerate(contains):
                selected = [p for i, p in enumerate(pairs) if selected_mask >> i & 1]
                actual = contraction_count(rows, selected)
                need(actual == expected, "Linear-forest contraction identity")
                if forest_blocks(n, selected) is None:
                    forbidden_support_count += 1
                    need(expected == 0, "Cycle or branching pair-set cannot occur in a path")
                response_count += 1
                local_responses += 1
            graph_count += 1
        by_order.append(dict(n=n, graphs=count, selected_pair_responses=local_responses))

    minimality = []
    for n in (2, 3):
        buckets = {}
        for code in range(2 ** (n * (n - 1) // 2 - 1)):
            rows = from_code(n, code, radix=2, fixed_pair=(0, 1))
            record = one_pair_record(rows)
            vector = tuple(r['H'] for r in record['resolutions'])
            buckets.setdefault(vector, set()).add(record['H'])
        need(all(len(v) == 1 for v in buckets.values()), "No smaller one-digon scalar collision")
        minimality.append(dict(n=n, indexed_resolution_vectors=[
            dict(vector=k, possible_union_H=sorted(v)) for k, v in sorted(buckets.items())]))

    collision = [one_pair_record(rows) for rows in ((6, 9, 10, 1), (6, 5, 8, 3))]
    need([r['H'] for r in collision] == [6, 8], "Collision union counts")
    need(all(r['strong'] for r in collision), "Collision graphs strong")
    need(all([q['H'] for q in r['resolutions']] == [5, 5] for r in collision), "Collision resolution vector")
    need(len({tuple(q['canonical']) for r in collision for q in r['resolutions']}) == 1,
         "All four tournament resolutions are isomorphic")
    padding_checks = 0
    for record in collision:
        for before in range(3):
            for after in range(3):
                need(h_dp(padded(record['rows'], before, after)) == record['H'], "Padding retains count")
                padding_checks += 1

    # A tournament path contraction can create either defect state already at n=3.
    tt3 = (6, 4, 0)
    c3 = (2, 4, 1)
    defects = []
    for original, blocks, expected in ((tt3, [(0, 2), (1,)], (0, 0)),
                                        (c3, [(0, 1), (2,)], (2, 1))):
        quotient = endpoint_quotient(original, blocks)
        need(quotient == expected, "Minimal contraction defect")
        defects.append(dict(tournament_rows=original, blocks=blocks,
                            quotient_rows=quotient, quotient_H=h_dp(quotient)))

    embedding_controls = 0
    for n in range(4):
        for code in range(4 ** (n * (n - 1) // 2)):
            rows = from_code(n, code)
            tournament, arcs = tournament_embedding(rows)
            constrained = sum(all((u, v) in set(zip(p, p[1:])) for u, v in arcs)
                              for p in words(tournament))
            need(constrained == h_dp(rows), "Literal embedded constrained words")
            need(forced_arcs_inclusion_exclusion(tournament, arcs) == constrained,
                 "Independent arc-deletion inclusion-exclusion")
            embedding_controls += 1

    targets = []
    for rows, expected in (((14, 12, 9, 1), 7), ((30, 12, 24, 17, 3), 21)):
        need(h_dp(rows) == expected, "Inherited target rows")
        tournament, arcs = tournament_embedding(rows)
        literal = forced_arcs_inclusion_exclusion(tournament, arcs)
        need(literal == expected, "Seven/twenty-one constrained tournament response")
        full = h_dp(tournament)
        need(full % 2 == 1 and full not in (7, 21), "Full count is not the constrained count")
        targets.append(dict(input_rows=rows, input_H=expected, tournament_rows=tournament,
                            tournament_matrix=matrix(tournament), matching_arcs=arcs,
                            endpoint_quotient_rows=endpoint_quotient(tournament, arcs),
                            constrained_H=literal, full_tournament_H=full))

    # Branch/cycle supports are rejected; reversed forced arcs are not inferred.
    need(forest_blocks(4, [(0, 1), (0, 2), (0, 3)]) is None, "Degree-three hostile")
    need(forest_blocks(3, [(0, 1), (1, 2), (0, 2)]) is None, "Cycle hostile")
    need(contraction_count((2, 0), [(0, 1)]) == 1, "Only legal block orientation contributes")
    sources = [
        '05-knowledge/results/session234_digraphs_20261109/PROOFS.md',
        '05-knowledge/results/session235_arc_pair_models_20261109/theory/PROOFS.md',
        '01-canon/theorems/THM-5041-matrix-tiling-pair-dictionary.md',
        '01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md',
    ]
    report = dict(status='PASS', scope='All loopless labelled digraphs of orders0..4 and every selected unordered-pair set; exact additional named witnesses. No spectrum completeness claim.',
                  graph_count=graph_count, response_count=response_count,
                  forbidden_support_count=forbidden_support_count, by_order=by_order,
                  one_digon_collision=collision, one_digon_minimality=minimality,
                  padding_checks=padding_checks, minimal_contraction_defects=defects,
                  embedding_literal_controls=embedding_controls, constrained_target_witnesses=targets,
                  source_sha256={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in sources},
                  checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
