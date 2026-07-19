#!/usr/bin/env python3
"""Dependency-free exact referee for THM-1217.

For a sum/difference beat denominator q and fast speeds b_i, put

    g_i = gcd(b_i,q),  Q_i = q/g_i,
    d0 = gcd(q,b_1,...,b_6),  L = q/d0 = lcm_i Q_i.

The strict danger mask for b_i has period Q_i and A(Q_i) residues, where
A(Q)=2*ceil(Q/14)-1.  Its lift to Z/L has (L/Q_i)*A(Q_i) residues.  A
defining sum/difference pair gives the same mask, leaving five labelled
obligations.  They all contain zero.

This referee checks four nested certificates:

* the common-zero union budget C0 = 1 + sum_i (|M_i|-1);
* the Hunter/spanning-tree budget
      C_T = sum_i |M_i| - sum_(ij in T) |M_i intersect M_j|;
* the exact longest cyclic dangerous run ell(U);
* direct strict-radius safety at the supplied beat point.

The escape statements require U to be a proper subset of Z/L.  The
cardinality-only supplier C0+1 additionally uses C0<L.  The script audits
these hypotheses explicitly instead of silently treating a full mask as a
finite cyclic run.

Tournament audit
----------------
The pairwise observable on mask obligations is intersection cardinality.
It is symmetric, so orienting by (period,label) is only a tie gauge.  On five
obligations the gauge is transitive: score histogram 0,1,2,3,4; zero directed
cycles; five singleton SCCs; one Hamiltonian path; and ten edge flips under
gauge reversal.  The useful object is instead the rooted quotient-clock
incidence tree L -> Q_i together with the master residues and block phase.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd


def require(condition: bool, message: object) -> None:
    """Always-on certificate check, including under ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def ceil_div(n: int, d: int) -> int:
    require(d > 0, (n, d))
    return -((-n) // d)


def lcm(x: int, y: int) -> int:
    require(x > 0 and y > 0, (x, y))
    return x // gcd(x, y) * y


def lcm_many(values: tuple[int, ...] | list[int]) -> int:
    answer = 1
    for value in values:
        answer = lcm(answer, value)
    return answer


def gcd_many(values: tuple[int, ...] | list[int]) -> int:
    answer = 0
    for value in values:
        answer = gcd(answer, value)
    return answer


def divisors(n: int) -> list[int]:
    require(n > 0, n)
    return [d for d in range(1, n + 1) if n % d == 0]


def window_count(Q: int) -> int:
    require(Q > 0, Q)
    return 2 * ceil_div(Q, 14) - 1


def reduced_mask(Q: int, unit: int) -> frozenset[int]:
    require(Q > 0 and gcd(unit, Q) == 1, (Q, unit))
    return frozenset(
        p for p in range(Q)
        if 14 * min((unit * p) % Q, (-(unit * p)) % Q) < Q
    )


def lifted_mask(L: int, Q: int, unit: int) -> frozenset[int]:
    require(L > 0 and Q > 0 and L % Q == 0, (L, Q))
    base = reduced_mask(Q, unit)
    result = frozenset(p for p in range(L) if p % Q in base)
    require(len(result) == (L // Q) * window_count(Q),
            (L, Q, unit, len(result)))
    return result


def direct_q_mask(q: int, speed: int) -> frozenset[int]:
    require(q > 0 and speed > 0, (q, speed))
    return frozenset(
        p for p in range(q)
        if 14 * min((speed * p) % q, (-speed * p) % q) < q
    )


def longest_proper_cyclic_run(values: set[int], L: int) -> int | None:
    """Longest run in a proper cyclic subset; ``None`` for the full clock."""
    require(L > 0 and values <= set(range(L)), (L, values))
    if len(values) == L:
        return None
    best = current = 0
    for index in range(2 * L):
        if index % L in values:
            current += 1
            best = max(best, current)
        else:
            current = 0
    require(best < L, (L, values, best))
    return best


def prufer_tree(code: tuple[int, ...], vertices: int = 5) -> tuple[tuple[int, int], ...]:
    """Decode a Prüfer word; all 5^3 labelled spanning trees are obtained."""
    require(vertices >= 2 and len(code) == vertices - 2, (vertices, code))
    degree = [1] * vertices
    for vertex in code:
        require(0 <= vertex < vertices, (vertices, code))
        degree[vertex] += 1
    edges: list[tuple[int, int]] = []
    for vertex in code:
        leaf = next(i for i, value in enumerate(degree) if value == 1)
        edges.append((leaf, vertex))
        degree[leaf] -= 1
        degree[vertex] -= 1
    remaining = [i for i, value in enumerate(degree) if value == 1]
    require(len(remaining) == 2, (code, degree))
    edges.append((remaining[0], remaining[1]))
    return tuple(edges)


ALL_TREES = tuple(prufer_tree(code) for code in product(range(5), repeat=3))
require(len(ALL_TREES) == 125 and len(set(ALL_TREES)) == 125, len(ALL_TREES))


def union_budgets(masks: tuple[frozenset[int], ...]) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    require(len(masks) == 5 and all(0 in mask for mask in masks), masks)
    sizes = [len(mask) for mask in masks]
    common_zero = 1 + sum(size - 1 for size in sizes)
    # Maximum-weight spanning tree on the complete five-vertex intersection
    # graph.  Kruskal is exact here and avoids scanning all 125 labelled trees
    # in every row of the exhaustive mixed-period census.
    weighted_edges = sorted(
        ((len(masks[i] & masks[j]), i, j) for i in range(5) for j in range(i + 1, 5)),
        reverse=True,
    )
    parent = list(range(5))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    chosen: list[tuple[int, int]] = []
    best_credit = 0
    for weight, left, right in weighted_edges:
        left_root, right_root = root(left), root(right)
        if left_root == right_root:
            continue
        parent[left_root] = right_root
        chosen.append((left, right))
        best_credit += weight
        if len(chosen) == 4:
            break
    best_tree = tuple(chosen)
    require(len(best_tree) == 4, (masks, best_tree))
    tree_bound = sum(sizes) - best_credit
    union_size = len(set().union(*masks))
    require(union_size <= tree_bound <= common_zero,
            (sizes, union_size, tree_bound, common_zero, best_tree))
    return common_zero, tree_bound, best_tree


def check_tree_optimizer() -> int:
    """Cross-check Kruskal against all 125 trees on deterministic packets."""
    rows = 0
    for L in range(1, 25):
        options = mask_options(L)
        for offset in range(min(25, len(options) ** 2)):
            indices = tuple((offset * (j + 2) + j * j) % len(options) for j in range(5))
            masks = tuple(options[index][2] for index in indices)
            _common, tree_bound, _tree = union_budgets(masks)
            exhaustive_credit = max(
                sum(len(masks[i] & masks[j]) for i, j in tree)
                for tree in ALL_TREES
            )
            require(tree_bound == sum(map(len, masks)) - exhaustive_credit,
                    (L, indices, tree_bound, exhaustive_credit))
            rows += 1
    return rows


def mask_options(L: int) -> list[tuple[int, int, frozenset[int]]]:
    """One option for each distinct unit-mask at every quotient Q|L."""
    options: list[tuple[int, int, frozenset[int]]] = []
    for Q in divisors(L):
        seen: set[frozenset[int]] = set()
        for unit in range(1, Q + 1):
            if gcd(unit, Q) != 1:
                continue
            mask = lifted_mask(L, Q, unit)
            if mask not in seen:
                seen.add(mask)
                options.append((Q, unit, mask))
    return options


def check_master_identity() -> int:
    """Deterministic mixed rows check L=q/d0=lcm_i Q_i."""
    state = 0x1217C0DE
    rows = 0
    for q in range(1, 501):
        for _ in range(40):
            speeds: list[int] = []
            for _index in range(6):
                state = (1664525 * state + 1013904223) & 0xFFFFFFFF
                speeds.append(1 + state % (4 * q + 17))
            d0 = gcd_many([q, *speeds])
            periods = [q // gcd(q, speed) for speed in speeds]
            require(q // d0 == lcm_many(periods),
                    (q, speeds, d0, periods))
            rows += 1
    return rows


def check_pair_coincidence(limit: int = 400) -> tuple[int, int]:
    """Exhaust sum and positive-difference defining pairs."""
    sums = differences = 0
    for q in range(1, limit + 1):
        for left in range(1, q):
            right = q - left
            require(direct_q_mask(q, left) == direct_q_mask(q, right),
                    ("sum", q, left, right))
            sums += 1
        for left in range(1, 2 * q + 1):
            right = left + q
            require(direct_q_mask(q, left) == direct_q_mask(q, right),
                    ("difference", q, left, right))
            differences += 1
    return sums, differences


def check_tournament_gauge() -> dict[str, object]:
    """Verify the deterministic tie gauge used in the tournament audit.

    Intersection cardinality is symmetric, so the orientation contains no
    mathematical information; this check only makes the reported gauge
    fingerprints executable rather than leaving them as prose constants.
    """
    names = ("M1", "M2", "M3", "M4", "M5")
    periods = (2, 4, 4, 2, 4)
    order = tuple(sorted(range(5), key=lambda i: (periods[i], i)))
    rank = {vertex: index for index, vertex in enumerate(order)}
    edges = frozenset(
        (left, right) if rank[left] < rank[right] else (right, left)
        for left, right in combinations(range(5), 2)
    )

    def has_edge(left: int, right: int) -> bool:
        return (left, right) in edges

    scores = tuple(sorted(sum(has_edge(i, j) for j in range(5)) for i in range(5)))
    cycles = sum(
        (has_edge(i, j) and has_edge(j, k) and has_edge(k, i))
        or (has_edge(i, k) and has_edge(k, j) and has_edge(j, i))
        for i, j, k in combinations(range(5), 3)
    )
    hamilton_paths = sum(
        all(has_edge(path[i], path[i + 1]) for i in range(4))
        for path in permutations(range(5))
    )
    reachability = [
        {j for j in range(5) if i == j or has_edge(i, j)} for i in range(5)
    ]
    for middle in range(5):
        for source in range(5):
            if middle in reachability[source]:
                reachability[source] |= reachability[middle]
    sccs = len({
        tuple(j for j in range(5)
              if j in reachability[i] and i in reachability[j])
        for i in range(5)
    })
    reverse_edges = frozenset((right, left) for left, right in edges)
    reverse_flips = sum((right, left) in reverse_edges for left, right in edges)
    require(
        (order, scores, cycles, sccs, hamilton_paths, reverse_flips)
        == ((0, 3, 1, 2, 4), (0, 1, 2, 3, 4), 0, 5, 1, 10),
        (order, scores, cycles, sccs, hamilton_paths, reverse_flips),
    )
    return {
        "path": "->".join(names[i] for i in order),
        "scores": scores,
        "cycles": cycles,
        "sccs": sccs,
        "hamilton_paths": hamilton_paths,
        "reverse_flips": reverse_flips,
    }


def check_mixed_mask_space(limit: int = 48) -> dict[str, object]:
    """Exhaust unordered five-mask packets on every master clock L<=limit."""
    configs = proper = full = common_active = tree_active = 0
    tree_strict = run_strict = 0
    first_full: tuple[object, ...] | None = None
    first_full_without_Q1: tuple[object, ...] | None = None
    largest_options = (0, 0)
    for L in range(1, limit + 1):
        options = mask_options(L)
        if len(options) > largest_options[1]:
            largest_options = (L, len(options))
        for chosen in combinations_with_replacement(range(len(options)), 5):
            rows = [options[index] for index in chosen]
            periods = tuple(row[0] for row in rows)
            if lcm_many(periods) != L:
                continue
            masks = tuple(row[2] for row in rows)
            require(all(0 in mask for mask in masks), (L, rows))
            for Q, _unit, mask in rows:
                require(len(mask) == (L // Q) * window_count(Q),
                        (L, Q, len(mask)))
            common_bound, tree_bound, _tree = union_budgets(masks)
            union = set().union(*masks)
            configs += 1
            if len(union) == L:
                full += 1
                if first_full is None:
                    first_full = (L, periods, tuple(len(mask) for mask in masks))
                if 1 not in periods and first_full_without_Q1 is None:
                    first_full_without_Q1 = (
                        L, periods, tuple(len(mask) for mask in masks)
                    )
                continue
            proper += 1
            ell = longest_proper_cyclic_run(union, L)
            require(ell is not None and ell <= len(union), (L, union, ell))
            for start in range(L):
                require(any((start + step) % L not in union
                            for step in range(ell + 1)),
                        (L, union, ell, start))
            if common_bound < L:
                common_active += 1
                for start in range(L):
                    require(any((start + step) % L not in union
                                for step in range(common_bound + 1)),
                            ("common", L, common_bound, start))
            if tree_bound < L:
                tree_active += 1
                for start in range(L):
                    require(any((start + step) % L not in union
                                for step in range(tree_bound + 1)),
                            ("tree", L, tree_bound, start))
            tree_strict += tree_bound < common_bound
            run_strict += ell < tree_bound
    require(first_full is not None, "expected the Q=1 obstruction")
    return {
        "master_limit": limit,
        "configurations": configs,
        "proper": proper,
        "full": full,
        "common_active": common_active,
        "tree_active": tree_active,
        "tree_strict": tree_strict,
        "run_strict": run_strict,
        "first_full": first_full,
        "first_full_without_Q1": first_full_without_Q1,
        "largest_options": largest_options,
    }


def gap_block(a: int, k: int, q: int) -> tuple[int, int]:
    denominator = 14 * a
    lo = ceil_div(q * (14 * k + 1), denominator)
    hi = q * (14 * k + 13) // denominator
    return lo, hi


def check_full_no_q1_guardrail() -> dict[str, object]:
    """Realize the first full mixed clock with no period-one obligation."""
    q = 16
    speeds = (17, 35, 53, 71, 88, 104)
    require(speeds[-1] - speeds[-2] == q, speeds)
    g_values = tuple(gcd(speed, q) for speed in speeds)
    periods = tuple(q // g for g in g_values)
    require(periods == (16, 16, 16, 16, 2, 2), periods)
    all_masks = tuple(direct_q_mask(q, speed) for speed in speeds)
    require(all_masks[4] == all_masks[5], all_masks)
    masks = all_masks[:5]
    expected = (
        frozenset({0, 1, 15}),
        frozenset({0, 5, 11}),
        frozenset({0, 3, 13}),
        frozenset({0, 7, 9}),
        frozenset(range(0, 16, 2)),
    )
    require(masks == expected, masks)
    union = set().union(*masks)
    require(union == set(range(q)), union)
    common_bound, tree_bound, tree = union_budgets(masks)
    require(longest_proper_cyclic_run(union, q) is None, union)
    return {
        "q": q,
        "speeds": speeds,
        "periods": periods,
        "masks": tuple(tuple(sorted(mask)) for mask in masks),
        "common_bound": common_bound,
        "tree_bound": tree_bound,
        "tree": tree,
        "union": "full Z/16Z",
    }


def check_residual_row() -> dict[str, object]:
    a = 79
    speeds = (140, 210, 350, 420, 490, 770)
    q = 280
    k = 11
    require(speeds[-1] - speeds[-2] == q, speeds)
    g_values = tuple(gcd(speed, q) for speed in speeds)
    periods = tuple(q // value for value in g_values)
    d0 = gcd_many([q, *speeds])
    L = q // d0
    require((d0, L, periods) == (70, 4, (2, 4, 4, 2, 4, 4)),
            (d0, L, periods))
    require(lcm_many(periods) == L, periods)

    all_masks = tuple(
        lifted_mask(L, Q, speed // g)
        for speed, g, Q in zip(speeds, g_values, periods)
    )
    require(all_masks[4] == all_masks[5], all_masks)
    masks = all_masks[:5]
    expected_masks = (
        frozenset({0, 2}), frozenset({0}), frozenset({0}),
        frozenset({0, 2}), frozenset({0}),
    )
    require(masks == expected_masks, masks)
    union = set().union(*masks)
    require(union == {0, 2}, union)
    common_bound, tree_bound, _best_tree = union_budgets(masks)
    chosen_tree = ((0, 3), (0, 1), (0, 2), (0, 4))
    tree_credit = sum(len(masks[i] & masks[j]) for i, j in chosen_tree)
    require((sum(map(len, masks)), tree_credit, tree_bound) == (7, 5, 2),
            (sum(map(len, masks)), tree_credit, tree_bound))
    ell = longest_proper_cyclic_run(union, L)
    require((common_bound, tree_bound, ell) == (3, 2, 1),
            (common_bound, tree_bound, ell))

    block = gap_block(a, k, q)
    numerators = tuple(range(block[0], block[1] + 1))
    witnesses = tuple(p for p in numerators if p % L not in union)
    require((block, numerators, witnesses) == ((40, 42), (40, 41, 42), (41,)),
            (block, numerators, witnesses))
    witness = witnesses[0]

    packet = tuple(range(1, 7)) + (a,) + speeds
    least_residues = tuple(
        min((speed * witness) % q, (-speed * witness) % q)
        for speed in packet
    )
    require(len(packet) == 13 and all(14 * residue >= q for residue in least_residues),
            (packet, least_residues))
    clearances = tuple(Fraction(residue, q) for residue in least_residues)

    harmonic = a * sum((Fraction(1, speed) for speed in speeds), Fraction(0))
    require(harmonic > 1, harmonic)
    require(Fraction(speeds[0], a) < Fraction(13, 6), speeds[0])
    require(Fraction(d0, a) < Fraction(397, 432), d0)

    intersection_matrix = tuple(
        tuple(len(left & right) for right in masks) for left in masks
    )
    return {
        "a": a,
        "speeds": speeds,
        "q": q,
        "k": k,
        "g_values": g_values,
        "periods": periods,
        "d0": d0,
        "L": L,
        "masks": tuple(tuple(sorted(mask)) for mask in masks),
        "union": tuple(sorted(union)),
        "common_bound": common_bound,
        "tree_edges": chosen_tree,
        "tree_credit": tree_credit,
        "tree_bound": tree_bound,
        "run": ell,
        "block": numerators,
        "witness": witness,
        "packet": packet,
        "least_residues": least_residues,
        "clearances": clearances,
        "harmonic": harmonic,
        "first_ratio": Fraction(speeds[0], a),
        "gcd_ratio": Fraction(d0, a),
        "intersection_matrix": intersection_matrix,
    }


def main() -> None:
    identity_rows = check_master_identity()
    sum_rows, difference_rows = check_pair_coincidence()
    tree_optimizer_rows = check_tree_optimizer()
    census = check_mixed_mask_space()
    full_guardrail = check_full_no_q1_guardrail()
    row = check_residual_row()
    tournament = check_tournament_gauge()

    print("THM-1217 MIXED-PERIOD BEAT-MASK TREE EXACT REFEREE")
    print("method=integer/Fraction only; always-on checks; no dependencies")
    print(f"master_identity_rows={identity_rows}")
    print(f"pair_mask_rows: sum={sum_rows} difference={difference_rows}")
    print(f"kruskal_vs_all_125_trees_rows={tree_optimizer_rows}")
    print(
        "mixed_census: "
        f"L<=${census['master_limit']} configs={census['configurations']} "
        f"proper={census['proper']} full={census['full']}"
    .replace("$", ""))
    print(
        "mixed_suppliers: "
        f"common_C<L={census['common_active']} "
        f"tree_CT<L={census['tree_active']} "
        f"tree_strict={census['tree_strict']} "
        f"run_strict={census['run_strict']}"
    )
    print(f"largest_option_clock={census['largest_options']}")
    print(f"first_full_union={census['first_full']}")
    print(f"first_full_without_Q1={census['first_full_without_Q1']}")
    print("properness_note=full U has no finite escape run; C<L is required for C+1")
    print("realized_full_no_Q1_guardrail:")
    for key in ("q", "speeds", "periods", "masks", "common_bound", "tree_bound", "tree", "union"):
        print(f"  {key}={full_guardrail[key]}")
    print()
    print("CORRECTED a=79 RESIDUAL ROW")
    for key in (
        "a", "speeds", "q", "k", "g_values", "periods", "d0", "L",
        "masks", "union", "intersection_matrix", "common_bound",
        "tree_edges", "tree_credit", "tree_bound", "run", "block",
        "witness", "packet", "least_residues", "clearances", "harmonic",
        "first_ratio", "gcd_ratio",
    ):
        print(f"{key}={row[key]}")
    print("thresholds: B_common=4 B_tree=3 B_run=2")
    print("witness_time=41/280; all 13 clearances >= 1/14")
    print("frontier_checks: b1/a<13/6; d0/a<397/432; a*sum(1/bi)>1")
    print()
    print("TOURNAMENT / ALTERNATE-VERTEX AUDIT")
    print("observable=|Mi intersect Mj| (symmetric)")
    print(f"gauge_path={tournament['path']} (period then label)")
    print(
        f"scores={str(tournament['scores']).replace(' ', '')}; "
        f"cycles={tournament['cycles']}; SCCs={tournament['sccs']}; "
        f"Hamilton_paths={tournament['hamilton_paths']}; "
        f"reverse_flips={tournament['reverse_flips']}"
    )
    print("faithful_vertices=master residues + quotient clocks + lifted mask obligations")
    print("destroyed_by_tournament=mask placement, quotient maps, block phase, run length")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
