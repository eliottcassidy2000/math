#!/usr/bin/env python3
"""Fraction-exact replay of THM-779/MISTAKE-147's raw-wall refuter.

The audit keeps four structures separate.

1. Metric base: P={1,2,11,12,13} is 1/14-safe on
   I=[25/182,27/182].
2. Prime-seven stalk: W0={1,4,5,6,8,9,10} has the fixed token
   permutation (0,5,4,1,6,3,2) throughout I.
3. Clock: A_m=182m+1 has exactly 2m walls in I, all covered by the
   persistent W0 stalk; the resulting 13-speed family is divisor-complete
   for every modulus 2,...,14.  A=1000 is checked separately.
4. Normalized algebra: states with offset multiset {0,0,1,...,6} form an
   A8 torsor.  The two unit anchor moves generate A8; the directed state
   graph has 20,160 vertices, 40,320 edges, one SCC, and 5,760
   monochromatic seven-cycles.
5. Visitor de-phasing: the incoming one-window estimate is refuted exactly
   by (f,c,c')=(11,6,8); four successive paired co-visits cross order.  The
   triangle-width proof requires two fast-window lengths, giving the corrected
   factor-two bound in THM-783.

Tournament Analysis guardrail: endpoint walls, rather than runners, are the
only vertices on which strict temporal precedence is naturally binary.  The
pairwise observable is exact Fraction order, the switch/gauge reverses the
marked interval orientation, and chronological order is the tie Hamiltonian
path.  This tournament is transitive and therefore forgets the persistent
seven-token stalk.  The normalized state graph retains that stalk dynamics
but is not a tournament.  This explicitly challenges the assumption that
runner vertices or a clean owner tournament carry the obstruction.

Only the Python standard library is used.  Every geometric comparison is a
Fraction comparison; no floating point enters a verdict.
"""

from __future__ import annotations

import hashlib
import itertools
from collections import Counter, deque
from fractions import Fraction
from math import comb, factorial, floor


PRIME = 7
DELTA = Fraction(1, 14)
P = (1, 2, 11, 12, 13)
W0 = (1, 4, 5, 6, 8, 9, 10)
I_LEFT = Fraction(25, 182)
I_RIGHT = Fraction(27, 182)
I_CENTER = Fraction(1, 7)
EXPECTED_STALK = (0, 5, 4, 1, 6, 3, 2)
EXPECTED_MULTISET = (0, 0, 1, 2, 3, 4, 5, 6)


def mod1(x: Fraction) -> Fraction:
    return x - floor(x)


def circle_norm(x: Fraction) -> Fraction:
    y = mod1(x)
    return min(y, 1 - y)


def ceil_fraction(x: Fraction) -> int:
    return -floor(-x)


def nearest_integer_off_wall(x: Fraction) -> int:
    """Nearest integer when x is not a half-integer."""
    assert x.denominator != 2 or x.numerator % 2 == 0
    return floor(x + Fraction(1, 2))


def token(w: int, x: Fraction) -> int:
    assert w % PRIME
    return (-pow(w, -1, PRIME) * nearest_integer_off_wall(w * x)) % PRIME


def bad_sheets(w: int, x: Fraction) -> tuple[int, ...]:
    return tuple(
        k
        for k in range(PRIME)
        if circle_norm(Fraction(w, PRIME) * (x + k)) < DELTA
    )


def covers_all_sheets(speeds: tuple[int, ...], x: Fraction) -> bool:
    return set().union(*(set(bad_sheets(w, x)) for w in speeds)) == set(range(PRIME))


def endpoint(w: int, j: int) -> Fraction:
    return Fraction(2 * j + 1, 2 * w)


def wall_index_bounds(w: int) -> tuple[int, int]:
    first = ceil_fraction(w * I_LEFT - Fraction(1, 2))
    last = floor(w * I_RIGHT - Fraction(1, 2))
    return max(0, first), min(w - 1, last)


def walls_in_interval(w: int) -> tuple[tuple[int, Fraction], ...]:
    first, last = wall_index_bounds(w)
    if first > last:
        return ()
    rows = tuple((j, endpoint(w, j)) for j in range(first, last + 1))
    assert all(I_LEFT <= x <= I_RIGHT for _, x in rows)
    if first:
        assert endpoint(w, first - 1) < I_LEFT
    if last + 1 < w:
        assert endpoint(w, last + 1) > I_RIGHT
    return rows


def interval_clearance(w: int) -> Fraction:
    """Exact minimum of ||w x|| on I.

    Distance to the integers has no interior local minimum except zero, so
    endpoints suffice unless an integer lies in wI.
    """
    left = w * I_LEFT
    right = w * I_RIGHT
    if ceil_fraction(left) <= floor(right):
        return Fraction(0)
    return min(circle_norm(left), circle_norm(right))


def stalk_audit() -> dict:
    clearances = {w: interval_clearance(w) for w in P}
    assert min(clearances.values()) == DELTA
    assert all(value >= DELTA for value in clearances.values())

    all_walls = sorted(
        (endpoint(w, j), w)
        for w in W0
        for j in range(w)
    )
    left_wall = max(x for x, _ in all_walls if x < I_CENTER)
    right_wall = min(x for x, _ in all_walls if x > I_CENTER)
    assert (left_wall, right_wall) == (Fraction(1, 8), Fraction(3, 20))
    assert left_wall < I_LEFT < I_RIGHT < right_wall
    assert all(not walls_in_interval(w) for w in W0)

    probes = (
        I_LEFT,
        (3 * I_LEFT + I_RIGHT) / 4,
        I_CENTER,
        (I_LEFT + 3 * I_RIGHT) / 4,
        I_RIGHT,
    )
    rows = []
    for x in probes:
        direct = tuple(bad_sheets(w, x) for w in W0)
        tokens = tuple(token(w, x) for w in W0)
        assert all(len(row) == 1 for row in direct)
        assert tuple(row[0] for row in direct) == tokens == EXPECTED_STALK
        assert set(tokens) == set(range(PRIME))
        rows.append((x, tokens))

    return {
        "clearances": clearances,
        "stalk_chamber": (left_wall, right_wall),
        "probe_rows": tuple(rows),
    }


def divisor_witnesses(family: tuple[int, ...]) -> dict[int, int]:
    answer = {}
    for modulus in range(2, 15):
        witnesses = tuple(w for w in family if w % modulus == 0)
        assert witnesses
        answer[modulus] = min(witnesses)
    return answer


def direct_family_audit(m: int) -> dict:
    a = 182 * m + 1
    walls = walls_in_interval(a)
    expected_indices = tuple(range(25 * m, 27 * m))
    assert tuple(j for j, _ in walls) == expected_indices
    assert len(walls) == 2 * m
    assert a % PRIME == 1

    for _, x in walls:
        stalk_rows = tuple(bad_sheets(w, x) for w in W0)
        assert tuple(row[0] for row in stalk_rows) == EXPECTED_STALK
        assert bad_sheets(a, x) == ()
        assert covers_all_sheets(W0 + (a,), x)

    cuts = (I_LEFT,) + tuple(x for _, x in walls) + (I_RIGHT,)
    chamber_probes = tuple((left + right) / 2 for left, right in zip(cuts, cuts[1:]))
    for x in chamber_probes:
        assert tuple(token(w, x) for w in W0) == EXPECTED_STALK
        assert covers_all_sheets(W0 + (a,), x)
    assert covers_all_sheets(W0 + (a,), I_LEFT)
    assert covers_all_sheets(W0 + (a,), I_RIGHT)

    family = tuple(sorted(tuple(PRIME * p for p in P) + W0 + (a,)))
    assert len(family) == 13 == len(set(family))
    assert family[0] == 1
    witnesses = divisor_witnesses(family)
    return {
        "m": m,
        "A": a,
        "first_wall": walls[0],
        "last_wall": walls[-1],
        "walls": len(walls),
        "chambers": len(chamber_probes),
        "family": family,
        "divisor_witnesses": witnesses,
    }


def count_family_audit(limit: int = 1000) -> dict:
    total = 0
    for m in range(1, limit + 1):
        a = 182 * m + 1
        first, last = wall_index_bounds(a)
        assert (first, last) == (25 * m, 27 * m - 1)
        total += last - first + 1
    assert total == limit * (limit + 1)
    return {"m_range": (1, limit), "families": limit, "total_walls": total}


def a1000_audit() -> dict:
    a = 1000
    assert a % PRIME
    walls = walls_in_interval(a)
    assert tuple(j for j, _ in walls) == tuple(range(137, 148))
    assert len(walls) == 11
    for _, x in walls:
        assert bad_sheets(a, x) == ()
        assert covers_all_sheets(W0 + (a,), x)
    cuts = (I_LEFT,) + tuple(x for _, x in walls) + (I_RIGHT,)
    assert all(covers_all_sheets(W0 + (a,), (x + y) / 2) for x, y in zip(cuts, cuts[1:]))
    family = tuple(sorted(tuple(PRIME * p for p in P) + W0 + (a,)))
    witnesses = divisor_witnesses(family)
    return {
        "A": a,
        "indices": (walls[0][0], walls[-1][0]),
        "walls": len(walls),
        "first_wall": walls[0][1],
        "last_wall": walls[-1][1],
        "divisor_witnesses": witnesses,
    }


def fast_period_index(fast: int, time: Fraction) -> int:
    """Index the open interval between consecutive fast-owner midpoint walls."""
    shifted = fast * time - Fraction(1, 2)
    assert shifted.denominator != 1
    return floor(shifted)


def dephase_factor_two_audit() -> dict:
    """Exact counterexample to the one-window de-phase estimate.

    Each listed c- and d-wall advances to the next wall, and both lie in the
    same f-period.  Their signed separation crosses zero, so the available
    separation interval has width 2/f rather than 1/f.
    """
    fast, c, d = 11, 6, 8
    pairs = ((1, 2), (2, 3), (3, 4), (4, 5))
    rows = []
    for m, n in pairs:
        tc = endpoint(c, m)
        td = endpoint(d, n)
        cell_c = fast_period_index(fast, tc)
        cell_d = fast_period_index(fast, td)
        assert cell_c == cell_d
        rows.append((m, n, tc, td, td - tc, cell_c))

    assert tuple(row[4] for row in rows) == (
        Fraction(1, 16), Fraction(1, 48),
        Fraction(-1, 48), Fraction(-1, 16),
    )
    assert (c + d) % PRIME == 0  # the pair is cluster-balanced mod 7
    length = len(rows)
    delta = d - c
    old_bound = Fraction(c * d, fast * delta) + 1
    corrected_bound = Fraction(2 * c * d, fast * delta) + 1
    assert length > old_bound
    assert length < corrected_bound
    return {
        "speeds": (fast, c, d),
        "paired_co_visits": tuple(rows),
        "length": length,
        "old_bound": old_bound,
        "corrected_bound": corrected_bound,
        "pair_balanced_mod7": True,
    }


def period_index_or_none(fast: int, time: Fraction) -> int | None:
    """Return the fast-period index, or None at a simultaneous fast wall."""
    shifted = fast * time - Fraction(1, 2)
    if shifted.denominator == 1:
        return None
    return floor(shifted)


def longest_circular_true(pattern: tuple[bool, ...]) -> int:
    current = best = 0
    for value in pattern + pattern:
        current = current + 1 if value else 0
        best = max(best, min(current, len(pattern)))
    return best


def fixed_companion_span_audit(limit: int = 50) -> dict:
    """Exhaust the corrected THM-786 fixed-companion span for small triples.

    A g-wall is served when some c-wall lies in the same open f-period.  The
    pattern is one-periodic; the g midpoint walls in [0,1) give its full word.
    """
    checked = 0
    old_factor_one_failures = []
    balanced_checked = 0
    balanced_old_failures = []
    largest_run = (0, None)
    for fast in range(3, limit + 1):
        for g in range(2, fast):
            for c in range(1, g):
                checked += 1
                c_cells = {
                    cell
                    for m in range(-1, c + 1)
                    if (cell := period_index_or_none(fast, endpoint(c, m))) is not None
                }
                pattern = tuple(
                    (cell := period_index_or_none(fast, endpoint(g, j))) is not None
                    and cell in c_cells
                    for j in range(g)
                )
                length = longest_circular_true(pattern)
                delta = g - c
                corrected = Fraction(2 * g * c, fast * delta) + 1
                assert length < corrected
                old_integer_bound = (g * c) // (fast * delta) + 1
                if length > old_integer_bound:
                    row = (fast, g, c, length, old_integer_bound, corrected)
                    old_factor_one_failures.append(row)
                if fast % PRIME and g % PRIME and c % PRIME and (g + c) % PRIME == 0:
                    balanced_checked += 1
                    if length > old_integer_bound:
                        balanced_old_failures.append(
                            (fast, g, c, length, old_integer_bound, corrected)
                        )
                if length > largest_run[0]:
                    largest_run = (length, (fast, g, c))

    assert checked == comb(limit, 3)
    assert old_factor_one_failures
    assert any(row[:4] == (11, 8, 6, 4) for row in old_factor_one_failures)
    assert balanced_old_failures
    assert any(row[:4] == (11, 8, 6, 4) for row in balanced_old_failures)
    return {
        "limit": limit,
        "triples": checked,
        "corrected_failures": 0,
        "old_factor_one_failures": len(old_factor_one_failures),
        "first_old_failure": old_factor_one_failures[0],
        "balanced_triples": balanced_checked,
        "balanced_old_failures": len(balanced_old_failures),
        "first_balanced_old_failure": balanced_old_failures[0],
        "largest_run": largest_run,
    }


Permutation = tuple[int, ...]
State = tuple[int, ...]


def identity_permutation() -> Permutation:
    return tuple(range(8))


def cycle_permutation(cycle: tuple[int, ...]) -> Permutation:
    result = list(range(8))
    for a, b in zip(cycle, cycle[1:] + cycle[:1]):
        result[a] = b
    return tuple(result)


def compose(left: Permutation, right: Permutation) -> Permutation:
    """left after right."""
    return tuple(left[right[i]] for i in range(8))


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * 8
    for i, image in enumerate(permutation):
        result[image] = i
    return tuple(result)


def permutation_order(permutation: Permutation) -> int:
    value = identity_permutation()
    for order in range(1, 100):
        value = compose(permutation, value)
        if value == identity_permutation():
            return order
    raise AssertionError("order bound exceeded")


def is_even(permutation: Permutation) -> bool:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(8)
        for j in range(i + 1, 8)
    )
    return inversions % 2 == 0


def generated_group(generators: tuple[Permutation, ...]) -> set[Permutation]:
    moves = generators + tuple(inverse(g) for g in generators)
    seen = {identity_permutation()}
    queue = deque(seen)
    while queue:
        value = queue.popleft()
        for move in moves:
            nxt = compose(move, value)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return seen


def act_on_state(permutation: Permutation, state: State) -> State:
    result = [0] * 8
    for owner, image in enumerate(permutation):
        result[image] = state[owner]
    return tuple(result)


def normalized_states() -> tuple[State, ...]:
    states = []
    for zero_pair in itertools.combinations(range(8), 2):
        nonzero_owners = tuple(a for a in range(8) if a not in zero_pair)
        for residues in itertools.permutations(range(1, 7)):
            state = [0] * 8
            for owner, residue in zip(nonzero_owners, residues):
                state[owner] = residue
            states.append(tuple(state))
    return tuple(states)


def anchor_move(state: State, owner: int) -> State:
    assert state[owner] == 0
    result = tuple(
        0 if b == owner else (residue + 1) % PRIME
        for b, residue in enumerate(state)
    )
    assert tuple(sorted(result)) == EXPECTED_MULTISET
    return result


def reachable_count(graph: list[list[int]], root: int) -> int:
    seen = {root}
    queue = deque((root,))
    while queue:
        vertex = queue.popleft()
        for nxt in graph[vertex]:
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return len(seen)


def a8_and_graph_audit(digest: "hashlib._Hash") -> dict:
    p = cycle_permutation((1, 2, 3, 4, 5, 6, 7))
    q = cycle_permutation((0, 2, 3, 4, 5, 6, 7))
    quotient = compose(p, inverse(q))
    assert quotient == cycle_permutation((0, 1, 2))
    assert permutation_order(p) == permutation_order(q) == 7
    assert is_even(p) and is_even(q) and is_even(quotient)

    group = generated_group((p, q))
    even_count = sum(is_even(permutation) for permutation in itertools.permutations(range(8)))
    assert len(group) == even_count == factorial(8) // 2 == 20160
    assert all(is_even(permutation) for permutation in group)

    states = normalized_states()
    state_set = set(states)
    assert len(states) == len(state_set) == comb(8, 2) * factorial(6) == 20160
    base = EXPECTED_MULTISET
    orbit = {act_on_state(permutation, base) for permutation in group}
    assert orbit == state_set

    # The inverses of the displayed seven-cycles are exactly the two anchor
    # moves enabled at the base state.
    assert act_on_state(inverse(p), base) == anchor_move(base, 0)
    assert act_on_state(inverse(q), base) == anchor_move(base, 1)

    index = {state: i for i, state in enumerate(states)}
    adjacency = [[] for _ in states]
    reverse = [[] for _ in states]
    indegrees = [0] * len(states)
    edges = []
    for source, state in enumerate(states):
        enabled = tuple(a for a, residue in enumerate(state) if residue == 0)
        assert len(enabled) == 2
        for owner in enabled:
            target = index[anchor_move(state, owner)]
            adjacency[source].append(target)
            reverse[target].append(source)
            indegrees[target] += 1
            edges.append((source, owner, target))
            digest.update(f"E|{state}|{owner}|{states[target]}\n".encode())

    assert len(edges) == 40320
    assert Counter(map(len, adjacency)) == {2: 20160}
    assert Counter(indegrees) == {2: 20160}
    root = index[base]
    assert reachable_count(adjacency, root) == len(states)
    assert reachable_count(reverse, root) == len(states)

    cycle_histogram = Counter()
    cycles_by_owner = []
    for owner in range(8):
        unseen = {i for i, state in enumerate(states) if state[owner] == 0}
        count = 0
        while unseen:
            start = min(unseen)
            current = start
            cycle = []
            while current not in cycle:
                cycle.append(current)
                current = index[anchor_move(states[current], owner)]
            assert current == start
            assert all(vertex in unseen for vertex in cycle)
            unseen.difference_update(cycle)
            cycle_histogram[len(cycle)] += 1
            count += 1
        assert count == 720
        cycles_by_owner.append(count)
    assert cycle_histogram == {7: 5760}

    for permutation in sorted(group):
        digest.update(f"G|{permutation}\n".encode())

    return {
        "generator_orders": (permutation_order(p), permutation_order(q)),
        "quotient": quotient,
        "generated_group": len(group),
        "even_permutations": even_count,
        "torsor_orbit": len(orbit),
        "states": len(states),
        "directed_edges": len(edges),
        "indegree_histogram": dict(Counter(indegrees)),
        "outdegree_histogram": dict(Counter(map(len, adjacency))),
        "forward_reach": reachable_count(adjacency, root),
        "reverse_reach": reachable_count(reverse, root),
        "scc_count": 1,
        "monochromatic_cycles_by_owner": tuple(cycles_by_owner),
        "cycle_length_histogram": dict(cycle_histogram),
    }


def tournament_guardrail(m: int = 5) -> dict:
    """Fingerprint the wall-event tournament, then state its loss."""
    walls = walls_in_interval(182 * m + 1)
    n = len(walls)
    assert n == 2 * m
    # Earlier event points to later event.  Reversing the marked interval cut
    # reverses every edge.
    scores = tuple(n - 1 - i for i in range(n))
    score_histogram = dict(sorted(Counter(scores).items()))
    assert score_histogram == {score: 1 for score in range(n)}
    assert all(walls[i][1] < walls[j][1] for i in range(n) for j in range(i + 1, n))
    return {
        "vertices": "A_m wall events in I (not runners)",
        "event_count": n,
        "pairwise_observable": "strict Fraction time precedence",
        "switch_gauge": "reverse marked interval orientation",
        "tie_hamiltonian_path": tuple(j for j, _ in walls),
        "score_histogram": score_histogram,
        "directed_3cycles": 0,
        "scc_sizes": (1,) * n,
        "edge_flips_under_switch": comb(n, 2),
        "hamiltonian_paths": 1,
        "owner_switches": 0,
        "guardrail": "owner quotient has one active vertex; tournament forgets the exact W0 stalk",
    }


def render_mapping(mapping: dict) -> str:
    return "{" + ", ".join(f"{key}: {value}" for key, value in mapping.items()) + "}"


def main() -> None:
    digest = hashlib.sha256()
    stalk = stalk_audit()
    count_audit = count_family_audit()
    detailed = tuple(direct_family_audit(m) for m in (1, 2, 3, 5, 10, 25, 100))
    thousand = a1000_audit()
    dephase = dephase_factor_two_audit()
    serving_span = fixed_companion_span_audit()

    for label, payload in (
        ("stalk", stalk),
        ("count", count_audit),
        ("families", detailed),
        ("A1000", thousand),
        ("dephase", dephase),
        ("serving_span", serving_span),
    ):
        digest.update(f"{label}|{payload}\n".encode())

    group_graph = a8_and_graph_audit(digest)
    tournament = tournament_guardrail()
    digest.update(f"tournament|{tournament}\n".encode())

    print("LRC14 r=8 RAW-WALL REFUTER / A8 STATE GRAPH — EXACT REPLAY")
    print("=" * 78)
    print("METHOD")
    print("  arithmetic: fractions.Fraction only; all wall and clearance tests exact")
    print("  geometry: direct strict bad-sheet sets, not only token formulas")
    print("  group: explicit permutation closure and full normalized-state graph")
    print()
    print("CORE-SAFE BASE AND PERSISTENT SEVEN-STALK")
    print(f"  P={P}  I=[{I_LEFT},{I_RIGHT}]  center={I_CENTER}  delta={DELTA}")
    print(f"  exact interval clearances: {render_mapping(stalk['clearances'])}")
    print(f"  minimum clearance: {min(stalk['clearances'].values())}  safe=True")
    print(f"  enclosing W0 chamber: {stalk['stalk_chamber']}")
    print(f"  W0={W0}")
    print(f"  token permutation at all {len(stalk['probe_rows'])} probes: {EXPECTED_STALK}")
    print("  W0 covers F7 throughout I; W0 walls in I: 0")
    print()
    print("INFINITE RAW-WALL FAMILY")
    print(
        f"  symbolic index audit m={count_audit['m_range'][0]}..{count_audit['m_range'][1]}: "
        f"families={count_audit['families']} total_walls={count_audit['total_walls']}"
    )
    print("  A_m=182m+1; exact indices 25m..27m-1; every wall/chamber covered")
    for row in detailed:
        print(
            f"  m={row['m']:3d} A={row['A']:5d} walls={row['walls']:3d} "
            f"j=[{row['first_wall'][0]},{row['last_wall'][0]}] "
            f"x=[{row['first_wall'][1]},{row['last_wall'][1]}] owners={len(row['family'])}"
        )
    witnesses = detailed[0]["divisor_witnesses"]
    print(f"  divisor witnesses 2..14 (independent of m): {render_mapping(witnesses)}")
    print("  primitive=True; distinct owners=13; divisor_complete_2_through_14=True")
    print()
    print("CONCRETE A=1000 ROW")
    print(
        f"  indices={thousand['indices']} walls={thousand['walls']} "
        f"x=[{thousand['first_wall']},{thousand['last_wall']}] all_covered=True"
    )
    print()
    print("DE-PHASE / SERVING FACTOR-TWO AUDIT (THM-783/786 SCOPE CORRECTION)")
    print(
        f"  (fast,c,c')={dephase['speeds']} paired_co_visits={dephase['length']} "
        f"pair_balanced_mod7={dephase['pair_balanced_mod7']}"
    )
    for m, n, tc, td, difference, cell in dephase["paired_co_visits"]:
        print(
            f"  (m,n)=({m},{n}) times=({tc},{td}) signed_difference={difference} "
            f"fast_period={cell}"
        )
    print(
        f"  old_one-window_bound={dephase['old_bound']} refuted; "
        f"corrected_two-window_bound={dephase['corrected_bound']}"
    )
    print(
        f"  exhaustive f<= {serving_span['limit']}: triples={serving_span['triples']} "
        f"corrected_factor_two_failures={serving_span['corrected_failures']} "
        f"old_factor_one_failures={serving_span['old_factor_one_failures']}"
    )
    first = serving_span["first_old_failure"]
    first_balanced = serving_span["first_balanced_old_failure"]
    print(
        f"  first_old_failure={(first[0], first[1], first[2], first[3])} "
        f"old_bound={first[4]} corrected_bound={first[5]}"
    )
    print(
        f"  lens-balanced triples={serving_span['balanced_triples']} "
        f"old_factor_one_failures={serving_span['balanced_old_failures']} "
        f"first={(first_balanced[0], first_balanced[1], first_balanced[2], first_balanced[3])} "
        f"old_bound={first_balanced[4]} corrected_bound={first_balanced[5]} "
        f"largest_small_run={serving_span['largest_run']}"
    )
    print()
    print("A8 TORSOR AND NORMALIZED COVERED-STATE GRAPH")
    print(
        "  generators: p=(1 2 3 4 5 6 7), q=(0 2 3 4 5 6 7); "
        f"orders={group_graph['generator_orders']}"
    )
    print(f"  p*q^-1={group_graph['quotient']} = cycle (0 1 2)")
    print(
        f"  generated group={group_graph['generated_group']} "
        f"all even permutations={group_graph['even_permutations']} torsor orbit={group_graph['torsor_orbit']}"
    )
    print(
        f"  states={group_graph['states']} directed_edges={group_graph['directed_edges']} "
        f"in={group_graph['indegree_histogram']} out={group_graph['outdegree_histogram']}"
    )
    print(
        f"  forward/reverse reach={group_graph['forward_reach']}/{group_graph['reverse_reach']} "
        f"SCCs={group_graph['scc_count']}"
    )
    print(
        f"  monochromatic cycles by owner={group_graph['monochromatic_cycles_by_owner']} "
        f"length histogram={group_graph['cycle_length_histogram']}"
    )
    print()
    print("TOURNAMENT ANALYSIS / ASSUMPTION CHALLENGE")
    print(f"  vertices: {tournament['vertices']}")
    print(f"  observable: {tournament['pairwise_observable']}")
    print(f"  switch/gauge: {tournament['switch_gauge']}")
    print(f"  tie Hamiltonian path: {tournament['tie_hamiltonian_path']}")
    print(
        f"  fingerprint: scores={tournament['score_histogram']} cycles3={tournament['directed_3cycles']} "
        f"SCCs={tournament['scc_sizes']} flips={tournament['edge_flips_under_switch']} "
        f"Hamiltonian_paths={tournament['hamiltonian_paths']}"
    )
    print(f"  no-clean-owner-switch guardrail: {tournament['guardrail']}")
    print("  challenged vertices: runners -> wall events -> normalized token states")
    print("  preserved by state graph: duplicate root, anchor choice, token translation")
    print("  destroyed by tournament quotient: stalk identity, residues, metric core interval")
    print()
    print(f"DETERMINISTIC AUDIT SHA256: {digest.hexdigest()}")
    print("VERDICT: exact PASS; raw covered-wall length is unbounded, while owner-switch complexity is zero")


if __name__ == "__main__":
    main()
