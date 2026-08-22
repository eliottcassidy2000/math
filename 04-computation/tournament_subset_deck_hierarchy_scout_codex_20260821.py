#!/usr/bin/env python3
"""Exact ordinary-tournament subset/deck hierarchy scout (Codex, 2026-08-21).

Universe
========
For every ``1 <= n <= 6`` this program enumerates all
``2**binom(n, 2)`` *labelled* orientations of ``K_n``.  It then forms the
unlabelled isomorphism quotient by an independent orbit sweep.  No random or
external data enter the load-bearing census.

The scout compares progressively richer, genuinely unlabelled summaries:

* induced k-subtournament type profiles;
* vertex-rooted triangle and four-set ("square") incidence profiles;
* k-profile deletion decks and the full vertex-deleted deck;
* ranks over Q and F_2 of vertex/type incidence and deletion matrices;
* card-boundary pairings, especially (deleted card, missing outdegree).

It also checks the exact double-counting identities behind these summaries.
Every truth-bearing gate uses ``require`` rather than ``assert``, so normal and
``python -O`` executions have identical logical content and identical output.

Reproduction
------------
  python 04-computation/tournament_subset_deck_hierarchy_scout_codex_20260821.py
  python -O 04-computation/tournament_subset_deck_hierarchy_scout_codex_20260821.py
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, factorial


MAX_N = 6
EXPECTED_CLASS_COUNTS = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56}


def require(condition: bool, message: str) -> None:
    """Optimization-safe validation gate."""
    if not condition:
        raise RuntimeError(message)


PAIR_LIST: dict[int, tuple[tuple[int, int], ...]] = {}
PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}
ORDERS: dict[int, tuple[tuple[int, ...], ...]] = {}
CANON: dict[int, dict[int, int]] = {}
CLASS_REPS: dict[int, tuple[int, ...]] = {}
ORBIT_SIZE: dict[int, dict[int, int]] = {}
UNROOTED_TYPES: dict[int, tuple[int, ...]] = {}
UNROOTED_INDEX: dict[int, dict[int, int]] = {}
ROOTED_TYPES: dict[int, tuple[int, ...]] = {}
ROOTED_INDEX: dict[int, dict[int, int]] = {}


def prepare_static(n: int) -> None:
    pairs = tuple(combinations(range(n), 2))
    PAIR_LIST[n] = pairs
    PAIR_INDEX[n] = {pair: i for i, pair in enumerate(pairs)}
    ORDERS[n] = tuple(permutations(range(n)))


def has_edge(bits: int, n: int, u: int, v: int) -> bool:
    """Whether u -> v; bit 1 means i -> j for i < j."""
    if u == v:
        return False
    if u < v:
        return bool((bits >> PAIR_INDEX[n][(u, v)]) & 1)
    return not bool((bits >> PAIR_INDEX[n][(v, u)]) & 1)


def relabel(bits: int, n: int, order: tuple[int, ...]) -> int:
    """New vertex i is old vertex order[i]."""
    out = 0
    for index, (i, j) in enumerate(PAIR_LIST[n]):
        if has_edge(bits, n, order[i], order[j]):
            out |= 1 << index
    return out


def prepare_isomorphism_classes(n: int) -> None:
    """Exhaustive labelled orbit sweep, producing a total canonical map."""
    total = 1 << len(PAIR_LIST[n])
    canon: dict[int, int] = {}
    reps: list[int] = []
    orbit_sizes: dict[int, int] = {}
    for bits in range(total):
        if bits in canon:
            continue
        orbit = {relabel(bits, n, order) for order in ORDERS[n]}
        rep = min(orbit)
        require(rep == bits, f"orbit sweep lost least representative n={n} bits={bits}")
        for image in orbit:
            old = canon.setdefault(image, rep)
            require(old == rep, f"overlapping isomorphism orbits n={n}")
        reps.append(rep)
        orbit_sizes[rep] = len(orbit)
    require(len(canon) == total, f"canonical map incomplete n={n}")
    require(len(reps) == EXPECTED_CLASS_COUNTS[n], f"class-count control failed n={n}")
    require(sum(orbit_sizes.values()) == total, f"orbit-size sum failed n={n}")
    for rep, size in orbit_sizes.items():
        require(factorial(n) % size == 0, f"bad orbit divisor n={n} rep={rep}")
    CANON[n] = canon
    CLASS_REPS[n] = tuple(reps)
    ORBIT_SIZE[n] = orbit_sizes
    UNROOTED_TYPES[n] = tuple(reps)
    UNROOTED_INDEX[n] = {rep: i for i, rep in enumerate(reps)}


@lru_cache(maxsize=None)
def rooted_canon(bits: int, n: int, root: int) -> int:
    """Canonical rooted code with the distinguished vertex placed at 0."""
    others = tuple(v for v in range(n) if v != root)
    return min(relabel(bits, n, (root, *tail)) for tail in permutations(others))


def prepare_rooted_types(n: int) -> None:
    types = {
        rooted_canon(rep, n, root)
        for rep in CLASS_REPS[n]
        for root in range(n)
    }
    ROOTED_TYPES[n] = tuple(sorted(types))
    ROOTED_INDEX[n] = {code: i for i, code in enumerate(ROOTED_TYPES[n])}


def prepare() -> None:
    # n=7 gets only a pair table for the explicitly declared targeted family;
    # no n=7 isomorphism quotient or exhaustive tournament census is built.
    for n in range(1, MAX_N + 2):
        prepare_static(n)
    for n in range(1, MAX_N + 1):
        prepare_isomorphism_classes(n)
    for n in range(1, MAX_N + 1):
        prepare_rooted_types(n)


def induced(bits: int, n: int, vertices: tuple[int, ...]) -> int:
    out = 0
    for index, (i, j) in enumerate(combinations(range(len(vertices)), 2)):
        if has_edge(bits, n, vertices[i], vertices[j]):
            out |= 1 << index
    return out


def delete_vertex(bits: int, n: int, vertex: int) -> int:
    return induced(bits, n, tuple(v for v in range(n) if v != vertex))


def converse(bits: int, n: int) -> int:
    return bits ^ ((1 << len(PAIR_LIST[n])) - 1)


def is_converse_pair(a: int, b: int, n: int) -> bool:
    return CANON[n][converse(a, n)] == b


@lru_cache(maxsize=None)
def outdegrees(bits: int, n: int) -> tuple[int, ...]:
    return tuple(sum(has_edge(bits, n, v, w) for w in range(n)) for v in range(n))


@lru_cache(maxsize=None)
def score_sequence(bits: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(outdegrees(bits, n)))


def is_cyclic_triangle(bits: int, n: int, triple: tuple[int, int, int]) -> bool:
    sub = induced(bits, n, triple)
    return score_sequence(sub, 3) == (1, 1, 1)


@lru_cache(maxsize=None)
def c3_count(bits: int, n: int) -> int:
    return sum(is_cyclic_triangle(bits, n, triple) for triple in combinations(range(n), 3))


def four_type_indices() -> dict[str, int]:
    by_score = {
        score_sequence(rep, 4): UNROOTED_INDEX[4][rep]
        for rep in UNROOTED_TYPES[4]
    }
    return {
        "transitive": by_score[(0, 1, 2, 3)],
        "source_c3": by_score[(1, 1, 1, 3)],
        "c4": by_score[(1, 1, 2, 2)],
        "sink_c3": by_score[(0, 2, 2, 2)],
    }


def four_profile_named(bits: int, n: int) -> tuple[int, int, int, int]:
    """Counts (q=0, source-C3, q=2/C4, sink-C3) among four-sets."""
    indices = four_type_indices()
    profile = induced_profile(bits, n, 4)
    return (
        profile[indices["transitive"]],
        profile[indices["source_c3"]],
        profile[indices["c4"]],
        profile[indices["sink_c3"]],
    )


def four_chirality(bits: int, n: int) -> int:
    _n0, source, _c4, sink = four_profile_named(bits, n)
    return source - sink


def four_chirality_from_scores(bits: int, n: int) -> int:
    return sum(comb(d, 3) - comb(n - 1 - d, 3) for d in outdegrees(bits, n))


def validate_four_type_qh() -> dict[tuple[int, int], int]:
    qh_to_type: dict[tuple[int, int], int] = {}
    for rep in UNROOTED_TYPES[4]:
        q = c3_count(rep, 4)
        degrees = outdegrees(rep, 4)
        h = degrees.count(3) - degrees.count(0)
        key = (q, h)
        require(key not in qh_to_type, f"(q,h) fails to classify 4-types at {key}")
        qh_to_type[key] = UNROOTED_INDEX[4][rep]
    require(set(qh_to_type) == {(0, 0), (1, 1), (2, 0), (1, -1)}, "bad (q,h) atlas")
    return qh_to_type


def maximal_flag_signature(bits: int, n: int) -> tuple[tuple[tuple[int, ...], int], ...]:
    """Distribution of proper induced-type sequences along all maximal flags.

    The final n-type is omitted: retaining it would identify the input by
    definition.  The signature records gluing/alignment of the separate
    k-profiles across one nested chain.
    """
    counts: Counter[tuple[int, ...]] = Counter()
    for order in ORDERS[n]:
        sequence = tuple(
            UNROOTED_INDEX[k][CANON[k][induced(bits, n, order[:k])]]
            for k in range(3, n)
        )
        counts[sequence] += 1
    return tuple(sorted(counts.items()))


@lru_cache(maxsize=None)
def induced_profile(bits: int, n: int, k: int) -> tuple[int, ...]:
    counts = [0] * len(UNROOTED_TYPES[k])
    for vertices in combinations(range(n), k):
        code = CANON[k][induced(bits, n, vertices)]
        counts[UNROOTED_INDEX[k][code]] += 1
    return tuple(counts)


@lru_cache(maxsize=None)
def unrooted_incidence_rows(bits: int, n: int, k: int) -> tuple[tuple[int, ...], ...]:
    rows = [[0] * len(UNROOTED_TYPES[k]) for _ in range(n)]
    for vertices in combinations(range(n), k):
        sub = induced(bits, n, vertices)
        column = UNROOTED_INDEX[k][CANON[k][sub]]
        for vertex in vertices:
            rows[vertex][column] += 1
    return tuple(tuple(row) for row in rows)


@lru_cache(maxsize=None)
def rooted_incidence_rows(bits: int, n: int, k: int) -> tuple[tuple[int, ...], ...]:
    rows = [[0] * len(ROOTED_TYPES[k]) for _ in range(n)]
    for vertices in combinations(range(n), k):
        sub = induced(bits, n, vertices)
        for local_root, vertex in enumerate(vertices):
            code = rooted_canon(sub, k, local_root)
            rows[vertex][ROOTED_INDEX[k][code]] += 1
    return tuple(tuple(row) for row in rows)


@lru_cache(maxsize=None)
def deletion_profile_rows(bits: int, n: int, k: int) -> tuple[tuple[int, ...], ...]:
    rows = []
    for vertex in range(n):
        card = delete_vertex(bits, n, vertex)
        rows.append(induced_profile(card, n - 1, k))
    return tuple(rows)


def triangle_role_rows(bits: int, n: int) -> tuple[tuple[int, int, int, int], ...]:
    """Rows are (source, middle, sink, cyclic) counts at each vertex."""
    rows = []
    for vertex in range(n):
        source = middle = sink = cyclic = 0
        others = [w for w in range(n) if w != vertex]
        for a, b in combinations(others, 2):
            wins = int(has_edge(bits, n, vertex, a)) + int(has_edge(bits, n, vertex, b))
            if is_cyclic_triangle(bits, n, tuple(sorted((vertex, a, b)))):
                cyclic += 1
            elif wins == 2:
                source += 1
            elif wins == 1:
                middle += 1
            else:
                sink += 1
        rows.append((source, middle, sink, cyclic))
    return tuple(rows)


def rank_q(matrix: tuple[tuple[int, ...], ...]) -> int:
    if not matrix or not matrix[0]:
        return 0
    a = [[Fraction(x) for x in row] for row in matrix]
    rows = len(a)
    columns = len(a[0])
    rank = 0
    for column in range(columns):
        pivot = next((r for r in range(rank, rows) if a[r][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        pivot_value = a[rank][column]
        a[rank] = [x / pivot_value for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][column]:
                factor = a[r][column]
                a[r] = [x - factor * y for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def rank_f2(matrix: tuple[tuple[int, ...], ...]) -> int:
    if not matrix or not matrix[0]:
        return 0
    bitrows = []
    for row in matrix:
        value = 0
        for column, entry in enumerate(row):
            if entry & 1:
                value |= 1 << column
        bitrows.append(value)
    rank = 0
    pivot = max((value.bit_length() for value in bitrows), default=0) - 1
    while pivot >= 0:
        hit = next((r for r in range(rank, len(bitrows)) if (bitrows[r] >> pivot) & 1), None)
        if hit is None:
            pivot -= 1
            continue
        bitrows[rank], bitrows[hit] = bitrows[hit], bitrows[rank]
        for r in range(len(bitrows)):
            if r != rank and ((bitrows[r] >> pivot) & 1):
                bitrows[r] ^= bitrows[rank]
        rank += 1
        pivot -= 1
    return rank


def matrix_add_rows(matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    if not matrix:
        return ()
    return tuple(sum(row[j] for row in matrix) for j in range(len(matrix[0])))


def tuple_scale(scale: int, row: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(scale * x for x in row)


def tuple_subtract(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x - y for x, y in zip(a, b))


def validate_all_labelled_identities() -> dict[str, int]:
    """Run cheap exact identities on every labelled tournament n <= 6."""
    tested = 0
    profile_checks = 0
    triangle_checks = 0
    four_card_checks = 0
    for n in range(1, MAX_N + 1):
        total = 1 << len(PAIR_LIST[n])
        for bits in range(total):
            tested += 1
            degrees = outdegrees(bits, n)
            if n >= 3:
                cyclic = c3_count(bits, n)
                score_formula = comb(n, 3) - sum(comb(d, 2) for d in degrees)
                require(cyclic == score_formula, f"C3 score identity failed n={n} bits={bits}")
                rows = triangle_role_rows(bits, n)
                require(sum(row[3] for row in rows) == 3 * cyclic, "cyclic incidence sum failed")
                for vertex, row in enumerate(rows):
                    d = degrees[vertex]
                    source, middle, sink, c_vertex = row
                    require(source == comb(d, 2), "rooted triangle source formula failed")
                    require(sink == comb(n - 1 - d, 2), "rooted triangle sink formula failed")
                    require(middle + c_vertex == d * (n - 1 - d), "rooted cross-pair formula failed")
                triangle_checks += 1
            if n >= 4:
                n0, n_plus, n_c4, n_minus = four_profile_named(bits, n)
                c3 = c3_count(bits, n)
                chirality = four_chirality_from_scores(bits, n)
                q_total = (n - 3) * c3
                require(n_plus - n_minus == chirality, "four-card chirality/score identity failed")
                require(n_plus + n_minus + 2 * n_c4 == q_total, "four-card q-face sum failed")
                require(
                    2 * n_plus == q_total - 2 * n_c4 + chirality,
                    "source-C3 decomposition failed",
                )
                require(
                    2 * n_minus == q_total - 2 * n_c4 - chirality,
                    "sink-C3 decomposition failed",
                )
                require(n0 == comb(n, 4) - q_total + n_c4, "q=0 decomposition failed")
                four_card_checks += 1
            # Test every proper k-profile identity.  These are deliberately on
            # all labels, not just isomorphism representatives.
            for k in range(2, n):
                profile = induced_profile(bits, n, k)
                incidence = unrooted_incidence_rows(bits, n, k)
                deleted = deletion_profile_rows(bits, n, k)
                require(matrix_add_rows(incidence) == tuple_scale(k, profile), "incidence sum failed")
                require(matrix_add_rows(deleted) == tuple_scale(n - k, profile), "Kelly deletion sum failed")
                for a_row, d_row in zip(incidence, deleted):
                    require(tuple_subtract(profile, a_row) == d_row, "incidence/deletion complement failed")
                profile_checks += 1
    return {
        "labelled_tournaments": tested,
        "triangle_checks": triangle_checks,
        "four_card_checks": four_card_checks,
        "profile_level_checks": profile_checks,
    }


def full_deck(bits: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(CANON[n - 1][delete_vertex(bits, n, v)] for v in range(n)))


def all_proper_profiles(bits: int, n: int) -> tuple[tuple[int, tuple[int, ...]], ...]:
    return tuple((k, induced_profile(bits, n, k)) for k in range(3, n))


def deleted_profile_deck(bits: int, n: int) -> tuple[tuple[int, tuple[tuple[int, ...], ...]], ...]:
    """Separate unlabelled multisets of deleted k-profiles, 3 <= k <= n-2."""
    return tuple(
        (k, tuple(sorted(deletion_profile_rows(bits, n, k))))
        for k in range(3, n - 1)
    )


def paired_deleted_profile_rows(bits: int, n: int) -> tuple[tuple[tuple[int, ...], ...], ...]:
    rows = []
    for vertex in range(n):
        rows.append(tuple(deletion_profile_rows(bits, n, k)[vertex] for k in range(3, n - 1)))
    return tuple(sorted(rows))


def incidence_rank_signature(bits: int, n: int) -> tuple[tuple[int, int, int, int, int, int], ...]:
    rows = []
    for k in (3, 4):
        if k >= n:
            continue
        incidence = unrooted_incidence_rows(bits, n, k)
        deleted = deletion_profile_rows(bits, n, k)
        rooted = rooted_incidence_rows(bits, n, k)
        rows.append(
            (
                k,
                rank_q(incidence),
                rank_q(deleted),
                rank_f2(incidence),
                rank_f2(deleted),
                rank_q(rooted) * 10 + rank_f2(rooted),
            )
        )
    return tuple(rows)


def record(bits: int, n: int) -> dict[str, object]:
    degrees = outdegrees(bits, n)
    triangles = triangle_role_rows(bits, n) if n >= 3 else tuple(() for _ in range(n))
    squares = rooted_incidence_rows(bits, n, 4) if n >= 4 else tuple(() for _ in range(n))
    deck = full_deck(bits, n) if n >= 2 else ()
    paired_card_degree = tuple(sorted((deck_card, degrees[v]) for v, deck_card in enumerate(
        CANON[n - 1][delete_vertex(bits, n, w)] for w in range(n)
    ))) if n >= 2 else ()
    paired_card_triangle = tuple(
        sorted(
            (CANON[n - 1][delete_vertex(bits, n, v)], degrees[v], triangles[v])
            for v in range(n)
        )
    ) if n >= 3 else paired_card_degree
    paired_card_local = tuple(
        sorted(
            (CANON[n - 1][delete_vertex(bits, n, v)], degrees[v], triangles[v], squares[v])
            for v in range(n)
        )
    ) if n >= 4 else paired_card_triangle
    return {
        "bits": bits,
        "score": tuple(sorted(degrees)),
        "c3": c3_count(bits, n) if n >= 3 else 0,
        "proper": all_proper_profiles(bits, n),
        "maximal_flags": maximal_flag_signature(bits, n),
        "tri_rows": tuple(sorted(triangles)),
        "square_rows": tuple(sorted(squares)),
        "paired_local_rows": tuple(sorted(zip(triangles, squares))),
        "ranks": incidence_rank_signature(bits, n),
        "deleted_profiles": deleted_profile_deck(bits, n),
        "paired_deleted_profiles": paired_deleted_profile_rows(bits, n),
        "deck": deck,
        "paired_card_degree": paired_card_degree,
        "paired_card_triangle": paired_card_triangle,
        "paired_card_local": paired_card_local,
        "four_chirality": four_chirality(bits, n) if n >= 4 else 0,
        "four_c4": four_profile_named(bits, n)[2] if n >= 4 else 0,
    }


def hierarchy_signatures(rec: dict[str, object]) -> dict[str, object]:
    score = rec["score"]
    proper = rec["proper"]
    maximal_flags = rec["maximal_flags"]
    tri = rec["tri_rows"]
    square = rec["square_rows"]
    paired_local = rec["paired_local_rows"]
    ranks = rec["ranks"]
    deleted_profiles = rec["deleted_profiles"]
    paired_deleted_profiles = rec["paired_deleted_profiles"]
    deck = rec["deck"]
    paired_card_degree = rec["paired_card_degree"]
    paired_card_triangle = rec["paired_card_triangle"]
    paired_card_local = rec["paired_card_local"]
    return {
        "G0_c3": rec["c3"],
        "G1_score": score,
        "G2_score_proper_profiles": (score, proper),
        "G2F_maximal_flag_type_sequences": maximal_flags,
        "G3_plus_incidence_ranks": (score, proper, ranks),
        "L1_plus_triangle_vertex_multiset": (score, proper, tri),
        "L2_plus_unpaired_square_multiset": (score, proper, tri, square),
        "L3_paired_triangle_square_rows": (score, proper, paired_local),
        "D1_deleted_profile_multisets": (score, proper, deleted_profiles),
        "D2_paired_deleted_profile_rows": (score, proper, paired_deleted_profiles),
        "D3_full_vertex_deck": deck,
        "D4_deck_plus_unpaired_score": (deck, score),
        "D5_paired_card_outdegree": paired_card_degree,
        "D6_paired_card_triangle": paired_card_triangle,
        "D7_paired_card_all_local": paired_card_local,
    }


HIERARCHY_ORDER = (
    "G0_c3",
    "G1_score",
    "G2_score_proper_profiles",
    "G2F_maximal_flag_type_sequences",
    "G3_plus_incidence_ranks",
    "L1_plus_triangle_vertex_multiset",
    "L2_plus_unpaired_square_multiset",
    "L3_paired_triangle_square_rows",
    "D1_deleted_profile_multisets",
    "D2_paired_deleted_profile_rows",
    "D3_full_vertex_deck",
    "D4_deck_plus_unpaired_score",
    "D5_paired_card_outdegree",
    "D6_paired_card_triangle",
    "D7_paired_card_all_local",
)


def collision_stats(
    reps: tuple[int, ...],
    signatures: dict[int, dict[str, object]],
    name: str,
    n: int,
) -> dict[str, object]:
    buckets: dict[object, list[int]] = defaultdict(list)
    for rep in reps:
        buckets[signatures[rep][name]].append(rep)
    collisions = [tuple(values) for values in buckets.values() if len(values) > 1]
    collisions.sort()
    witness = collisions[0][:2] if collisions else ()
    return {
        "buckets": len(buckets),
        "max": max(map(len, buckets.values())),
        "collision_buckets": len(collisions),
        "collision_classes": sum(map(len, collisions)),
        "witness": witness,
        "witness_converse": bool(witness and is_converse_pair(witness[0], witness[1], n)),
        "all_collisions": collisions,
    }


def incidence_rank_atlas() -> tuple[dict[tuple[int, int], Counter[tuple[int, ...]]], list[tuple[int, int, int, int, int]]]:
    distributions: dict[tuple[int, int], Counter[tuple[int, ...]]] = {}
    f2_witnesses: list[tuple[int, int, int, int, int]] = []
    for n in range(3, MAX_N + 1):
        for k in (3, 4):
            if k >= n:
                continue
            counter: Counter[tuple[int, ...]] = Counter()
            for rep in CLASS_REPS[n]:
                incidence = unrooted_incidence_rows(rep, n, k)
                deleted = deletion_profile_rows(rep, n, k)
                rq_i, rq_d = rank_q(incidence), rank_q(deleted)
                rf_i, rf_d = rank_f2(incidence), rank_f2(deleted)
                rr_q, rr_f = rank_q(rooted_incidence_rows(rep, n, k)), rank_f2(
                    rooted_incidence_rows(rep, n, k)
                )
                require(rq_i == rq_d, f"Q row-space/rank identity failed n={n} k={k}")
                require(abs(rf_i - rf_d) <= 1, f"F2 rank gap exceeds one n={n} k={k}")
                if n % 2 == 0 and k % 2 == 1:
                    require(rf_i == rf_d, f"odd/odd F2 equality failed n={n} k={k}")
                if rf_i != rf_d:
                    f2_witnesses.append((n, k, rep, rf_i, rf_d))
                counter[(rq_i, rf_i, rf_d, rr_q, rr_f)] += 1
            distributions[(n, k)] = counter
    return distributions, f2_witnesses


def full_deck_collision_anatomy(
    n: int,
    stats: dict[str, dict[str, object]],
    records: dict[int, dict[str, object]],
) -> list[dict[str, object]]:
    rows = []
    for bucket in stats["D3_full_vertex_deck"]["all_collisions"]:
        for a, b in combinations(bucket, 2):
            rows.append(
                {
                    "pair": (hex(a), hex(b)),
                    "converse": is_converse_pair(a, b, n),
                    "score_equal": records[a]["score"] == records[b]["score"],
                    "proper_profiles_equal": records[a]["proper"] == records[b]["proper"],
                    "triangle_rows_equal": records[a]["tri_rows"] == records[b]["tri_rows"],
                    "square_rows_equal": records[a]["square_rows"] == records[b]["square_rows"],
                    "paired_local_equal": records[a]["paired_local_rows"] == records[b]["paired_local_rows"],
                    "paired_card_degree_separates": records[a]["paired_card_degree"]
                    != records[b]["paired_card_degree"],
                }
            )
    return rows


def type_legend(k: int) -> list[str]:
    rows = []
    for index, rep in enumerate(UNROOTED_TYPES[k]):
        rows.append(f"u{index}=0x{rep:x}/score{score_sequence(rep, k)}")
    return rows


def rooted_type_legend(k: int) -> list[str]:
    rows = []
    for index, code in enumerate(ROOTED_TYPES[k]):
        root_degree = sum(has_edge(code, k, 0, w) for w in range(1, k))
        unrooted = UNROOTED_INDEX[k][CANON[k][code]]
        rows.append(f"r{index}=0x{code:x}/u{unrooted}/root_out{root_degree}")
    return rows


def subset_masks(n: int, k: int) -> tuple[int, ...]:
    return tuple(sum(1 << i for i in subset) for subset in combinations(range(n), k))


def rank_f2_bitrows(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def kneser_inclusion_inverse_atlas() -> list[dict[str, object]]:
    """Check the proposed W_{r,r+1}^{-1} on intersection classes, r<=6.

    Complementing an (r+1)-column turns inclusion into the symmetric
    disjointness adjacency of KG(2r+1,r).  For fixed S,R with j=|S cap R|,
    the sparse product AB has only r+1 summands, and depends only on j.
    Checking all j is therefore an exact full-matrix check, not sampling.
    """
    atlas = []
    for r in range(1, 7):
        n = 2 * r + 1
        masks = subset_masks(n, r)
        index = {mask: i for i, mask in enumerate(masks)}

        def inverse_value(intersection: int) -> Fraction:
            return Fraction((-1) ** intersection, (r + 1) * comb(r, intersection))

        products = []
        for j in range(r + 1):
            a = r - j
            value = Fraction(0)
            if a:
                value += a * inverse_value(a - 1)
            value += (j + 1) * inverse_value(a)
            products.append(value)
            require(value == (1 if j == r else 0), f"Kneser inverse failed r={r} j={j}")

        # Independent exact F2 row construction.  A row has the r+1
        # r-subsets obtained by omitting one point from S^c.
        universe = (1 << n) - 1
        bitrows: list[int] = []
        for mask in masks:
            complement = universe ^ mask
            row = 0
            point = complement
            while point:
                low = point & -point
                disjoint_r_set = complement ^ low
                row |= 1 << index[disjoint_r_set]
                point ^= low
            bitrows.append(row)
        f2_rank = rank_f2_bitrows(bitrows)
        require(f2_rank == comb(2 * r, r), f"unexpected finite F2 rank r={r}")

        multiplicities = tuple(
            comb(n, i) - (comb(n, i - 1) if i else 0)
            for i in range(r + 1)
        )
        eigenvalues = tuple(((-1) ** i) * (r + 1 - i) for i in range(r + 1))
        require(sum(multiplicities) == len(masks), f"Kneser spectrum multiplicity failed r={r}")
        require(all(value for value in eigenvalues), f"singular Kneser spectrum r={r}")
        atlas.append(
            {
                "r": r,
                "n": n,
                "size": len(masks),
                "AB_by_intersection_j=0..r": tuple(products),
                "inverse_values_j=0..r": tuple(inverse_value(j) for j in range(r + 1)),
                "Q_rank": len(masks),
                "F2_rank": f2_rank,
                "F2_nullity": len(masks) - f2_rank,
                "spectrum": tuple(zip(eigenvalues, multiplicities)),
            }
        )
    return atlas


def extend_by_observer(base_bits: int, base_n: int, incident_word: int) -> int:
    """Add vertex base_n; bit v of word means v -> new vertex."""
    n = base_n + 1
    out = 0
    for i, j in PAIR_LIST[base_n]:
        if has_edge(base_bits, base_n, i, j):
            out |= 1 << PAIR_INDEX[n][(i, j)]
    for vertex in range(base_n):
        if (incident_word >> vertex) & 1:
            out |= 1 << PAIR_INDEX[n][(vertex, base_n)]
    return out


def targeted_n7_four_card_audit(base_reps: tuple[int, ...]) -> dict[str, object]:
    """All 64 observer extensions of the eight n=6 full-deck colliders."""
    chirality_hist: Counter[int] = Counter()
    c4_hist: Counter[int] = Counter()
    score_to_chirality: dict[tuple[int, ...], int] = {}
    total = 0
    for base in base_reps:
        for word in range(1 << 6):
            bits = extend_by_observer(base, 6, word)
            total += 1
            n0, n_plus, n_c4, n_minus = four_profile_named(bits, 7)
            chirality = four_chirality_from_scores(bits, 7)
            q_total = 4 * c3_count(bits, 7)
            require(n_plus - n_minus == chirality, "target n=7 chirality failure")
            require(n_plus + n_minus + 2 * n_c4 == q_total, "target n=7 q sum failure")
            require(n0 == comb(7, 4) - q_total + n_c4, "target n=7 n0 failure")
            score = score_sequence(bits, 7)
            previous = score_to_chirality.setdefault(score, chirality)
            require(previous == chirality, "chirality not score-determined in target n=7 family")
            chirality_hist[chirality] += 1
            c4_hist[n_c4] += 1
    return {
        "family": "all 64 one-observer extensions of the 8 n=6 full-deck-collision representatives",
        "objects": total,
        "distinct_scores": len(score_to_chirality),
        "chirality_hist": dict(sorted(chirality_hist.items())),
        "c4_hist": dict(sorted(c4_hist.items())),
    }


def permutation_parity(order: tuple[int, ...]) -> int:
    return sum(order[i] > order[j] for i in range(len(order)) for j in range(i + 1, len(order))) & 1


def canonical_parity_map(d: int) -> dict[int, int]:
    """Parity of the isomorphism from a labelled d-tournament to its rep.

    ``relabel(rep, order)`` produces a labelled image; the inverse order has
    the same parity.  Duplicate images must have the same parity exactly when
    every tournament automorphism is even.  Iterating every representative
    and every permutation validates both well-definedness and totality.
    """
    parity_map: dict[int, int] = {}
    for rep in CLASS_REPS[d]:
        automorphism_parities = []
        for order in ORDERS[d]:
            image = relabel(rep, d, order)
            parity = permutation_parity(order)
            if image == rep:
                automorphism_parities.append(parity)
            if image in parity_map:
                require(parity_map[image] == parity, f"canonical parity ambiguity d={d} image={image}")
            else:
                parity_map[image] = parity
        require(automorphism_parities and not any(automorphism_parities), f"odd automorphism d={d}")
    require(len(parity_map) == 1 << len(PAIR_LIST[d]), f"canonical parity map incomplete d={d}")
    return parity_map


def directed_boundary_descriptor(
    edge_bits: int,
    d: int,
    parity_map: dict[int, int],
) -> tuple[tuple[int, int], ...]:
    """(unlabelled face type, parity xor boundary sign) for d+1 faces."""
    rows = []
    for omitted in range(d + 1):
        vertices = tuple(v for v in range(d + 1) if v != omitted)
        face = induced(edge_bits, d + 1, vertices)
        type_index = UNROOTED_INDEX[d][CANON[d][face]]
        rows.append((type_index, parity_map[face] ^ (omitted & 1)))
    return tuple(rows)


def descriptor_satisfied(descriptor: tuple[tuple[int, int], ...], rule: int) -> bool:
    first_type, first_constant = descriptor[0]
    target = ((rule >> first_type) & 1) ^ first_constant
    return all((((rule >> face_type) & 1) ^ constant) == target for face_type, constant in descriptor[1:])


def rule_orientation_sign(bits: int, d: int, rule: int, parity_map: dict[int, int]) -> int:
    type_index = UNROOTED_INDEX[d][CANON[d][bits]]
    oriented_bit = ((rule >> type_index) & 1) ^ parity_map[bits]
    return 1 if oriented_bit else -1


def switch_tournament(bits: int, n: int, vertex_mask: int) -> int:
    out = bits
    for index, (i, j) in enumerate(PAIR_LIST[n]):
        if ((vertex_mask >> i) & 1) ^ ((vertex_mask >> j) & 1):
            out ^= 1 << index
    return out


def local_rule_optimization(d: int) -> tuple[dict[str, object], set[int], dict[int, int]]:
    """Optimize every equivariant alternating d-face rule exactly."""
    parity_map = canonical_parity_map(d)
    expected_reps = {
        4: (0, 2, 4, 5),
        5: (0, 2, 4, 5, 8, 9, 10, 11, 12, 40, 41, 76),
    }[d]
    require(UNROOTED_TYPES[d] == expected_reps, f"unexpected rule-bit ordering d={d}")
    total_edges = 1 << len(PAIR_LIST[d + 1])
    descriptors = [directed_boundary_descriptor(bits, d, parity_map) for bits in range(total_edges)]
    descriptor_hist = Counter(descriptors)
    type_count = len(UNROOTED_TYPES[d])
    full_rule_mask = (1 << type_count) - 1
    counts = [0] * (1 << type_count)
    inconsistent_edge_tournaments = 0

    # A descriptor imposes rule[type_i] = y xor constant_i for one global
    # y.  Enumerate the two y choices and every free rule bit, rather than the
    # 4096 rules separately for every one of the 32768 edge tournaments.
    for descriptor, multiplicity in descriptor_hist.items():
        fixed_mask = 0
        base_assignment = 0
        inconsistent = False
        for face_type, constant in descriptor:
            bit = 1 << face_type
            if fixed_mask & bit:
                if ((base_assignment >> face_type) & 1) != constant:
                    inconsistent = True
                    break
            else:
                fixed_mask |= bit
                if constant:
                    base_assignment |= bit
        if inconsistent:
            inconsistent_edge_tournaments += multiplicity
            continue
        free_mask = full_rule_mask ^ fixed_mask
        for assignment in (base_assignment, base_assignment ^ fixed_mask):
            free = free_mask
            while True:
                counts[assignment | free] += multiplicity
                if free == 0:
                    break
                free = (free - 1) & free_mask

    for rule in range(1 << type_count):
        require(counts[rule] == counts[rule ^ full_rule_mask], f"global-negative mismatch d={d}")
    best_count = max(counts)
    best_rules = tuple(i for i, count in enumerate(counts) if count == best_count)
    expected_best = {4: 144, 5: 4608}[d]
    expected_rules = {4: (7, 8), 5: (1415, 2680)}[d]
    require(best_count == expected_best, f"local-rule optimum mismatch d={d}")
    require(best_rules == expected_rules, f"best-rule IDs mismatch d={d}")
    best_successes = {
        bits
        for bits, descriptor in enumerate(descriptors)
        if descriptor_satisfied(descriptor, best_rules[0])
    }
    require(len(best_successes) == best_count, f"success-set size mismatch d={d}")
    return (
        {
            "d": d,
            "rule_reps": UNROOTED_TYPES[d],
            "rules": 1 << type_count,
            "edge_tournaments": total_edges,
            "descriptor_types": len(descriptor_hist),
            "intrinsically_inconsistent_edges": inconsistent_edge_tournaments,
            "count_histogram": dict(sorted(Counter(counts).items())),
            "best_count": best_count,
            "best_probability": Fraction(best_count, total_edges),
            "best_rules": best_rules,
        },
        best_successes,
        parity_map,
    )


def oriented_two_graph_gate(best_successes: set[int]) -> dict[str, object]:
    """Identify the d=5 success set with the seed-10 switching/relabel orbit."""
    seed = 10
    require(min(best_successes) == seed, "unexpected first d=5 success")
    # S and its complement induce the same cut, so fix vertex 0 outside S.
    switch_masks = tuple(range(0, 1 << 6, 2))
    orbit = {
        switch_tournament(relabel(seed, 6, order), 6, switch_mask)
        for order in ORDERS[6]
        for switch_mask in switch_masks
    }
    require(len(orbit) == 4608, "seed-10 switching/relabel orbit size mismatch")
    require(orbit == best_successes, "best-rule successes differ from oriented two-graph orbit")
    stabilizer_order = sum(
        switch_tournament(relabel(seed, 6, order), 6, switch_mask) == seed
        for order in ORDERS[6]
        for switch_mask in switch_masks
    )
    require(stabilizer_order == 5, "triangle-flux stabilizer order mismatch")
    switching_reps = tuple(
        sorted(CANON[6][switch_tournament(seed, 6, switch_mask)] for switch_mask in switch_masks)
    )
    switching_reps = tuple(sorted(set(switching_reps)))
    expected_reps = (10, 21, 81, 140, 147, 148, 313, 314)
    require(switching_reps == expected_reps, "switching-class ordinary types mismatch")
    automorphism_orders = tuple(factorial(6) // ORBIT_SIZE[6][rep] for rep in switching_reps)
    require(automorphism_orders == (1, 1, 1, 5, 1, 1, 1, 5), "switching type aut orders mismatch")
    return {
        "seed": seed,
        "switching_relabel_orbit": len(orbit),
        "equals_best_success_set": True,
        "triangle_flux_stabilizer_order": stabilizer_order,
        "ordinary_iso_reps": switching_reps,
        "ordinary_aut_orders": automorphism_orders,
    }


def edge_sign(bits: int, n: int, i: int, j: int) -> int:
    return 1 if has_edge(bits, n, i, j) else -1


def triangle_flux(bits: int, n: int, v: int, a: int, b: int) -> int:
    return edge_sign(bits, n, v, a) * edge_sign(bits, n, a, b) * edge_sign(bits, n, b, v)


def flux_pfaffian(bits: int, omitted_vertex: int) -> int:
    """Pfaffian of H^(v) in the ascending order of the other four vertices."""
    a, b, c, d = [vertex for vertex in range(5) if vertex != omitted_vertex]

    def h(i: int, j: int) -> int:
        return triangle_flux(bits, 5, omitted_vertex, i, j)

    return h(a, b) * h(c, d) - h(a, c) * h(b, d) + h(a, d) * h(b, c)


def optimal_d5_pfaffian_formula_gate(parity_map: dict[int, int]) -> dict[str, object]:
    """Identify rule 1415 with its 16-term exact edge-sign formula."""
    rule = 1415
    values: list[int] = []
    rhs_hist: Counter[int] = Counter()
    for bits in range(1 << len(PAIR_LIST[5])):
        # Orientation sign convention: bit 1 is +1 and bit 0 is -1.
        f_value = rule_orientation_sign(bits, 5, rule, parity_map)
        edge_product = 1
        for i, j in PAIR_LIST[5]:
            edge_product *= edge_sign(bits, 5, i, j)
        rhs = sum(((-1) ** v) * flux_pfaffian(bits, v) for v in range(5)) - edge_product
        require(rhs in (-4, 4), f"Pfaffian RHS escaped +/-4 at bits={bits}")
        require(4 * f_value == rhs, f"Pfaffian local-rule formula failed bits={bits}")
        values.append(f_value)
        rhs_hist[rhs] += 1

    # Exact Walsh transform.  Standard bit characters use (-1)^bit, while
    # x_e=+1 on bit 1 and -1 on bit 0; multiply by (-1)^degree to convert.
    walsh = values[:]
    stride = 1
    while stride < len(walsh):
        for start in range(0, len(walsh), 2 * stride):
            for offset in range(stride):
                x = walsh[start + offset]
                y = walsh[start + offset + stride]
                walsh[start + offset] = x + y
                walsh[start + offset + stride] = x - y
        stride *= 2
    fourier_numerators = {
        mask: coefficient * ((-1) ** mask.bit_count())
        for mask, coefficient in enumerate(walsh)
        if coefficient
    }
    degree_hist = Counter(mask.bit_count() for mask in fourier_numerators)
    require(degree_hist == Counter({6: 15, 10: 1}), "unexpected Pfaffian-rule Fourier support")
    require(set(map(abs, fourier_numerators.values())) == {256}, "unexpected Fourier coefficient")
    full_mask = (1 << len(PAIR_LIST[5])) - 1
    require(fourier_numerators[full_mask] == -256, "full edge-parity coefficient mismatch")
    bowtie_masks = [mask for mask in fourier_numerators if mask.bit_count() == 6]
    for mask in bowtie_masks:
        degrees = [0] * 5
        for edge_index, (i, j) in enumerate(PAIR_LIST[5]):
            if (mask >> edge_index) & 1:
                degrees[i] += 1
                degrees[j] += 1
        require(sorted(degrees) == [2, 2, 2, 2, 4], f"non-bowtie degree-six support {mask}")
    support = tuple(
        (hex(mask), numerator, mask.bit_count())
        for mask, numerator in sorted(fourier_numerators.items())
    )
    return {
        "checked_labelled_5_tournaments": len(values),
        "rhs_hist": dict(sorted(rhs_hist.items())),
        "formula": "4f=sum_v (-1)^v Pf(H^(v))-product_edges(x_ij)",
        "fourier_support": support,
        "fourier_coefficients": "numerator/1024, hence each displayed +/-256 is +/-1/4",
        "support_summary": "15 degree-6 bowties plus degree-10 full edge parity",
    }


def apex_normalized_card(bits: int, apex: int) -> int:
    """Switch a labelled 5-tournament so apex dominates, then delete apex."""
    signs = [edge_sign(bits, 5, apex, vertex) if vertex != apex else 1 for vertex in range(5)]
    remaining = [vertex for vertex in range(5) if vertex != apex]
    out = 0
    for index, (i, j) in enumerate(PAIR_LIST[4]):
        old_i, old_j = remaining[i], remaining[j]
        switched_sign = signs[old_i] * edge_sign(bits, 5, old_i, old_j) * signs[old_j]
        if switched_sign == 1:
            out |= 1 << index
    # Explicit gauge check: x'_{apex,u}=1 for every remaining vertex.
    require(
        all(signs[apex] * edge_sign(bits, 5, apex, u) * signs[u] == 1 for u in remaining),
        "apex normalization failed",
    )
    return out


def apex_lift_recursion_gate(parity4: dict[int, int], parity5: dict[int, int]) -> dict[str, object]:
    ratios = [set() for _ in range(5)]
    checked = 0
    for bits in range(1 << len(PAIR_LIST[5])):
        f5 = rule_orientation_sign(bits, 5, 1415, parity5)
        for apex in range(5):
            card = apex_normalized_card(bits, apex)
            f4 = rule_orientation_sign(card, 4, 8, parity4)
            ratio = f5 * f4
            ratios[apex].add(ratio)
            require(f5 == ((-1) ** (apex + 1)) * f4, f"apex lift failed bits={bits} v={apex}")
            checked += 1
    ratio_tuple = tuple(next(iter(values)) for values in ratios)
    require(ratio_tuple == (-1, 1, -1, 1, -1), "apex ratio word mismatch")
    return {
        "checked_(T,v)_pairs": checked,
        "f4_rule": 8,
        "f4_rep_signs": (-1, -1, -1, 1),
        "f5_rule": 1415,
        "ratio_by_apex_0..4": ratio_tuple,
        "identity": "f5(T)=(-1)^(v+1) f4(apex-normalized T-v)",
    }


def recursive_density_gate(parity4: dict[int, int], parity5: dict[int, int]) -> dict[str, object]:
    """Gate the exact 9/64 density inheritance and its stopping obstruction.

    Add a sixth vertex which dominates a labelled five-tournament U.  The
    f5-oriented boundary condition on the resulting six-tournament is then
    equivalent, object by object, to the f4-oriented boundary condition on U.
    We also exhaustively check that f5 is invariant under every vertex switch.
    """
    directed_f4: set[int] = set()
    directed_f5_apex_lifts: set[int] = set()
    for bits in range(1 << len(PAIR_LIST[5])):
        apex_lift = extend_by_observer(bits, 5, 0)  # new vertex 5 dominates all old vertices
        require(induced(apex_lift, 6, tuple(range(5))) == bits, "apex lift changed base card")
        is_f4_directed = descriptor_satisfied(
            directed_boundary_descriptor(bits, 4, parity4),
            8,
        )
        is_f5_directed = descriptor_satisfied(
            directed_boundary_descriptor(apex_lift, 5, parity5),
            1415,
        )
        require(is_f5_directed == is_f4_directed, f"recursive density mismatch U={bits}")
        if is_f4_directed:
            directed_f4.add(bits)
        if is_f5_directed:
            directed_f5_apex_lifts.add(bits)
    require(directed_f4 == directed_f5_apex_lifts, "recursive success sets differ")
    require(len(directed_f4) == 144, "recursive 9/64 success count mismatch")

    switching_pairs = 0
    switching_failures = 0
    # A cut and its complement agree, so fix vertex 0 outside the switching set.
    for bits in range(1 << len(PAIR_LIST[5])):
        original = rule_orientation_sign(bits, 5, 1415, parity5)
        for switch_mask in range(0, 1 << 5, 2):
            switched = switch_tournament(bits, 5, switch_mask)
            if rule_orientation_sign(switched, 5, 1415, parity5) != original:
                switching_failures += 1
            switching_pairs += 1
    require(switching_failures == 0, "f5 is not switching-invariant")

    return {
        "normalized_U_checked": 1 << len(PAIR_LIST[5]),
        "directed_f4_U": len(directed_f4),
        "directed_f5_apex_lifts": len(directed_f5_apex_lifts),
        "success_sets_equal": True,
        "density": Fraction(len(directed_f4), 1 << len(PAIR_LIST[5])),
        "f5_switching_pairs_checked": switching_pairs,
        "f5_switching_failures": switching_failures,
        "universal_next_apex_lift_possible": len(directed_f4) == (1 << len(PAIR_LIST[5])),
    }


def normalized_flat_edge_sign(mask: int, n: int, i: int, j: int) -> int:
    """Skew edge field with x_{0j}=+1; other upper entries come from mask."""
    if i > j:
        return -normalized_flat_edge_sign(mask, n, j, i)
    require(i < j, "flat edge sign requested on diagonal")
    if i == 0:
        return 1
    free_index = PAIR_INDEX[n - 1][(i - 1, j - 1)]
    return 1 if (mask >> free_index) & 1 else -1


def normalized_flat_chi(mask: int, n: int, i: int, j: int, k: int) -> int:
    return (
        normalized_flat_edge_sign(mask, n, i, j)
        * normalized_flat_edge_sign(mask, n, j, k)
        * normalized_flat_edge_sign(mask, n, k, i)
    )


def uniform_rank3_gp_holds(mask: int, n: int) -> bool:
    """Uniform rank-3 Grassmann-Pluecker sign axiom for a flat triangle field."""
    for five_set in combinations(range(n), 5):
        for pivot in five_set:
            b, c, d, e = [vertex for vertex in five_set if vertex != pivot]
            terms = (
                normalized_flat_chi(mask, n, pivot, b, c)
                * normalized_flat_chi(mask, n, pivot, d, e),
                -normalized_flat_chi(mask, n, pivot, b, d)
                * normalized_flat_chi(mask, n, pivot, c, e),
                normalized_flat_chi(mask, n, pivot, b, e)
                * normalized_flat_chi(mask, n, pivot, c, d),
            )
            if terms[0] == terms[1] == terms[2]:
                return False
    return True


def arbitrary_alternating_chi(
    ascending_values: dict[tuple[int, int, int], int],
    i: int,
    j: int,
    k: int,
) -> int:
    ordered = (i, j, k)
    ascending = tuple(sorted(ordered))
    local_order = tuple(ascending.index(vertex) for vertex in ordered)
    return ascending_values[ascending] * ((-1) ** permutation_parity(local_order))


def hostile_gp_field_holds() -> tuple[bool, tuple[tuple[int, tuple[int, int, int]], ...]]:
    """All ascending chi=+1 except chi_{135}=-1 (one-based labels)."""
    values = {triple: 1 for triple in combinations(range(5), 3)}
    values[(0, 2, 4)] = -1
    failures = []
    five_set = tuple(range(5))
    for pivot in five_set:
        b, c, d, e = [vertex for vertex in five_set if vertex != pivot]
        terms = (
            arbitrary_alternating_chi(values, pivot, b, c)
            * arbitrary_alternating_chi(values, pivot, d, e),
            -arbitrary_alternating_chi(values, pivot, b, d)
            * arbitrary_alternating_chi(values, pivot, c, e),
            arbitrary_alternating_chi(values, pivot, b, e)
            * arbitrary_alternating_chi(values, pivot, c, d),
        )
        if terms[0] == terms[1] == terms[2]:
            failures.append((pivot, terms))
    return not failures, tuple(failures)


def flat_gp_atlas() -> list[dict[str, object]]:
    atlas = []
    for n in (5, 6, 7):
        free_edges = comb(n - 1, 2)
        total = 1 << free_edges
        good = sum(uniform_rank3_gp_holds(mask, n) for mask in range(total))
        require(good == factorial(n - 1), f"flat GP count mismatch n={n}")
        # Positive flatness control: on every increasing four-set, the product
        # of its four increasing triangle fluxes is +1.
        for mask in range(total):
            for a, b, c, d in combinations(range(n), 4):
                product = (
                    normalized_flat_chi(mask, n, a, b, c)
                    * normalized_flat_chi(mask, n, a, b, d)
                    * normalized_flat_chi(mask, n, a, c, d)
                    * normalized_flat_chi(mask, n, b, c, d)
                )
                require(product == 1, f"triangle field not flat n={n} mask={mask}")
        atlas.append(
            {
                "n": n,
                "flat_fields": total,
                "uniform_rank3_GP": good,
                "factorial_(n-1)": factorial(n - 1),
            }
        )
    hostile_ok, hostile_failures = hostile_gp_field_holds()
    require(not hostile_ok, "single-flipped chi135 hostile unexpectedly passed GP")
    atlas.append(
        {
            "hostile": "n=5 ascending chi all +1 except chi_135=-1",
            "GP_holds": hostile_ok,
            "failure_pivots_zero_based": hostile_failures,
        }
    )
    return atlas


def main() -> None:
    prepare()
    four_qh_atlas = validate_four_type_qh()
    identity_counts = validate_all_labelled_identities()
    rank_distributions, f2_witnesses = incidence_rank_atlas()
    kneser_atlas = kneser_inclusion_inverse_atlas()
    d4_rule_atlas, _d4_successes, _d4_parity = local_rule_optimization(4)
    d5_rule_atlas, d5_successes, _d5_parity = local_rule_optimization(5)
    two_graph_gate = oriented_two_graph_gate(d5_successes)
    d5_pfaffian_gate = optimal_d5_pfaffian_formula_gate(_d5_parity)
    apex_recursion_gate = apex_lift_recursion_gate(_d4_parity, _d5_parity)
    density_recursion_gate = recursive_density_gate(_d4_parity, _d5_parity)
    gp_atlas = flat_gp_atlas()

    all_records: dict[int, dict[int, dict[str, object]]] = {}
    all_signatures: dict[int, dict[int, dict[str, object]]] = {}
    all_stats: dict[int, dict[str, dict[str, object]]] = {}
    for n in range(3, MAX_N + 1):
        records = {rep: record(rep, n) for rep in CLASS_REPS[n]}
        signatures = {rep: hierarchy_signatures(rec) for rep, rec in records.items()}
        stats = {
            name: collision_stats(CLASS_REPS[n], signatures, name, n)
            for name in HIERARCHY_ORDER
        }
        all_records[n] = records
        all_signatures[n] = signatures
        all_stats[n] = stats
        # The proper maximal-flag distribution and the full deck carry exactly
        # the same unmarked information.  Check equality of their fibres.
        for a, b in combinations(CLASS_REPS[n], 2):
            same_flag = records[a]["maximal_flags"] == records[b]["maximal_flags"]
            same_deck = records[a]["deck"] == records[b]["deck"]
            require(same_flag == same_deck, f"flag/deck fibre mismatch n={n} a={a} b={b}")

    n6_deck_colliders = tuple(
        sorted(
            {
                rep
                for bucket in all_stats[6]["D3_full_vertex_deck"]["all_collisions"]
                for rep in bucket
            }
        )
    )
    require(len(n6_deck_colliders) == 8, "unexpected n=6 full-deck collider count")
    targeted_n7 = targeted_n7_four_card_audit(n6_deck_colliders)

    print("TOURNAMENT SUBSET / DECK HIERARCHY EXACT SCOUT")
    print("=" * 78)
    print("status=FINITE-EXACT universe=all labelled ordinary tournaments n<=6")
    print("truth_gates=require() active_in_normal_and_python_-O")

    print("\nA. Labelled universe and isomorphism quotient")
    for n in range(1, MAX_N + 1):
        orbit_hist = Counter(ORBIT_SIZE[n].values())
        aut_hist = Counter(factorial(n) // size for size in ORBIT_SIZE[n].values())
        print(
            f"  n={n} labelled={1 << len(PAIR_LIST[n])} classes={len(CLASS_REPS[n])} "
            f"orbit_sizes={dict(sorted(orbit_hist.items()))} aut_sizes={dict(sorted(aut_hist.items()))}"
        )
    print(f"  identity_control_counts={identity_counts}")

    print("\nB. Type-column legends")
    for k in (3, 4):
        print(f"  k={k} unrooted: {'; '.join(type_legend(k))}")
        print(f"  k={k} rooted:   {'; '.join(rooted_type_legend(k))}")

    print("\nC. Reconstruction collisions under richer unlabelled summaries")
    print("  fields: buckets/max/collision_buckets/collision_classes/first_witness/converse")
    for n in range(3, MAX_N + 1):
        print(f"  n={n} classes={len(CLASS_REPS[n])}")
        for name in HIERARCHY_ORDER:
            row = all_stats[n][name]
            witness = tuple(hex(x) for x in row["witness"])
            print(
                f"    {name:38s} {row['buckets']:3d}/{row['max']:2d}/"
                f"{row['collision_buckets']:2d}/{row['collision_classes']:2d} "
                f"witness={witness} converse={row['witness_converse']}"
            )

    print("\nC2. Separate proper-k profiles versus their joint and flag-glued forms")
    for n in range(4, MAX_N + 1):
        print(f"  n={n}")
        for k in range(3, n):
            buckets: dict[tuple[int, ...], list[int]] = defaultdict(list)
            for rep in CLASS_REPS[n]:
                buckets[induced_profile(rep, n, k)].append(rep)
            collisions = [members for members in buckets.values() if len(members) > 1]
            witness = tuple(hex(x) for x in collisions[0][:2]) if collisions else ()
            print(
                f"    separate_k={k} buckets={len(buckets)} max={max(map(len, buckets.values()))} "
                f"collision_buckets={len(collisions)} witness={witness}"
            )
        joint_buckets: dict[object, list[int]] = defaultdict(list)
        flag_buckets: dict[object, list[int]] = defaultdict(list)
        for rep in CLASS_REPS[n]:
            joint_buckets[all_records[n][rep]["proper"]].append(rep)
            flag_buckets[all_records[n][rep]["maximal_flags"]].append(rep)
        joint_collisions = sum(len(values) > 1 for values in joint_buckets.values())
        flag_collisions = sum(len(values) > 1 for values in flag_buckets.values())
        print(
            f"    joint_all_proper buckets={len(joint_buckets)} collision_buckets={joint_collisions} "
            f"flag_glued buckets={len(flag_buckets)} collision_buckets={flag_collisions}"
        )

    print("\nD. Full-deck collision anatomy")
    for n in range(3, MAX_N + 1):
        rows = full_deck_collision_anatomy(n, all_stats[n], all_records[n])
        print(f"  n={n} collision_pairs={len(rows)}")
        for row in rows:
            print(f"    {row}")

    print("\nE. Incidence/deletion rank atlas")
    print("  distribution key=(rank_Q,rank_F2_inc,rank_F2_del,rank_Q_rooted,rank_F2_rooted)")
    for (n, k), counter in sorted(rank_distributions.items()):
        print(f"  n={n} k={k} distribution={dict(sorted(counter.items()))}")
    f2_witness_hist = Counter((n, k, ri, rd) for n, k, _rep, ri, rd in f2_witnesses)
    print(f"  F2_rank_difference_hist={dict(sorted(f2_witness_hist.items()))}")
    print(
        "  first_F2_rank_difference_witness="
        f"{(f2_witnesses[0][0], f2_witnesses[0][1], hex(f2_witnesses[0][2]), f2_witnesses[0][3], f2_witnesses[0][4])}"
    )

    print("\nE2. Odd inclusion / Kneser exact inverse and ranks")
    print("  W_{r,r+1}, after complementing columns, is A[S,R]=1[S cap R=empty].")
    print("  Proposed B[S,R]=(-1)^|S cap R|/((r+1) C(r,|S cap R|)).")
    for row in kneser_atlas:
        print(f"  {row}")
    print("  n=7,r=3 specialization: B intersection values=(1/4,-1/12,1/12,-1/4),")
    print("  spectrum=((4,1),(-3,6),(2,14),(-1,14)), Q-rank=35, F2-rank=20/nullity=15.")

    print("\nE3. Four-card q/chirality decomposition")
    print(f"  four_type_classifier_(q,h)->type_index={dict(sorted(four_qh_atlas.items()))}")
    for n in range(4, MAX_N + 1):
        chirality_hist = Counter(int(all_records[n][rep]["four_chirality"]) for rep in CLASS_REPS[n])
        c4_hist = Counter(int(all_records[n][rep]["four_c4"]) for rep in CLASS_REPS[n])
        score_buckets: dict[tuple[int, ...], set[int]] = defaultdict(set)
        for rep in CLASS_REPS[n]:
            score_buckets[all_records[n][rep]["score"]].add(int(all_records[n][rep]["four_chirality"]))
        require(all(len(values) == 1 for values in score_buckets.values()), "chirality/score collision")
        print(
            f"  n={n} chirality_hist={dict(sorted(chirality_hist.items()))} "
            f"c4_hist={dict(sorted(c4_hist.items()))} score_determines_chirality=True"
        )
    print(f"  targeted_n7={targeted_n7}")

    print("\nE4. Exact equivariant alternating local-rule optimization")
    print("  Rule bit order is the listed canonical d-tournament representatives; canonicalizing")
    print("  permutation parity is well-defined because every checked tournament automorphism is even.")
    print(f"  d=4 {d4_rule_atlas}")
    print(f"  d=5 {d5_rule_atlas}")
    print(f"  d=5_oriented_two_graph_gate={two_graph_gate}")
    print(f"  d=5_pfaffian_formula_gate={d5_pfaffian_gate}")
    print(f"  apex_lift_recursion_gate={apex_recursion_gate}")
    print(f"  recursive_density_gate={density_recursion_gate}")
    print("  Both optima are exactly 9/64 and have exactly two globally-negative rules.")
    print("  The d=5 success set is exactly the seed-10 vertex-switch/S6 orbit, giving an explicit")
    print("  oriented-two-graph bridge.  Any external c5 interpretation or novelty claim remains unaudited.")

    print("\nE5. Flat oriented-two-graph fields and uniform rank-3 GP")
    for row in gp_atlas:
        print(f"  {row}")
    print("  Counts 24,120,720=(n-1)! suggest circular-order fields; that classification is not claimed.")

    print("\nF. Exact identities and mechanisms")
    print("  I1 c3(T)=C(n,3)-sum_v C(d+(v),2); global triangle count adds no score data.")
    print("  I2 rooted triangle row v=(C(d,2), d(n-1-d)-c_v, C(n-1-d,2), c_v).")
    print("  I3 for each induced k-type vector p and vertex v: A_v + D_v = p,")
    print("     sum_v A_v=k p, sum_v D_v=(n-k)p (A=containing incidence, D=deleted profile).")
    print("  I4 over Q, rowspace(A)=rowspace(D), hence ranks agree exactly for 0<k<n.")
    print("     Over F2 ranks differ by at most one; if n even and k odd they agree.")
    print("  I5 the full deck determines every proper induced-type profile by Kelly double counting.")
    print("     It also determines paired (card,c3-loss), so c3-loss is not a genuine boundary repair.")
    print("  I6 card + missing-vertex outdegree is a genuine marked-boundary coordinate; it resolves")
    print("     every isomorphism collision checked through n=6, while deck + unpaired scores does not.")
    print("  I7 every four-card is classified by (q,h): q=#cyclic triangle faces in {0,1,2},")
    print("     h=#universal sources-#universal sinks in {-1,0,1}.  With D=sum_v[")
    print("     C(d+(v),3)-C(d-(v),3)] and C4=#q=2 cards,")
    print("     N_plus=((n-3)C3-2C4+D)/2 and N_minus=((n-3)C3-2C4-D)/2,")
    print("     and N0=C(n,4)-(n-3)C3+C4.  Thus chirality is variable but score-determined;")
    print("     the directed Hamiltonian C4 count is the only new scalar in the four-profile.")
    print("  I8 complementing W_{r,r+1} gives KG(2r+1,r), and the displayed intersection-kernel")
    print("     is its exact rational inverse.  The check is exhaustive over all intersection classes")
    print("     for r=1..6; the checked finite F2 ranks equal C(2r,r).")
    print("  I9 the unmarked maximal proper-flag type-sequence distribution is exactly deck-equivalent:")
    print("     partition flags by their omitted last vertex to derive it from the full deck; conversely")
    print("     the terminal (n-1)-type in each sequence recovers every card multiplicity.  Therefore")
    print("     unmarked flag gluing cannot repair a deck collision; a boundary/card mark is essential.")
    print("  I10 exhaustive local-rule optimization gives 144/1024=4608/32768=9/64 for d=4,5.")
    print("      At d=5 the 4608 successful edge tournaments are exactly one switching/relabel orbit")
    print("      with seed 10 and stabilizer order 5; this is finite-exact, not a literature-priority claim.")
    print("  I11 rule 1415 is not an opaque table: on every labelled five-tournament its sign is")
    print("      one quarter of an alternating sum of five triangle-flux Pfaffians minus full edge parity.")
    print("      Its exact Fourier support is 15 degree-six bowties plus the degree-ten parity monomial.")
    print("  I12 the d=5 rule is the unique alternating apex lift of d=4 rule 8 in the checked convention:")
    print("      after switching any apex v to dominate, f5(T)=(-1)^(v+1) f4(T-v).")
    print("  I13 among normalized flat triangle fields, the uniform rank-3 GP counts at n=5,6,7 are")
    print("      24,120,720=(n-1)!; the one-flipped chi_135 hostile fails at pivots 1 and 5.")
    print("  I14 after switching a labelled six-tournament so vertex 5 dominates, its f5-directed")
    print("      boundary condition is equivalent, for each of all 1024 remaining U, to U being")
    print("      f4-directed.  The common success set has size 144, so 9/64 inherits exactly.")
    print("      Rule f5 is invariant under all 16384 distinct labelled cut checks.  A universal next")
    print("      apex lift is obstructed in this recursion because 144<1024; no general nonexistence")
    print("      theorem for other d=6 rules is claimed.")

    print("\nG. Scope and hostile boundary")
    print("  The transitive T3 and cyclic C3 share the same ordinary two-vertex deletion deck:")
    print("  reconstruction failure begins at n=3.  At n=6 the four full-deck collision pairs")
    print("  are converse pairs and survive score, all proper subset profiles, card-internal")
    print("  profile losses, and unpaired boundary multisets; pairing card with boundary")
    print("  outdegree separates them.  This is finite evidence only: no all-n reconstruction")
    print("  theorem is claimed, and no n=7 exhaustive census is imported into this result.")
    print("\nALL REQUIRE GATES PASSED")


if __name__ == "__main__":
    main()
