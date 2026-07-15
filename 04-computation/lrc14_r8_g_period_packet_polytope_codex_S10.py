#!/usr/bin/env python3
"""Exact marginal-obstruction certificates for the dense THM-786 regime.

There are two complementary one-dimensional quotients of a run containing
consecutive second-fastest-owner (g) walls.

* An f-period containing g has a balanced companion cluster B.  A fractional
  transversal of the balanced-cluster hypergraph converts companion wall-grid
  capacity into an extent bound.
* A complete g-period has q or q+1 f-walls and at most one wall of each slower
  companion.  The period-sum law restricts its binary packet type.  Long runs
  force the normalized clock-rate vector to approach the convex hull of these
  allowed packet types.

This script uses only exact integer/Fraction arithmetic.  It certifies a dense
profile caught by the g-period polytope although its cluster-cover cost exceeds
g, then certifies a second profile on which both marginal tests (and the fixed-
companion span test) admit an infinite-frequency relaxation.  The second
profile is a no-go for those *marginal proof methods*, not a blocking-run or
LRC counterexample; its exact event walk has maximum covered run two.

Tournament Analysis / assumption challenge
-------------------------------------------
Putting packet types in any deterministic total order gives only a transitive
tournament (one Hamiltonian path, no directed cycles, singleton SCCs).  That
quotient forgets cluster-speed incidence, convex coefficients, and, crucially,
which f-cell and g-cell contain the same physical wall.  The theorem-bearing
object is instead a labelled two-clock incidence transport: balanced f-cell
hyperedges and allowed g-cell packets coupled by the exact Beatty event word.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256


MODULUS = 7
EXPECTED_DIGEST = "0a552ab80738d3089b8967a6f83caf842fe35f467e19d1af5a4791d02de103dc"


def inverse_mod7(value: int) -> int:
    assert value % MODULUS
    return pow(value % MODULUS, -1, MODULUS)


def selected(mask: int, companions: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        speed for index, speed in enumerate(companions) if (mask >> index) & 1
    )


def mask_of(speeds: tuple[int, ...], companions: tuple[int, ...]) -> int:
    return sum(1 << companions.index(speed) for speed in speeds)


def balanced_clusters(g: int, companions: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        mask
        for mask in range(1, 1 << len(companions))
        if (
            inverse_mod7(g)
            + sum(inverse_mod7(c) for c in selected(mask, companions))
        )
        % MODULUS
        == 0
    )


def packet_types(
    f: int, g: int, companions: tuple[int, ...]
) -> tuple[tuple[int, int], ...]:
    """Return (epsilon,mask), where an allowed g-period has q+epsilon f-walls."""
    q = f // g
    return tuple(
        (epsilon, mask)
        for epsilon in (0, 1)
        for mask in range(1 << len(companions))
        if (
            (q + epsilon) * inverse_mod7(f)
            + sum(inverse_mod7(c) for c in selected(mask, companions))
        )
        % MODULUS
        == 0
    )


def cover_cost(weights: tuple[F, ...], companions: tuple[int, ...]) -> F:
    return sum((weight * speed for weight, speed in zip(weights, companions)), F(0))


def assert_fractional_cover(
    weights: tuple[F, ...], clusters: tuple[int, ...]
) -> None:
    for cluster in clusters:
        assert sum(
            (weight for index, weight in enumerate(weights) if cluster >> index & 1),
            F(0),
        ) >= 1


def assert_fractional_packing(
    packing: tuple[tuple[int, F], ...], companions: tuple[int, ...]
) -> F:
    """Check the dual constraints sum_{B contains c} y_B <= c."""
    for index, speed in enumerate(companions):
        load = sum((weight for mask, weight in packing if mask >> index & 1), F(0))
        assert load <= speed
    return sum((weight for _, weight in packing), F(0))


def cyclic_run_length(word: tuple[int, ...], owner_index: int) -> int:
    hits = tuple(bool(mask >> owner_index & 1) for mask in word)
    if all(hits):
        return len(hits)
    best = current = 0
    for hit in hits + hits:
        current = current + 1 if hit else 0
        best = max(best, current)
    return min(best, len(hits))


def token(speed: int, x: F) -> int | None:
    y = speed * x
    quotient, remainder = divmod(y.numerator, y.denominator)
    if 2 * remainder == y.denominator:
        return None
    nearest = quotient + int(2 * remainder > y.denominator)
    return (-inverse_mod7(speed) * nearest) % MODULUS


def exact_covered_runs(speeds: tuple[int, ...]) -> tuple[tuple[tuple[F, int], ...], ...]:
    events: dict[F, list[int]] = defaultdict(list)
    for speed in speeds:
        for wall_index in range(speed):
            events[F(2 * wall_index + 1, 2 * speed)].append(speed)

    runs: list[tuple[tuple[F, int], ...]] = []
    current: list[tuple[F, int]] = []
    for x in sorted(events):
        owners = events[x]
        covered = False
        if len(owners) == 1:
            owner = owners[0]
            other_tokens = [token(speed, x) for speed in speeds if speed != owner]
            covered = None not in other_tokens and set(other_tokens) == set(range(7))
        if covered:
            current.append((x, owners[0]))
        else:
            if current:
                runs.append(tuple(current))
            current = []
    if current:
        runs.append(tuple(current))
    return tuple(sorted(runs, key=lambda run: (-len(run), run)))


def main() -> None:
    # ------------------------------------------------------------------
    # Profile A: cluster marginals are inconclusive, but g-period packets
    # have an exact cardinality separator.
    # ------------------------------------------------------------------
    f_a, g_a = 65, 64
    companions_a = (26, 33, 40, 47, 54, 61)
    clusters_a = balanced_clusters(g_a, companions_a)
    types_a = packet_types(f_a, g_a, companions_a)

    assert inverse_mod7(g_a) == 1
    assert inverse_mod7(f_a) == 4
    assert tuple(inverse_mod7(c) for c in companions_a) == (3,) * 6
    assert len(clusters_a) == 15
    assert all(mask.bit_count() == 2 for mask in clusters_a)
    assert len(types_a) == 21
    assert all(mask.bit_count() == 1 + epsilon for epsilon, mask in types_a)

    # Primal cover lambda_c=1/2.  The exact fractional-matching packing below
    # has incident load c at every vertex, so weak duality proves optimality.
    cover_a = (F(1, 2),) * 6
    assert_fractional_cover(cover_a, clusters_a)
    cover_cost_a = cover_cost(cover_a, companions_a)
    matching_edges = (
        (0, 3, F(31, 2)),
        (0, 4, F(21, 2)),
        (1, 4, F(33)),
        (2, 4, F(21, 2)),
        (2, 5, F(59, 2)),
        (3, 5, F(63, 2)),
    )
    packing_a = tuple((1 << i | 1 << j, weight) for i, j, weight in matching_edges)
    packing_value_a = assert_fractional_packing(packing_a, companions_a)
    assert cover_cost_a == packing_value_a == F(261, 2) > g_a

    # Every allowed packet has |D|-epsilon=1.  The target clock-rate vector has
    # value sum(c/g)-theta=65/16, so eta=49/16.  Coordinate lattice discrepancy
    # is <1/T over T complete g-periods; the separator l1 norm is seven.
    theta_a = F(f_a % g_a, g_a)
    target_separator_a = sum((F(c, g_a) for c in companions_a), F(0)) - theta_a
    allowed_separator_a = F(1)
    eta_a = target_separator_a - allowed_separator_a
    separator_l1_a = 7
    period_bound_a = F(separator_l1_a, 1) / eta_a
    assert theta_a == F(1, 64)
    assert target_separator_a == F(65, 16)
    assert eta_a == F(49, 16)
    assert period_bound_a == F(16, 7)

    # ------------------------------------------------------------------
    # Profile B: an exact no-go for both separate marginal quotients.
    # ------------------------------------------------------------------
    f_b, g_b = 69, 29
    companions_b = (4, 5, 12, 13, 16, 27)
    clusters_b = balanced_clusters(g_b, companions_b)
    types_b = packet_types(f_b, g_b, companions_b)
    cluster_sets_b = tuple(selected(mask, companions_b) for mask in clusters_b)

    assert tuple(inverse_mod7(c) for c in companions_b) == (2, 3, 3, 6, 4, 6)
    assert inverse_mod7(f_b) == 6 and inverse_mod7(g_b) == 1
    assert len(clusters_b) == 9
    assert len(types_b) == 18

    # Exact primal/dual certificate for fractional cluster-cover optimum 49.
    cover_b = tuple(F(int(c in (4, 5, 13, 27))) for c in companions_b)
    assert_fractional_cover(cover_b, clusters_b)
    cover_cost_b = cover_cost(cover_b, companions_b)
    dual_b = (
        (mask_of((13,), companions_b), F(13)),
        (mask_of((27,), companions_b), F(27)),
        (mask_of((5, 12), companions_b), F(5)),
        (mask_of((4, 16), companions_b), F(4)),
    )
    assert all(mask in clusters_b for mask, _ in dual_b)
    dual_value_b = assert_fractional_packing(dual_b, companions_b)
    assert cover_cost_b == dual_value_b == 49 > g_b

    # A 29-packet cyclic balanced-cluster word respects all c/g capacities and
    # every corrected fixed-companion span.  Five triple packets separate
    # singleton-27 blocks of lengths 5,5,5,5,4.
    singleton_27 = mask_of((27,), companions_b)
    triple_5_13_16 = mask_of((5, 13, 16), companions_b)
    assert singleton_27 in clusters_b and triple_5_13_16 in clusters_b
    cluster_word: list[int] = []
    for block_length in (5, 5, 5, 5, 4):
        cluster_word.append(triple_5_13_16)
        cluster_word.extend([singleton_27] * block_length)
    cluster_word_b = tuple(cluster_word)
    assert len(cluster_word_b) == g_b

    span_data_b = []
    for index, companion in enumerate(companions_b):
        incidence = sum(bool(mask >> index & 1) for mask in cluster_word_b)
        assert F(incidence, len(cluster_word_b)) <= F(companion, g_b)
        max_run = cyclic_run_length(cluster_word_b, index)
        span_bound = F(1) + F(2 * g_b * companion, f_b * (g_b - companion))
        if incidence:
            assert max_run < span_bound
        span_data_b.append((companion, incidence, max_run, span_bound))

    # The signed visitor-set difference law is automatic for every adjacent
    # pair because both full packets {g}+B have inverse sum zero.
    for left, right in zip(cluster_word_b, cluster_word_b[1:] + cluster_word_b[:1]):
        entrants = right & ~left
        leavers = left & ~right
        entrant_sum = sum(
            inverse_mod7(c)
            for index, c in enumerate(companions_b)
            if entrants >> index & 1
        ) % MODULUS
        leaver_sum = sum(
            inverse_mod7(c)
            for index, c in enumerate(companions_b)
            if leavers >> index & 1
        ) % MODULUS
        assert entrant_sum == leaver_sum

    # Exact convex decomposition of the target g-period rate vector.  The
    # integer coefficient is the numerator over the common denominator g=29.
    decomposition_b = (
        (0, (5, 13), 2),
        (0, (5, 27), 1),
        (0, (12, 27), 6),
        (0, (5, 12, 16, 27), 2),
        (0, (13, 16, 27), 7),
        (1, (4, 12, 13, 27), 4),
        (1, (16, 27), 7),
    )
    allowed_b = set(types_b)
    assert sum(coefficient for _, _, coefficient in decomposition_b) == g_b
    assert sum(epsilon * coefficient for epsilon, _, coefficient in decomposition_b) == 11
    for epsilon, packet, _ in decomposition_b:
        assert (epsilon, mask_of(packet, companions_b)) in allowed_b
    for companion in companions_b:
        marginal_numerator = sum(
            coefficient
            for _, packet, coefficient in decomposition_b
            if companion in packet
        )
        assert marginal_numerator == companion

    # This is not a run certificate: exact event chronology destroys the
    # relaxed construction almost immediately.
    speeds_b = companions_b + (g_b, f_b)
    runs_b = exact_covered_runs(speeds_b)
    max_exact_run_b = len(runs_b[0]) if runs_b else 0
    longest_b = tuple(run for run in runs_b if len(run) == max_exact_run_b)
    assert max_exact_run_b == 2

    digest = sha256()
    for label, rows in (
        ("CA", tuple((mask, selected(mask, companions_a)) for mask in clusters_a)),
        ("PA", types_a),
        ("CB", tuple((mask, selected(mask, companions_b)) for mask in clusters_b)),
        ("PB", types_b),
        ("DB", decomposition_b),
        ("WB", cluster_word_b),
        ("RB", longest_b),
    ):
        digest.update(f"{label}|{rows!r}\n".encode())
    hexdigest = digest.hexdigest()

    print("THM-786 dense marginal packet-polytope certificate")
    print("profile A: g-period separator beyond cluster density")
    print(f"f={f_a} g={g_a} companions={companions_a} sum={sum(companions_a)}")
    print(f"inverse_residues={(inverse_mod7(f_a), inverse_mod7(g_a)) + tuple(inverse_mod7(c) for c in companions_a)} order=(f,g,C)")
    print(f"balanced_cluster_count={len(clusters_a)} characterization=all_2_subsets")
    print(f"fractional_cluster_cover_optimum={cover_cost_a} dual_value={packing_value_a} > g={g_a}")
    print(f"g_packet_type_count={len(types_a)} epsilon_counts={(sum(e == 0 for e, _ in types_a), sum(e == 1 for e, _ in types_a))}")
    print("allowed_packet_characterization=|D|=1+epsilon")
    print(f"separator=|D|-epsilon<=1 target={target_separator_a} eta={eta_a} l1={separator_l1_a}")
    print(f"complete_g_period_bound=T<{period_bound_a} hence integer_T<=2")
    print()
    print("profile B: exact no-go for separate marginal arguments")
    print(f"f={f_b} g={g_b} companions={companions_b}")
    print(f"inverse_residues={(inverse_mod7(f_b), inverse_mod7(g_b)) + tuple(inverse_mod7(c) for c in companions_b)} order=(f,g,C)")
    print(f"balanced_cluster_count={len(clusters_b)} clusters={cluster_sets_b}")
    print(f"fractional_cluster_cover_optimum={cover_cost_b} dual_value={dual_value_b} > g={g_b}")
    print("cluster_frequency_relaxation={27}:24/29, {5,13,16}:5/29")
    print("fixed_span_checks=" + str(tuple(span_data_b)))
    print("signed_difference_law=PASS on all 29 cyclic transitions")
    print(f"g_packet_type_count={len(types_b)} epsilon_counts={(sum(e == 0 for e, _ in types_b), sum(e == 1 for e, _ in types_b))}")
    print("exact_target_decomposition_numerators_over_29=" + str(decomposition_b))
    print("target_marginals=(epsilon:11/29,C:4/29,5/29,12/29,13/29,16/29,27/29)")
    print(f"exact_event_walk_max_covered_run={max_exact_run_b}")
    print("longest_exact_runs=" + str(tuple(tuple((str(x), owner) for x, owner in run) for run in longest_b)))
    print("interpretation=marginal no-go only; exact f/g cell incidence and event order remain essential")
    print(f"certificate_digest={hexdigest}")
    print()
    print("Tournament Analysis / assumption challenge:")
    print("  vertices: allowed g-period packet types (21 in A, 18 in B)")
    print("  gauge: (epsilon, packet mask), yielding a transitive total-order tournament")
    print("  A fingerprint: scores=0..20 cycles=0 sccs=21 hp_count=1")
    print("  B fingerprint: scores=0..17 cycles=0 sccs=18 hp_count=1")
    print("  information lost: hypergraph incidence, convex weights, f-cell/g-cell wall identity")
    print("  exact carrier: labelled two-clock incidence transport + Beatty event word")

    if EXPECTED_DIGEST != "TO_BE_FILLED":
        assert hexdigest == EXPECTED_DIGEST
    print("PASS")


if __name__ == "__main__":
    main()
