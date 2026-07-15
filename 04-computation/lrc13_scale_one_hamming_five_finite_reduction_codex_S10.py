#!/usr/bin/env python3
"""Exact replay for THM-820, the scale-one Hamming-five finite reduction.

The theorem has an infinite symbolic part and three finite exact audits.

* Seven retained speeds have a 1/8 lonely point.  The resulting safe interval
  first forces x <= 416; the exact minimal-cycle bank sharpens this to 388.
  The collar graph then gives either a recursive doubling box or an induced
  exceptional four-cycle.  Strict reciprocal mass sharpens the exceptional
  anchor to 228.  The lower bounds from both the seven- and eight-speed safe
  intervals then have an exact upper envelope v <= 1986 for the second finite
  box.  ``finite_reduction_constants`` checks every rational constant.
* ``height_one_censuses`` evaluates all C(12,5)=792 radius-five rows and all
  C(12,6)=924 radius-six rows.  Maxima are exact: a maximum of a minimum of
  triangular waves occurs at a self-cusp (denominator 2u) or a pair crossing
  (denominator u+v or |u-v|).  Every comparison below is integer-only.
* ``five_cycle_audit`` enumerates the complete integer-centre bank for a
  spanning collar five-cycle and exhibits chordless realizations of every
  feasible band type.

Tournament Analysis / assumption challenge
------------------------------------------
The local pair observable is oriented left-handoff at an owner exit.  Silent
pairs are completed by increasing or decreasing label to obtain two gauges
and a tie Hamiltonian path.  This shadow is telemetry: two examples below
have the same labelled five-cycle and tournament fingerprints but closing
bands 9 and 22 and different maximin values.  The theorem-bearing carrier is
the typed incidence

    (owner exit, provider tooth, z, integer centre k, speed ratio, half-open flag)

followed by the dynamic strict-safe-component / remaining-tooth incidence.
That component carrier preserves exact emptiness of the final safe set.  Bare
runner, residue, event, and tournament vertices lose width or scale.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
from math import comb
from struct import pack
from typing import Iterable, Sequence


P = 13
BASE = tuple(range(1, P))
DELTA = F(1, P)


def fmt(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def exact_maximin(speeds: Sequence[int]) -> tuple[F, tuple[F, ...]]:
    """Return M(speeds) and every maximizer, using exact integer crossings."""

    denominators = {2 * u for u in speeds}
    for u, v in combinations(speeds, 2):
        denominators.add(u + v)
        denominators.add(abs(u - v))
    denominators.discard(0)

    best_num, best_den = 0, 1
    witnesses: set[F] = set()
    for denominator in sorted(denominators):
        for numerator in range(1, denominator):
            value_num = min(
                min((u * numerator) % denominator,
                    denominator - (u * numerator) % denominator)
                for u in speeds
            )
            comparison = value_num * best_den - best_num * denominator
            if comparison > 0:
                best_num, best_den = value_num, denominator
                witnesses = {F(numerator, denominator)}
            elif comparison == 0 and value_num:
                witnesses.add(F(numerator, denominator))
    return F(best_num, best_den), tuple(sorted(witnesses))


def digest_u64(digest: object, value: int) -> None:
    digest.update(pack(">Q", value))


def height_one_census(radius: int) -> tuple[
    list[tuple[F, tuple[int, ...], tuple[int, ...], tuple[F, ...]]], str
]:
    rows = []
    digest = sha256()
    for missing in combinations(BASE, radius):
        packet = tuple(r + P if r in missing else r for r in BASE)
        maximum, witnesses = exact_maximin(packet)
        rows.append((maximum, missing, packet, witnesses))
        for value in missing + (maximum.numerator, maximum.denominator, len(witnesses)):
            digest_u64(digest, value)
        for witness in witnesses:
            digest_u64(digest, witness.numerator)
            digest_u64(digest, witness.denominator)
    return rows, digest.hexdigest()


def height_one_censuses() -> None:
    rows5, digest5 = height_one_census(5)
    minimum5 = min(row[0] for row in rows5)
    minimizers5 = [row for row in rows5 if row[0] == minimum5]
    expected_missing5 = (1, 3, 5, 7, 9)
    expected_packet5 = (14, 2, 16, 4, 18, 6, 20, 8, 22, 10, 11, 12)
    expected_witnesses5 = tuple(F(a, 24) for a in (1, 5, 7, 17, 19, 23))
    assert len(rows5) == comb(12, 5) == 792
    assert minimum5 == F(1, 12)
    assert minimizers5 == [
        (F(1, 12), expected_missing5, expected_packet5, expected_witnesses5)
    ]
    assert digest5 == "323913c07710f8321d29d07c9dfdfb4d13107549063587584b4ee0ceab923dc9"

    value_histogram5 = Counter(row[0] for row in rows5)
    witness_histogram5 = Counter(len(row[3]) for row in rows5)
    threshold_counts = {
        threshold: sum(count for value, count in value_histogram5.items()
                       if value <= threshold)
        for threshold in (F(1, 12), F(1, 11), F(3, 31), F(1, 10), F(1, 9))
    }
    assert len(value_histogram5) == 66
    assert threshold_counts == {
        F(1, 12): 1, F(1, 11): 3, F(3, 31): 8, F(1, 10): 11, F(1, 9): 65
    }
    assert witness_histogram5 == {2: 589, 4: 100, 6: 80, 8: 18, 10: 4, 18: 1}

    # Structural interpretation of the unique minimum.
    sorted_minimum_packet = tuple(sorted(expected_packet5))
    assert sorted_minimum_packet == tuple(sorted(tuple(range(2, 23, 2)) + (11,)))
    ap11_witnesses = tuple(F(a, 24) for a in (1, 5, 7, 11, 13, 17, 19, 23))
    assert exact_maximin(tuple(range(2, 23, 2)))[0] == F(1, 12)
    assert all(min((11 * t) % 1, 1 - (11 * t) % 1) >= F(1, 12)
               for t in expected_witnesses5)
    assert tuple(t for t in ap11_witnesses
                 if min((11 * t) % 1, 1 - (11 * t) % 1) >= F(1, 12)) \
        == expected_witnesses5

    rows6, digest6 = height_one_census(6)
    minimum6 = min(row[0] for row in rows6)
    minimizers6 = [row for row in rows6 if row[0] == minimum6]
    expected_missing6 = (1, 3, 5, 7, 9, 11)
    expected_packet6 = (14, 2, 16, 4, 18, 6, 20, 8, 22, 10, 24, 12)
    expected_witnesses6 = tuple(F(a, 26) for a in range(1, 26) if a != 13)
    assert len(rows6) == comb(12, 6) == 924
    assert minimum6 == DELTA
    assert minimizers6 == [
        (DELTA, expected_missing6, expected_packet6, expected_witnesses6)
    ]
    assert tuple(sorted(expected_packet6)) == tuple(range(2, 25, 2))
    assert digest6 == "c41767013f5403a00004609ba31b165896a9a94416e0693393aecd9915737ccc"

    print("HEIGHT_ONE_HAMMING_FIVE_CENSUS")
    print(f"rows={len(rows5)} value_count={len(value_histogram5)} minimum={fmt(minimum5)} unique=1")
    print("missing_labels=(1,3,5,7,9) packet=2*[11]+{11}")
    print("maximizers={1,5,7,17,19,23}/24")
    print("lowest_value_counts={1/12:1,1/11:2,3/31:5,1/10:3}")
    print("cumulative_counts_through={1/12:1,1/11:3,3/31:8,1/10:11,1/9:65}")
    print(f"witness_count_histogram={dict(sorted(witness_histogram5.items()))}")
    print(f"certificate_digest={digest5}")
    print()
    print("HEIGHT_ONE_RADIUS_SIX_METHOD_BOUNDARY")
    print(f"rows={len(rows6)} minimum={fmt(minimum6)} unique=1")
    print("missing_labels=(1,3,5,7,9,11) packet=2*[12]")
    print("maximizers={a/26:1<=a<=25,a!=13}")
    print(f"certificate_digest={digest6}")
    print("interpretation=exact_AP_scale_change_not_a_radius_six_closure")
    print()


def finite_reduction_constants() -> None:
    # Seven retained speeds: M >= 1/8.  With m=max(P), the delta-safe interval
    # has length 2*(1/8-1/13)/m = 5/(52m).
    seven_gap = F(1, 8) - DELTA
    seven_length_coefficient = 2 * seven_gap
    assert seven_gap == F(5, 104)
    assert seven_length_coefficient == F(5, 52)

    # Five danger combs contribute at most 10L/13 + (2/13)sum 1/u.
    # Coverage therefore requires sum 1/u >= (3/2)L = 15/(104m).
    five_bulk = F(10, 13)
    five_reciprocal_coefficient = (1 - five_bulk) * F(13, 2)
    required_five_reciprocal = five_reciprocal_coefficient * seven_length_coefficient
    assert five_reciprocal_coefficient == F(3, 2)
    assert required_five_reciprocal == F(15, 104)
    least_multiplier = F(5, 1) / required_five_reciprocal
    assert least_multiplier == F(104, 3)
    assert (104 * 12) // 3 == 416

    # If x>24, every owner has positive collar indegree.  A minimal cycle has
    # length four or five.  A four-cycle gives reciprocal sum <=9/(2x); the
    # exact five-cycle bank has a centre at least four, hence one speed at
    # least 3x/2 and reciprocal sum <=14/(3x).  The latter is the weaker cap.
    all_large_reciprocal_upper = F(14, 3)
    sharpened_least_multiplier = all_large_reciprocal_upper / required_five_reciprocal
    assert sharpened_least_multiplier == F(1456, 45)
    sharpened_x_cap = (1456 * 12) // 45
    assert sharpened_x_cap == 388

    # Exceptional branch: P union {x} has eight speeds and M >= 1/9.  Since
    # x>=14>max(P), its maximum is x and its safe length is 8/(117x).
    eight_gap = F(1, 9) - DELTA
    eight_length_coefficient = 2 * eight_gap
    assert eight_gap == F(4, 117)
    assert eight_length_coefficient == F(8, 117)

    # Four top combs leave coefficient 5/13, hence require
    # sum_top 1/u >= (5/2)L = 20/(117x).  The exceptional k=5 edge makes one
    # speed at least twice another, so with v=min(top), sum_top 1/u<=7/(2v).
    four_bulk = F(8, 13)
    four_reciprocal_coefficient = (1 - four_bulk) * F(13, 2)
    required_four_reciprocal = four_reciprocal_coefficient * eight_length_coefficient
    assert four_reciprocal_coefficient == F(5, 2)
    assert required_four_reciprocal == F(20, 117)
    exceptional_multiplier = F(7, 2) / required_four_reciprocal
    assert exceptional_multiplier == F(819, 40)

    # In the v>2x branch, sum_all=1/x+sum_top and
    # sum_top<=7/(2v)<7/(4x), so sum_all<11/(4x).  The inequality is strict.
    exceptional_all_upper = F(11, 4)
    exceptional_x_multiplier = exceptional_all_upper / required_five_reciprocal
    assert exceptional_x_multiplier == F(286, 15)
    exceptional_x_cap = (286 * 12 - 1) // 15
    assert exceptional_x_cap == 228
    exceptional_v_cap = (819 * exceptional_x_cap) // 40
    assert exceptional_v_cap == 4668
    assert 4 * exceptional_v_cap == 18672

    # The seven-core mass also leaves a top-only lower bound
    # 15/(104m)-1/x after subtracting the anchor reciprocal.  When positive,
    # combine it with sum_top<=7/(2v).  Audit the exact integer envelope of
    # that cap and the eight-core cap over the full uniform superset box.
    envelope = -1
    envelope_rows = []
    for m in range(7, 13):
        for x in range(14, exceptional_x_cap + 1):
            eight_core_cap = (819 * x) // 40
            residual_lower = F(15, 104 * m) - F(1, x)
            if residual_lower > 0:
                ratio_cap_fraction = F(7, 2) / residual_lower
                ratio_cap = ratio_cap_fraction.numerator // ratio_cap_fraction.denominator
                cap = min(eight_core_cap, ratio_cap)
            else:
                ratio_cap = None
                cap = eight_core_cap
            row = (m, x, eight_core_cap, ratio_cap, residual_lower)
            if cap > envelope:
                envelope = cap
                envelope_rows = [row]
            elif cap == envelope:
                envelope_rows.append(row)
    assert envelope == 1986
    assert envelope_rows == [(12, 97, 1986, 2046, F(69, 40352))]

    crossover = F(4384, 45)
    eight_at_cross = F(819, 40) * crossover
    residual_at_cross = F(15, 104 * 12) - 1 / crossover
    assert eight_at_cross == F(7, 2) / residual_at_cross
    assert 97 < crossover < 98

    residual_98 = F(15, 104 * 12) - F(1, 98)
    cap_98_fraction = F(7, 2) / residual_98
    cap_98 = cap_98_fraction.numerator // cap_98_fraction.denominator
    assert residual_98 == F(37, 20384) and cap_98 == 1928

    print("UNIFORM_TWO_BOX_CONSTANTS")
    print("seven_core_threshold=1/8 gap=5/104 safe_length=5/(52m)")
    print("five_comb_cover_requires=sum(1/u)>=15/(104m)")
    print("preliminary_least_bound=x<=floor(104m/3)<=416")
    print("all_large_minimal_cycle=length4_sum<=9/(2x)_or_length5_sum<=14/(3x)")
    print("sharpened_least_bound=x<=floor(1456m/45)<=388")
    print("doubling_box=v<=2x,w<=2v,y<=2w,z<=2y; numeric_caps=(388,776,1552,3104,6208)")
    print("exceptional_trigger=v>2x => induced_top_four_positive_indegree")
    print("exceptional_cycle=THM815_a*{1,2,4,8} spread<=4")
    print("exceptional_global_mass=sum_all<11/(4x) => x<286m/15 => integer_x<=228")
    print("eight_core_threshold=1/9 gap=4/117 safe_length=8/(117x)")
    print("four_comb_cover_requires=sum_top(1/u)>=20/(117x)")
    print("cycle_reciprocal_upper=sum_top(1/u)<=7/(2v)")
    print("preliminary_exceptional_box=v<=floor(819x/40)<=4668")
    print("seven_core_top_residual=sum_top>=15/(104m)-1/x")
    print("x_dependent_cap=v<=min(floor(819x/40),floor((7/2)/(15/(104m)-1/x)))_when_positive")
    print("exact_envelope=m12_x97_caps=(1986,2046)_min1986; x98_second_cap=1928")
    print("exceptional_box=v<=1986 max_top<=4v<=7944")
    print()


def inverse_mod13(x: int) -> int:
    return pow(x % P, -1, P)


def left_handoff(provider: int, owner: int) -> tuple[bool, int | None]:
    """Return eligibility and its integer centre k=z-13m."""

    z = (provider % P) * inverse_mod13(owner) % P
    numerator = z * owner - 2 * provider
    quotient = numerator // (P * owner)
    for m in range(quotient - 2, quotient + 3):
        scaled = z * owner - 2 * provider - P * m * owner
        if -owner < scaled <= owner:
            return True, z - P * m
    return False, None


def canonical_rotation(word: tuple[int, ...]) -> tuple[int, ...]:
    return min(word[i:] + word[:i] for i in range(len(word)))


def cycle_labels(word: tuple[int, ...]) -> tuple[int, ...]:
    labels = [1]
    for centre in word[:-1]:
        labels.append(labels[-1] * inverse_mod13(centre) % P)
    return tuple(labels)


def label_orbit_representative(labels: Iterable[int]) -> tuple[int, ...]:
    labels = tuple(labels)
    return min(tuple(sorted(a * x % P for x in labels)) for a in BASE)


def strongly_connected_component_sizes(edge: Sequence[Sequence[bool]]) -> tuple[int, ...]:
    n = len(edge)
    reach = [list(row) for row in edge]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    used: set[int] = set()
    sizes = []
    for i in range(n):
        if i not in used:
            component = {j for j in range(n) if reach[i][j] and reach[j][i]}
            used.update(component)
            sizes.append(len(component))
    return tuple(sorted(sizes, reverse=True))


def tournament_fingerprint(edge: Sequence[Sequence[bool]]) -> tuple[
    tuple[tuple[int, int], ...], int, tuple[int, ...], int, tuple[int, ...]
]:
    n = len(edge)
    score = tuple(sorted(Counter(sum(row) for row in edge).items()))
    triangles = 0
    for i, j, k in combinations(range(n), 3):
        triangles += bool(
            (edge[i][j] and edge[j][k] and edge[k][i])
            or (edge[j][i] and edge[i][k] and edge[k][j])
        )
    hamiltonian_paths = []
    for path in permutations(range(n)):
        if all(edge[path[i]][path[i + 1]] for i in range(n - 1)):
            hamiltonian_paths.append(path)
    return (
        score,
        triangles,
        strongly_connected_component_sizes(edge),
        len(hamiltonian_paths),
        hamiltonian_paths[0],
    )


def complete_tournament(
    live: Sequence[Sequence[bool]], labels: Sequence[int], reverse_ties: bool
) -> tuple[list[list[bool]], int]:
    edge = [list(row) for row in live]
    silent = 0
    for i, j in combinations(range(len(labels)), 2):
        assert not (edge[i][j] and edge[j][i])
        if not edge[i][j] and not edge[j][i]:
            silent += 1
            forward = (labels[i] < labels[j]) ^ reverse_ties
            edge[i][j], edge[j][i] = forward, not forward
    return edge, silent


def live_graph(speeds: Sequence[int]) -> list[list[bool]]:
    n = len(speeds)
    edge = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                edge[i][j] = left_handoff(speeds[i], speeds[j])[0]
    return edge


def five_cycle_audit() -> None:
    # In a five-cycle every lambda>=1/2 and product lambda=1, so lambda<=16.
    # Distinct labels exclude centres congruent to 0 or 1 modulo 13.
    centres = tuple(k for k in range(2, 34) if k % P not in (0, 1))
    survivors = []
    rows = 0
    for word in combinations_with_replacement(centres, 5):
        rows += 1
        lower = upper = residue_product = 1
        for k in word:
            lower *= k - 1
            upper *= k + 1
            residue_product = residue_product * (k % P) % P
        if lower <= 32 < upper and residue_product == 1:
            survivors.append(word)
    numerical_survivors = [
        (2, 2, 2, 2, 9),
        (2, 2, 2, 2, 22),
        (2, 2, 2, 3, 6),
        (2, 2, 2, 4, 11),
        (2, 2, 3, 3, 4),
        (2, 2, 3, 5, 5),
    ]
    assert len(centres) == 28
    assert rows == comb(32, 5) == 201_376
    assert survivors == numerical_survivors

    # The last numerical word has product of lower endpoints exactly 32.
    # Therefore every edge is at its lower endpoint; its k=3 edge has
    # lambda=1, impossible for two distinct residue labels.
    impossible = numerical_survivors[-1]
    assert __import__("math").prod(k - 1 for k in impossible) == 32
    feasible_types = numerical_survivors[:-1]
    assert all(max(word) >= 4 for word in feasible_types)

    cyclic_words: set[tuple[int, ...]] = set()
    label_orbits: set[tuple[int, ...]] = set()
    for multiset in feasible_types:
        for word in set(permutations(multiset)):
            word = canonical_rotation(word)
            labels = cycle_labels(word)
            assert len(set(labels)) == 5
            cyclic_words.add(word)
            label_orbits.add(label_orbit_representative(labels))
    assert len(cyclic_words) == 16
    assert len(label_orbits) == 15

    # One chordless all-large realization for each feasible band multiset.
    examples = [
        ((2, 2, 2, 2, 9), (1, 7, 10, 5, 9), (105, 137, 270, 447, 438), F(21, 109)),
        ((2, 2, 2, 2, 22), (1, 7, 10, 5, 9), (105, 189, 374, 681, 1114), F(107, 559)),
        ((2, 2, 2, 3, 6), (1, 7, 10, 5, 6), (105, 163, 205, 317, 266), F(27, 137)),
        ((2, 2, 2, 4, 11), (1, 7, 10, 5, 11), (105, 202, 400, 798, 531), F(21, 113)),
        ((2, 2, 3, 3, 4), (1, 7, 10, 12, 4), (105, 176, 270, 246, 173), F(9, 62)),
    ]
    fingerprints = []
    for word, labels, speeds, expected_maximum in examples:
        assert all(u > 24 and u % P == label for u, label in zip(speeds, labels))
        live = live_graph(speeds)
        expected_edges = {(i, (i + 1) % 5) for i in range(5)}
        actual_edges = {(i, j) for i in range(5) for j in range(5) if live[i][j]}
        assert actual_edges == expected_edges
        centres_seen = tuple(left_handoff(speeds[i], speeds[(i + 1) % 5])[1]
                             for i in range(5))
        assert centres_seen == word
        core = tuple(x for x in BASE if x not in labels)
        maximum, _ = exact_maximin(core + speeds)
        assert maximum == expected_maximum
        forward, silent_f = complete_tournament(live, labels, False)
        reverse, silent_r = complete_tournament(live, labels, True)
        assert silent_f == silent_r == 5
        fingerprints.append((tournament_fingerprint(forward),
                             tournament_fingerprint(reverse)))

    # Same labelled live cycle and tournament telemetry; radically different
    # final band centre and different exact LRC value.
    assert examples[0][1] == examples[1][1]
    assert fingerprints[0] == fingerprints[1]
    assert examples[0][0][-1] == 9 and examples[1][0][-1] == 22
    assert examples[0][3] != examples[1][3]

    print("SPANNING_FIVE_CYCLE_BAND_AUDIT")
    print(f"allowed_centres={centres} sorted_multisets={rows}")
    print("numerical_survivors=(22229,2222-22,22236,2224-11,22334,22355)")
    print("22355_status=impossible_distinct_labels_lower_product_32_forces_k3_ratio_1")
    print("feasible_band_types=5 cyclic_words_up_to_rotation=16 multiplicative_label_orbits=15")
    for word, labels, speeds, maximum in examples:
        print(f"word={word} labels={labels} speeds={speeds} live_edges=5 M={fmt(maximum)}")
    f0, r0 = fingerprints[0]
    print("band9_vs_band22=same_labels_same_live_cycle_same_tournament_different_metric")
    print(f"forward_fingerprint=score{dict(f0[0])},triangles{f0[1]},scc{f0[2]},HP{f0[3]}")
    print(f"reverse_fingerprint=score{dict(r0[0])},triangles{r0[1]},scc{r0[2]},HP{r0[3]}")
    print("shared_tie_Hamiltonian_path=(1,7,10,5,9)")
    print()


def persistent_four_cycle_family() -> None:
    labels = (1, 2, 3, 4, 8)
    top_base = {1: 79, 2: 54, 4: 30, 8: 34}
    expected_top_edges = {(1, 8), (8, 4), (4, 2), (2, 1)}
    expected_anchor_hits = (
        (1,), (1, 8), (4,), (), (), (8,), (4,), ()
    )
    silent_counts = []
    for n in range(16):
        c = 1 + P * n
        speeds = tuple(16 if label == 3 else c * top_base[label] for label in labels)
        live = live_graph(speeds)
        top_edges = {
            (labels[i], labels[j])
            for i in range(5) for j in range(5)
            if labels[i] != 3 and labels[j] != 3 and live[i][j]
        }
        assert top_edges == expected_top_edges
        assert strongly_connected_component_sizes(live) == (4, 1)
        hits = tuple(label for i, label in enumerate(labels) if label != 3 and live[i][2])
        assert hits == expected_anchor_hits[n % 8]
        assert not any(live[2][j] for j in range(5) if j != 2)

        forward, silent_f = complete_tournament(live, labels, False)
        reverse, silent_r = complete_tournament(live, labels, True)
        expected_fp = (((1, 2), (2, 1), (3, 2)), 3, (5,), 9)
        fp_f = tournament_fingerprint(forward)
        fp_r = tournament_fingerprint(reverse)
        assert fp_f[:4] == expected_fp and fp_r[:4] == expected_fp
        if n < 8:
            silent_counts.append(silent_f)
        assert silent_f == silent_r

    packet0 = (5, 6, 7, 9, 10, 11, 12, 79, 54, 16, 30, 34)
    assert exact_maximin(packet0)[0] == F(5, 21)
    assert silent_counts == [5, 4, 5, 6, 6, 5, 5, 6]

    print("PERSISTENT_TOP_FOUR_SCC_FAMILY")
    print("labels=(1,2,3,4,8) anchor_u3=16 top=(79,54,30,34)*(1+13N)")
    print("top_live_cycle=1->8->4->2->1 full_live_SCC_sizes=(4,1)")
    print("top_to_anchor_hits_by_N_mod8=((1),(1,8),(4),(),(),(8),(4),())")
    print(f"silent_pair_counts_by_N_mod8={tuple(silent_counts)}")
    print("both_completed_gauges=score{1:2,2:1,3:2},triangles3,SCC(5),HP9")
    print("N0_packet_M=5/21")
    print("interpretation=collar_graph_alone_allows_unbounded_top_scale; reciprocal_component_carrier_restores_finiteness")
    print()


def main() -> None:
    print("LRC13 SCALE-ONE HAMMING-FIVE FINITE REDUCTION - EXACT REPLAY")
    print("theorem=THM-820 scope=finite_reduction_not_full_Hamming-five_closure")
    print("renumber_note=THM-819_was_live-assigned_to_primitive_harmonic_law")
    print("arithmetic=integer_maximin_and_endpoint_bands_plus_exact_Fraction_constants")
    print()
    finite_reduction_constants()
    height_one_censuses()
    five_cycle_audit()
    persistent_four_cycle_family()
    print("RESULT")
    print("height_one_H5_all_loose=792/792 unique_minimum=1/12")
    print("uniform_scale_one_H5_search_is_finite_in_two_explicit_boxes")
    print("full_arbitrary_height_H5_closure=OPEN")


if __name__ == "__main__":
    main()
