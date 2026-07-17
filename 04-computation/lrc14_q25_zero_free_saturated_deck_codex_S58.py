#!/usr/bin/env python3
"""Exact audit of structural repairs to the false LRC(14) q<=25 claim.

THM-762/764 says that, for a covering family and 15 <= q <= 28, a
q-witness exists exactly when q has no zero owner and the signed unit-pair
deck misses a card.  The older counterexamples left open the tempting repair
"exclude zero loading and coherent/common-prime packets".  This file refutes
that repair with the exact row

    S0 = {43,55,61,70,73,79,83,99,103,104,109,113,156}.

Every q=15,...,25 is zero-free and has a complete signed deck.  The row is
primitive and covering, passes the uncapped S312 band-residual predicate,
has diameter 113, ratio 156/43 < 4, no four-term arithmetic progression,
all leave-one-out gcds one, and maximum common-prime packet three.  Three is
best possible for *any* zero-free covering row: owners of 8, 12, and 10 must
be distinct, since sharing would create a 24-, 20-, or 15-multiple.

There is also an infinite exact obstruction orbit.  Put

    L = lcm(2,...,25),       S_k = {s+kL : s in S0}.

For every k>=0, S_k has the same covering and signed-deck incidence through
q=25, is primitive, has diameter 113, and has maximum prime packet exactly
three.  Its speed ratio tends to one, while the explicit time

    t_k = 1 / ((43+kL)+(156+kL))

has clearance (43+kL)/(199+2kL), tending to 1/2.  Thus bounded-q failure can
coexist with arbitrarily strong metric looseness.  The mechanism is a
CRT-compatible irredundant signed-card covering code, not divisor loading.

The strongest honest local replacement is the collision-surplus form of the
signed-deck theorem.  If

    N_q = # {s : gcd(s,q)=1},
    C_q = N_q - |B_q(S)|,
    h_q = phi(q)/2,

then, for covering S and 15<=q<=28 with no zero owner, q is a witness iff

    C_q > N_q - h_q.

In particular N_q<h_q is a residue-blind sufficient condition.  For the
prime moduli 17,19,23 and a 13-set, the exact collision thresholds are
respectively 5,4,2; S0 attains equality at all three, showing sharpness.

Tournament Analysis is deliberately telemetry.  Its vertices are runners,
the pairwise observable is lexicographic role dominance, and the switch is
deck-first versus covering-first.  The fixed tie Hamiltonian path is the
increasing speed order.  The exact carrier remains the bipartite incidence
hypergraph between runners and signed-card/covering obligations; the common
translation orbit proves that any runner tournament forgets essential metric
geometry.

All witness, interval, maximum, gcd, and incidence checks are exact.
Provenance: codex-S58.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, lcm

from lrc14_certificates import M_exact, good_intervals, is_covering


S0 = (43, 55, 61, 70, 73, 79, 83, 99, 103, 104, 109, 113, 156)
DEEP_WELL = tuple(range(1, 13)) + (182,)
SMALL_COVER_MODULI = tuple(range(2, 15))
TARGET_MODULI = tuple(range(15, 26))
PAIR_DECK_MODULI = tuple(range(15, 29))
TRANSLATION_PERIOD = lcm(*range(2, 26))
PI_UPPER = F(22, 7)


def dist_num(speed: int, a: int, q: int) -> int:
    residue = speed * a % q
    return min(residue, q - residue)


def is_q_witness(speeds: tuple[int, ...], a: int, q: int) -> bool:
    return all(14 * dist_num(speed, a, q) >= q for speed in speeds)


def clearance(speeds: tuple[int, ...], t: F) -> F:
    return min(min((speed * t) % 1, 1 - (speed * t) % 1) for speed in speeds)


def first_witness(speeds: tuple[int, ...], qmax: int = 100):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if is_q_witness(speeds, a, q):
                return q, a, F(min(dist_num(speed, a, q) for speed in speeds), q)
    return None


def unit_sign_pairs(q: int) -> tuple[int, ...]:
    return tuple(a for a in range(1, q // 2 + 1) if gcd(a, q) == 1)


def signed_card(speed: int, q: int) -> int:
    residue = speed % q
    return min(residue, q - residue)


def deck_data(speeds: tuple[int, ...], q: int):
    units = unit_sign_pairs(q)
    unit_speeds = tuple(speed for speed in speeds if gcd(speed, q) == 1)
    owners = {
        card: tuple(speed for speed in unit_speeds if signed_card(speed, q) == card)
        for card in units
    }
    blocked = tuple(card for card in units if owners[card])
    zeros = tuple(speed for speed in speeds if speed % q == 0)
    good_as = tuple(a for a in range(1, q) if is_q_witness(speeds, a, q))
    return {
        "units": units,
        "unit_speeds": unit_speeds,
        "owners": owners,
        "blocked": blocked,
        "zeros": zeros,
        "good_as": good_as,
    }


def prime_factors(number: int) -> tuple[int, ...]:
    factors = []
    candidate = 2
    while candidate * candidate <= number:
        if number % candidate == 0:
            factors.append(candidate)
            while number % candidate == 0:
                number //= candidate
        candidate += 1
    if number > 1:
        factors.append(number)
    return tuple(factors)


def primes_up_to(bound: int) -> tuple[int, ...]:
    return tuple(
        candidate
        for candidate in range(2, bound + 1)
        if all(candidate % divisor for divisor in range(2, int(candidate**0.5) + 1))
    )


def prime_packet_counts(speeds: tuple[int, ...]):
    return tuple(
        (prime, sum(speed % prime == 0 for speed in speeds))
        for prime in primes_up_to(max(speeds))
        if sum(speed % prime == 0 for speed in speeds) >= 2
    )


def longest_arithmetic_progression(speeds: tuple[int, ...]) -> tuple[int, ...]:
    speed_set = set(speeds)
    best: tuple[int, ...] = ()
    for left, right in combinations(speeds, 2):
        difference = right - left
        progression = []
        value = left
        while value in speed_set:
            progression.append(value)
            value += difference
        if len(progression) > len(best):
            best = tuple(progression)
    return best


def private_cards(speeds: tuple[int, ...]):
    data = {q: deck_data(speeds, q) for q in TARGET_MODULI}
    return {
        speed: tuple(
            (q, card)
            for q in TARGET_MODULI
            for card, owners in data[q]["owners"].items()
            if owners == (speed,)
        )
        for speed in speeds
    }


def band_residual_rows(speeds: tuple[int, ...]):
    rows = []
    for removed in speeds:
        core = tuple(speed for speed in speeds if speed != removed)
        intervals = good_intervals(core)
        measure = sum((right - left for left, right in intervals), F(0))
        components = len(intervals)
        margin = F(components) - PI_UPPER * removed * measure
        assert measure > 0 and margin > 0
        rows.append((removed, components, measure, removed * measure, margin))
    return tuple(rows)


def translated(k: int) -> tuple[int, ...]:
    assert k >= 0
    return tuple(speed + k * TRANSLATION_PERIOD for speed in S0)


def verify_uniform_packet_three() -> tuple[int, Counter[int]]:
    """Certify that no prime divides four members of any translated S_k.

    If p divided four translated speeds, it would divide their base-speed
    differences.  Every nontrivial four-subset difference gcd below has only
    factors 2,3,5 (hence factors of L), while the common base residue is a
    unit for those factors.  Translation by kL cannot turn it into zero.
    """

    nontrivial = 0
    gcd_hist: Counter[int] = Counter()
    for quad in combinations(S0, 4):
        difference_gcd = reduce(gcd, (abs(speed - quad[0]) for speed in quad[1:]))
        if difference_gcd == 1:
            continue
        nontrivial += 1
        gcd_hist[difference_gcd] += 1
        factors = prime_factors(difference_gcd)
        assert set(factors) <= {2, 3, 5}
        assert all(TRANSLATION_PERIOD % prime == 0 for prime in factors)
        assert all(quad[0] % prime for prime in factors)
    assert sum(speed % 2 == 0 for speed in S0) == 3
    return nontrivial, gcd_hist


def make_tournament(vertices, keys):
    edge = set()
    for left, right in combinations(vertices, 2):
        if keys[left] > keys[right] or (keys[left] == keys[right] and left < right):
            edge.add((left, right))
        else:
            edge.add((right, left))
    return edge


def tournament_fingerprint(vertices, edge):
    scores = {
        vertex: sum((vertex, other) in edge for other in vertices if other != vertex)
        for vertex in vertices
    }
    score_hist = dict(sorted(Counter(scores.values()).items()))
    cycles = sum(
        ((a, b) in edge and (b, c) in edge and (c, a) in edge)
        or ((b, a) in edge and (c, b) in edge and (a, c) in edge)
        for a, b, c in combinations(vertices, 3)
    )

    reach = {(u, v): u == v or (u, v) in edge for u in vertices for v in vertices}
    for middle in vertices:
        for left in vertices:
            for right in vertices:
                reach[left, right] = reach[left, right] or (
                    reach[left, middle] and reach[middle, right]
                )
    unseen = set(vertices)
    sccs = []
    while unseen:
        vertex = min(unseen)
        component = tuple(
            sorted(other for other in unseen if reach[vertex, other] and reach[other, vertex])
        )
        sccs.append(component)
        unseen.difference_update(component)

    index = {vertex: i for i, vertex in enumerate(vertices)}
    dp = [[0] * len(vertices) for _ in range(1 << len(vertices))]
    for vertex in vertices:
        dp[1 << index[vertex]][index[vertex]] = 1
    for mask in range(1 << len(vertices)):
        for last in vertices:
            count = dp[mask][index[last]]
            if not count:
                continue
            for nxt in vertices:
                bit = 1 << index[nxt]
                if not mask & bit and (last, nxt) in edge:
                    dp[mask | bit][index[nxt]] += count
    hamiltonian_paths = sum(dp[-1])
    path = tuple(sorted(vertices, key=lambda vertex: (-scores[vertex], vertex)))
    assert all((path[i], path[i + 1]) in edge for i in range(len(path) - 1))
    return score_hist, cycles, tuple(sccs), hamiltonian_paths, path


def tournament_report(speeds: tuple[int, ...], private):
    deck_first = {
        speed: (
            len(private[speed]),
            sum(gcd(speed, q) == 1 for q in TARGET_MODULI),
            sum(speed % divisor == 0 for divisor in SMALL_COVER_MODULI),
        )
        for speed in speeds
    }
    covering_first = {
        speed: (
            sum(speed % divisor == 0 for divisor in SMALL_COVER_MODULI),
            len(private[speed]),
            sum(gcd(speed, q) == 1 for q in TARGET_MODULI),
        )
        for speed in speeds
    }
    deck_edge = make_tournament(speeds, deck_first)
    covering_edge = make_tournament(speeds, covering_first)
    deck_fp = tournament_fingerprint(speeds, deck_edge)
    covering_fp = tournament_fingerprint(speeds, covering_edge)
    flips = sum(
        ((left, right) in deck_edge) != ((left, right) in covering_edge)
        for left, right in combinations(speeds, 2)
    )

    print("TOURNAMENT ANALYSIS (lossy telemetry)")
    print("  vertices=runners; exact carrier=runner x (signed-card or covering-obligation) incidence hypergraph")
    print("  alternatives_considered=moduli,signed_cards,zero_events,gaps,fixed_sections,section_boundaries,wall_events,residues,cover_arcs,Fourier_modes,matroid_circuits,proof_obligations")
    print("  pairwise_observable=lexicographic role dominance")
    print("  switch/gauges=deck-private-first -> small-divisor-covering-first")
    print("  tie_Hamiltonian_path=" + "->".join(map(str, speeds)))
    print("  preserved=ranked private-card/unit-incidence/covering counts")
    print("  destroyed=which card is owned, cross-q compatibility, metric clearance, and higher-denominator behavior")
    labels = ("score_hist", "directed_3cycles", "SCCs", "Hamiltonian_paths", "Hamiltonian_path")
    print("  deck-first " + " ".join(f"{key}={value}" for key, value in zip(labels, deck_fp)))
    print("  covering-first " + " ".join(f"{key}={value}" for key, value in zip(labels, covering_fp)))
    print(f"  edge_flips={flips}/{len(speeds) * (len(speeds) - 1) // 2}")
    print("  verdict=the bipartite incidence code carries q<=25; its translation orbit proves the tournament cannot carry metric truth")


def main() -> None:
    assert TRANSLATION_PERIOD == 26_771_144_400
    assert len(S0) == len(set(S0)) == 13
    assert reduce(gcd, S0) == 1
    assert gcd(70 - 43, 83 - 43) == 1
    assert is_covering(S0)
    assert min(S0) > 14 and max(S0) - min(S0) == 113

    print("LRC14 q<=25 ZERO-FREE SATURATED-DECK AUDIT -- exact certificate")
    print("source=codex-S58")
    print(f"S0={S0}")
    print("primitive=True covering=True all_speeds_above_14=True")
    print(f"diameter={max(S0)-min(S0)} ratio={F(max(S0), min(S0))} < 4")
    print("covering_owners=" + str(tuple(
        (divisor, tuple(speed for speed in S0 if speed % divisor == 0))
        for divisor in SMALL_COVER_MODULI
    )))

    deck = {q: deck_data(S0, q) for q in PAIR_DECK_MODULI}
    print("q15_25_zero_free_complete_decks:")
    for q in TARGET_MODULI:
        row = deck[q]
        assert not row["zeros"]
        assert row["blocked"] == row["units"]
        assert not row["good_as"]
        unit_count = len(row["unit_speeds"])
        collision_excess = unit_count - len(row["blocked"])
        saturation_threshold = unit_count - len(row["units"])
        assert collision_excess == saturation_threshold
        multiplicities = tuple(sorted((len(row["owners"][card]) for card in row["units"]), reverse=True))
        print(
            f"  q={q} cards={row['units']} zeros=() unit_count={unit_count} "
            f"owner_multiplicities={multiplicities} collision_excess={collision_excess} "
            f"full_deck_threshold={saturation_threshold} good_multipliers=()"
        )

    assert tuple(
        len(deck[prime]["unit_speeds"]) - len(deck[prime]["blocked"])
        for prime in (17, 19, 23)
    ) == (5, 4, 2)
    print("prime_collision_thresholds_sharp: q17=5 q19=4 q23=2 (witness iff strict excess)")
    print("residue_blind_positive_repair: zero-free q with unit_count < phi(q)/2 implies a q-witness")

    private = private_cards(S0)
    assert all(private[speed] for speed in S0)
    print("private_signed_cards (every runner is deletion-critical):")
    deletion_witnesses = []
    for speed in S0:
        print(f"  speed={speed} private={private[speed]}")
        removed = tuple(other for other in S0 if other != speed)
        assert any(
            any(is_q_witness(removed, a, q) for a in range(1, q))
            for q, _card in private[speed]
        )
        q, card = private[speed][0]
        multiplier = pow(card, -1, q)
        assert is_q_witness(removed, multiplier, q)
        assert dist_num(speed, multiplier, q) == 1
        remaining_clearance = F(
            min(dist_num(other, multiplier, q) for other in removed), q
        )
        assert remaining_clearance > F(1, 13)
        deletion_witnesses.append((speed, q, multiplier, remaining_clearance, card))
    print("deletion_witness_circuit=(removed_speed,q,a,12_speed_clearance,reattached_private_card):")
    print("  " + str(tuple(deletion_witnesses)))
    print("  every deletion is strictly above the 1/13 twelve-speed threshold, while reattachment blocks it at signed distance 1")

    witness = first_witness(S0)
    assert witness == (27, 2, F(2, 27))
    assert deck[26]["zeros"] == (104, 156)
    assert deck[27]["good_as"] == (2, 8, 19, 25)
    exact_m = M_exact(S0)
    assert exact_m == F(43, 199)
    assert clearance(S0, F(1, 199)) == exact_m
    print(f"first_rational_witness={witness}")
    print(f"M_exact={exact_m} at t=1/199 and 198/199")

    packets = prime_packet_counts(S0)
    assert max(count for _prime, count in packets) == 3
    leave_gcds = tuple(
        reduce(gcd, (other for other in S0 if other != speed)) for speed in S0
    )
    assert leave_gcds == (1,) * 13
    progression = longest_arithmetic_progression(S0)
    assert len(progression) == 3
    print(f"prime_packet_counts_ge2={packets}")
    print("max_common_prime_packet=3 (sharp under covering + q15_25 zero-free)")
    print("sharpness_reason=owners of 8,12,10 are distinct: lcm(8,12)=24, lcm(8,10)=40 is 20-divisible, lcm(12,10)=60 is 15-divisible")
    assert lcm(8, 12) == 24 and lcm(8, 10) % 20 == 0 and lcm(12, 10) % 15 == 0
    print(f"leave_one_out_gcds={leave_gcds}")
    print(f"longest_arithmetic_progression={progression} length=3")

    band_rows = band_residual_rows(S0)
    components = tuple(row[1] for row in band_rows)
    weakest = min(band_rows, key=lambda row: row[4])
    print("S312_exact_band_residual=True")
    print(f"leave_component_vector={components}")
    print(
        f"weakest_22_over_7_certificate=removed:{weakest[0]} "
        f"r-(22/7)w|G|:{weakest[4]} > 0"
    )

    nontrivial_quads, gcd_hist = verify_uniform_packet_three()
    print("translation_orbit:")
    print(f"  L=lcm(2,...,25)={TRANSLATION_PERIOD}")
    print("  S_k=S0+kL preserves every residue/covering/deck obligation through q=25")
    print("  primitive_for_all_k=True because gcd(70-43,83-43)=gcd(27,40)=1")
    print("  diameter_for_all_k=113")
    print("  ratio=(156+kL)/(43+kL) -> 1")
    print("  explicit_clearance=(43+kL)/(199+2kL) at t=1/(199+2kL) -> 1/2")
    print(
        f"  no_four_common_prime_packet_certificate=nontrivial_difference_gcd_quads:"
        f"{nontrivial_quads}/715 gcd_hist:{dict(sorted(gcd_hist.items()))}"
    )
    print("  every_nontrivial_quad_gcd_has_only_2,3,5_factors_and_nonzero_common_residue")
    print("  max_common_prime_packet_for_all_k=3")
    for k in (0, 1, 10):
        speeds = translated(k)
        assert reduce(gcd, speeds) == 1
        assert is_covering(speeds)
        assert max(speeds) - min(speeds) == 113
        assert all(not deck_data(speeds, q)["zeros"] for q in TARGET_MODULI)
        assert all(deck_data(speeds, q)["blocked"] == unit_sign_pairs(q) for q in TARGET_MODULI)
        assert all(not deck_data(speeds, q)["good_as"] for q in TARGET_MODULI)
        denominator = speeds[0] + speeds[-1]
        expected_clearance = F(speeds[0], denominator)
        assert clearance(speeds, F(1, denominator)) == expected_clearance
        print(
            f"  sample_k={k} min={speeds[0]} max={speeds[-1]} "
            f"ratio={F(speeds[-1], speeds[0])} explicit_clearance={expected_clearance}"
        )

    tournament_report(S0, private)

    assert first_witness(DEEP_WELL) == (27, 2, F(2, 27))
    assert M_exact(DEEP_WELL) == F(14, 183)

    print("FINAL VERDICT")
    print("zero_owner_exclusion_repairs_q25=False")
    print("primitive_or_covering_or_ratio_or_bounded_diameter_or_small_prime_packet_repairs_q25=False")
    print("coherent_sheet_exclusion_repairs_q25=False")
    print("near_tightness_note={1,...,12,182} already has M=14/183 and first witness 2/27")
    print("strongest_correct_replacement=zero-owner plus collision-surplus/signed-deck criterion, not a raw structural scalar")


if __name__ == "__main__":
    main()
