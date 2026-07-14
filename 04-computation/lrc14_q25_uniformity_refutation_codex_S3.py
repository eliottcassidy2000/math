#!/usr/bin/env python3
"""Exact refutation of the proposed uniform q<=25 LRC(14) witness bound.

The S312 observation was:

    every sampled covering 13-set in the uncapped, diameter-<=339 residual
    has a rational lonely witness a/q with 15 <= q <= 25.

This script gives two exact counterexamples.

* V26 = 26*{1,...,12} union {339} is the transparent coherent example.
  In fact, c*{1,...,12} blocks every denominator 15..25 for every c.
* SSTAR is a stronger gcd-incoherent example.  No prime divides seven of its
  speeds, so it has no common-factor pack large enough for the d<=6 form of
  THM-737.  Nevertheless it blocks every q in 15..25.

For both families the script verifies, using the repository's exact Fraction
interval engine, every predicate in S312's ``is_band_residual``:

* divisor covering for 2..14;
* at least four speeds above 14;
* diameter at most 339;
* for every leave-one-out core P=S\{w}, w <= r_P/(pi*|G'_P|).

The last item is certified without floating point: pi < 22/7 and the exact
inequality (22/7)*w*|G'_P| < r_P is checked for every w.

The valid residue theorem left behind by the refutation is also checked.  For
a covering family and 15<=q<=25, a q-witness exists exactly when q divides no
speed and the signed unit-pair deck modulo q misses a pair.  The two examples
have a complete blocker at every q, but explicit witnesses immediately above
25 (q=27 for V26 and q=26 for SSTAR).

Tournament Analysis uses moduli, not runners, as vertices.  Its exact blocker
deck is the mandatory sidecar; the tournament is only telemetry.  The pair-
first gauge ranks moduli by retained signed-pair structure, while the
compression-first gauge ranks them by the cost of a shortest blocker
certificate.  The fixed tie Hamiltonian path is 15->16->...->25.  The script
reports score histograms, directed 3-cycles, SCCs, edge flips, and exact
Hamiltonian-path counts.

All mathematical comparisons are integral or Fraction-exact.  The only use of
pi is the rigorous classical upper bound pi < 22/7.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd

from lrc14_certificates import M_exact, good_intervals, is_covering


V26 = tuple([26 * i for i in range(1, 13)] + [339])
SSTAR = (81, 91, 131, 151, 157, 196, 258, 274, 313, 328, 330, 339, 348)
MODULI = tuple(range(15, 26))
PI_UPPER = F(22, 7)


def dist_num(s: int, a: int, q: int) -> int:
    """Numerator of ||s*a/q|| in [0,q/2]."""
    r = (s * a) % q
    return min(r, q - r)


def clearance_at_rational(speeds: tuple[int, ...], a: int, q: int) -> F:
    return F(min(dist_num(s, a, q) for s in speeds), q)


def is_q_witness(speeds: tuple[int, ...], a: int, q: int) -> bool:
    """Closed LRC(14) test, written without division."""
    return all(14 * dist_num(s, a, q) >= q for s in speeds)


def first_witness(speeds: tuple[int, ...], qmax: int = 500):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if is_q_witness(speeds, a, q):
                return q, a, clearance_at_rational(speeds, a, q)
    return None


def unit_sign_pairs(q: int) -> tuple[int, ...]:
    """Least positive representatives for (Z/qZ)^*/{+1,-1}."""
    return tuple(b for b in range(1, q // 2 + 1) if gcd(b, q) == 1)


def blocked_unit_pairs(speeds: tuple[int, ...], q: int) -> tuple[int, ...]:
    pairs = {
        min(s % q, q - (s % q))
        for s in speeds
        if gcd(s, q) == 1
    }
    return tuple(sorted(pairs))


def zero_owners(speeds: tuple[int, ...], q: int) -> tuple[int, ...]:
    return tuple(s for s in speeds if s % q == 0)


def blocker_data(speeds: tuple[int, ...]):
    out = {}
    for q in MODULI:
        units = unit_sign_pairs(q)
        blocked = blocked_unit_pairs(speeds, q)
        zeros = zero_owners(speeds, q)
        good_as = tuple(a for a in range(1, q) if is_q_witness(speeds, a, q))
        # For a covering family, nonunit a reduces to q/gcd(a,q)<=12 and is
        # killed by the small-denominator divisibility covering.  A unit a is
        # killed iff q has a zero owner or its inverse sign-pair is blocked.
        residue_blocked = bool(zeros) or blocked == units
        assert residue_blocked == (not good_as)
        out[q] = {
            "units": units,
            "blocked": blocked,
            "zeros": zeros,
            "good_as": good_as,
        }
    return out


def leave_one_out_stats(speeds: tuple[int, ...]):
    rows = []
    for w in speeds:
        core = tuple(s for s in speeds if s != w)
        intervals = good_intervals(core)
        measure = sum((b - a for a, b in intervals), F(0))
        components = len(intervals)
        weighted_measure = w * measure
        margin_22_7 = F(components) - PI_UPPER * weighted_measure
        assert measure > 0
        assert margin_22_7 > 0
        rows.append((w, components, measure, weighted_measure, margin_22_7))
    return rows


def primes_up_to(n: int) -> tuple[int, ...]:
    primes = []
    for candidate in range(2, n + 1):
        if all(candidate % p for p in primes if p * p <= candidate):
            primes.append(candidate)
    return tuple(primes)


def pack_stats(speeds: tuple[int, ...]):
    counts = tuple(
        (p, sum(s % p == 0 for s in speeds))
        for p in primes_up_to(max(speeds))
        if sum(s % p == 0 for s in speeds) >= 2
    )
    max_prime_packet = max((count for _, count in counts), default=1)
    leave_gcds = tuple(
        reduce(gcd, (s for s in speeds if s != w))
        for w in speeds
    )
    return counts, max_prime_packet, leave_gcds


def covering_owners(speeds: tuple[int, ...]):
    return tuple((q, next(s for s in speeds if s % q == 0)) for q in range(2, 15))


def make_tournament(vertices, keys, tie_path):
    """Orient by a lexicographic gauge; exact tie orientation follows tie_path."""
    tie_rank = {v: i for i, v in enumerate(tie_path)}
    edge = set()
    for u, v in combinations(vertices, 2):
        ku, kv = keys[u], keys[v]
        if ku > kv or (ku == kv and tie_rank[u] < tie_rank[v]):
            edge.add((u, v))
        else:
            edge.add((v, u))
    return edge


def tournament_fingerprint(vertices, edge):
    scores = {v: sum((v, w) in edge for w in vertices if w != v) for v in vertices}
    score_hist = dict(sorted(Counter(scores.values()).items()))
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        if (
            ((a, b) in edge and (b, c) in edge and (c, a) in edge)
            or ((b, a) in edge and (c, b) in edge and (a, c) in edge)
        ):
            cycles += 1

    # Reachability gives SCCs without an external graph dependency.
    reach = {(u, v): (u == v or (u, v) in edge) for u in vertices for v in vertices}
    for k in vertices:
        for i in vertices:
            for j in vertices:
                reach[i, j] = reach[i, j] or (reach[i, k] and reach[k, j])
    unseen = set(vertices)
    sccs = []
    while unseen:
        v = min(unseen)
        component = tuple(sorted(w for w in unseen if reach[v, w] and reach[w, v]))
        sccs.append(component)
        unseen.difference_update(component)

    # Exact count of directed Hamiltonian paths by subset DP.
    index = {v: i for i, v in enumerate(vertices)}
    size = 1 << len(vertices)
    dp = [[0] * len(vertices) for _ in range(size)]
    for v in vertices:
        dp[1 << index[v]][index[v]] = 1
    for mask in range(size):
        for last in vertices:
            count = dp[mask][index[last]]
            if not count:
                continue
            for nxt in vertices:
                bit = 1 << index[nxt]
                if not (mask & bit) and (last, nxt) in edge:
                    dp[mask | bit][index[nxt]] += count
    hamiltonian_count = sum(dp[-1])

    # In these lexicographic gauges the tournament is transitive, so score
    # order is the unique Hamiltonian path; still verify rather than assume.
    path = tuple(sorted(vertices, key=lambda v: (-scores[v], v)))
    assert all((path[i], path[i + 1]) in edge for i in range(len(path) - 1))
    return score_hist, cycles, tuple(sccs), hamiltonian_count, path


def tournament_report(name: str, data):
    vertices = MODULI
    tie_path = MODULI
    pair_keys = {
        q: (len(data[q]["blocked"]), int(not data[q]["zeros"]))
        for q in vertices
    }
    compression_keys = {}
    for q in vertices:
        cost = 1 if data[q]["zeros"] else len(data[q]["units"])
        compression_keys[q] = (-cost, len(data[q]["zeros"]))

    pair_edge = make_tournament(vertices, pair_keys, tie_path)
    compression_edge = make_tournament(vertices, compression_keys, tie_path)
    pair_fp = tournament_fingerprint(vertices, pair_edge)
    compression_fp = tournament_fingerprint(vertices, compression_edge)
    flips = sum(
        ((u, v) in pair_edge) != ((u, v) in compression_edge)
        for u, v in combinations(vertices, 2)
    )

    print(f"TOURNAMENT {name}")
    print("  vertices=moduli q=15..25 (runner vertices challenged and rejected)")
    print("  alternatives_considered=runners,gaps,fixed_sections,section_boundaries,wall_events,residues,cover_arcs,Fourier_modes,matroid_circuits,proof_obligations")
    print("  pairwise_observable=lex advantage in (signed-pairs retained, pair-only blocker) vs (negative shortest blocker cost, zero-owner multiplicity)")
    print("  switch/gauges=pair-first -> compression-first")
    print("  tie_Hamiltonian_path=" + "->".join(map(str, tie_path)))
    print("  preserved=exact q-blocked predicate plus zero/pair-deck sidecar")
    print("  destroyed=metric clearance M, endpoint order, scale, and cross-q compatibility")
    labels = ("score_hist", "directed_3cycles", "SCCs", "Hamiltonian_paths", "Hamiltonian_path")
    print("  pair-first " + " ".join(f"{k}={v}" for k, v in zip(labels, pair_fp)))
    print("  compression-first " + " ".join(f"{k}={v}" for k, v in zip(labels, compression_fp)))
    print(f"  edge_flips={flips}/{len(vertices) * (len(vertices) - 1) // 2}")
    print("  verdict=tournament telemetry is transitive here; the blocker deck, not the orientation, carries the proof")


def verify_dilated_block_lemma(c: int):
    """Finite replay of the elementary all-c proof on the requested c=26 row."""
    block = tuple(c * i for i in range(1, 13))
    for q in MODULI:
        for a in range(1, q):
            assert not is_q_witness(block, a, q)


def family_report(name: str, speeds: tuple[int, ...], expected_first_witness):
    print("=" * 88)
    print(name)
    print("S=" + str(speeds))
    assert len(speeds) == 13 and len(set(speeds)) == 13
    assert is_covering(speeds)
    assert sum(s > 14 for s in speeds) >= 4
    assert max(speeds) - min(speeds) <= 339
    print(f"primitive_gcd={reduce(gcd, speeds)} covering=True f_above_14={sum(s > 14 for s in speeds)} diameter={max(speeds)-min(speeds)}")
    print("covering_owners=" + str(covering_owners(speeds)))

    rows = leave_one_out_stats(speeds)
    r_vector = tuple(row[1] for row in rows)
    max_row = max(rows, key=lambda row: row[3])
    min_r = min(r_vector)
    assert PI_UPPER * max_row[3] < min_r
    print("S312_exact_residual=True")
    print("leave_component_vector=" + str(r_vector))
    print(f"max_w_times_G=w:{max_row[0]} value:{max_row[3]}")
    print(f"uniform_no_cap_chain: pi*w*G < (22/7)*{max_row[3]} = {PI_UPPER*max_row[3]} < min_r={min_r} <= r_w")

    data = blocker_data(speeds)
    print("q15_25_blocker_decks:")
    for q in MODULI:
        row = data[q]
        print(
            f"  q={q} units={row['units']} blocked={row['blocked']} "
            f"zeros={row['zeros']} good_multipliers={row['good_as']}"
        )
    assert all(not data[q]["good_as"] for q in MODULI)
    print("no_q15_25_witness=True")

    witness = first_witness(speeds)
    assert witness == expected_first_witness
    print(f"first_rational_witness={witness}")
    exact_m = M_exact(speeds)
    print(f"M_exact={exact_m}")

    prime_counts, max_prime_packet, leave_gcds = pack_stats(speeds)
    print(f"prime_packet_counts_ge2={prime_counts}")
    print(f"max_common_prime_packet={max_prime_packet}")
    print(f"leave_one_out_gcds={leave_gcds}")
    print(f"THM737_d_le_6_pack_excluded={max_prime_packet <= 6}")
    tournament_report(name, data)


def main():
    print("LRC14 q<=25 UNIFORMITY REFUTATION -- exact certificate")
    print("valid_pair_deck_criterion: for covering S and 15<=q<=25, a witness exists iff q has no zero owner and B_q(S) misses a unit sign-pair")
    print("nonunit_reason: reducing a/q gives q0<=12, and covering supplies a q0-divisible speed")
    print("unit_reason: q/14 is strictly between 1 and 2, so only residues 0,+1,-1 are unsafe")
    print("elementary_dilated_block_reason: if gcd(c,q)>1 choose i=q/gcd(c,q)<=12; otherwise choose the signed least representative of (a*c)^(-1), also <=12")

    verify_dilated_block_lemma(26)
    family_report("V26 coherent near-dilate", V26, (27, 2, F(2, 27)))
    assert clearance_at_rational(V26, 27, 338) == F(1, 13)
    print("V26_explicit_THM757_point=t=27/338 clearance=1/13")

    family_report("SSTAR gcd-incoherent residual", SSTAR, (26, 3, F(1, 13)))
    assert clearance_at_rational(SSTAR, 167, 470) == F(101, 470)
    print("SSTAR_exact_maximizer=t=167/470 clearance=101/470")

    print("=" * 88)
    print("FINAL VERDICT")
    print("uniform_q_le_25_bound=FALSE")
    print("loose_iff_q_le_25=FALSE")
    print("gcd_incoherent_rescue=FALSE")
    print("narrow_valid_statement=the pair-deck equivalence and per-family certificates only")
    print("research_redirect=classify simultaneous zero/pair blocker decks and allow adaptive denominators or another exact certificate route")


if __name__ == "__main__":
    main()
