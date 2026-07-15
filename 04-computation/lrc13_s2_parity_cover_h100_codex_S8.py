#!/usr/bin/env python3
"""Finite-exact THM-776 s=2 parity-cover exclusion through height 100.

For

    A = 2U union {x,y},  |U|=10,  x,y odd,

Using THM-769, A is tight at 1/13 exactly when the loose set

    G_U = {tau : min_(u in U) ||u tau|| > 1/13}

is contained in the closed opposite-parity eligibility locus H_{x,y}.  For a
fixed odd pair, put B_{x,y}=(R/Z)\H_{x,y}.  A core U works exactly when its
closed 1/13-danger teeth cover B_{x,y}.

This program atomizes that statement over exact rational endpoints.  On each
open atom it records the candidate quotient speeds whose danger teeth contain
the atom.  The existence of a core is then a finite hitting-set CNF.  We test
all odd x<y<=99 against every quotient speed u<=50.

Exactness conventions:

* G_U is open and H_{x,y} is closed, so any failure G_U intersect B_{x,y}
  is open and contains an atom midpoint.  No isolated endpoint can be missed.
* atom endpoints include ||u tau||=1/13 for all u<=50 and
  ||w tau||=2/13 for every odd w<=99;
* all geometric predicates use fractions.Fraction;
* inclusion-minimal coverer masks are a lossless monotone-CNF compression;
* CaDiCaL proves UNSAT with at most 11 quotient speeds, while Glucose gives an
  independent UNSAT check at the needed bound 10;
* {1,...,12} covers the whole circle, providing the exact upper bound 12.

Tournament Analysis.  On a good atom the two odd runners form a two-vertex
owner tournament: x->y when x owns sheet 0, and y->x otherwise.  The gauge
tau->1-tau reverses the edge.  Every such tournament is transitive, with score
histogram (0,1), no directed cycles, two singleton SCCs, and one Hamiltonian
path.  The natural x<y order is the tie Hamiltonian path.  This is telemetry,
not the proof carrier: runner tournaments destroy simultaneous coverage.  The
predicate-preserving vertices are bad atoms/proof obligations, with quotient
speeds incident to the atoms they kill.

Requires python-sat with the bundled cadical195 and glucose4 solvers.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations
from math import gcd
from multiprocessing import get_context
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


DELTA = F(1, 13)
STATE: dict[str, object] = {}


def circle_distance(z: F) -> F:
    residue = z.numerator % z.denominator
    return F(min(residue, z.denominator - residue), z.denominator)


def nearest_integer(z: F) -> int:
    """Unique nearest integer on all eligibility atoms used below."""
    return (2 * z.numerator + z.denominator) // (2 * z.denominator)


def boundary_grid(height: int) -> tuple[list[F], list[F], list[F]]:
    umax = height // 2
    odds = range(1, height + 1, 2)
    points = {F(0), F(1)}

    # ||u tau||=1/13.
    for u in range(1, umax + 1):
        for k in range(u + 1):
            for sign in (-1, 1):
                z = F(13 * k + sign, 13 * u)
                if 0 <= z <= 1:
                    points.add(z)

    # ||w tau||=2/13.
    for w in odds:
        for k in range(w + 1):
            for sign in (-2, 2):
                z = F(13 * k + sign, 13 * w)
                if 0 <= z <= 1:
                    points.add(z)

    endpoints = sorted(points)
    widths = [b - a for a, b in zip(endpoints, endpoints[1:])]
    midpoints = [(a + b) / 2 for a, b in zip(endpoints, endpoints[1:])]
    return endpoints, widths, midpoints


def make_state(height: int) -> dict[str, object]:
    assert height == 100, "the stored THM-776 certificate is for height 100"
    umax = height // 2
    odds = tuple(range(1, height + 1, 2))
    endpoints, widths, midpoints = boundary_grid(height)

    coverers = []
    for tau in midpoints:
        mask = 0
        for u in range(1, umax + 1):
            if circle_distance(u * tau) <= DELTA:
                mask |= 1 << (u - 1)
        coverers.append(mask)

    cardinality_clauses = {
        bound: CardEnc.atmost(
            lits=list(range(1, umax + 1)),
            bound=bound,
            top_id=umax,
            encoding=EncType.seqcounter,
        ).clauses
        for bound in (10, 11)
    }

    return {
        "height": height,
        "umax": umax,
        "odds": odds,
        "endpoints": endpoints,
        "widths": widths,
        "midpoints": midpoints,
        "coverers": coverers,
        "cardinality_clauses": cardinality_clauses,
    }


def pair_record(pair: tuple[int, int]):
    """Build the exact obstruction hypergraph and solve its two CNFs."""
    x, y = pair
    umax = STATE["umax"]
    midpoints = STATE["midpoints"]
    widths = STATE["widths"]
    coverers = STATE["coverers"]
    cardinality_clauses = STATE["cardinality_clauses"]
    assert isinstance(umax, int)

    raw_constraints: set[int] = set()
    labels: list[int] = []
    # 0=x->y, 1=y->x, 2=same-sheet tie, 3=exactly one eligible, 4=neither.
    measures = [F(0) for _ in range(5)]

    for tau, width, coverer_mask in zip(midpoints, widths, coverers):
        eligible_x = circle_distance(x * tau) <= 2 * DELTA
        eligible_y = circle_distance(y * tau) <= 2 * DELTA
        if eligible_x and eligible_y:
            parity_x = nearest_integer(x * tau) & 1
            parity_y = nearest_integer(y * tau) & 1
            if parity_x != parity_y:
                label = parity_x
            else:
                label = 2
        elif eligible_x or eligible_y:
            label = 3
        else:
            label = 4

        labels.append(label)
        measures[label] += width
        if label >= 2:
            raw_constraints.add(coverer_mask)

    # If C1 is contained in C2, hitting C1 automatically hits C2.
    ordered = sorted(raw_constraints, key=lambda c: (c.bit_count(), c))
    minimal: list[int] = []
    for constraint in ordered:
        if not any((smaller & constraint) == smaller for smaller in minimal):
            minimal.append(constraint)

    clauses = [
        [u + 1 for u in range(umax) if (constraint >> u) & 1]
        for constraint in minimal
    ]
    assert all(clauses), "an empty clause would be an immediate geometric witness"

    answers: dict[int, bool] = {}
    statistics: dict[int, tuple[int, int, int, int]] = {}
    for bound, solver_name in ((11, "cadical195"), (10, "glucose4")):
        with Solver(
            name=solver_name,
            bootstrap_with=clauses + cardinality_clauses[bound],
        ) as solver:
            satisfiable = solver.solve()
            stats = solver.accum_stats()
        answers[bound] = satisfiable
        statistics[bound] = (
            stats.get("restarts", 0),
            stats.get("conflicts", 0),
            stats.get("decisions", 0),
            stats.get("propagations", 0),
        )

    # Near tau=0 the pair is a same-sheet tie, so linear component counting is
    # already cyclic component counting.  A direct good-to-good orientation
    # flip cannot occur because 2/13<1/2 leaves a gap between parity teeth.
    assert labels[0] >= 2 and labels[-1] >= 2
    assert all(
        not (labels[i] < 2 and labels[i - 1] < 2 and labels[i] != labels[i - 1])
        for i in range(1, len(labels))
    )
    component_orientations = []
    previous = 2
    for label in labels:
        if label < 2 and previous >= 2:
            component_orientations.append(label)
        previous = label

    digest = sha256()
    digest.update(f"{x},{y}\n".encode())
    for constraint in minimal:
        digest.update(f"{constraint:x}\n".encode())

    return {
        "pair": pair,
        "answers": answers,
        "statistics": statistics,
        "constraint_count": len(minimal),
        "constraint_cardinalities": tuple(
            sorted(Counter(c.bit_count() for c in minimal).items())
        ),
        "measures": tuple(measures),
        "component_orientations": tuple(Counter(component_orientations).items()),
        "digest": digest.hexdigest(),
    }


def exact_maximin(speeds: tuple[int, ...]) -> tuple[F, F]:
    """Exact M(S), using all pair crossings and self-cusps."""
    denominators = {2 * v for v in speeds}
    for a, b in combinations(speeds, 2):
        denominators.add(a + b)
        denominators.add(b - a)
    denominators.discard(0)

    best = F(0)
    witness = F(0)
    for q in denominators:
        for p in range(1, q):
            numerator = min(
                min((v * p) % q, q - (v * p) % q) for v in speeds
            )
            value = F(numerator, q)
            if value > best:
                best = value
                witness = F(p, q)
    return best, witness


def direct_height24_crosscheck():
    rows = 0
    tight = 0
    closest = None
    for core in combinations(range(1, 13), 10):
        for x, y in combinations(range(1, 24, 2), 2):
            speeds = tuple(sorted((*(2 * u for u in core), x, y)))
            rows += 1
            value, witness = exact_maximin(speeds)
            if value == DELTA:
                tight += 1
            record = (value, witness, speeds, core, x, y)
            if closest is None or record < closest:
                closest = record
    return rows, tight, closest


def divisor_complete_core_count(umax: int):
    """Quantify the THM-772-filtered rows avoided by the hypergraph method."""
    moduli = tuple(range(2, 13))
    full = (1 << len(moduli)) - 1
    dp = {(0, 0, 0): 1}
    first_max = None
    for u in range(1, umax + 1):
        divisor_mask = sum(
            1 << i for i, modulus in enumerate(moduli) if u % modulus == 0
        )
        next_dp = dict(dp)
        created_full = False
        for (size, mask, common_gcd), count in dp.items():
            if size == 10:
                continue
            key = (size + 1, mask | divisor_mask, gcd(common_gcd, u))
            next_dp[key] = next_dp.get(key, 0) + count
            if size + 1 == 10 and key[1] == full:
                created_full = True
        if created_full and first_max is None:
            first_max = u
        dp = next_dp

    all_rows = sum(
        count for (size, mask, _), count in dp.items()
        if size == 10 and mask == full
    )
    primitive_rows = sum(
        count for (size, mask, common_gcd), count in dp.items()
        if size == 10 and mask == full and common_gcd == 1
    )
    return all_rows, primitive_rows, first_max


def fmt_fraction(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=100)
    parser.add_argument(
        "--workers", type=int, default=min(8, os.cpu_count() or 2)
    )
    parser.add_argument(
        "--quiet-progress", action="store_true",
        help="suppress progress messages on stderr",
    )
    args = parser.parse_args()

    global STATE
    STATE = make_state(args.height)
    height = STATE["height"]
    umax = STATE["umax"]
    odds = STATE["odds"]
    endpoints = STATE["endpoints"]
    assert isinstance(height, int) and isinstance(umax, int)

    pairs = [
        (x, y) for i, x in enumerate(odds) for y in odds[i + 1:]
    ]
    if args.workers == 1:
        iterator = map(pair_record, pairs)
        pool = None
    else:
        pool = get_context("fork").Pool(args.workers)
        iterator = pool.imap_unordered(pair_record, pairs, chunksize=1)

    records = []
    try:
        for index, record in enumerate(iterator, 1):
            records.append(record)
            if not args.quiet_progress and (index % 100 == 0 or index == len(pairs)):
                print(f"completed {index}/{len(pairs)} pairs", file=sys.stderr, flush=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    records.sort(key=lambda r: r["pair"])

    assert all(not r["answers"][11] for r in records)
    assert all(not r["answers"][10] for r in records)
    universal_cover = (1 << 12) - 1
    # Reconstruct every minimal mask from the per-pair digest is impossible,
    # so the upper bound is also checked geometrically on every atom here.
    assert all(mask & universal_cover for mask in STATE["coverers"])

    grid_digest = sha256(
        "\n".join(fmt_fraction(z) for z in endpoints).encode()
    ).hexdigest()
    atlas_digest = sha256()
    for record in records:
        x, y = record["pair"]
        atlas_digest.update(f"{x},{y}:{record['digest']}\n".encode())

    constraint_count_hist = Counter(r["constraint_count"] for r in records)
    cardinality_total = Counter()
    for record in records:
        cardinality_total.update(dict(record["constraint_cardinalities"]))

    solver_totals = {}
    for bound in (11, 10):
        solver_totals[bound] = tuple(
            sum(record["statistics"][bound][i] for record in records)
            for i in range(4)
        )

    measure_totals = tuple(
        sum((record["measures"][i] for record in records), F(0))
        for i in range(5)
    )
    assert measure_totals[0] == measure_totals[1]
    assert sum(measure_totals) == len(pairs)

    component_hist = Counter(
        sum(dict(record["component_orientations"]).values())
        for record in records
    )
    orientation_components = tuple(
        sum(dict(record["component_orientations"]).get(i, 0) for record in records)
        for i in (0, 1)
    )
    assert orientation_components[0] == orientation_components[1]
    zero_good_pairs = tuple(
        record["pair"] for record in records
        if not dict(record["component_orientations"])
    )

    direct_rows, direct_tight, closest = direct_height24_crosscheck()
    assert direct_rows == 4356 and direct_tight == 0
    assert closest[0] == F(2, 23)

    divisor_all, divisor_primitive, first_divisor_max = divisor_complete_core_count(umax)
    assert (divisor_all, divisor_primitive, first_divisor_max) == (
        885_427_750, 884_640_190, 12
    )

    script_digest = sha256(Path(__file__).read_bytes()).hexdigest()

    print("THM-776 finite-exact s=2 parity-cover certificate")
    print("strict core-safe / closed exception-danger endpoint convention")
    print(f"script_sha256={script_digest}")
    print(f"height={height} quotient_umax={umax} odd_pair_max={height-1}")
    print(
        f"endpoints={len(endpoints)} open_atoms={len(STATE['midpoints'])} "
        f"odd_pairs={len(pairs)}"
    )
    print(f"grid_sha256={grid_digest}")
    print(f"constraint_atlas_sha256={atlas_digest.hexdigest()}")
    print()

    print("HITTING-SET VERDICT")
    print("sat_at_most_11=0/1225 solver=cadical195")
    print("sat_at_most_10=0/1225 solver=glucose4")
    print("universal_upper_cover=(1,2,...,12)")
    print("transversal_number_hist=((12,1225),)")
    print(f"minimal_constraint_count_hist={tuple(sorted(constraint_count_hist.items()))}")
    print(f"minimal_constraint_cardinality_total={tuple(sorted(cardinality_total.items()))}")
    print(
        "solver_stats_at_most_11(restarts,conflicts,decisions,propagations)="
        f"{solver_totals[11]}"
    )
    print(
        "solver_stats_at_most_10(restarts,conflicts,decisions,propagations)="
        f"{solver_totals[10]}"
    )
    print()

    print("THM-772 FILTER QUANTIFICATION")
    print(f"divisor_complete_10cores_upto50={divisor_all}")
    print(f"primitive_divisor_complete_10cores_upto50={divisor_primitive}")
    print(f"first_possible_core_max={first_divisor_max}")
    print(f"primitive_filtered_core_pair_packets={divisor_primitive*len(pairs)}")
    print("all are subsumed by the ambient per-pair transversal exclusion")
    print()

    print("DIRECT EXACT-M CROSSCHECK")
    print(f"height24_rows={direct_rows} tight_rows={direct_tight}")
    print(
        f"closest_M={fmt_fraction(closest[0])} witness_t={fmt_fraction(closest[1])} "
        f"speeds={closest[2]}"
    )
    print()

    print("TOURNAMENT / INCIDENCE AUDIT")
    print("pairwise_observable=which odd runner owns sheet 0 on a good atom")
    print("gauge=tau->1-tau (reverses every owner edge)")
    print("tie_hamiltonian_path=natural odd-runner order x<y")
    print("each_good_tournament: score_hist=(0,1) c3=0 scc=(1,1) hp=1")
    print(f"good_component_orientation_totals={orientation_components}")
    print(f"good_component_hist={tuple(sorted(component_hist.items()))}")
    print(f"zero_good_pairs={zero_good_pairs}")
    label_names = ("x_to_y", "y_to_x", "same_sheet", "exactly_one", "neither")
    for name, total in zip(label_names, measure_totals):
        average = total / len(pairs)
        print(
            f"average_measure_{name}={fmt_fraction(average)} "
            f"decimal={float(average):.12f}"
        )
    print("runner_tournament_loses=simultaneous cover and transversal number")
    print("predicate_preserving_carrier=bad-atom/proof-obligation incidence hypergraph")
    print()
    print("FINAL: no ten-even/two-odd tight 12-set has max(A)<=100")
    print("SCOPE: FINITE-EXACT, not a uniform exclusion")


if __name__ == "__main__":
    main()
