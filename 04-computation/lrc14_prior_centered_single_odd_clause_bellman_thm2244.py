#!/usr/bin/env python3
"""Exact prior-centered single-odd-clause Bellman certificate for THM-2244.

This companion keeps every available even checkpoint in the positive carrier
and selects odd checkpoint 1 in the negative carrier.  For each labelled core
j it centers the selected atom at a possibly different earlier base time s_j:

    q = |O| - sum_(k in O,j)
        (X_(j,lambda_j-k) - rho^(lambda_j-k-s_j) X_(j,s_j)),

where rho=-1/13.  When all relative exponents are odd, every correction
coefficient is negative and q<=0 on the selected negative carrier.

The theorem uses the earliest prior time of the required parity.  The Bellman
relaxation eliminates one factor of thirteen at a time.  At each
step the three current danger bits may have arbitrary dependence subject only
to their exact conditional marginals (2-next_bit)/13.  The automaton retains
the even-clause mask, the number of selected raw odd atoms, and all three
labelled base bits.  ``--scout`` uses floating point vertices only;
``--exact-profile`` replays one named configuration, and ``--exact-ledger``
runs the complete Fraction/primal-dual certificate.

The wider scout modes and the 223-row low-first reference census are retained
as overlap and failure-boundary controls.  Relative to THM-2239's proved
194-row ledger, the theorem's novel content is the exclusion of 28 of the 29
first-depth-three profiles.
"""

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product


PRIME = 13
RHO = Fraction(-1, PRIME)
TARGET = Fraction(961, 6930)
STATES = tuple(range(8))
STATE_BITS = tuple(
    tuple((state >> core) & 1 for core in range(3))
    for state in STATES
)
COUPLING_COLUMNS = tuple(
    (Fraction(1),) + tuple(Fraction(bit) for bit in bits)
    for bits in STATE_BITS
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def solve_square(matrix, rhs):
    """Solve a square rational system, returning None when singular."""
    size = len(rhs)
    augmented = [list(matrix[row]) + [rhs[row]] for row in range(size)]
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if augmented[row][column]
            ),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = augmented[column][column]
        augmented[column] = [
            value / scale for value in augmented[column]
        ]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    augmented[row][entry]
                    - scale * augmented[column][entry]
                    for entry in range(size + 1)
                ]
    return tuple(augmented[row][-1] for row in range(size))


INVERTIBLE_BASES = tuple(
    basis
    for basis in combinations(STATES, 4)
    if solve_square(
        tuple(
            tuple(
                COUPLING_COLUMNS[basis[column]][row]
                for column in range(4)
            )
            for row in range(4)
        ),
        (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    )
    is not None
)


@lru_cache(maxsize=None)
def coupling_vertex_witnesses(marginals):
    """Map every exact vertex to all bases that generate it."""
    rhs = (Fraction(1),) + tuple(marginals)
    witnesses = {}
    for basis in INVERTIBLE_BASES:
        matrix = tuple(
            tuple(
                COUPLING_COLUMNS[basis[column]][row]
                for column in range(4)
            )
            for row in range(4)
        )
        basic = solve_square(matrix, rhs)
        require(basic is not None, "invertible coupling basis became singular")
        if any(value < 0 for value in basic):
            continue
        vertex = [Fraction(0)] * 8
        for column, state in enumerate(basis):
            vertex[state] = basic[column]
        vertex = tuple(vertex)
        witnesses.setdefault(vertex, []).append(basis)
    require(witnesses, "coupling polytope has no vertex")
    return tuple(
        (vertex, tuple(bases))
        for vertex, bases in sorted(witnesses.items())
    )


@lru_cache(maxsize=None)
def coupling_vertices(marginals):
    return tuple(
        vertex for vertex, _ in coupling_vertex_witnesses(tuple(marginals))
    )


TRANSITION_VERTICES = tuple(
    coupling_vertices(
        tuple(
            Fraction(2 - STATE_BITS[parent][core], PRIME)
            for core in range(3)
        )
    )
    for parent in STATES
)
TERMINAL_VERTICES = coupling_vertices((Fraction(1, 7),) * 3)
TRANSITION_VERTICES_FLOAT = tuple(
    tuple(tuple(float(entry) for entry in vertex) for vertex in vertices)
    for vertices in TRANSITION_VERTICES
)
TERMINAL_VERTICES_FLOAT = tuple(
    tuple(float(entry) for entry in vertex)
    for vertex in TERMINAL_VERTICES
)


def maximize_float(values, vertices):
    return max(
        sum(vertex[state] * values[state] for state in STATES)
        for vertex in vertices
    )


EXACT_LP_CHECKS = 0


@lru_cache(maxsize=None)
def maximize_exact(values, marginals):
    """Exact primal vertex maximum with an independently checked dual."""
    global EXACT_LP_CHECKS
    EXACT_LP_CHECKS += 1
    witnessed_vertices = coupling_vertex_witnesses(tuple(marginals))
    vertex_values = tuple(
        (
            sum(vertex[state] * values[state] for state in STATES),
            vertex,
            bases,
        )
        for vertex, bases in witnessed_vertices
    )
    primal = max(row[0] for row in vertex_values)
    candidate_bases = []
    for value, _, bases in vertex_values:
        if value == primal:
            candidate_bases.extend(bases)
    candidate_bases.extend(INVERTIBLE_BASES)
    checked_bases = set()
    for basis in candidate_bases:
        if basis in checked_bases:
            continue
        checked_bases.add(basis)
        solution = solve_square(
            tuple(COUPLING_COLUMNS[state] for state in basis),
            tuple(values[state] for state in basis),
        )
        require(solution is not None, "dual basis became singular")
        alpha, *beta = solution
        if any(
            alpha
            + sum(
                beta[core] * STATE_BITS[state][core]
                for core in range(3)
            )
            < values[state]
            for state in STATES
        ):
            continue
        candidate = alpha + sum(
            beta[core] * marginals[core] for core in range(3)
        )
        require(candidate >= primal, "weak duality failure")
        if candidate == primal:
            return primal
    raise RuntimeError("no exact dual certificate matches the primal optimum")


def lowfirst_reference_ledger():
    """The 223 low-first profiles in THM-2233's residue."""
    depth_one = tuple(
        (1, second, third)
        for second in range(1, 20)
        for third in range(second + 1, 20)
        if third >= 5
    )
    depth_two = (
        ((2, 3, 5),)
        + tuple((2, 4, third) for third in range(5, 20))
        + tuple((2, second, second + 2) for second in range(5, 18))
    )
    depth_three = (
        ((3, 3, 5), (3, 4, 5), (3, 4, 6))
        + tuple((3, 5, third) for third in range(6, 20))
        + tuple((3, second, second + 2) for second in range(6, 18))
    )
    ledger = tuple(dict.fromkeys(depth_one + depth_two + depth_three))
    require(len(depth_one) == 165, "depth-one ledger drift")
    require(len(depth_two) == 29, "depth-two ledger drift")
    require(len(depth_three) == 29, "depth-three ledger drift")
    require(len(ledger) == 223, "combined ledger drift")
    return ledger


def current_scalar_ledger():
    """THM-2239's 194-row ledger: depth one plus depth three."""
    ledger = tuple(
        profile
        for profile in lowfirst_reference_ledger()
        if profile[0] in (1, 3)
    )
    require(len(ledger) == 194, "THM-2239 current ledger drift")
    return ledger


def nonempty_subsets(values):
    values = tuple(values)
    return tuple(
        subset
        for size in range(1, len(values) + 1)
        for subset in combinations(values, size)
    )


def bases_for(profile, odd_checkpoints, mode):
    """Choose per-core centers.

    ``prior_latest`` and ``prior_earliest`` require s_j to precede every
    selected atom.  ``parity_earliest`` permits a later center and chooses
    the least nonnegative integer of the required parity; this is useful as
    a hostile extension at the first-depth boundary.
    """
    bases = []
    for valuation in profile:
        atom_times = tuple(valuation - k for k in odd_checkpoints)
        parity = valuation % 2
        if mode.startswith("prior_"):
            choices = tuple(
                time
                for time in range(min(atom_times))
                if time % 2 == parity
            )
            if not choices:
                return None
            bases.append(
                choices[-1] if mode == "prior_latest" else choices[0]
            )
        elif mode == "parity_earliest":
            bases.append(parity)
        else:
            raise ValueError(f"unknown base mode: {mode}")
    return tuple(bases)


def make_schedule(profile, odd_checkpoints, base_times):
    even_checkpoints = tuple(range(0, profile[0] + 1, 2))
    horizon = max(profile + base_times)
    even_at_time = {time: [] for time in range(horizon + 1)}
    odd_at_time = {time: [] for time in range(horizon + 1)}
    base_at_time = {time: [] for time in range(horizon + 1)}
    for clause, checkpoint in enumerate(even_checkpoints):
        for core, valuation in enumerate(profile):
            even_at_time[valuation - checkpoint].append((clause, core))
    for checkpoint in odd_checkpoints:
        for core, valuation in enumerate(profile):
            odd_at_time[valuation - checkpoint].append(core)
    for core, time in enumerate(base_times):
        base_at_time[time].append(core)
    correction = tuple(
        sum(
            (
                RHO ** (profile[core] - checkpoint - base_times[core])
                for checkpoint in odd_checkpoints
            ),
            Fraction(0),
        )
        for core in range(3)
    )
    require(all(value < 0 for value in correction), "correction sign drift")
    return {
        "profile": tuple(profile),
        "even_checkpoints": even_checkpoints,
        "odd_checkpoints": tuple(odd_checkpoints),
        "base_times": tuple(base_times),
        "horizon": horizon,
        "even_at_time": even_at_time,
        "odd_at_time": odd_at_time,
        "base_at_time": base_at_time,
        "correction": correction,
        "full_even_mask": (1 << len(even_checkpoints)) - 1,
    }


def update_summary(summary, time, state, schedule):
    even_mask, odd_count, base_mask = summary
    bits = STATE_BITS[state]
    for clause, core in schedule["even_at_time"][time]:
        if bits[core]:
            even_mask |= 1 << clause
    odd_count += sum(
        bits[core] for core in schedule["odd_at_time"][time]
    )
    for core in schedule["base_at_time"][time]:
        if bits[core]:
            base_mask |= 1 << core
    return even_mask, odd_count, base_mask


def reachable_future_summaries(schedule):
    horizon = schedule["horizon"]
    future = {horizon: {(0, 0, 0)}}
    current = {(0, 0, 0)}
    for time in range(horizon, 0, -1):
        current = {
            update_summary(summary, time, state, schedule)
            for summary in current
            for state in STATES
        }
        future[time - 1] = current
    return future


def payoff(summary, schedule, exact):
    even_mask, odd_count, base_mask = summary
    zero = Fraction(0) if exact else 0.0
    if even_mask != schedule["full_even_mask"]:
        return zero
    one_minus_q = (
        Fraction(1 - len(schedule["odd_checkpoints"]) + odd_count)
        if exact
        else float(1 - len(schedule["odd_checkpoints"]) + odd_count)
    )
    for core, correction in enumerate(schedule["correction"]):
        if (base_mask >> core) & 1:
            one_minus_q -= correction if exact else float(correction)
    return max(zero, one_minus_q)


def bellman_bound(schedule, exact=False):
    """Return the arbitrary-coupling Bellman bound."""
    future = reachable_future_summaries(schedule)
    horizon = schedule["horizon"]
    value = {}

    for parent in STATES:
        if exact:
            marginals = tuple(
                Fraction(2 - STATE_BITS[parent][core], PRIME)
                for core in range(3)
            )
        for summary in sorted(future[0]):
            child_values = tuple(
                payoff(update_summary(summary, 0, child, schedule), schedule, exact)
                for child in STATES
            )
            if exact:
                value[parent, summary] = maximize_exact(
                    child_values,
                    marginals,
                )
            else:
                value[parent, summary] = maximize_float(
                    child_values,
                    TRANSITION_VERTICES_FLOAT[parent],
                )

    for time in range(1, horizon):
        next_value = {}
        for parent in STATES:
            if exact:
                marginals = tuple(
                    Fraction(2 - STATE_BITS[parent][core], PRIME)
                    for core in range(3)
                )
            for summary in sorted(future[time]):
                child_values = tuple(
                    value[
                        child,
                        update_summary(summary, time, child, schedule),
                    ]
                    for child in STATES
                )
                if exact:
                    next_value[parent, summary] = maximize_exact(
                        child_values,
                        marginals,
                    )
                else:
                    next_value[parent, summary] = maximize_float(
                        child_values,
                        TRANSITION_VERTICES_FLOAT[parent],
                    )
        value = next_value

    terminal_values = tuple(
        value[
            state,
            update_summary((0, 0, 0), horizon, state, schedule),
        ]
        for state in STATES
    )
    if exact:
        return maximize_exact(
            terminal_values,
            (Fraction(1, 7),) * 3,
        ), tuple(len(future[time]) for time in range(horizon))
    return maximize_float(
        terminal_values,
        TERMINAL_VERTICES_FLOAT,
    ), tuple(len(future[time]) for time in range(horizon))


def audit_negative_carrier(schedule):
    """Hostile Boolean-box check of q<=0 on the selected carrier."""
    atoms = tuple(
        (checkpoint, core)
        for checkpoint in schedule["odd_checkpoints"]
        for core in range(3)
    )
    maximum = None
    checks = 0
    for atom_bits in product((0, 1), repeat=len(atoms)):
        if not all(
            any(
                atom_bits[index]
                for index, (owner, _) in enumerate(atoms)
                if owner == checkpoint
            )
            for checkpoint in schedule["odd_checkpoints"]
        ):
            continue
        for base_bits in product((0, 1), repeat=3):
            q_value = (
                Fraction(
                    len(schedule["odd_checkpoints"]) - sum(atom_bits)
                )
                + sum(
                    schedule["correction"][core] * base_bits[core]
                    for core in range(3)
                )
            )
            maximum = q_value if maximum is None else max(maximum, q_value)
            require(q_value <= 0, "q positive on hostile negative-carrier box")
            checks += 1
    require(maximum == 0, "negative-carrier maximum drift")
    return checks


def signed_configurations(profile, base_mode):
    candidate_odd = tuple(range(1, profile[0] + 1, 2))
    configurations = []
    for odd_checkpoints in nonempty_subsets(candidate_odd):
        base_times = bases_for(profile, odd_checkpoints, base_mode)
        if base_times is None:
            continue
        configurations.append(
            make_schedule(profile, odd_checkpoints, base_times)
        )
    return tuple(configurations)


def scout(modes):
    ledger = lowfirst_reference_ledger()
    print(
        f"lowfirst_reference_ledger_size={len(ledger)} "
        f"target={TARGET} decimal={float(TARGET):.15f}"
    )
    all_results = {}
    for mode in modes:
        results = {}
        inadmissible = []
        for index, profile in enumerate(ledger, 1):
            configurations = signed_configurations(profile, mode)
            if not configurations:
                inadmissible.append(profile)
                continue
            rows = []
            for schedule in configurations:
                bound, state_counts = bellman_bound(schedule, exact=False)
                rows.append(
                    (
                        bound,
                        schedule["odd_checkpoints"],
                        schedule["base_times"],
                        schedule["correction"],
                        max(state_counts),
                    )
                )
            results[profile] = min(rows)
            if index % 25 == 0:
                print(f"progress mode={mode} profiles={index}/{len(ledger)}")
        passing = tuple(
            profile
            for profile, row in results.items()
            if row[0] < float(TARGET) - 1e-13
        )
        near = tuple(
            sorted(results.items(), key=lambda item: item[1][0])[:12]
        )
        print(
            f"mode={mode} admissible={len(results)} "
            f"inadmissible={len(inadmissible)} passing={len(passing)}"
        )
        print(f"mode={mode} passing_profiles={passing}")
        print(f"mode={mode} best_rows={near}")
        all_results[mode] = (results, inadmissible, passing)
    return all_results


def exact_replay(profile, odd_checkpoints, base_mode):
    base_times = bases_for(profile, odd_checkpoints, base_mode)
    require(base_times is not None, "requested exact configuration inadmissible")
    schedule = make_schedule(profile, odd_checkpoints, base_times)
    checks = audit_negative_carrier(schedule)
    bound, state_counts = bellman_bound(schedule, exact=True)
    print(f"exact_profile={profile}")
    print(f"exact_odd_checkpoints={odd_checkpoints}")
    print(f"exact_base_mode={base_mode} base_times={base_times}")
    print(f"exact_corrections={schedule['correction']}")
    print(f"exact_negative_carrier_box_checks={checks}")
    print(f"exact_future_state_counts={state_counts}")
    print(f"exact_bound={bound} decimal={float(bound):.15f}")
    print(f"exact_target_margin={TARGET-bound}")
    print(f"exact_pass={bound < TARGET}")
    print(f"exact_distinct_primal_dual_lp_checks={EXACT_LP_CHECKS}")


def exact_current_ledger():
    """Exact audit of every low-first profile admitting an earlier center."""
    reference_ledger = lowfirst_reference_ledger()
    current_ledger = current_scalar_ledger()
    admissible_rows = []
    inadmissible = []
    for profile in reference_ledger:
        configurations = signed_configurations(profile, "prior_earliest")
        if not configurations:
            inadmissible.append(profile)
            continue
        require(
            len(configurations) == 1,
            "low-first prior-center census unexpectedly has multiple choices",
        )
        schedule = configurations[0]
        require(
            schedule["odd_checkpoints"] == (1,),
            "low-first admissible odd checkpoint drift",
        )
        require(
            all(
                base < valuation - 1
                and (valuation - 1 - base) % 2 == 1
                for valuation, base in zip(
                    profile,
                    schedule["base_times"],
                )
            ),
            "strict prior-center parity drift",
        )
        sign_checks = audit_negative_carrier(schedule)
        bound, state_counts = bellman_bound(schedule, exact=True)
        admissible_rows.append(
            (
                profile,
                schedule["base_times"],
                schedule["correction"],
                bound,
                TARGET - bound,
                sign_checks,
                max(state_counts),
            )
        )

    passing = tuple(row[0] for row in admissible_rows if row[3] < TARGET)
    failing = tuple(row[0] for row in admissible_rows if row[3] >= TARGET)
    expected_passing = tuple(
        profile
        for profile in reference_ledger
        if profile[0] == 2 or (profile[0] == 3 and profile != (3, 4, 5))
    )
    require(tuple(inadmissible) == tuple(
        profile for profile in reference_ledger if profile[0] == 1
    ), "inadmissible depth-one classification drift")
    require(passing == expected_passing, "exact passing classification drift")
    require(failing == ((3, 4, 5),), "exact failure classification drift")
    novel_passing = tuple(
        profile for profile in passing if profile in current_ledger
    )
    expected_novel_passing = tuple(
        profile
        for profile in current_ledger
        if profile[0] == 3 and profile != (3, 4, 5)
    )
    require(
        novel_passing == expected_novel_passing
        and len(novel_passing) == 28,
        "novel depth-three passing classification drift",
    )

    # The exceptional prior-centered row has four possible per-core centers.
    # Audit all of them exactly rather than assuming that smaller correction
    # coefficients automatically give a smaller relaxed Bellman value.
    exceptional = (3, 4, 5)
    exceptional_base_options = ((1,), (0, 2), (1, 3))
    exceptional_rows = []
    for base_times in product(*exceptional_base_options):
        schedule = make_schedule(exceptional, (1,), base_times)
        bound, _ = bellman_bound(schedule, exact=True)
        exceptional_rows.append((base_times, schedule["correction"], bound))
    exceptional_rows = tuple(sorted(exceptional_rows, key=lambda row: row[2]))
    require(
        exceptional_rows[0][0] == (1, 0, 1)
        and exceptional_rows[0][2] == Fraction(1393650030, 5710115047),
        "exceptional best prior-center bound drift",
    )
    require(
        all(row[2] > TARGET for row in exceptional_rows),
        "an exceptional prior center unexpectedly passes",
    )

    # Positive control: the guard-free three-bit weakening of THM-2239 still
    # reproduces the earlier signed (4,6,8) closure value.
    control_schedule = make_schedule((4, 6, 8), (1, 3), (0, 0, 0))
    control_bound, _ = bellman_bound(control_schedule, exact=True)
    require(
        control_bound
        == Fraction(17322925655936326, 358301251098635299),
        "weakened geometric positive-control bound drift",
    )

    row_transcript = "\n".join(
        (
            f"{profile}|{base_times}|"
            f"{','.join(map(str, correction))}|{bound}|{margin}"
        )
        for (
            profile,
            base_times,
            correction,
            bound,
            margin,
            _,
            _,
        ) in admissible_rows
    )
    row_digest = sha256(row_transcript.encode()).hexdigest()
    closest_passing = min(
        (row for row in admissible_rows if row[3] < TARGET),
        key=lambda row: row[4],
    )
    strongest_passing = min(
        (row for row in admissible_rows if row[3] < TARGET),
        key=lambda row: row[3],
    )
    closest_novel_passing = min(
        (
            row
            for row in admissible_rows
            if row[0] in novel_passing
        ),
        key=lambda row: row[4],
    )

    print("theorem=THM-2244")
    print("thm2233_lowfirst_reference_ledger=223")
    print("thm2239_current_ledger=194")
    print("exact_prior_center_inadmissible_depth_one=165")
    print("exact_prior_center_admissible_depth_two_three=58")
    print("selected_odd_checkpoints=(1,)")
    print("retained_even_checkpoints=all")
    print(f"exact_row_digest={row_digest}")
    print(f"exact_raw_passing_count={len(passing)}")
    print(f"exact_raw_passing_profiles={passing}")
    print(f"exact_novel_passing_count={len(novel_passing)}")
    print(f"exact_novel_passing_profiles={novel_passing}")
    print(f"exact_failing_profiles={failing}")
    print(f"closest_raw_passing_row={closest_passing}")
    print(f"closest_novel_passing_row={closest_novel_passing}")
    print(f"strongest_passing_row={strongest_passing}")
    print(f"exceptional_345_all_prior_centers={exceptional_rows}")
    print(f"weakened_geometric_positive_control={control_bound}")
    print("exact_primal_dual_parity=PASS")
    print(f"exact_distinct_primal_dual_lp_checks={EXACT_LP_CHECKS}")
    print("raw_exclusions=57")
    print("novel_exclusions_after_thm2239=28")
    print("resulting_ledger=166")
    print("resulting_ledger_classification=165_depth_one_plus_(3,4,5)")


def parse_profile(text):
    values = tuple(int(part) for part in text.split(","))
    if len(values) != 3:
        raise ValueError("profile must have three comma-separated integers")
    return values


def parse_checkpoints(text):
    return tuple(int(part) for part in text.split(",") if part)


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--scout",
        action="store_true",
        help="float-scout the complete current ledger",
    )
    parser.add_argument(
        "--modes",
        default="prior_latest,prior_earliest,parity_earliest",
        help="comma-separated base modes for --scout",
    )
    parser.add_argument(
        "--exact-profile",
        type=parse_profile,
        help="exact replay profile, e.g. 2,3,5",
    )
    parser.add_argument(
        "--exact-ledger",
        action="store_true",
        help="exactly audit every current row with a prior center",
    )
    parser.add_argument(
        "--odd",
        type=parse_checkpoints,
        default=(1,),
        help="odd checkpoint subset for exact replay",
    )
    parser.add_argument(
        "--base-mode",
        default="prior_earliest",
        help="base mode for exact replay",
    )
    args = parser.parse_args()

    require(len(INVERTIBLE_BASES) == 58, "invertible basis census drift")
    if args.scout:
        scout(tuple(args.modes.split(",")))
    if args.exact_ledger:
        exact_current_ledger()
    if args.exact_profile is not None:
        exact_replay(args.exact_profile, args.odd, args.base_mode)
    if not args.scout and not args.exact_ledger and args.exact_profile is None:
        parser.error("choose --scout, --exact-ledger, and/or --exact-profile")


if __name__ == "__main__":
    main()
