#!/usr/bin/env python3
"""Exact unit-time positive-set Bellman referee for THM-2229.

For a scalar valuation profile (lambda_1,lambda_2,lambda_3), this script
keeps the positive unit residual A_+ as one clause and every available even
divided-blocker union U_k as another clause.  It resolves time one factor of
13 at a time.  The only relaxation is arbitrary coupling subject to the
pointwise root-fibre marginal caps

    A_+ : 10/13,
    one unit danger atom : 2/13.

Every marginal LP is solved twice, by exact primal basic-feasible-solution
enumeration and by independent exact dual-vertex enumeration.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations


P = 13
TARGET = Fraction(961, 6930)
UNIT_CAP = Fraction(2, P)
POSITIVE_SET_CAP = Fraction(10, P)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def solve_square(matrix, rhs):
    """Solve a square rational system, returning None when singular."""
    n = len(rhs)
    aug = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if aug[r][col]), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for row in range(n):
            if row == col:
                continue
            scale = aug[row][col]
            if scale:
                aug[row] = [
                    aug[row][j] - scale * aug[col][j]
                    for j in range(n + 1)
                ]
    return tuple(aug[i][-1] for i in range(n))


def minimize_marginal_dual(values, active_caps):
    """Independently enumerate vertices of the exact marginal-LP dual."""
    m = len(active_caps)
    if m == 0:
        return values[0]
    patterns = 1 << m
    rows = []
    rhs = []
    for pattern in range(patterns):
        rows.append(
            (Fraction(1),)
            + tuple(Fraction((pattern >> i) & 1) for i in range(m))
        )
        rhs.append(values[pattern])
    for i in range(m):
        rows.append(
            (Fraction(0),)
            + tuple(Fraction(int(i == j)) for j in range(m))
        )
        rhs.append(Fraction(0))

    best = None
    rank = m + 1
    for active in combinations(range(len(rows)), rank):
        matrix = tuple(rows[i] for i in active)
        solution = solve_square(matrix, tuple(rhs[i] for i in active))
        if solution is None:
            continue
        alpha, *beta = solution
        if any(value < 0 for value in beta):
            continue
        if any(
            alpha
            + sum(beta[i] for i in range(m) if (pattern >> i) & 1)
            < values[pattern]
            for pattern in range(patterns)
        ):
            continue
        candidate = alpha + sum(
            active_caps[i] * beta[i] for i in range(m)
        )
        if best is None or candidate < best:
            best = candidate
    require(best is not None, "marginal dual has no enumerated vertex")
    return best


def maximize_marginal_lp(values, active_caps):
    """Maximize an expectation under upper bounds on active marginals."""
    m = len(active_caps)
    if m == 0:
        return values[0]
    patterns = 1 << m
    columns = []
    objective = []

    for pattern in range(patterns):
        columns.append(
            (Fraction(1),)
            + tuple(Fraction((pattern >> i) & 1) for i in range(m))
        )
        objective.append(values[pattern])
    for i in range(m):
        columns.append(
            (Fraction(0),)
            + tuple(Fraction(int(i == j)) for j in range(m))
        )
        objective.append(Fraction(0))

    rhs = (Fraction(1),) + tuple(active_caps)
    rank = m + 1
    best = None
    for basis in combinations(range(len(columns)), rank):
        matrix = tuple(
            tuple(columns[basis[col]][row] for col in range(rank))
            for row in range(rank)
        )
        solution = solve_square(matrix, rhs)
        if solution is None or any(value < 0 for value in solution):
            continue
        candidate = sum(
            objective[basis[i]] * solution[i] for i in range(rank)
        )
        if best is None or candidate > best:
            best = candidate

    require(best is not None, "marginal primal has no basic feasible solution")
    dual = minimize_marginal_dual(values, active_caps)
    require(best == dual, f"primal/dual drift: {best} != {dual}")
    return best


def eliminate_time(value, clause_caps, clause_count):
    """Eliminate one root digit, allowing arbitrary active-clause coupling."""
    active = tuple(sorted(clause_caps))
    local_clause_bits = tuple(1 << clause for clause in active)
    output = []
    for later in range(1 << clause_count):
        local_values = []
        for pattern in range(1 << len(active)):
            newly_satisfied = 0
            for i, bit in enumerate(local_clause_bits):
                if (pattern >> i) & 1:
                    newly_satisfied |= bit
            local_values.append(value[later | newly_satisfied])
        output.append(
            maximize_marginal_lp(
                tuple(local_values),
                tuple(clause_caps[clause] for clause in active),
            )
        )
    return tuple(output)


def profile_bound(profile):
    """Return the exact positive-set/even-checkpoint Bellman bound."""
    lambda_1, _, lambda_3 = profile
    require(2 <= lambda_1 <= 5, f"profile outside theorem census: {profile}")
    even_checkpoints = tuple(range(0, lambda_1 + 1, 2))
    clause_count = 1 + len(even_checkpoints)
    all_clauses = (1 << clause_count) - 1

    caps_by_time = {0: {0: POSITIVE_SET_CAP}}
    for clause, checkpoint in enumerate(even_checkpoints, start=1):
        for valuation in profile:
            time = valuation - checkpoint
            require(time >= 0, "negative unit-time offset")
            caps = caps_by_time.setdefault(time, {})
            caps[clause] = caps.get(clause, Fraction(0)) + UNIT_CAP

    value = tuple(
        Fraction(int(state == all_clauses))
        for state in range(1 << clause_count)
    )
    for time in range(max(caps_by_time) + 1):
        value = eliminate_time(
            value,
            caps_by_time.get(time, {}),
            clause_count,
        )
    return value[0]


def all_profiles():
    return tuple(
        (lambda_1, lambda_2, lambda_3)
        for lambda_3 in range(2, 20)
        for lambda_2 in range(1, lambda_3)
        for lambda_1 in range(1, lambda_2 + 1)
    )


def profiles_with_first(first):
    return tuple(
        (first, lambda_2, lambda_3)
        for lambda_2 in range(first, 19)
        for lambda_3 in range(lambda_2 + 1, 20)
    )


def digest_profiles(profiles):
    payload = ";".join(",".join(map(str, profile)) for profile in profiles)
    return sha256(payload.encode("ascii")).hexdigest()


def expected_closed_first_two():
    return {
        profile
        for profile in profiles_with_first(2)
        if (
            (profile[1] == 3 and profile[2] >= 6)
            or (profile[1] >= 5 and profile[2] != profile[1] + 2)
        )
    }


def expected_closed_first_three():
    return {
        profile
        for profile in profiles_with_first(3)
        if (
            (
                profile[1] == 3
                and (profile[2] == 4 or profile[2] >= 6)
            )
            or (profile[1] == 4 and profile[2] >= 7)
            or (profile[1] >= 6 and profile[2] != profile[1] + 2)
        )
    }


def main():
    universe = set(all_profiles())
    require(len(universe) == 1140, "legal profile universe drift")

    low_exclusions = {
        (1, 1, 2),
        (1, 1, 3),
        (1, 2, 3),
        (2, 2, 3),
        (1, 1, 4),
        (1, 2, 4),
        (1, 3, 4),
        (2, 2, 4),
        (2, 3, 4),
        (3, 3, 4),
    }
    thm2226_residue = {
        (4, 6, 8),
        (4, 6, 9),
        (4, 7, 8),
        (4, 7, 9),
        (5, 7, 9),
        (5, 7, 10),
        (5, 8, 9),
        (5, 8, 10),
    }
    first_four_five = set(profiles_with_first(4)) | set(
        profiles_with_first(5)
    )
    thm2226_closed = first_four_five - thm2226_residue
    thm2224_closed = {
        profile for profile in universe if profile[0] >= 6
    }
    thm2227_closed = {
        (4, 6, 9),
        (4, 7, 8),
        (4, 7, 9),
        (5, 7, 10),
        (5, 8, 9),
        (5, 8, 10),
    }

    require(len(low_exclusions) == 10, "low-exclusion ledger drift")
    require(len(thm2226_closed) == 217, "THM-2226 census drift")
    require(len(thm2224_closed) == 455, "THM-2224 census drift")
    prior_closed = low_exclusions | thm2226_closed | thm2224_closed
    require(len(prior_closed) == 682, "prior combined ledger drift")
    require(len(universe - prior_closed) == 458, "prior residue drift")
    current_prior_closed = prior_closed | thm2227_closed
    require(len(current_prior_closed) == 688, "current prior ledger drift")
    require(
        len(universe - current_prior_closed) == 452,
        "post-THM-2227 residue drift",
    )

    bounds = {}
    raw_closed = set()
    for first in range(2, 6):
        for profile in profiles_with_first(first):
            bound = profile_bound(profile)
            bounds[profile] = bound
            if bound < TARGET:
                raw_closed.add(profile)

    by_first = {
        first: {
            profile
            for profile in raw_closed
            if profile[0] == first
        }
        for first in range(2, 6)
    }
    require(
        tuple(len(by_first[first]) for first in range(2, 6))
        == (106, 107, 119, 104),
        "per-first-depth closure drift",
    )
    require(
        by_first[2] == expected_closed_first_two(),
        "lambda_1=2 closed-profile classification drift",
    )
    require(
        by_first[3] == expected_closed_first_three(),
        "lambda_1=3 closed-profile classification drift",
    )
    require(
        set(profiles_with_first(4)) - by_first[4] == {(4, 6, 8)},
        "lambda_1=4 survivor drift",
    )
    require(
        set(profiles_with_first(5)) - by_first[5] == {(5, 7, 9)},
        "lambda_1=5 survivor drift",
    )
    require(len(raw_closed) == 436, "raw closure count drift")

    expected_passing_bounds = {
        2: Fraction(625252, 4826809),
        3: Fraction(49000, 371293),
        4: Fraction(2896, 28561),
        5: Fraction(74441360, 815730721),
    }
    expected_failing_bounds = {
        2: Fraction(4480, 28561),
        3: Fraction(824660, 4826809),
        4: Fraction(60792, 371293),
        5: Fraction(59340, 371293),
    }
    for first in range(2, 6):
        row_profiles = set(profiles_with_first(first))
        passing = tuple(
            bounds[profile] for profile in row_profiles if profile in raw_closed
        )
        failing = tuple(
            bounds[profile] for profile in row_profiles if profile not in raw_closed
        )
        require(
            max(passing) == expected_passing_bounds[first],
            f"first-depth {first} worst passing bound drift",
        )
        require(
            min(failing) == expected_failing_bounds[first],
            f"first-depth {first} least failing bound drift",
        )

    low_overlap = raw_closed & low_exclusions
    sieve_overlap = raw_closed & thm2226_closed
    high_overlap = raw_closed & thm2224_closed
    hidden_parity_overlap = raw_closed & thm2227_closed
    require(low_overlap == {(3, 3, 4)}, "low overlap drift")
    require(len(sieve_overlap) == 217, "THM-2226 overlap drift")
    require(not high_overlap, "THM-2224 overlap should be empty")
    require(
        hidden_parity_overlap == thm2227_closed,
        "THM-2227 overlap drift",
    )

    novel_pre_2227 = raw_closed - prior_closed
    novel_current = raw_closed - current_prior_closed
    expected_new_sieve_residue = {
        (4, 6, 9),
        (4, 7, 8),
        (4, 7, 9),
        (5, 7, 10),
        (5, 8, 9),
        (5, 8, 10),
    }
    require(
        len(novel_current & set(profiles_with_first(2))) == 106,
        "novel lambda_1=2 count drift",
    )
    require(
        len(novel_current & set(profiles_with_first(3))) == 106,
        "novel lambda_1=3 count drift",
    )
    require(
        novel_pre_2227 & thm2226_residue == expected_new_sieve_residue,
        "new THM-2226-residue exclusions drift",
    )
    require(len(novel_pre_2227) == 218, "pre-THM-2227 closure count drift")
    require(len(novel_current) == 212, "current additive closure count drift")

    combined_closed = current_prior_closed | raw_closed
    remaining = tuple(sorted(universe - combined_closed))
    novel_pre_2227_sorted = tuple(sorted(novel_pre_2227))
    novel_current_sorted = tuple(sorted(novel_current))
    require(len(combined_closed) == 900, "combined closure count drift")
    require(len(remaining) == 240, "combined remaining-ledger drift")
    require(
        sum(profile[0] <= 3 for profile in remaining) == 238,
        "low-first remaining count drift",
    )
    require(
        {profile for profile in remaining if profile[0] >= 4}
        == {(4, 6, 8), (5, 7, 9)},
        "high-first remaining pair drift",
    )

    print(f"target={TARGET}")
    print(f"positive_set_root_cap={POSITIVE_SET_CAP}")
    print(f"unit_danger_root_cap={UNIT_CAP}")
    for first in range(2, 6):
        row_profiles = set(profiles_with_first(first))
        passing = tuple(
            bounds[profile] for profile in row_profiles if profile in raw_closed
        )
        failing = tuple(
            bounds[profile] for profile in row_profiles if profile not in raw_closed
        )
        print(
            f"lambda1={first} "
            f"closed={len(by_first[first])}/{len(row_profiles)} "
            f"worst_passing={max(passing)} "
            f"target_minus_worst={TARGET-max(passing)} "
            f"least_failing={min(failing)}"
        )
    print(f"raw_closed={len(raw_closed)}")
    print(f"overlap_low_ten={len(low_overlap)} profiles={tuple(sorted(low_overlap))}")
    print(f"overlap_thm2226={len(sieve_overlap)}")
    print(f"overlap_thm2224={len(high_overlap)}")
    print(f"overlap_thm2227={len(hidden_parity_overlap)}")
    print(
        "new_thm2226_residue_exclusions="
        f"{tuple(sorted(expected_new_sieve_residue))}"
    )
    print(f"independent_new_from_458_ledger={len(novel_pre_2227_sorted)}")
    print(f"additive_new_after_thm2227={len(novel_current_sorted)}")
    print(
        "novel_lambda1_2_rule="
        "(lambda2=3 and lambda3>=6) or "
        "(lambda2>=5 and lambda3!=lambda2+2)"
    )
    print(
        "novel_lambda1_3_rule="
        "(lambda2=3 and lambda3>=6) or "
        "(lambda2=4 and lambda3>=7) or "
        "(lambda2>=6 and lambda3!=lambda2+2)"
    )
    print(
        "independent_new_digest="
        f"{digest_profiles(novel_pre_2227_sorted)}"
    )
    print(f"additive_new_digest={digest_profiles(novel_current_sorted)}")
    print(f"combined_remaining={len(remaining)}")
    print("remaining_lambda1_1=all 165 profiles with lambda3>=5")
    print(
        "remaining_lambda1_2="
        "(2,2,c) for 5<=c<=19; (2,3,5); "
        "(2,4,c) for 5<=c<=19; (2,b,b+2) for 5<=b<=17"
    )
    print(
        "remaining_lambda1_3="
        "(3,3,5); (3,4,5); (3,4,6); "
        "(3,5,c) for 6<=c<=19; (3,b,b+2) for 6<=b<=17"
    )
    print(f"remaining_digest={digest_profiles(remaining)}")
    print("remaining_lambda1_le_3=238")
    print("remaining_high_first=((4, 6, 8), (5, 7, 9))")
    print("exact_primal_dual_vertex_parity=PASS")


if __name__ == "__main__":
    main()
