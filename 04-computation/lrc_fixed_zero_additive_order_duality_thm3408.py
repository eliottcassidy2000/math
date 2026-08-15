#!/usr/bin/env python3
"""Exact companion for THM-3408's fixed-zero additive-order duality.

The program is self-contained and standard-library only.  It checks the
literal stratum-density formula, performs two complete arithmetic cutoffs for
the exceptional quotient orders, verifies the six-core prime-loss table and
the q=21,22,102 hostiles, and solves the finite q<=20000 fractional games by
an exact rational simplex.  Every reported game value is accompanied and
checked by both a sheet-stratum distribution and an owner-order distribution,
so the finite scan does not depend on floating point or an unverified solver
status.  Runtime gates survive ``python -O``.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, isqrt, lcm
from pathlib import Path


BASES = (15, 16, 18, 20, 24)
B_EXCEPTIONS = (
    21, 22, 26, 28, 33, 35, 39, 42, 44, 46, 50, 52, 56, 57, 63, 65,
    70, 74, 76, 77, 78, 84, 91, 99, 102, 104, 114, 117, 130, 143,
    156, 186, 190,
)
H_CORRIDOR = (22, 33, 44, 46, 50, 102)
LOSS_PRIMES = (7, 11, 13, 19, 31, 37)
Q102_SHEET_WEIGHTS = {
    6: Fraction(1, 25),
    34: Fraction(4, 25),
    51: Fraction(4, 25),
    102: Fraction(16, 25),
}
Q102_OWNER_WEIGHTS = {
    2: Fraction(7, 150),
    3: Fraction(7, 150),
    17: Fraction(4, 25),
    102: Fraction(56, 75),
}
EXPECTED_SEMANTIC_DIGEST = "ab976baf418e2610aee1558eb249a160998426812fa96dddbc83b5a4a528541e"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


def prime_factorization(n):
    require(n >= 1, n)
    answer = []
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            exponent = 0
            while n % divisor == 0:
                n //= divisor
                exponent += 1
            answer.append((divisor, exponent))
        divisor += 1
    if n > 1:
        answer.append((n, 1))
    return tuple(answer)


def prime_divisors(n):
    return tuple(prime for prime, _ in prime_factorization(n))


def phi(n):
    answer = n
    for prime in prime_divisors(n):
        answer = answer // prime * (prime - 1)
    return answer


def divisors_from_factorization(factors):
    answer = [1]
    for prime, exponent in factors:
        old = tuple(answer)
        power = 1
        for _ in range(exponent):
            power *= prime
            answer.extend(divisor * power for divisor in old)
    return tuple(sorted(answer))


def divisors(n):
    return divisors_from_factorization(prime_factorization(n))


def is_prime(n):
    if n < 2:
        return False
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            return False
        divisor += 1
    return True


def coprime_prefix_count(n, bound):
    """Count 1<=r<=bound coprime to n by exact inclusion-exclusion."""
    primes = prime_divisors(n)
    answer = 0
    for mask in range(1 << len(primes)):
        divisor = 1
        sign = 1
        for index, prime in enumerate(primes):
            if mask & (1 << index):
                divisor *= prime
                sign = -sign
        answer += sign * (bound // divisor)
    return answer


def alpha(order):
    if order == 1:
        return Fraction(1)
    cutoff = (order - 1) // 14
    return Fraction(2 * coprime_prefix_count(order, cutoff), phi(order))


def base_free(n):
    return not any(n % base == 0 for base in BASES)


def score(n):
    """Return phi(n)/2^omega(n)."""
    return Fraction(phi(n), 1 << len(prime_divisors(n)))


def primes_below(bound):
    return tuple(number for number in range(2, bound) if is_prime(number))


def enumerate_score_below(bound):
    """Enumerate all n with phi(n)/2^omega(n)<bound, without an n cutoff.

    For n=prod p^a the score is prod p^(a-1)(p-1)/2.  Once the
    primes already chosen have score rho, a new prime must satisfy
    rho*(p-1)/2<bound.  The only factor below one is the first power of 2,
    so all possible primes occur below 4*bound+1.  Increasing-prime DFS and
    increasing exponents therefore exhaust every factorization exactly once.
    """
    prime_limit = (4 * bound).numerator // (4 * bound).denominator + 2
    primes = primes_below(prime_limit)
    found = set()

    def visit(start, number, rho):
        require(number not in found, ("duplicate factorization", number))
        found.add(number)
        for index in range(start, len(primes)):
            prime = primes[index]
            new_rho = rho * Fraction(prime - 1, 2)
            if new_rho >= bound:
                break
            power = prime
            while new_rho < bound:
                visit(index + 1, number * power, new_rho)
                power *= prime
                new_rho *= prime

    visit(0, 1, Fraction(1))
    result = tuple(sorted(found))
    for number in result:
        require(score(number) < bound, (number, score(number), bound))
    return result


def cyclic_distance_numerator(residue, modulus):
    residue %= modulus
    return min(residue, modulus - residue)


def additive_order(sheet, modulus):
    return modulus // gcd(sheet, modulus)


def danger_set(q, speed):
    require(q >= 2 and speed % q != 0, (q, speed))
    return frozenset(
        sheet
        for sheet in range(q)
        if 14 * cyclic_distance_numerator(speed * sheet, q) < q
    )


def stratum(q, order):
    require(q % order == 0, (q, order))
    return frozenset(sheet for sheet in range(q) if additive_order(sheet, q) == order)


def image_order(q, owner_order, sheet_order):
    require(q % owner_order == 0 and q % sheet_order == 0, (q, owner_order, sheet_order))
    return owner_order // gcd(owner_order, q // sheet_order)


def units(modulus):
    return tuple(residue for residue in range(1, modulus) if gcd(residue, modulus) == 1)


def sign_representatives(modulus):
    return tuple(residue for residue in units(modulus) if residue <= (-residue) % modulus)


def mask_of(values):
    answer = 0
    for value in values:
        answer |= 1 << value
    return answer


def literal_coverage_types(q):
    by_mask = {}
    for speed in range(1, q):
        mask = mask_of(danger_set(q, speed))
        if mask not in by_mask or speed < by_mask[mask]:
            by_mask[mask] = speed
    return tuple(sorted((speed, mask) for mask, speed in by_mask.items()))


def minimum_literal_cover(q):
    coverage_types = literal_coverage_types(q)
    full = (1 << q) - 1
    tested = 0
    for rank in range(1, len(coverage_types) + 1):
        for chosen in combinations(coverage_types, rank):
            tested += 1
            union = 0
            for _, mask in chosen:
                union |= mask
            if union == full:
                return rank, tuple(speed for speed, _ in chosen), len(coverage_types), tested
    return None, (), len(coverage_types), tested


class ExactLinearProgram:
    """Exact two-phase simplex for max c*x subject to A*x<=b and x>=0.

    This is used only to discover certificates.  The caller separately checks
    primal feasibility, dual feasibility, and equality of their values.
    """

    def __init__(self, matrix, bounds, objective):
        self.row_count = len(bounds)
        self.variable_count = len(objective)
        rows = self.row_count
        columns = self.variable_count
        require(rows >= 1 and columns >= 1, (rows, columns))
        require(all(len(row) == columns for row in matrix), "ragged LP matrix")
        self.basic = [columns + row for row in range(rows)]
        self.nonbasic = list(range(columns)) + [-1]
        self.tableau = [
            [Fraction(0) for _ in range(columns + 2)]
            for _ in range(rows + 2)
        ]
        for row in range(rows):
            for column in range(columns):
                self.tableau[row][column] = Fraction(matrix[row][column])
            self.tableau[row][columns] = Fraction(-1)
            self.tableau[row][columns + 1] = Fraction(bounds[row])
        for column in range(columns):
            self.tableau[rows][column] = -Fraction(objective[column])
        self.tableau[rows + 1][columns] = Fraction(1)

    def pivot(self, leaving_row, entering_column):
        rows = self.row_count
        columns = self.variable_count
        tableau = self.tableau
        pivot_value = tableau[leaving_row][entering_column]
        require(pivot_value != 0, ("zero pivot", leaving_row, entering_column))
        for row in range(rows + 2):
            if row == leaving_row:
                continue
            multiplier = tableau[row][entering_column]
            if multiplier == 0:
                continue
            for column in range(columns + 2):
                if column != entering_column:
                    tableau[row][column] -= (
                        tableau[leaving_row][column] * multiplier / pivot_value
                    )
        for column in range(columns + 2):
            if column != entering_column:
                tableau[leaving_row][column] /= pivot_value
        for row in range(rows + 2):
            if row != leaving_row:
                tableau[row][entering_column] /= -pivot_value
        tableau[leaving_row][entering_column] = Fraction(1, 1) / pivot_value
        self.basic[leaving_row], self.nonbasic[entering_column] = (
            self.nonbasic[entering_column],
            self.basic[leaving_row],
        )

    def simplex(self, phase):
        rows = self.row_count
        columns = self.variable_count
        objective_row = rows + 1 if phase == 1 else rows
        while True:
            eligible = tuple(
                column
                for column in range(columns + 1)
                if not (phase == 2 and self.nonbasic[column] == -1)
            )
            entering = min(
                eligible,
                key=lambda column: (
                    self.tableau[objective_row][column],
                    self.nonbasic[column],
                ),
            )
            if self.tableau[objective_row][entering] >= 0:
                return True
            possible = tuple(
                row
                for row in range(rows)
                if self.tableau[row][entering] > 0
            )
            if not possible:
                return False
            leaving = min(
                possible,
                key=lambda row: (
                    self.tableau[row][columns + 1]
                    / self.tableau[row][entering],
                    self.basic[row],
                ),
            )
            self.pivot(leaving, entering)

    def solve(self):
        rows = self.row_count
        columns = self.variable_count
        first_leaving = min(
            range(rows), key=lambda row: self.tableau[row][columns + 1]
        )
        if self.tableau[first_leaving][columns + 1] < 0:
            self.pivot(first_leaving, columns)
            require(self.simplex(1), "infeasible phase-I ray")
            require(self.tableau[rows + 1][columns + 1] == 0, "infeasible LP")
            for row in range(rows):
                if self.basic[row] == -1:
                    candidates = tuple(
                        column
                        for column in range(columns + 1)
                        if self.nonbasic[column] != -1
                        and self.tableau[row][column] != 0
                    )
                    if candidates:
                        entering = max(
                            candidates,
                            key=lambda column: (
                                abs(self.tableau[row][column]),
                                -self.nonbasic[column],
                            ),
                        )
                        self.pivot(row, entering)
                    break
        require(self.simplex(2), "unbounded LP")
        solution = [Fraction(0) for _ in range(columns)]
        for row, label in enumerate(self.basic):
            if label < columns:
                solution[label] = self.tableau[row][columns + 1]
        dual = []
        for constraint in range(rows):
            slack_label = columns + constraint
            positions = tuple(
                column
                for column, label in enumerate(self.nonbasic)
                if label == slack_label
            )
            dual.append(
                self.tableau[rows][positions[0]] if positions else Fraction(0)
            )
        return self.tableau[rows][columns + 1], tuple(solution), tuple(dual)


def fractional_game(q):
    orders = divisors(q)[1:]
    dimension = len(orders)
    payoff = tuple(
        tuple(alpha(image_order(q, owner, sheet)) for sheet in orders)
        for owner in orders
    )
    matrix = tuple(
        row + (Fraction(-1),)
        for row in payoff
    ) + (
        (Fraction(1),) * dimension + (Fraction(0),),
        (Fraction(-1),) * dimension + (Fraction(0),),
    )
    bounds = (Fraction(0),) * dimension + (Fraction(1), Fraction(-1))
    objective = (Fraction(0),) * dimension + (Fraction(-1),)
    objective_value, solution, dual = ExactLinearProgram(
        matrix, bounds, objective
    ).solve()
    sheet_weights = solution[:-1]
    value = solution[-1]
    owner_weights = dual[:dimension]
    owner_loads = tuple(
        sum(row[index] * sheet_weights[index] for index in range(dimension))
        for row in payoff
    )
    sheet_loads = tuple(
        sum(owner_weights[index] * payoff[index][column] for index in range(dimension))
        for column in range(dimension)
    )
    require(all(weight >= 0 for weight in sheet_weights), (q, "sheet weights"))
    require(all(weight >= 0 for weight in owner_weights), (q, "owner weights"))
    require(sum(sheet_weights) == 1, (q, "sheet normalization", sum(sheet_weights)))
    require(sum(owner_weights) == 1, (q, "owner normalization", sum(owner_weights)))
    require(max(owner_loads) == value, (q, "primal", max(owner_loads), value))
    require(min(sheet_loads) == value, (q, "dual", min(sheet_loads), value))
    require(objective_value == -value, (q, objective_value, value))
    require(all(multiplier >= 0 for multiplier in dual), (q, "LP dual sign"))
    sheet_support = tuple(
        (order, weight)
        for order, weight in zip(orders, sheet_weights)
        if weight
    )
    owner_support = tuple(
        (order, weight)
        for order, weight in zip(orders, owner_weights)
        if weight
    )
    active_owners = tuple(
        order for order, load in zip(orders, owner_loads) if load == value
    )
    tight_sheets = tuple(
        order for order, load in zip(orders, sheet_loads) if load == value
    )
    return value, sheet_support, owner_support, active_owners, tight_sheets


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    semantic = ExactDigest()

    density_checks = 0
    fibre_checks = 0
    pullback_checks = 0
    for q in range(2, 121):
        for owner_order in divisors(q)[1:]:
            for mode in units(owner_order):
                speed = (q // owner_order) * mode
                literal = danger_set(q, speed)
                reduced = frozenset(
                    sheet
                    for sheet in range(q)
                    if 14
                    * cyclic_distance_numerator(mode * sheet, owner_order)
                    < owner_order
                )
                require(literal == reduced, ("mode pullback", q, owner_order, mode))
                pullback_checks += 1
                for sheet_order in divisors(q)[1:]:
                    sheets = stratum(q, sheet_order)
                    image = image_order(q, owner_order, sheet_order)
                    covered = literal & sheets
                    require(
                        Fraction(len(covered), len(sheets)) == alpha(image),
                        (
                            "density",
                            q,
                            owner_order,
                            mode,
                            sheet_order,
                            len(covered),
                            len(sheets),
                            image,
                            alpha(image),
                        ),
                    )
                    density_checks += 1
            for sheet_order in divisors(q)[1:]:
                image = image_order(q, owner_order, sheet_order)
                reductions = tuple(
                    ((q // sheet_order) * unit) % owner_order
                    for unit in units(sheet_order)
                )
                multiplicities = tuple(
                    reductions.count(residue)
                    for residue in sorted(set(reductions))
                )
                require(
                    len(set(reductions)) == phi(image),
                    ("image support", q, owner_order, sheet_order, image),
                )
                require(
                    len(set(multiplicities)) == 1
                    and multiplicities[0] == phi(sheet_order) // phi(image),
                    ("equal fibres", q, owner_order, sheet_order, image, multiplicities),
                )
                fibre_checks += 1

    score_4_25 = enumerate_score_below(Fraction(350, 3))
    score_7_44 = enumerate_score_below(Fraction(616, 5))
    require(len(score_4_25) == 1710 and max(score_4_25) == 30030, "4/25 cutoff")
    require(len(score_7_44) == 1829 and max(score_7_44) == 39270, "7/44 cutoff")
    over_4_25 = tuple(
        number
        for number in score_4_25
        if number >= 2 and base_free(number) and alpha(number) > Fraction(4, 25)
    )
    require(over_4_25 == B_EXCEPTIONS, ("B classification", over_4_25))
    over_7_44_outside_b = tuple(
        number
        for number in score_7_44
        if number >= 2
        and base_free(number)
        and number not in B_EXCEPTIONS
        and alpha(number) > Fraction(7, 44)
    )
    require(not over_7_44_outside_b, over_7_44_outside_b)
    equality_7_44 = tuple(
        number
        for number in score_7_44
        if number >= 2
        and base_free(number)
        and number not in B_EXCEPTIONS
        and alpha(number) == Fraction(7, 44)
    )

    e_set = tuple(number for number in B_EXCEPTIONS if number not in H_CORRIDOR)
    require(all(alpha(number) == Fraction(1, 6) for number in e_set), "E density")
    require(all(alpha(number) > Fraction(1, 6) for number in H_CORRIDOR), "H density")
    loss_table = []
    for prime in LOSS_PRIMES:
        bucket = tuple(number for number in e_set if number % prime == 0)
        require(
            all(2 <= number // prime <= 14 for number in bucket),
            ("loss bucket", prime, bucket),
        )
        loss_table.append((prime, bucket))
    require(
        all(any(number in bucket for _, bucket in loss_table) for number in e_set),
        "loss buckets do not cover E",
    )

    lcm_checks = 0
    for q in range(2, 81):
        mode_rows = tuple(
            (owner_order, mode)
            for owner_order in divisors(q)[1:]
            for mode in sign_representatives(owner_order)
        )
        for size in (1, 2, 3):
            for family in combinations(mode_rows, size):
                modulus = lcm(*(owner_order for owner_order, _ in family))
                if modulus > 30:
                    continue
                q_union = set()
                modulus_union = set()
                for owner_order, mode in family:
                    q_union.update(danger_set(q, (q // owner_order) * mode))
                    modulus_union.update(
                        danger_set(modulus, (modulus // owner_order) * mode)
                    )
                pulled = {
                    sheet for sheet in range(q) if sheet % modulus in modulus_union
                }
                require(q_union == pulled, ("lcm pullback", q, family, modulus))
                require((len(q_union) == q) == (len(modulus_union) == modulus), (
                    "lcm cover equivalence", q, family, modulus
                ))
                lcm_checks += 1

    q21_rank = minimum_literal_cover(21)
    q22_rank = minimum_literal_cover(22)
    require(q21_rank[:3] == (8, (1, 2, 3, 4, 5, 7, 8, 10), 8), q21_rank)
    require(q22_rank[:3] == (7, (1, 2, 3, 5, 7, 9, 11), 7), q22_rank)
    q21_unit_speeds = tuple(speed for speed in range(1, 11) if gcd(speed, 21) == 1)
    q21_unit_union = frozenset().union(*(danger_set(21, speed) for speed in q21_unit_speeds))
    q21_missing = tuple(sorted(set(range(21)) - q21_unit_union))
    require(
        q21_missing == (3, 6, 7, 9, 12, 14, 15, 18),
        ("q21 prime-sheet loss", q21_missing),
    )
    require(not (q21_unit_union & stratum(21, 3)), "q21 order-three stratum")
    require(not (q21_unit_union & stratum(21, 7)), "q21 order-seven stratum")

    q102_value, q102_sheets, q102_owners, q102_active, q102_tight = fractional_game(102)
    require(q102_value == Fraction(4, 25), q102_value)
    require(dict(q102_sheets) == Q102_SHEET_WEIGHTS, q102_sheets)
    require(dict(q102_owners) == Q102_OWNER_WEIGHTS, q102_owners)
    require(q102_active == (2, 3, 17, 102), q102_active)
    require(q102_tight == (6, 34, 51, 102), q102_tight)

    scan_digest = ExactDigest()
    scanned = 0
    equality_moduli = []
    largest_dimension = 0
    support_total = 0
    dual_support_total = 0
    for q in range(15, 20001):
        if not base_free(q):
            continue
        value, sheet_support, owner_support, active_owners, tight_sheets = fractional_game(q)
        require(value <= Fraction(4, 25), ("scan ceiling", q, value))
        require(
            (value == Fraction(4, 25)) == (q % 102 == 0),
            ("scan equality", q, value),
        )
        if value == Fraction(4, 25):
            equality_moduli.append(q)
        largest_dimension = max(largest_dimension, len(divisors(q)) - 1)
        support_total += len(sheet_support)
        dual_support_total += len(owner_support)
        scan_digest.add((
            q,
            value,
            sheet_support,
            owner_support,
            active_owners,
            tight_sheets,
        ))
        scanned += 1
    require(scanned == 15985, scanned)
    require(len(equality_moduli) == 78, len(equality_moduli))
    for q in equality_moduli:
        orders = divisors(q)[1:]
        projected_dual_loads = tuple(
            sum(
                weight * alpha(image_order(q, owner_order, sheet_order))
                for owner_order, weight in Q102_OWNER_WEIGHTS.items()
            )
            for sheet_order in orders
        )
        require(
            min(projected_dual_loads) == Fraction(4, 25),
            ("q102 lower-certificate lift", q, min(projected_dual_loads)),
        )

    semantic.add(("density", density_checks, fibre_checks, pullback_checks))
    semantic.add(("cutoff_4_25", score_4_25))
    semantic.add(("cutoff_7_44", score_7_44))
    semantic.add(("B", B_EXCEPTIONS, tuple((m, alpha(m)) for m in B_EXCEPTIONS)))
    semantic.add(("E", e_set))
    semantic.add(("H", H_CORRIDOR, tuple((m, alpha(m)) for m in H_CORRIDOR)))
    semantic.add(("loss", tuple(loss_table)))
    semantic.add(("lcm", lcm_checks))
    semantic.add(("q21", q21_rank, q21_unit_speeds, q21_missing))
    semantic.add(("q22", q22_rank))
    semantic.add(("q102", q102_value, q102_sheets, q102_owners, q102_active, q102_tight))
    semantic.add(("scan", scan_digest.hexdigest(), tuple(equality_moduli)))
    semantic_digest = semantic.hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (
            semantic_digest, EXPECTED_SEMANTIC_DIGEST
        ))

    print("THM-3408 FIXED-ZERO ADDITIVE-ORDER DUALITY EXACT COMPANION")
    print(f"source_sha256_lf={lf_hash(source)}")
    print("status=PROVED structural laws plus FINITE-EXACT arithmetic and q<=20000 game census")
    print("scope=fixed source centre zero; transverse cyclic-sheet owners; no mobile-centre or LRC14 closure")
    print(
        f"literal_density_checks={density_checks};equal_fibre_checks={fibre_checks};"
        f"mode_pullback_checks={pullback_checks};lcm_family_checks={lcm_checks}"
    )
    print(
        f"score_cutoff_4_25=count_{len(score_4_25)},max_{max(score_4_25)};"
        f"score_cutoff_7_44=count_{len(score_7_44)},max_{max(score_7_44)}"
    )
    print(f"B_over_4_25={B_EXCEPTIONS}")
    print(f"H_strict_over_1_6={tuple((m, alpha(m)) for m in H_CORRIDOR)}")
    print(f"E_equal_1_6={e_set}")
    print(f"outside_B_alpha_max=7/44;equality_moduli={equality_7_44}")
    print(f"prime_loss_buckets={tuple(loss_table)}")
    print(
        f"q21_rank={q21_rank[0]};search_subsets={q21_rank[3]};"
        f"unit_modes={q21_unit_speeds};missing={q21_missing}"
    )
    print(
        f"q22_rank={q22_rank[0]};search_subsets={q22_rank[3]};"
        f"witness={q22_rank[1]};kernel_orders=(11,2)"
    )
    print(
        f"q102_game_value={q102_value};sheet_weights={q102_sheets};"
        f"owner_weights={q102_owners}"
    )
    print(f"q102_active_owners={q102_active};q102_tight_sheet_orders={q102_tight}")
    print(
        f"finite_scan=q15..20000_base_free;moduli={scanned};"
        f"largest_game_dimension={largest_dimension};max_value=4/25;"
        f"six_owner_margin=1/150"
    )
    print(
        f"finite_scan_equality_count={len(equality_moduli)};"
        f"equality_iff_102_divides_q=YES;"
        f"equality_first_last=({equality_moduli[0]},{equality_moduli[-1]})"
    )
    print(
        f"finite_scan_sheet_support_total={support_total};"
        f"owner_support_total={dual_support_total};"
        f"certificate_digest={scan_digest.hexdigest()}"
    )
    print("corridor=any_base_free_fixed_zero_cover_by_at_most_six_owners_uses_order_in_(22,33,44,46,50,102)")
    print(f"semantic_sha256={semantic_digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
