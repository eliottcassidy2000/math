#!/usr/bin/env python3
"""Independent exact verifier for the reserved THM-4382 candidate.

The script deliberately does not import or inspect the producer artifact in
.scratch/lrc_nonresonant_next_20260903.  It checks the signed-(1,4,1)
classification, the owner/determinant gate, two independent exact measure
implementations, the complete pq < 289 table, the global ceiling, and finite
hostiles including disjointness from the signed-(1,2,1) sector.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd


F = Fraction
R = F(3, 14)
ONE_FOURTEENTH = F(1, 14)
SHARP = F(12, 301)
CLASSIFICATION_LIMIT = 199

check_count = 0


def check(condition: bool, message: str) -> None:
    global check_count
    check_count += 1
    if not condition:
        raise RuntimeError(message)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def mod_one(x: F) -> F:
    return x - floor_fraction(x)


def nearest_integer(x: F) -> int:
    """Nearest integer for x >= 0; half-integer choice is immaterial here."""
    return floor_fraction(x + F(1, 2))


def circle_distance(x: F) -> F:
    u = mod_one(x)
    return min(u, 1 - u)


def is_odd_three_unit(w: int) -> bool:
    return w > 0 and w % 2 == 1 and w % 3 != 0


def branch_from_congruence(p: int, q: int) -> str | None:
    """Return the unique primitive endpoint branch, before distinctness."""
    branches: list[str] = []
    if (q - 19 * p) % 24 == 0:
        branches.append("mean")
    if (q - 5 * p) % 24 == 0:
        branches.append("difference")
    check(len(branches) <= 1, f"branch overlap at endpoints {(p, q)}")
    return branches[0] if branches else None


def middle_speed(p: int, q: int, branch: str) -> int:
    numerator = p + q if branch == "mean" else q - p
    check(numerator > 0 and numerator % 4 == 0,
          f"nonintegral middle at {(p, q, branch)}")
    return numerator // 4


def branch_sign(branch: str) -> int:
    # Uniform identity 4b = q + s p.
    return 1 if branch == "mean" else -1


def signed_presentations(triple: tuple[int, int, int], coefficient: int):
    """All choices of the coefficient-carrying speed and positive branch."""
    found = []
    for i, middle in enumerate(triple):
        endpoints = [triple[j] for j in range(3) if j != i]
        p, q = sorted(endpoints)
        if coefficient * middle == p + q:
            found.append((middle, p, q, "mean"))
        if coefficient * middle == q - p:
            found.append((middle, p, q, "difference"))
    return found


def positive_part(x: F) -> F:
    return max(F(0), x)


def thresholds(p: int, q: int) -> tuple[F, F]:
    return F(3 * (q - p), 56), F(3 * (q + p), 56)


def component_length(p: int, q: int, k: int) -> F:
    """Length for determinant K=4k, from direct interval intersection."""
    separation = F(4 * abs(k), p * q)
    return max(F(0), min(F(2) * R / q, R / p + R / q - separation))


def f_length(p: int, q: int, t: int | F) -> F:
    A, B = thresholds(p, q)
    t = F(t)
    return F(4, p * q) * (positive_part(B - t) - positive_part(A - t))


def measure_series(p: int, q: int) -> F:
    _, B = thresholds(p, q)
    total = F(0)
    k = 1
    while F(k) < B:
        if k % 3 != 0:
            direct = component_length(p, q, k)
            folded = f_length(p, q, k)
            check(direct == folded,
                  f"component/f mismatch at {(p, q, k, direct, folded)}")
            total += folded
        k += 1
    return 2 * total


def E(t: F) -> F:
    rho = t - 3 * floor_fraction(t / 3)
    check(F(0) <= rho < F(3), f"bad E residue {rho}")
    if rho <= 1:
        value = -rho * rho / 3
    elif rho <= 2:
        value = rho - 1 - rho * rho / 3
    else:
        value = 2 * rho - 3 - rho * rho / 3
    check(F(-1, 3) <= value <= 0, f"E out of range at {t}: {value}")
    return value


def measure_quadrature(p: int, q: int) -> F:
    A, B = thresholds(p, q)
    return F(3, 98) + F(8, p * q) * (E(B) - E(A))


def measure_ceiling(p: int, q: int) -> F:
    return F(3, 98) + F(8, 3 * p * q)


def danger_x(w: int, sheet: int, x: F) -> bool:
    return circle_distance(w * (x + F(sheet, 3))) < ONE_FOURTEENTH


def x_boundaries(triple: tuple[int, int, int]) -> list[F]:
    points = {F(0), F(1)}
    for w in triple:
        radius = F(1, 14 * w)
        for sheet in range(3):
            shift = F(sheet, 3)
            for n in range(w):
                center = F(n, w) - shift
                points.add(mod_one(center - radius))
                points.add(mod_one(center + radius))
    return sorted(points)


def measure_direct_x(triple: tuple[int, int, int]) -> F:
    """Definition-level circle-cell union over all six sheet assignments."""
    points = x_boundaries(triple)
    total = F(0)
    for left, right in zip(points, points[1:]):
        if left == right:
            continue
        x = (left + right) / 2
        active = 0
        for sheet_assignment in permutations((0, 1, 2)):
            if all(danger_x(w, sheet, x)
                   for w, sheet in zip(triple, sheet_assignment)):
                active += 1
        check(active <= 1,
              f"sheet-assignment overlap at {triple}, x={x}: {active}")
        if active:
            total += right - left
    return total


def y_boundaries(triple: tuple[int, int, int]) -> list[F]:
    points = {F(0), F(1)}
    for w in triple:
        radius = R / w
        for n in range(w):
            center = F(n, w)
            points.add(mod_one(center - radius))
            points.add(mod_one(center + radius))
    return sorted(points)


def audit_owner_cells(p: int, b: int, q: int, branch: str) -> tuple[F, int]:
    """Check the full owner gate on every exact y-cell for one triple."""
    s = branch_sign(branch)
    triple = (p, b, q)
    points = y_boundaries(triple)
    total = F(0)
    cells = 0
    for left, right in zip(points, points[1:]):
        if left == right:
            continue
        y = (left + right) / 2
        nearest = {w: nearest_integer(w * y) for w in triple}
        eligible = {w: abs(w * y - nearest[w]) < R for w in triple}
        owners = {
            w: (-pow(w, -1, 3) * nearest[w]) % 3
            for w in triple
        }
        member = all(eligible.values()) and len(set(owners.values())) == 3

        n_p, n_b, n_q = nearest[p], nearest[b], nearest[q]
        delta = n_q + s * n_p - 4 * n_b
        K = q * n_p - p * n_q
        endpoint_gate = eligible[p] and eligible[q] and K % 4 == 0
        if endpoint_gate:
            k = K // 4
            endpoint_gate = k % 3 != 0

        if all(eligible.values()):
            # |delta| < 6r = 9/7 forces this three-value range.
            check(delta in (-1, 0, 1),
                  f"analytic defect range failed at {(triple, y, delta)}")
        if member:
            check(delta == 0,
                  f"distinct owners did not force delta=0 at {(triple, y)}")
            check(K % 4 == 0,
                  f"owner member has K not divisible by four at {(triple, y)}")
            check((K // 4) % 3 != 0,
                  f"owner member has 3|k at {(triple, y)}")
        check(member == endpoint_gate,
              f"owner/determinant equivalence failed at {(triple, y)}")
        if member:
            total += right - left
        cells += 1
    return total, cells


def admissible_endpoint_pairs(product_limit: int) -> list[tuple[int, int, int, str]]:
    rows = []
    for p in range(1, product_limit):
        if not is_odd_three_unit(p):
            continue
        for q in range(p + 1, product_limit):
            if p * q >= product_limit:
                break
            if not is_odd_three_unit(q) or gcd(p, q) != 1:
                continue
            branch = branch_from_congruence(p, q)
            if branch is None:
                continue
            b = middle_speed(p, q, branch)
            if not is_odd_three_unit(b) or len({p, b, q}) != 3:
                continue
            rows.append((p, b, q, branch))
    return rows


def main() -> None:
    semantic_records: list[str] = []

    # Exhaustive classification and overlap hostile on primitive odd 3-units.
    pool = [w for w in range(1, CLASSIFICATION_LIMIT + 1)
            if is_odd_three_unit(w)]
    primitive_triples = 0
    count_141 = 0
    count_121 = 0
    overlap = 0
    for triple in combinations(pool, 3):
        if gcd(gcd(triple[0], triple[1]), triple[2]) != 1:
            continue
        primitive_triples += 1
        rel141 = signed_presentations(triple, 4)
        rel121 = signed_presentations(triple, 2)
        check(len(rel141) <= 1,
              f"multiple signed-(1,4,1) presentations at {triple}: {rel141}")
        if rel141:
            count_141 += 1
            b, p, q, branch = rel141[0]
            check(gcd(p, q) == 1,
                  f"primitive triple did not give coprime endpoints: {triple}")
            check(branch_from_congruence(p, q) == branch,
                  f"congruence classification failed at {triple}: {rel141}")
            check(middle_speed(p, q, branch) == b,
                  f"middle reconstruction failed at {triple}: {rel141}")
        if rel121:
            count_121 += 1
        if rel141 and rel121:
            overlap += 1
    check(overlap == 0,
          f"signed-(1,4,1)/(1,2,1) overlap found through {CLASSIFICATION_LIMIT}")

    # Check both endpoint congruences as iff statements on a wider pair box.
    pair_congruence_checks = 0
    for p in range(1, 300):
        if not is_odd_three_unit(p):
            continue
        for q in range(p + 1, 300):
            if not is_odd_three_unit(q):
                continue
            mean_integral_unit = (p + q) % 4 == 0 and is_odd_three_unit((p + q) // 4)
            difference_integral_unit = ((q - p) % 4 == 0
                                        and is_odd_three_unit((q - p) // 4))
            check(mean_integral_unit == ((q - 19 * p) % 24 == 0),
                  f"mean congruence iff failed at {(p, q)}")
            check(difference_integral_unit == ((q - 5 * p) % 24 == 0),
                  f"difference congruence iff failed at {(p, q)}")
            pair_congruence_checks += 2

    small_rows = admissible_endpoint_pairs(289)
    check(small_rows, "small endpoint table unexpectedly empty")

    # Definition-level and owner-level exact checks on every small case.
    measured_rows = []
    owner_cell_count = 0
    direct_case_count = 0
    for p, b, q, branch in small_rows:
        series = measure_series(p, q)
        quadrature = measure_quadrature(p, q)
        owner_measure, cells = audit_owner_cells(p, b, q, branch)
        direct = measure_direct_x((p, b, q))
        check(series == quadrature,
              f"series/quadrature mismatch at {(p, b, q)}")
        check(series == owner_measure,
              f"series/owner-cell mismatch at {(p, b, q)}")
        check(series == direct,
              f"series/definition-cell mismatch at {(p, b, q)}")
        check(series <= measure_ceiling(p, q),
              f"period-three ceiling failed at {(p, b, q)}")
        measured_rows.append((p, b, q, branch, series))
        owner_cell_count += cells
        direct_case_count += 1

    # Additional larger controls chosen without reference to the producer.
    larger_controls = []
    for p in range(1, 40):
        if not is_odd_three_unit(p):
            continue
        for q in range(p + 1, 340):
            if not is_odd_three_unit(q) or gcd(p, q) != 1 or p * q < 289:
                continue
            branch = branch_from_congruence(p, q)
            if branch is None:
                continue
            b = middle_speed(p, q, branch)
            if is_odd_three_unit(b) and len({p, b, q}) == 3:
                larger_controls.append((p, b, q, branch))
            if len(larger_controls) == 8:
                break
        if len(larger_controls) == 8:
            break
    check(len(larger_controls) == 8, "failed to select eight larger controls")
    for p, b, q, branch in larger_controls:
        series = measure_series(p, q)
        quadrature = measure_quadrature(p, q)
        owner_measure, cells = audit_owner_cells(p, b, q, branch)
        direct = measure_direct_x((p, b, q))
        check(series == quadrature == owner_measure == direct,
              f"larger exact paths disagree at {(p, b, q)}")
        owner_cell_count += cells
        direct_case_count += 1

    maximum = max(row[4] for row in measured_rows)
    maximizers = [row for row in measured_rows if row[4] == maximum]
    check(maximum == SHARP,
          f"small maximum {maximum} differs from proposed sharp value {SHARP}")
    check(maximizers == [(1, 11, 43, "mean", SHARP)],
          f"small maximizer not unique/expected: {maximizers}")

    large_at_threshold = F(3, 98) + F(8, 3 * 289)
    check(large_at_threshold < SHARP,
          f"pq>=289 ceiling does not separate: {large_at_threshold} >= {SHARP}")

    semantic_records.extend([
        f"classification_limit={CLASSIFICATION_LIMIT}",
        f"primitive_triples={primitive_triples}",
        f"signed141={count_141}",
        f"signed121={count_121}",
        f"overlap={overlap}",
        f"pair_congruence_checks={pair_congruence_checks}",
        f"large_threshold={large_at_threshold}",
    ])
    for p, b, q, branch, measure in measured_rows:
        semantic_records.append(f"{p},{b},{q},{branch},{measure}")
    for p, b, q, branch in larger_controls:
        semantic_records.append(f"control:{p},{b},{q},{branch},{measure_series(p,q)}")
    semantic_hash = sha256(("\n".join(semantic_records) + "\n").encode("ascii")).hexdigest()

    print("THM-4382 clean-room exact verifier")
    print("scope=primitive distinct positive odd 3-units with a signed-(1,4,1) relation")
    print(f"classification_limit={CLASSIFICATION_LIMIT}")
    print(f"primitive_unit_triples={primitive_triples}")
    print(f"signed_141_triples={count_141}")
    print(f"signed_121_triples={count_121}")
    print(f"signed_sector_overlap={overlap}")
    print(f"pair_congruence_checks={pair_congruence_checks}")
    print(f"owner_cells_checked={owner_cell_count}")
    print(f"definition_measure_cases={direct_case_count}")
    print("large_product_regime=pq>=289")
    print(f"ceiling_at_pq_289={large_at_threshold}")
    print(f"sharp_target={SHARP}")
    print(f"small_product_regime=pq<289 cases={len(measured_rows)}")
    print("p b q branch measure ceiling")
    for p, b, q, branch, measure in measured_rows:
        print(f"{p} {b} {q} {branch} {measure} {measure_ceiling(p,q)}")
    print(f"unique_maximizer={(1, 11, 43)} measure={maximum}")
    print("larger_controls=" + ",".join(
        f"({p},{b},{q},{branch})" for p, b, q, branch in larger_controls
    ))
    print(f"semantic_sha256={semantic_hash}")
    print(f"checks={check_count}")
    print("PASS")


if __name__ == "__main__":
    main()
