"""Exact first-reset Collatz cylinders and spatial-density controls.

No convergence or exhaustive cycle claim. Normal and -O use identical
explicit checks. The JSON declares every finite universe.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def v2(n):
    check(n > 0, "valuation domain")
    return (n & -n).bit_length() - 1


def step(n):
    check(n > 0 and n % 2 == 1, "odd positive domain")
    z = 3 * n + 1
    k = v2(z)
    return z >> k, k


def follow(n, length):
    nodes = [n]
    exponents = []
    for _ in range(length):
        n, k = step(n)
        nodes.append(n)
        exponents.append(k)
    return nodes, exponents


def block(n):
    r = v2(n + 1) - 1
    u = (n + 1) >> (r + 1)
    k = 1 + v2(3 ** (r + 1) * u - 1)
    y = (3 ** (r + 1) * u - 1) >> (k - 1)
    return r, u, k, y


def kappa(r):
    k = 2
    while 2 ** (r + k) <= 3 ** (r + 1):
        k += 1
    return k


def minimum_reset_for_inequality(n, r):
    numerator = 3 ** (r + 1) * n + 3 ** (r + 1) - 2 ** (r + 1)
    denominator = 2 ** r * n
    k = 0
    while 2 ** k * denominator <= numerator:
        k += 1
    return k


def cylinder(r, k):
    check(r >= 0 and k >= 2, "first-reset word")
    u = (pow(3 ** (r + 1), -1, 2 ** k) * (1 + 2 ** (k - 1))) % 2 ** k
    residue = 2 ** (r + 1) * u - 1
    modulus = 2 ** (r + k + 1)
    carry = 3 ** (r + 1) - 2 ** (r + 1)
    independent = ((2 ** (r + k) - carry) *
                   pow(3 ** (r + 1), -1, modulus)) % modulus
    check(residue == independent, "independent cylinder formulas")
    return residue, modulus


def rational(q):
    return {"numerator": q.numerator, "denominator": q.denominator}


def main():
    source_count = cylinder_case_count = full_residue_checks = 0
    returns = []
    for n in range(1, 2 ** 17, 2):
        r, u, k, y = block(n)
        nodes, ks = follow(n, r + 1)
        check(ks == [1] * r + [k] and nodes[-1] == y, "first-reset formula")
        check(u % 2 == 1 and k >= 2, "run/reset boundary")
        d = 2 ** (r + k) - 3 ** (r + 1)
        check(d != 0, "exponent gap nonzero")
        check(Q(y - n) == 1 - Q(d * u + 1, 2 ** (k - 1)), "carry trichotomy")
        if d < 0:
            check(y > n, "growth side")
        else:
            check(y <= n, "nonincrease side")
        check((y == n) == (d * u == 2 ** (k - 1) - 1), "return boundary")
        check((y < n) == (k >= minimum_reset_for_inequality(n, r)),
              "exact strict reset threshold")
        if y == n:
            returns.append(n)
        source_count += 1

    for r in range(11):
        for k in range(2, 13):
            s, mod = cylinder(r, k)
            for j in range(3):
                n = s + j * mod
                rr, _, kk, _ = block(n)
                check((rr, kk) == (r, k), "cylinder lift")
            if r + k <= 10:
                hits = []
                for n in range(1, mod, 2):
                    rr, _, kk, _ = block(n)
                    if (rr, kk) == (r, k):
                        hits.append(n)
                    full_residue_checks += 1
                check(hits == [s], "complete dyadic cylinder")
            cylinder_case_count += 1

    check(cylinder(2, 5) == (23, 256), "23 cylinder")
    valve_hits = []
    for n in range(5, 2305, 18):
        nodes, ks = follow(n, 3)
        if ks == [1, 1, 5]:
            valve_hits.append(n)
            check(nodes[-1] < n and 128 * nodes[-1] == 27 * n + 19,
                  "valve descent/carry")
    check(valve_hits == [23], "one among 128 odd5mod9 sources")
    for q in range(100):
        n = 23 + 2304 * q
        nodes, ks = follow(n, 3)
        check(ks == [1, 1, 5] and nodes[-1] == 5 + 486 * q,
              "row5 exact progression")
    check(follow(95, 3) == ([95, 143, 215, 323], [1, 1, 1]), "all-growth hostile")
    check(all(follow(n, 3)[1] != [1, 1, 1] for n in range(5, 95, 18)),
          "minimal all-growth row5 hostile")

    family = []
    for a in range(1, 201):
        n = 3 * 2 ** (2 * a + 1) - 1
        r, u, k, y = block(n)
        check(n % 9 == 5 and (r, u) == (2 * a, 3), "family growth run")
        check(k == 4 + v2(a + 1), "LTE reset exponent")
        check(y == (3 ** (2 * a + 2) - 1) // 2 ** (3 + v2(a + 1)),
              "family reset endpoint")
        nodes, ks = follow(n, 2 * a + 1)
        check(ks == [1] * (2 * a) + [k] and nodes[-1] == y,
              "independent family trajectory")
        check((y < n) == (a <= 3), "family first-reset boundary in finite control")
        if a <= 8:
            family.append({"a": a, "n": n, "growth_length": r,
                           "reset_exponent": k, "reset_endpoint": y,
                           "descends": y < n})

    transition_rows = []
    for a in range(9):
        targets = [(3 * a + 1) * pow(2, -k, 9) % 9 for k in range(1, 7)]
        check(set(targets) == {1, 2, 4, 5, 7, 8}, "mod9 full support")
        for k, target in enumerate(targets, 1):
            dyadic = 2 ** (k + 1)
            residue = (pow(3, -1, dyadic) * (2 ** k - 1)) % dyadic
            t = (a - residue) * pow(dyadic, -1, 9) % 9
            n = residue + dyadic * t
            check(n > 0 and n % 9 == a, "CRT source")
            nxt, actual_k = step(n)
            check(actual_k == k and nxt % 9 == target, "CRT target")
        transition_rows.append({"source_mod9": a, "targets_for_k1_to6": targets})

    n_terms = 64
    floors = [(3 ** ell).bit_length() - 1 for ell in range(1, n_terms + 1)]
    lower = sum((Q(1, 2 ** exponent) for exponent in floors), Q(0))
    upper = lower + Q(1, 3 ** n_terms)
    check(floors == sorted(set(floors)), "distinct Beatty digit positions")
    check(all(r + kappa(r) - 1 == floors[r] for r in range(n_terms)),
          "mechanical threshold equals density digit")
    with localcontext() as ctx:
        ctx.prec = 70
        decimal_lower = str(Decimal(lower.numerator) / Decimal(lower.denominator))
        decimal_upper = str(Decimal(upper.numerator) / Decimal(upper.denominator))

    print(json.dumps({
        "status": "FINITE-EXACT replay; density theorem is spatial, not orbital",
        "universes": {"all_odd_sources_below": 2 ** 17, "source_count": source_count,
                      "cylinder_r_range": [0, 10], "cylinder_k_range": [2, 12],
                      "cylinder_pairs": cylinder_case_count,
                      "complete_cylinder_enumeration_r_plus_k_at_most": 10,
                      "complete_cylinder_source_checks": full_residue_checks,
                      "family_a_range": [1, 200]},
        "returns_in_finite_source_universe": returns,
        "valve": {"word": [1, 1, 5], "source_residue": 23, "modulus": 256,
                  "odd_row5_source": "23+2304q", "endpoint": "5+486q",
                  "relative_odd_row_density": "1/128", "period_source_count": 128,
                  "period_hits": valve_hits},
        "family_first_eight": family,
        "mod9_transition_rows": transition_rows,
        "spatial_weights_for_k1_to6": [str(Q(2 ** (6 - k), 63)) for k in range(1, 7)],
        "first_reset_descent_density": {"N": n_terms, "partial_sum": rational(lower),
                                        "strict_upper_bound": rational(upper),
                                        "tail_bound": "3^-64",
                                        "decimal_lower": decimal_lower,
                                        "decimal_upper": decimal_upper},
        "all_checks_passed": True,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
