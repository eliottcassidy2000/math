#!/usr/bin/env python3
"""THM-3181 exact audit: tournament half-grid reciprocity and join tails.

The primary engine uses direct permutation enumeration.  It checks parity,
the positive binomial coordinate, rational generating function, Worpitzky
integration, half-grid inversion, SCC order loss, and repeated-join tails.
"""

from itertools import permutations
from math import comb, factorial


def require(condition, payload):
    if not condition:
        raise AssertionError(payload)


def tournament_from_code(n, code):
    beats = [[False] * n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> bit) & 1:
                beats[i][j] = True
            else:
                beats[j][i] = True
            bit += 1
    return tuple(tuple(row) for row in beats)


def transitive(n):
    return tuple(tuple(i < j for j in range(n)) for i in range(n))


def directed_triangle():
    return (
        (False, True, False),
        (False, False, True),
        (True, False, False),
    )


def order_join(left, right):
    a, b = len(left), len(right)
    out = [[False] * (a + b) for _ in range(a + b)]
    for i in range(a):
        for j in range(a):
            out[i][j] = left[i][j]
    for i in range(b):
        for j in range(b):
            out[a + i][a + j] = right[i][j]
    for i in range(a):
        for j in range(b):
            out[i][a + j] = True
    return tuple(tuple(row) for row in out)


def repeated_join(base, r):
    require(r >= 1, r)
    out = base
    for _ in range(r - 1):
        out = order_join(out, base)
    return out


def backward_distribution(tournament):
    n = len(tournament)
    counts = [0] * n
    for order in permutations(range(n)):
        backward = sum(
            tournament[order[i + 1]][order[i]] for i in range(n - 1)
        )
        counts[backward] += 1
    return tuple(counts)


def ordered_profile_from_backward(backward):
    n = len(backward)
    profile = [0] * (n + 1)
    for d in range(1, n + 1):
        profile[d] = sum(
            comb(n - 1 - b, d - 1 - b) * backward[b]
            for b in range(d)
        )
    return tuple(profile)


def falling(t, d):
    value = 1
    for j in range(d):
        value *= t - j
    return value


def q_value_from_ordered_profile(ordered_profile, t):
    return sum(
        (ordered_profile[d] // factorial(d)) * falling(t, d)
        for d in range(1, len(ordered_profile))
    )


def q_value_from_backward(backward, m):
    n = len(backward)
    return sum(backward[b] * comb(m + b, n) for b in range(n))


def half_grid_recover(n, q_values):
    """Recover the palindromic A-vector from Q(0),...,Q(ceil(n/2))."""
    half = (n + 1) // 2
    recovered = [None] * n
    for k in range(1, half + 1):
        value = sum(
            (-1) ** (k - m) * comb(n + 1, k - m) * q_values[m]
            for m in range(k + 1)
        )
        recovered[k - 1] = value
        recovered[n - k] = value
    require(all(value is not None for value in recovered), (n, recovered))
    return tuple(recovered)


def worpitzky_value(backward, m):
    """THM-087 degree-(n-1) Worpitzky value a_m."""
    n = len(backward)
    return sum(
        backward[b] * comb(m + n - 1 - b, n - 1)
        for b in range(n)
    )


def gf_coefficient(backward, m):
    """Coefficient of x^m in x B_T(x)/(1-x)^(n+1)."""
    n = len(backward)
    return sum(
        backward[b] * comb(n + m - b - 1, n)
        for b in range(min(n, m))
    )


def profile_q_values(tournament, upper):
    backward = backward_distribution(tournament)
    ordered = ordered_profile_from_backward(backward)
    return tuple(
        q_value_from_ordered_profile(ordered, m) for m in range(upper + 1)
    )


def fixed_tail_formula(base, r, k):
    order = len(base)
    q_values = profile_q_values(base, k)
    return sum(
        (-1) ** (k - m)
        * comb(order * r + 1, k - m)
        * q_values[m] ** r
        for m in range(1, k + 1)
    )


def apply_shift_minus(sequence, root, multiplicity):
    values = list(sequence)
    for _ in range(multiplicity):
        values = [
            values[i + 1] - root * values[i] for i in range(len(values) - 1)
        ]
    return values


def score_sequence(tournament):
    return tuple(sorted((sum(row) for row in tournament), reverse=True))


def main():
    labelled = 0
    parity_checks = 0
    binomial_checks = 0
    gf_checks = 0
    worpitzky_checks = 0
    half_recoveries = 0

    # Explicit universe: every labelled tournament of orders 1 through 5.
    for n in range(1, 6):
        for code in range(1 << comb(n, 2)):
            tournament = tournament_from_code(n, code)
            backward = backward_distribution(tournament)
            require(sum(backward) == factorial(n), (n, code, backward))
            require(
                backward == tuple(reversed(backward)),
                ("palindrome", n, code),
            )
            ordered = ordered_profile_from_backward(backward)
            require(
                all(ordered[d] % factorial(d) == 0 for d in range(1, n + 1)),
                ("profile integrality", n, code, ordered),
            )

            q_values = [
                q_value_from_ordered_profile(ordered, m) for m in range(9)
            ]
            require(q_values[0] == 0, ("zero", n, code, q_values[0]))
            for m in range(1, 9):
                require(
                    q_value_from_ordered_profile(ordered, -m)
                    == (-1) ** n * q_values[m],
                    ("Q parity", n, code, m),
                )
                require(
                    q_values[m] == q_value_from_backward(backward, m),
                    ("positive binomial", n, code, m),
                )
                require(
                    q_values[m] == gf_coefficient(backward, m),
                    ("rational GF", n, code, m),
                )
                require(
                    q_values[m] - q_values[m - 1]
                    == worpitzky_value(backward, m - 1),
                    ("Worpitzky discrete integral", n, code, m),
                )
                parity_checks += 1
                binomial_checks += 1
                gf_checks += 1
                worpitzky_checks += 1

            half = (n + 1) // 2
            require(
                half_grid_recover(n, q_values[: half + 1]) == backward,
                ("half-grid recovery", n, code),
            )
            half_recoveries += 1
            labelled += 1

    require(labelled == 1099, labelled)

    # SCC/order-loss hostile and positive product control.
    k1 = transitive(1)
    c3 = directed_triangle()
    left = order_join(k1, c3)
    right = order_join(c3, k1)
    left_backward = backward_distribution(left)
    right_backward = backward_distribution(right)
    require(
        score_sequence(left) != score_sequence(right),
        (
            "hostiles accidentally isomorphic by score",
            score_sequence(left),
            score_sequence(right),
        ),
    )
    require(
        left_backward == right_backward,
        ("SCC order should be lost", left_backward, right_backward),
    )

    c3c3 = order_join(c3, c3)
    direct = backward_distribution(c3c3)
    n = len(c3c3)
    half = (n + 1) // 2
    q_c3 = profile_q_values(c3, half)
    product_q = tuple(value * value for value in q_c3)
    require(
        half_grid_recover(n, product_q) == direct,
        ("SCC product plus half inversion", direct, product_q),
    )

    # Repeated-join exact tails, checked against direct permutations to order 8.
    repeated_direct_checks = 0
    for base in (
        k1,
        c3,
        tournament_from_code(2, 1),
        tournament_from_code(3, 1),
    ):
        base_order = len(base)
        for r in range(1, 1 + 8 // base_order):
            joined = repeated_join(base, r)
            direct_backward = backward_distribution(joined)
            upper = min((len(joined) + 1) // 2, 4)
            for k in range(1, upper + 1):
                formula = fixed_tail_formula(base, r, k)
                require(
                    formula
                    == direct_backward[k - 1]
                    == direct_backward[len(joined) - k],
                    (
                        "repeated tail",
                        base_order,
                        r,
                        k,
                        formula,
                        direct_backward,
                    ),
                )
                repeated_direct_checks += 1

    # The C3, k=3 predicted annihilator has roots Q_C3(1..3).
    base = c3
    k = 3
    roots = profile_q_values(base, k)[1:]
    require(
        all(roots[i] < roots[i + 1] for i in range(len(roots) - 1)),
        roots,
    )
    recurrence_order = k * (k + 1) // 2
    sequence = [
        fixed_tail_formula(base, r, k)
        for r in range(1, 1 + recurrence_order + 12)
    ]
    residual = sequence
    for m, root in enumerate(roots, start=1):
        residual = apply_shift_minus(residual, root, k - m + 1)
    require(
        residual and all(value == 0 for value in residual),
        ("recurrence annihilator", roots, residual),
    )

    # Transitive boundary: K1 repeated joins recover ordinary Eulerian numbers.
    eulerian_tail = [fixed_tail_formula(k1, r, 3) for r in range(3, 11)]
    expected_eulerian_tail = [
        backward_distribution(transitive(r))[2] for r in range(3, 11)
    ]
    require(
        eulerian_tail == expected_eulerian_tail,
        ("Eulerian boundary", eulerian_tail, expected_eulerian_tail),
    )

    print("TOURNAMENT HALF-GRID RECIPROCITY -- THM-3181 EXACT COMPANION")
    print(
        "status=PROVED+VERIFIED-EXACT;"
        "no_canonical_runtime_imports=true;direct-permutation-engine"
    )
    print("universe=all labelled tournaments n=1..5;count=1099")
    print(
        f"Q_parity_checks={parity_checks};"
        "Q(-m)=(-1)^n Q(m);m=1..8;PASS"
    )
    print(
        f"positive_binomial_checks={binomial_checks};"
        "Q(m)=sum_b A(b)binom(m+b,n);PASS"
    )
    print(
        f"rational_GF_checks={gf_checks};"
        "sum_m Q(m)x^m=x B_T(x)/(1-x)^(n+1);PASS"
    )
    print(
        f"worpitzky_integration_checks={worpitzky_checks};"
        "Delta Q(m)=a_(m-1);PASS"
    )
    print(
        f"half_grid_recoveries={half_recoveries};"
        "ceil(n/2)_values_recover_full_A;PASS"
    )
    print(
        "SCC_order_hostile=K1_join_C3_vs_C3_join_K1;"
        "score_sequences_differ;A_distributions_equal;PASS"
    )
    print(
        "SCC_positive_control=C3_join_C3;"
        "Q_jet_product_plus_half_inversion=direct_A;PASS"
    )
    print(
        f"repeated_join_direct_tail_checks={repeated_direct_checks};"
        "total_order<=8;k<=4;PASS"
    )
    print(
        f"C3_k3_roots={roots};annihilator_multiplicities=(3,2,1);"
        f"order={recurrence_order};PASS"
    )
    print(
        f"K1_k3_Eulerian_boundary_r3_to_r10={eulerian_tail};PASS"
    )


if __name__ == "__main__":
    main()
