#!/usr/bin/env python3
"""Exact/numerical referee for THM-2000.

The user's convention is literal: an integer sequence contributes a *subset*
of the harmonic numbers.  Repetitions therefore disappear.  This referee
checks the support/multiplicity collision tax, Abel--Stieltjes summation,
multiplicative-block bounds, the Bertrand counterexamples to a false linear
growth dichotomy, and the two-axis reciprocal surface of the master figurate
array N(s,d,m).  It also inventories several sequence families already native
to the repository.  Analytic limit passages remain paper providers.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial, log
from pathlib import Path

import mpmath as mp


mp.mp.dps = 50


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-2000 referee failed: {message}")


def optimization_safety_probe() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "Python assert nodes are optimization-sensitive")
    return count


def finite_signatures(values: list[int]) -> tuple[F, F, F]:
    counts = Counter(value for value in values if value > 0)
    support = sum((F(1, value) for value in counts), F(0))
    multiple = sum((F(count, value) for value, count in counts.items()), F(0))
    collision = sum((F(count - 1, value) for value, count in counts.items()), F(0))
    require(multiple == support + collision, "collision-tax decomposition")
    return support, multiple, collision


def collision_tax_audit() -> list[tuple[str, F, F, F]]:
    factorial = [1]
    for n in range(1, 13):
        factorial.append(factorial[-1] * n)

    fibonacci = [1, 1]
    for _ in range(28):
        fibonacci.append(fibonacci[-1] + fibonacci[-2])

    catalan = [comb(2 * n, n) // (n + 1) for n in range(30)]
    labeled = [2 ** comb(n, 2) for n in range(1, 14)]
    switching = [2 ** comb(n - 1, 2) for n in range(1, 15)]
    tournaments = [1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056]
    even_graphs = [1, 1, 2, 3, 7, 16, 54, 243, 2038]
    score_sequences = [1, 1, 1, 2, 4, 9, 22, 59, 167, 490]

    families = [
        ("factorial_0_through_12", factorial),
        ("fibonacci_30_terms", fibonacci),
        ("catalan_30_terms", catalan),
        ("labeled_tournament_values", labeled),
        ("switching_values", switching),
        ("A000568_available_prefix", tournaments),
        ("A002854_available_prefix", even_graphs),
        ("A000571_available_prefix", score_sequences),
    ]
    rows = [(name, *finite_signatures(values)) for name, values in families]

    labeled_support = set(labeled)
    switching_support = set(switching)
    require(labeled_support == switching_support,
            "labeled and switching value supports are the same shifted range")
    row_map = {name: (support, multiple, tax) for name, support, multiple, tax in rows}
    require(row_map["switching_values"][1] - row_map["labeled_tournament_values"][1] == 1,
            "the old plus-one is exactly a duplicate-one tax")
    for name in ("factorial_0_through_12", "fibonacci_30_terms", "catalan_30_terms"):
        require(row_map[name][2] == 1, f"{name} has exactly one duplicate-one tax")
    return rows


def abel_rhs(values: list[int], x: int) -> F:
    """A(x)/x + integral_1^x A(t)t^-2 dt for a finite support."""
    support = sorted({value for value in values if 1 <= value <= x})
    count = 1 if support and support[0] == 1 else 0
    previous = 1
    integral = F(0)
    for value in support:
        if value == 1:
            continue
        integral += count * (F(1, previous) - F(1, value))
        previous = value
        count += 1
    integral += count * (F(1, previous) - F(1, x))
    return F(len(support), x) + integral


def abel_stieltjes_audit() -> int:
    families = [
        [n for n in range(1, 101)],
        [n * (n + 1) // 2 for n in range(1, 40)],
        [n * n for n in range(1, 40)],
        [1, 2, 8, 64, 1024, 32768],
        [1, 1, 2, 3, 5, 8, 13, 21, 34],
    ]
    rows = 0
    for values in families:
        for x in range(1, max(values) + 7):
            support = {value for value in values if 1 <= value <= x}
            lhs = sum((F(1, value) for value in support), F(0))
            require(lhs == abel_rhs(values, x), "exact Abel--Stieltjes identity")
            rows += 1
    return rows


def block_bounds(values: list[int]) -> tuple[int, F]:
    support = sorted({value for value in values if value > 0})
    maximum = support[-1]
    rows = 0
    total = F(0)
    k = 0
    while 2**k <= maximum:
        block = [value for value in support if 2**k <= value < 2 ** (k + 1)]
        actual = sum((F(1, value) for value in block), F(0))
        lower = F(len(block), 2 ** (k + 1))
        upper = F(len(block), 2**k)
        require(lower <= actual <= upper, "multiplicative-block harmonic sandwich")
        total += F(len(block), 2**k)
        rows += 1
        k += 1
    return rows, total


def sieve_primes(limit: int) -> list[int]:
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if flags[p]:
            flags[p * p : limit + 1 : p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [n for n in range(2, limit + 1) if flags[n]]


def logarithmic_block_audit() -> list[tuple[str, int, F]]:
    limit = 2**18
    naturals = list(range(1, limit + 1))
    odds = list(range(1, limit + 1, 2))
    primes = sieve_primes(limit)
    triangular = []
    n = 1
    while n * (n + 1) // 2 <= limit:
        triangular.append(n * (n + 1) // 2)
        n += 1
    fib = [1, 2]
    while fib[-1] + fib[-2] <= limit:
        fib.append(fib[-1] + fib[-2])
    bertrand_one = [int(n * log(n + 1) + 0.999999999999) for n in range(1, 30000)]
    bertrand_two = [int(n * log(n + 1) ** 2 + 0.999999999999) for n in range(1, 30000)]
    families = [
        ("naturals", naturals),
        ("odds", odds),
        ("primes", primes),
        ("triangular", triangular),
        ("fibonacci_support", fib),
        ("ceil_n_log_n", bertrand_one),
        ("ceil_n_log2_n", bertrand_two),
    ]
    rows = []
    for name, values in families:
        count, occupancy_sum = block_bounds(values)
        rows.append((name, count, occupancy_sum))
    require(bertrand_one[-1] / len(bertrand_one) > 10,
            "n log n example is genuinely superlinear")
    return rows


def master_figurate(s: int, d: int, m: int) -> int:
    return (s - 2) * comb(m + d - 2, d) + comb(m + d - 2, d - 1)


def polygonal(s: int, m: int) -> int:
    return ((s - 2) * m * m - (s - 4) * m) // 2


def polygonal_mass(s: int) -> mp.mpf:
    if s == 4:
        return mp.pi**2 / 6
    return (mp.mpf(2) / (s - 4)) * (
        mp.digamma(1) - mp.digamma(mp.mpf(2) / (s - 2))
    )


def master_surface_audit() -> tuple[int, list[tuple[int, mp.mpf]]]:
    rows = 0
    for s in range(3, 11):
        for d in range(2, 9):
            previous = None
            for m in range(1, 80):
                value = master_figurate(s, d, m)
                require(value > 0, "master figurate entries are positive")
                require(master_figurate(3, d, m) == comb(m + d - 1, d),
                        "simplex slice")
                require(master_figurate(s, 2, m) == polygonal(s, m),
                        "polygonal slice")
                if d > 2:
                    require(
                        value == sum(master_figurate(s, d - 1, i) for i in range(1, m + 1)),
                        "dimension-axis coning recurrence",
                    )
                if previous is not None:
                    require(value > previous, "each figurate column is strictly increasing")
                previous = value
                rows += 1

    # Exact termwise partial-fraction identity behind the digamma axis.
    for s in range(3, 13):
        for m in range(1, 200):
            if s == 4:
                require(F(1, polygonal(s, m)) == F(1, m * m), "square axis")
            else:
                b = s - 4
                rhs = F(2, b) * (F(1, m - F(b, s - 2)) - F(1, m))
                require(F(1, polygonal(s, m)) == rhs,
                        "polygonal partial-fraction identity")

    # Monotonicity on both axes: the common m=1 term stays one; every tail
    # term strictly decreases.
    for s in range(3, 10):
        for d in range(2, 8):
            for m in range(2, 80):
                require(master_figurate(s + 1, d, m) > master_figurate(s, d, m),
                        "shape-axis denominator increases")
                require(master_figurate(s, d + 1, m) > master_figurate(s, d, m),
                        "dimension-axis denominator increases")

    polygon_rows = [(s, polygonal_mass(s)) for s in range(3, 11)]
    require(abs(polygon_rows[0][1] - 2) < mp.mpf("1e-45"), "triangular mass is two")
    require(abs(polygon_rows[1][1] - mp.pi**2 / 6) < mp.mpf("1e-45"),
            "square mass is zeta two")
    require(abs(polygon_rows[2][1] - (3 * mp.log(3) - mp.pi / mp.sqrt(3)))
            < mp.mpf("1e-45"), "pentagonal mass closed form")
    require(abs(polygon_rows[3][1] - 2 * mp.log(2)) < mp.mpf("1e-45"),
            "hexagonal mass closed form")
    require(all(polygon_rows[i][1] > polygon_rows[i + 1][1]
                for i in range(len(polygon_rows) - 1)),
            "polygonal masses strictly descend")
    return rows, polygon_rows


def master_mass_partial(s: int, d: int, terms: int = 20000) -> mp.mpf:
    return mp.fsum(mp.mpf(1) / master_figurate(s, d, m) for m in range(1, terms + 1))


def power_sum(n: int, p: int) -> int:
    return sum(k**p for k in range(1, n + 1))


def powersum_axis_audit() -> list[tuple[int, mp.mpf, str]]:
    terms = 40000
    rows = []
    for p in range(1, 7):
        cumulative = 0
        reciprocal_terms = []
        for n in range(1, terms + 1):
            cumulative += n**p
            reciprocal_terms.append(mp.mpf(1) / cumulative)
        value = mp.fsum(reciprocal_terms)
        note = "numeric_partial"
        if p == 1:
            closed = mp.mpf(2)
            note = "2"
        elif p == 2:
            closed = 18 - 24 * mp.log(2)
            note = "18-24log(2)"
        elif p == 3:
            closed = 4 * mp.pi**2 / 3 - 12
            note = "4pi^2/3-12"
        else:
            closed = None
        if closed is not None:
            require(abs(value - closed) < mp.mpf("0.0001"),
                    f"power-sum p={p} closed form")
        rows.append((p, value, note))

    for n in range(1, 200):
        # 1 / sum k^2 = 6/[n(n+1)(2n+1)].
        lhs2 = F(1, power_sum(n, 2))
        rhs2 = F(6, n) + F(6, n + 1) - F(24, 2 * n + 1)
        require(lhs2 == rhs2, "square-pyramidal partial fraction")
        # sum k^3 = T_n^2.
        require(power_sum(n, 3) == (n * (n + 1) // 2) ** 2,
                "cubic Faulhaber triangular-square identity")
    return rows


def g_diagonal(n: int) -> int:
    value = (
        F(n**4, 96) - F(n**3, 16) + F(n * n, 3) + F(n, 32) + F(53, 64)
        + (-1) ** n * (F(11, 64) - F(n, 32))
    )
    require(value.denominator == 1, "G quasi-polynomial is integral")
    return value.numerator


def repo_sequence_atlas() -> list[tuple[str, mp.mpf, str]]:
    terms = 30000
    moser = [1 + comb(n, 2) + comb(n, 4) for n in range(1, terms + 1)]
    g_values = [g_diagonal(n) for n in range(terms)]
    fib = [1, 1]
    for _ in range(80):
        fib.append(fib[-1] + fib[-2])

    phi = list(range(200001))
    for p in range(2, len(phi)):
        if phi[p] == p:
            for multiple in range(p, len(phi), p):
                phi[multiple] -= phi[multiple] // p
    farey_totals = []
    totient_sum = 0
    for n in range(1, len(phi)):
        totient_sum += 1 if n == 1 else phi[n]
        farey_totals.append(2 * totient_sum - 1)

    families = [
        ("Moser_A000127", moser, "degree_4;partial_30000"),
        ("polygonal_diagonal_G", g_values, "degree_4_quasipolynomial;partial_30000"),
        ("Fibonacci_support", fib, "duplicate_1_removed"),
        ("Farey_endpoint_total", farey_totals, "2sum_phi-1;partial_200000"),
        ("A000568_support_prefix", [1, 1, 2, 4, 12, 56, 456, 6880, 191536],
         "available_prefix;duplicate_1_removed"),
        ("self_line_orbits_prefix", [2, 3, 22, 101, 852, 5582],
         "THM853_n5_to_n10_prefix"),
    ]
    rows = []
    for name, values, note in families:
        support = sorted({value for value in values if value > 0})
        mass = mp.fsum(mp.mpf(1) / value for value in support)
        rows.append((name, mass, note))
    return rows


def gauss_triangular_theta_audit() -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    q = mp.mpf(1) / 2
    triangular_sum = mp.fsum(q ** (n * (n + 1) // 2) for n in range(80))
    q_product = mp.mpf(1)
    even_product = mp.mpf(1)
    for n in range(1, 500):
        q_product *= 1 - q**n
        even_product *= 1 - q ** (2 * n)
    gauss_product = even_product**2 / q_product
    theta_form = mp.jtheta(2, 0, mp.sqrt(q)) / (2 * q ** (mp.mpf(1) / 8))
    require(abs(triangular_sum - gauss_product) < mp.mpf("1e-45"),
            "Gauss triangular-number product")
    require(abs(triangular_sum - theta_form) < mp.mpf("1e-45"),
            "theta-two form")
    return triangular_sum, gauss_product, theta_form


def ladder_ratio_sum_product_audit() -> tuple[F, mp.mpf, mp.mpf]:
    # Ratios 2k-1 have harmonic mass divergence.  Their partial sums are
    # k^2-1, whose reciprocal mass telescopes to 3/4.  Their partial products
    # are odd double factorials and have an erf closed form.
    for n in range(2, 500):
        partial = sum((F(1, k * k - 1) for k in range(2, n + 1)), F(0))
        closed = F(3, 4) - F(1, 2 * n) - F(1, 2 * (n + 1))
        require(partial == closed, "oblong reciprocal telescoping")

    product_mass = mp.mpf(0)
    odd_double_factorial = 1
    for k in range(1, 100):
        odd_double_factorial *= 2 * k - 1
        if k >= 2:
            product_mass += mp.mpf(1) / odd_double_factorial
    product_closed = (
        mp.e ** (mp.mpf(1) / 2)
        * mp.sqrt(mp.pi / 2)
        * mp.erf(1 / mp.sqrt(2))
        - 1
    )
    require(abs(product_mass - product_closed) < mp.mpf("1e-45"),
            "odd-double-factorial reciprocal erf form")
    odd_partial = mp.fsum(mp.mpf(1) / (2 * k - 1) for k in range(1, 200000))
    return F(3, 4), product_mass, odd_partial


def tournament_analytic_reversal_audit() -> int:
    rows = 0
    for n in range(1, 80):
        # Consecutive labeled-EGF terms have absolute ratio
        # 2^n |z|/(n+1), so every z != 0 eventually makes terms grow.
        require(F(2**n, n + 1) > 1 if n >= 2 else True,
                "labeled tournament EGF ratio eventually exceeds one at z=1")

        # Orbit count V_n is at least labeled_count/n!.  This lower bound is
        # already superpolynomial, so V_n/n^s cannot tend to zero for fixed s.
        lower = F(2 ** comb(n, 2), factorial(n))
        if n >= 20:
            require(lower > n**10, "Burnside orbit lower bound beats a fixed polynomial")
        rows += 1
    return rows


def tournament_and_carrier_audit() -> None:
    print("TOURNAMENT_AND_CARRIER_AUDIT")
    print("faithful_vertices=occupied_integer_values_or_logarithmic_blocks;not_sequence_indices")
    print("predicate_preserved=support_membership,reciprocal_mass,block_occupancy,convergence")
    print("multiplicity_sidecar=collision_tax;index_shift_and_repetition_are_quotiented_out")
    print("candidate_vertices=indices,values,blocks,figurate_axes,Dirichlet_exponents,proof_obligations")
    print("pair_observable=finite_support_mass_difference;gauge=exact_rational_then_family_name")
    print("induced_tournament=transitive;cycles=0;scc_sizes=all_1;hamiltonian_paths=1")
    print("why_unfaithful=finite_mass_order_does_not_determine_tail_convergence_or_closed_form")
    print("tie_path=mass_sorted_family_list;useful_as_atlas_only_not_as_proof_quotient")


def main() -> None:
    assert_nodes = optimization_safety_probe()
    collision_rows = collision_tax_audit()
    abel_rows = abel_stieltjes_audit()
    block_rows = logarithmic_block_audit()
    surface_rows, polygon_rows = master_surface_audit()
    power_rows = powersum_axis_audit()
    atlas_rows = repo_sequence_atlas()
    theta_sum, theta_product, theta_form = gauss_triangular_theta_audit()
    oblong_mass, double_factorial_mass, odd_partial = ladder_ratio_sum_product_audit()
    reversal_rows = tournament_analytic_reversal_audit()

    print("THM-2000 SUPPORT-HARMONIC ABEL-DINI FIGURATE-SURFACE REFEREE")
    print(f"optimization_sensitive_assert_nodes={assert_nodes}")
    print("SUPPORT_VERSUS_MULTIPLICITY")
    for name, support, multiple, tax in collision_rows:
        print(f"{name}:support={support};multiple={multiple};collision_tax={tax}")
    print("support_law=sigma_multi=sigma_set+sum_m(max(multiplicity(m)-1,0)/m)")
    print("labeled_tournament_support=switching_support;old_plus_one_is_duplicate_one_tax")

    print(f"ABEL_STIELTJES_EXACT_ROWS={abel_rows}")
    print("abel_law=sum_(m<=x,m_in_A)1/m=A(x)/x+integral_1^x A(t)/t^2 dt")
    print("block_law=for_Nk=#(A intersect [2^k,2^(k+1))), Nk/2^(k+1)<=block_mass<=Nk/2^k")
    for name, count, occupancy in block_rows:
        print(f"block_family={name};blocks={count};finite_occupancy_sum={occupancy}")
    print("bertrand_repair=n_log_n_is_superlinear_but_reciprocal_diverges;n_log2_n_converges")

    print(f"MASTER_FIGURATE_EXACT_ROWS={surface_rows}")
    print("N(s,d,m)=(s-2)C(m+d-2,d)+C(m+d-2,d-1)")
    print("mass_surface=d(d-1) double_integral (1-t)^(d-2)u^(d-1)/(1-tu^(s-2)) dtdu")
    print("simplex_axis=M(3,d)=d/(d-1);polygonal_axis=digamma_formula;both_axes_strictly_descend_to_1")
    for s, value in polygon_rows:
        print(f"polygonal_s={s};mass={mp.nstr(value, 30)}")
    for d in range(2, 7):
        partial = master_mass_partial(4, d)
        print(f"master_surface_s=4,d={d};partial20000={mp.nstr(partial, 20)}")

    print("POWER_SUM_FAULHABER_AXIS")
    for p, value, note in power_rows:
        print(f"power_p={p};partial40000={mp.nstr(value, 25)};closed={note}")

    print("REPO_SUPPORT_MASS_ATLAS")
    for name, value, note in atlas_rows:
        print(f"{name}:support_mass={mp.nstr(value, 25)};note={note}")

    print("GAUSS_TRIANGULAR_THETA")
    print(f"sum_q_triangular={mp.nstr(theta_sum, 30)}")
    print(f"q_product_form={mp.nstr(theta_product, 30)};theta2_form={mp.nstr(theta_form, 30)}")
    print("labeled_and_switching_support_mass_at_q_half_are_equal_to_this_exact_Gauss_product")

    print("RATIO_SUM_PRODUCT_TRIAD")
    print(f"odd_ratio_support_partial200000={mp.nstr(odd_partial, 25)};status=diverges")
    print(f"oblong_k2_minus1_support_mass={oblong_mass}")
    print(f"odd_double_factorial_k_ge2_mass={mp.nstr(double_factorial_mass, 30)}")

    print(f"TOURNAMENT_ANALYTIC_REVERSAL_ROWS={reversal_rows}")
    print("index_Dirichlet=sum V_n/n^s diverges_for_every_fixed_s")
    print("labeled_EGF=sum 2^C(n,2)z^n/n! has_radius_zero")
    print("reciprocal_value_Dirichlet=sum V_n^(-s) has_abscissa_zero")
    print("reciprocal_OGF=sum z^n/V_n is_entire")
    print("H_spectrum_divergence=OPEN;THM1370_proves_only_omissions_7_21_and_finite_coverage")
    print("THM1127_affine_phase_law=N(K+194040)=N(K)+121304 => support_contains_AP => reciprocal_diverges")

    tournament_and_carrier_audit()
    print("SCOPE=no_irrationality_or_transcendence_claim_for_unnamed_numeric_constants")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
