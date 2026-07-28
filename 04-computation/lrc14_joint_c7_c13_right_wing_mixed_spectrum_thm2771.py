#!/usr/bin/env python3
"""Exact THM-2771 physical-present-clock x target-label right-wing probe.

This reconstructs the corrected THM-2751 carrier for
all seven *physical present* C7 clock sheets e and all thirteen target labels
t.  The delayed-clock coefficient is still evaluated exactly on its seven
marked Q_(3,{1,2}) prefixes; its constant-tail amplitude is the entry of the
joint e x t table.  No semantic-arm or endpoint attachment is inferred from
Fourier survival.
"""

from fractions import Fraction
from functools import reduce
from math import gcd
from math import comb
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751 as wing


P = 13
Q = 7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank_q(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    columns = len(work[0])
    for column in range(columns):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [entry / value for entry in work[rank]]
        for i in range(len(work)):
            if i == rank or not work[i][column]:
                continue
            value = work[i][column]
            work[i] = [a - value * b for a, b in zip(work[i], work[rank])]
        rank += 1
    return rank


def rank_mod(matrix, prime):
    work = [[value % prime for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(inverse * value) % prime for value in work[rank]]
        for i in range(len(work)):
            if i == rank or not work[i][column]:
                continue
            value = work[i][column]
            work[i] = [
                (a - value * b) % prime
                for a, b in zip(work[i], work[rank])
            ]
        rank += 1
    return rank


def crt_index(e, t):
    """n in [0,91) with n=e mod7 and n=t mod13."""
    n = e + Q * ((pow(Q, -1, P) * (t - e)) % P)
    require(n % Q == e and n % P == t, "CRT index failed")
    return n


def cyclotomic_audit(table):
    """Exact modular sector audits for the C91 polynomial.

    Degree-zero gcd modulo 2 proves the characteristic-zero cyclotomic gcd
    is one (the cyclotomic factor is monic).  Degree-zero gcds modulo 7 and
    13 prove that sector resultant is a 91-unit without constructing its
    enormous integer representative.
    """
    import sympy as sp

    x = sp.symbols("x")
    clock = tuple(sum(row) for row in table)
    target = tuple(sum(table[e][t] for e in range(Q)) for t in range(P))
    crt = [0] * (P * Q)
    for e in range(Q):
        for t in range(P):
            crt[crt_index(e, t)] = table[e][t]

    polynomials = {
        7: sp.Poly(sum(value * x**i for i, value in enumerate(clock)), x),
        13: sp.Poly(sum(value * x**i for i, value in enumerate(target)), x),
        91: sp.Poly(sum(value * x**i for i, value in enumerate(crt)), x),
    }
    rows = []
    exact_small_resultants = []
    for order in (7, 13, 91):
        phi = sp.Poly(sp.cyclotomic_poly(order, x), x)
        poly = polynomials[order]
        gcd_degrees = tuple(
            sp.gcd(
                sp.Poly(poly, x, modulus=prime),
                sp.Poly(phi, x, modulus=prime),
            ).degree()
            for prime in (2, 7, 13)
        )
        rows.append((order, gcd_degrees))
        if order in (7, 13):
            exact_small_resultants.append((
                order,
                int(poly.resultant(phi)),
            ))
    augmentation = sum(sum(row) for row in table)
    return (
        clock, target, tuple(crt), augmentation, tuple(rows),
        tuple(exact_small_resultants),
    )


def target_uniformizer_audit(target):
    """Write S(u)=epsilon V(epsilon), invert V, and certify SK=epsilon.

    Here epsilon=u-1 and F_13[C_13]=F_13[epsilon]/(epsilon^13).
    The returned K is given both in epsilon and ordinary u coordinates.
    """
    target_mod13 = tuple(value % P for value in target)
    epsilon_coefficients = tuple(
        sum(target_mod13[degree] * comb(degree, power)
            for degree in range(power, P)) % P
        for power in range(P)
    )
    require(epsilon_coefficients[0] == 0
            and epsilon_coefficients[1] == 7,
            "target polynomial lost its simple augmentation zero")
    v_epsilon = tuple(epsilon_coefficients[i + 1] if i + 1 < P else 0
                      for i in range(P))
    require(v_epsilon[0] == 7, "uniformizer factor ceased to be a unit")

    inverse = [0] * P
    inverse[0] = pow(v_epsilon[0], -1, P)
    # Only modulo epsilon^12 is needed after multiplication by S=epsilon V;
    # choose the canonical degree-12 coefficient zero.
    for degree in range(1, P - 1):
        tail = sum(
            v_epsilon[i] * inverse[degree - i]
            for i in range(1, degree + 1)
        )
        inverse[degree] = -inverse[0] * tail % P
    inverse = tuple(inverse)

    product_epsilon = tuple(
        sum(epsilon_coefficients[i] * inverse[degree - i]
            for i in range(degree + 1)) % P
        for degree in range(P)
    )
    require(product_epsilon == (0, 1) + (0,) * 11,
            "epsilon-adic inverse failed")

    inverse_u = tuple(
        sum(
            inverse[power] * comb(power, degree) * (-1) ** (power - degree)
            for power in range(degree, P)
        ) % P
        for degree in range(P)
    )
    cyclic_product = tuple(
        sum(
            target_mod13[i] * inverse_u[(degree - i) % P]
            for i in range(P)
        ) % P
        for degree in range(P)
    )
    require(cyclic_product == (12, 1) + (0,) * 11,
            "group-ring identity S*K=u-1 failed")
    derivative = sum(
        degree * target_mod13[degree] for degree in range(P)
    ) % P
    require(derivative == 7, "target tangent/Bockstein coefficient changed")
    require(target_mod13 == (0, 0, 9, 0, 0, 0, 0, 0, 2, 2, 2, 2, 9),
            "target mod-13 profile changed")
    require(epsilon_coefficients == (0, 7, 8, 9, 12, 2, 4, 4, 7, 6, 7, 6, 9),
            "epsilon expansion changed")
    require(v_epsilon == (7, 8, 9, 12, 2, 4, 4, 7, 6, 7, 6, 9, 0),
            "uniformizer factor changed")
    require(inverse == (2, 7, 8, 10, 1, 10, 11, 4, 8, 3, 2, 4, 0),
            "epsilon-basis decoder changed")
    require(inverse_u == (7, 3, 3, 3, 12, 2, 9, 10, 9, 8, 10, 4, 0),
            "u-basis decoder changed")
    return (
        target_mod13, epsilon_coefficients, v_epsilon,
        inverse, inverse_u, cyclic_product, derivative,
    )


def cyclic_convolution_mod13(left, right):
    return tuple(
        sum(left[i] * right[(degree - i) % P] for i in range(P)) % P
        for degree in range(P)
    )


def main():
    module, _prefixes, _a, _b, rails, present, _starts = (
        marked.lift.m.core.build_carrier_data()
    )
    delayed = marked.marked_prefixes(
        module,
        marked.private.build_pair_prefixes(module),
        marked.two.deepest_fork(module),
    )
    source = marked.two.exclusive_source(module, 3)
    clock_comb = tuple(
        module.make_comb(module.C1, 182, 26 * e - 13, 26 * e + 13)
        for e in range(Q)
    )
    source_weight, target_weight, rail_common = marked.rail_data(
        rails, marked.RAIL
    )

    tables = {
        name: [[0 for _t in range(P)] for _e in range(Q)]
        for name in ("A", "B", "Ms", "Mt", "L", "R")
    }
    reduced = {
        name: [[None for _t in range(P)] for _e in range(Q)]
        for name in tables
    }
    source_clock_carriers = []
    for e in range(Q):
        raw_static = tuple(marked.two.intersect_sorted(source, clock_comb[e]))
        raw_static = tuple(module.subtract_comb(
            raw_static, module.W[1], 182, -13, 13
        ))
        raw_static = tuple(module.subtract_comb(
            raw_static, module.C2, 182, -13, 13
        ))
        source_clock_carriers.append(raw_static)
        for t in range(P):
            raw_a = tuple(module.subtract_comb(
                raw_static, module.W[2], 182,
                -14 * t - 13, -14 * t + 13,
            ))
            raw_a = tuple(module.subtract_comb(
                raw_a, module.C3, 182,
                14 * t - 13, 14 * t + 13,
            ))
            A = marked.intersect(rail_common, raw_a)
            B = marked.intersect(
                rail_common, marked.shift_union(raw_a, marked.SHIFT)
            )
            M = marked.intersect(A, B)
            L = wing.difference(A, B)
            right = wing.difference(B, A)
            require(marked.merge(M + L) == marked.merge(A)
                    and marked.merge(M + right) == marked.merge(B),
                    f"physical wing decomposition failed at {(e,t)}")
            packed = wing.multiplexed_vectors(
                module, delayed, present, source_weight, target_weight,
                A, B, M,
            )
            vectors = {
                "A": packed["A"],
                "B": packed["B"],
                "Ms": packed["M"],
                "Mt": packed["M"],
            }
            vectors["L"] = tuple(
                vectors["A"][i] - vectors["Ms"][i] for i in range(Q)
            )
            vectors["R"] = tuple(
                vectors["B"][i] - vectors["Mt"][i] for i in range(Q)
            )
            for name, vector in vectors.items():
                tables[name][e][t] = wing.scalar_amplitude(vector)
            for name, root in (
                ("A", 12), ("B", 1), ("Ms", 12), ("Mt", 1),
                ("L", 12), ("R", 1),
            ):
                reduced[name][e][t] = marked.unit_reduction(
                    vectors[name], root
                )

    # Physical present-clock sheets are disjoint in one source coordinate.
    require(all(
        not marked.intersect(source_clock_carriers[e], source_clock_carriers[f])
        for e in range(Q) for f in range(e + 1, Q)
    ), "distinct physical present clocks acquired a cross-sheet intersection")

    # Replay the exact fixed-t=3 full-clock boundary from the hostile audit.
    # marked.unit_reduction uses the fixed-chart root convention.  The
    # full-clock hostile audit uses the common full-row normalization, which
    # differs by the global source sign 12.
    left_mod13 = tuple(
        12 * reduced["L"][e][3][0][0] % P for e in range(Q)
    )
    right_mod13 = tuple(
        12 * reduced["R"][e][3][0][0] % P for e in range(Q)
    )
    require(left_mod13 == (0, 0, 0, 12, 0, 0, 0)
            and right_mod13 == (0, 9, 2, 2, 0, 0, 0),
            f"full-clock t=3 wing boundary changed: {left_mod13=} {right_mod13=}")

    raw_right = tuple(tuple(row) for row in tables["R"])
    global_content = reduce(gcd, (abs(value) for row in raw_right for value in row))
    require(global_content > 0, "joint right wing vanished")
    primitive = tuple(tuple(value // global_content for value in row)
                      for row in raw_right)
    support = tuple(
        (e, t) for e in range(Q) for t in range(P) if primitive[e][t]
    )
    ranks = (
        rank_q(primitive), rank_mod(primitive, 7), rank_mod(primitive, 13)
    )
    clock, target, crt, augmentation, sectors, small_resultants = (
        cyclotomic_audit(primitive)
    )
    uniformizer = target_uniformizer_audit(target)
    sector_units = tuple(
        (order, degrees[1] == degrees[2] == 0)
        for order, degrees in sectors
    )
    full_c91_unit = (
        gcd(augmentation, 91) == 1
        and all(is_unit for _order, is_unit in sector_units)
    )
    augmentation_quotient_unit = all(
        is_unit for _order, is_unit in sector_units
    )

    expected_primitive = (
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 966606, 966606, 966606, 966606, 966606, 966606,
         966606, 966606, 966606, 58479663),
        (0, 0, 966574, 966574, 966574, 966574, 966574, 966574,
         966574, 966574, 0, 0, 128071055),
        (0, 0, 966574, 966574, 966574, 966574, 966574, 966574,
         0, 0, 966574, 966574, 122754898),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    )
    expected_support = (
        (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9),
        (1, 10), (1, 11), (1, 12),
        (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8),
        (2, 9), (2, 12),
        (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7),
        (3, 10), (3, 11), (3, 12),
    )
    expected_target = (
        0, 0, 1933148, 2899754, 2899754, 2899754, 2899754,
        2899754, 1933180, 1933180, 1933180, 1933180, 309305616,
    )
    expected_resultants = (
        (7, 2164833805035694204114534905543204737802313970779),
        (13, 708762831620208775530213399682753305758900865474119127375721370505898041395766813306566259706602082304),
    )
    require(global_content == 5905329039529920,
            "joint raw-table content changed")
    require(primitive == expected_primitive,
            "primitive joint table changed")
    require(support == expected_support and ranks == (3, 1, 3),
            "joint support or rank census changed")
    require(clock == (0, 67179117, 135803647, 130487490, 0, 0, 0)
            and target == expected_target and augmentation == 333470254,
            "clock/target collapse or augmentation changed")
    require(sectors == (
        (7, (0, 0, 0)), (13, (12, 0, 1)), (91, (0, 0, 0))
    ), "cyclotomic gcd-degree certificate changed")
    require(small_resultants == expected_resultants,
            "small exact resultants changed")
    require(small_resultants[0][1] % 91 == 50
            and small_resultants[1][1] % 91 == 78
            and small_resultants[1][1] % 13 == 0
            and small_resultants[1][1] % (13 * 13) != 0,
            "small resultant residue/valuation boundary changed")
    require(sector_units == ((7, True), (13, False), (91, True))
            and not full_c91_unit and not augmentation_quotient_unit,
            "sector-unit classification changed")
    require(tuple(i for i, value in enumerate(crt) if value) == (
        2, 3, 8, 9, 10, 16, 17, 22, 24, 29, 30, 31, 36, 38,
        43, 44, 45, 50, 51, 57, 58, 59, 64, 71, 72, 80, 85, 86,
    ), "CRT support changed")

    # Intrinsic coefficient Bockstein.  The raw table has exact 13-valuation
    # one in its common content, so raw/13 mod13 is canonical; the remaining
    # common content is the unit 11.  Rescale K by 11^{-1}=6.
    require(global_content % 13 == 0 and global_content % (13 * 13) != 0,
            "raw common content lost exact 13-valuation one")
    bockstein_unit = (global_content // 13) % 13
    require(bockstein_unit == 11, "raw Bockstein unit changed")
    bockstein_rows = tuple(
        tuple((raw_right[e][t] // 13) % 13 for t in range(P))
        for e in range(Q)
    )
    require(bockstein_rows == tuple(
        tuple(bockstein_unit * primitive[e][t] % 13 for t in range(P))
        for e in range(Q)
    ), "raw/13 Bockstein disagrees with primitive shape")
    k_bockstein = tuple(6 * value % 13 for value in uniformizer[4])
    require(k_bockstein == (3, 5, 5, 5, 7, 12, 2, 8, 2, 9, 8, 11, 0),
            "intrinsic Bockstein decoder changed")

    decoded_rows = tuple(
        cyclic_convolution_mod13(
            tuple(value % P for value in row), uniformizer[4]
        )
        for row in primitive
    )
    require(decoded_rows == (
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (2, 6, 7, 7, 9, 2, 8, 12, 6, 5, 7, 11, 6),
        (0, 11, 1, 7, 6, 9, 9, 4, 12, 8, 8, 5, 8),
        (10, 10, 5, 12, 11, 2, 9, 10, 8, 0, 11, 10, 12),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    ), "rowwise target decoder changed")
    decoded_chart_cochain = tuple(row[0] for row in decoded_rows)
    decoded_target_augmentation = tuple(
        sum(decoded_rows[e][t] for e in range(Q)) % P for t in range(P)
    )
    require(decoded_chart_cochain == (0, 2, 0, 10, 0, 0, 0)
            and sum(decoded_chart_cochain) % P == 12,
            "sparse chart cochain or its augmentation changed")
    require(decoded_target_augmentation == (12, 1) + (0,) * 11,
            "decoded target boundary changed")
    decoded_bockstein_rows = tuple(
        cyclic_convolution_mod13(row, k_bockstein)
        for row in bockstein_rows
    )
    require(decoded_bockstein_rows == decoded_rows,
            "intrinsic Bockstein decoder changed the sparse cochain")
    bockstein_target = tuple(
        sum(bockstein_rows[e][t] for e in range(Q)) % 13
        for t in range(P)
    )
    require(cyclic_convolution_mod13(bockstein_target, k_bockstein)
            == (12, 1) + (0,) * 11,
            "intrinsic Bockstein target decoder failed")
    uniform_chart_convolution = (sum(decoded_chart_cochain) % 13,) * Q
    require(uniform_chart_convolution == (12,) * Q,
            "uniform chart averaging ceased to be pointwise -1")

    print("JOINT C7xC13 CORRECTED RIGHT-WING/COFIBER SPECTRUM")
    print("status=THM-2771 FINITE-EXACT candidate; no LRC conclusion")
    print(f"raw_right_table={raw_right}")
    print(f"global_content={global_content}")
    print(f"primitive_right_table={primitive}")
    print(f"support_count={len(support)} support={support}")
    print(f"matrix_ranks_Q_F7_F13={ranks}")
    print(f"clock_collapse={clock}")
    print(f"target_collapse={target}")
    print(f"augmentation={augmentation} mod7={augmentation%7} mod13={augmentation%13}")
    print(f"target_mod13={uniformizer[0]}")
    print(f"target_epsilon_coefficients={uniformizer[1]}")
    print(f"target_uniformizer_V_epsilon={uniformizer[2]}")
    print(f"target_uniformizer_inverse_K_epsilon={uniformizer[3]}")
    print(f"target_uniformizer_inverse_K_u={uniformizer[4]}")
    print(f"target_group_ring_product_S_times_K={uniformizer[5]}")
    print(f"target_tangent_class_I_mod_I2={uniformizer[6]}*(u-1)")
    print(f"decoded_primitive_rows_mod13={decoded_rows}")
    print(f"decoded_t0_chart_cochain={decoded_chart_cochain} chart_sum={sum(decoded_chart_cochain)%P}")
    print(f"decoded_target_column_sums={decoded_target_augmentation}")
    print("scaled_holotopy_candidate=7*a*decoded_t0 has chart sum -7*a mod13")
    print(f"raw_content_v13=1 bockstein_unit_(g0/13)_mod13={bockstein_unit}")
    print(f"bockstein_decoder_K_u={k_bockstein}")
    print(f"decoded_bockstein_rows_mod13={decoded_bockstein_rows}")
    print(f"decoded_t0_convolve_uniform_C7={uniform_chart_convolution}")
    print("pointwise_holotopy_candidate=a*(C*N7)=(-a,...,-a), cancelling a constant (+a)-chart obstruction after a granted chart identification")
    print(f"cyclotomic_sectors_(order,gcd_degrees_mod2_mod7_mod13)={sectors}")
    print(f"exact_order7_order13_resultants={small_resultants}")
    print(f"sector_units_mod91={sector_units}")
    print(f"full_C91_group_ring_unit={full_c91_unit} augmentation_quotient_unit={augmentation_quotient_unit} fully_mixed_order91_unit={dict(sector_units)[91]}")
    print(f"crt_support={tuple(i for i,value in enumerate(crt) if value)}")
    print(f"fixed_t3_mod13_left={left_mod13} right={right_mod13}")
    print("SCOPE: a nonzero/invertible coefficient spectrum is not a semantic arm, endpoint map, or noncommuting transition")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
