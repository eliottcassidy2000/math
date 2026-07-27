#!/usr/bin/env python3
"""Exact companion for THM-2581: the same-fibre b/r=5 owner-clock host.

This script reuses the independently route-doubled interval/profile engine in
``lrc14_base_only_bridge_opus_20260728.py``.  It changes the mathematical
input, not the integration rules:

    E = E_1,
    Q = Q_{1,{b}},
    f = 1_Q P_{13^2} 1_E,
    d = 13^5,
    U = P_d f,
    V = P_d 1_E.

The root charts and owner-clock cells are exactly those of THM-2471 and
THM-2449:

    A(y,u)=U((y+u)/13),       F(y,u)=V((y+u)/13),
    c_s(y)=sum_u A(y,u+s)F(y,u),
    C_ell(s)=integral_{cell_ell} c_s(y)dy,
    W(k,ell)=1/169 sum_s C_ell(s) zeta_13^(-ks).

Every owner-cell restriction remains inside the same physical y-integral.
No separately marginalized controls are multiplied.  The script proves the
two exact integration routes agree, checks the stored return and I_5 anchors,
checks both marginals and Hermitian symmetry, computes every cyclotomic 2x2
minor, and gives an independent two-prime certificate for their
nonvanishing.  The sigma={a}, r=3 host is recomputed as a hostile control: its
central reflected columns coincide and precisely 90 minors vanish.

Scope: one canonical typed row, word {b}, packet clock K=2.  The result is a
host-array theorem only.  It is not the THM-2512 interaction defect, does not
perform the theta=t-2u contraction, and implies no physical current, row
exclusion, or LRC(14) conclusion.
"""

from fractions import Fraction
from itertools import combinations
from math import lcm

import lrc14_base_only_bridge_opus_20260728 as base


W = base.W
T_DEN = base.T_DEN
P = 13
QMOD = 7
RPACK = 13**2
DEPTH_B = 13**5
DEPTH_A = 13**3
I5_B = Fraction(48602521488933856, 337437093630814766589)
RHO_B = Fraction(35505957232, 16132831966251)

PAT_QB = {
    0: "gout",
    1: "out",
    2: "out",
    3: "out",
    4: "out",
    5: "out",
    6: "out",
    7: "out",
    8: "in",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def buckets(C):
    """Integer bucket vectors for W before its common denominator."""
    return [
        [
            [
                sum(C[s][ell] for s in range(P) if (-k * s) % P == j)
                for j in range(P)
            ]
            for ell in range(QMOD)
        ]
        for k in range(P)
    ]


def marginals(C):
    Cglob = [sum(row) for row in C]
    J = []
    for k in range(P):
        b = [0] * P
        for s in range(P):
            b[(-k * s) % P] += Cglob[s]
        J.append(b)
    return Cglob, J


def minor_data(WI):
    records = []
    for k1, k2 in combinations(range(P), 2):
        for ell1, ell2 in combinations(range(QMOD), 2):
            value = base.bsub(
                base.bmul(WI[k1][ell1], WI[k2][ell2]),
                base.bmul(WI[k1][ell2], WI[k2][ell1]),
            )
            records.append((k1, k2, ell1, ell2, value))
    require(len(records) == 1638, "minor universe changed")
    return records


def eval_bucket(bucket, modulus, root):
    return sum((coefficient % modulus) * pow(root, j, modulus)
               for j, coefficient in enumerate(bucket)) % modulus


def rank_mod(matrix, modulus):
    work = [[entry % modulus for entry in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, modulus)
        work[rank] = [(entry * inverse) % modulus for entry in work[rank]]
        for row in range(len(work)):
            if row != rank and work[row][column]:
                factor = work[row][column]
                work[row] = [
                    (left - factor * right) % modulus
                    for left, right in zip(work[row], work[rank])
                ]
        rank += 1
    return rank


def modular_minor_cover(WI, records):
    """Independent sufficient witnesses that all integer minor numerators
    are nonzero in Q(zeta_13).

    If a bucket polynomial represented zero in Q(zeta_13), it would be a
    multiple of Phi_13 and hence vanish at every primitive 13th root in every
    characteristic away from 13.  A single nonzero reduction is therefore a
    valid nonvanishing certificate.  The two fixed primes, neither of which
    divides the common denominator, cover all minors.
    """
    witnesses = ((79, 18), (131, 107))
    remaining = set(range(len(records)))
    cover = []
    for modulus, root in witnesses:
        require(pow(root, 13, modulus) == 1 and root != 1,
                "listed modular root is not primitive of order 13")
        values = [
            [eval_bucket(WI[k][ell], modulus, root)
             for ell in range(QMOD)]
            for k in range(P)
        ]
        caught = set()
        for index in remaining:
            k1, k2, ell1, ell2, _ = records[index]
            determinant = (
                values[k1][ell1] * values[k2][ell2]
                - values[k1][ell2] * values[k2][ell1]
            ) % modulus
            if determinant:
                caught.add(index)
        remaining.difference_update(caught)
        cover.append((modulus, root, len(caught), len(remaining)))
    require(not remaining, "two-prime minor certificate is incomplete")
    return cover


def centred_nonzero_count(WI, J):
    count = 0
    for k in range(P):
        for ell in range(QMOD):
            centred = base.bsub(base.bscale(WI[k][ell], QMOD), J[k])
            if not base.is_zero_b(centred):
                count += 1
    return count


def centred_denominator_data(WI, J, common_denominator):
    universal = 1
    minimum = None
    for k in range(P):
        for ell in range(QMOD):
            centred = base.bsub(base.bscale(WI[k][ell], QMOD), J[k])
            for coordinate in base.canon(centred):
                if coordinate:
                    value = Fraction(coordinate, QMOD * common_denominator)
                    universal = lcm(universal, value.denominator)
                    absolute = abs(value)
                    if minimum is None or absolute < minimum:
                        minimum = absolute
    return universal, minimum


def colour_saturation(WI):
    return [
        sum(not base.is_zero_b(WI[k][ell]) for k in range(1, P))
        for ell in range(QMOD)
    ]


def exact_host(E, Q, packet_clock, depth):
    result = base.bridge_tables(
        E, Q, P, QMOD, packet_clock, depth, T_DEN, verbose=False
    )
    C = result["C"]
    Cglob, J = marginals(C)
    require(Cglob == result["Cglob"], "marginal reconstruction mismatch")
    WI = buckets(C)
    return result, C, Cglob, WI, J


def main():
    print("== THM-2581: the b-word depth-five owner-clock host ==")
    print("row:", W)
    print("owner j=1; sigma={b}; packet clock K=2; r=5; d=13^5")
    print("all owner-clock restrictions occur inside one common y-integral")

    print("\n== hostile engine control ==")
    base.toy_test()

    E = base.build_set(base.PAT_E, base.ZELL)
    QB = base.build_set(PAT_QB, base.ZELL)
    lenE = base.check_intervals(E, T_DEN)
    lenQB = base.check_intervals(QB, T_DEN)
    require((len(E), Fraction(lenE, T_DEN))
            == (57072, Fraction(1882176, 28589561)), "E anchor changed")
    require((len(QB), Fraction(lenQB, T_DEN))
            == (131762, Fraction(4839079319, 190921088358)),
            "Q_b anchor changed")
    print("\n== exact typed atoms ==")
    print(f"E_1: intervals={len(E)} mass={Fraction(lenE, T_DEN)}")
    print(f"Q_1,{{b}}: intervals={len(QB)} mass={Fraction(lenQB, T_DEN)}")

    print("\n== b/r=5 table; route-doubled integration ==")
    result, C, Cglob, WI, J = exact_host(E, QB, RPACK, DEPTH_B)
    DENC = result["DEN"]
    DENW = P * P * DENC
    require(DENC == RPACK * DEPTH_B**2 * T_DEN,
            "common C denominator changed")
    require(Fraction(result["int_f_num"], RPACK * T_DEN) == RHO_B,
            "stored b-word return changed")
    require(all(value == 0 for value in C[0]),
            "first-collision disjointness failed cellwise")
    require(Fraction(sum(Cglob), DENC) == P * P * I5_B,
            "sum_s C(s) does not reproduce 169 I_5")
    print(f"rho^2({{b}})={RHO_B}")
    print(f"I_5({{b}},K=2)={I5_B}")
    print(f"DENC={DENC}={base.factor_str(DENC)}")
    print(f"U pieces={result['nU']} V pieces={result['nV']}")
    print("two integration routes agree on all 91 cells and 13 globals: PASS")

    for k in range(P):
        row_sum = [0] * P
        for ell in range(QMOD):
            row_sum = base.badd(row_sum, WI[k][ell])
        require(base.canon(row_sum) == base.canon(J[k]),
                "owner-clock marginal failed")
    for ell in range(QMOD):
        column_sum = [0] * P
        for k in range(P):
            column_sum = base.badd(column_sum, WI[k][ell])
        require(base.is_zero_b(column_sum), "root-colour marginal failed")
    require(Fraction(sum(Cglob), DENW) == I5_B, "J(0) != I_5")
    require(all(not base.is_zero_b(J[k]) for k in range(1, P)),
            "a nonzero global root colour vanished")
    for k in range(P):
        for ell in range(QMOD):
            require(
                base.canon(WI[(-k) % P][ell])
                == base.canon(base.bconj(WI[k][ell])),
                "Hermitian symmetry failed",
            )
    print("sum_ell W(k,ell)=J(k), J(0)=I_5: PASS")
    print("sum_k W(k,ell)=0 for every ell: PASS")
    print("all 12 global nonzero root colours survive: PASS")
    print("Hermitian symmetry: PASS")

    print("\n== reflection law and oriented H-drift ==")
    for ell in range(QMOD):
        for s in range(P):
            require(C[s][ell] == C[(-s) % P][(-ell) % QMOD],
                    "reflection law C_ell(s)=C_-ell(-s) failed")
    delta = [C[s][4] - C[s][3] for s in range(P)]
    require(delta[0] == 0 and all(delta[(-s) % P] == -delta[s]
                                  for s in range(P)),
            "central reflection drift is not odd")
    require(any(delta), "central reflection drift unexpectedly vanished")
    delta_hat = []
    for k in range(P):
        b = [0] * P
        for s in range(P):
            b[(-k * s) % P] += delta[s]
        delta_hat.append(b)
    require(base.is_zero_b(delta_hat[0]), "odd drift has nonzero mean")
    require(all(not base.is_zero_b(delta_hat[k]) for k in range(1, P)),
            "odd drift misses a nonzero root colour")
    print("C_ell(s)=C_{-ell}(-s): PASS")
    print("Delta(s)=C_4(s)-C_3(s) is odd, nonzero, and mean-zero")
    print("Delta numerator vector over DENC:", delta)
    print("Delta(1)=", Fraction(delta[1], DENC))
    print("all 12 nonzero Fourier colours of Delta survive: PASS")

    print("\n== exact mixedness ==")
    records = minor_data(WI)
    zero_minors = [record[:4] for record in records
                   if base.is_zero_b(record[4])]
    require(not zero_minors, "a b/r=5 2x2 minor vanished")
    print("cyclotomic reduction: 1638/1638 2x2 minors nonzero")
    cover = modular_minor_cover(WI, records)
    for modulus, root, caught, remaining in cover:
        print(
            f"mod {modulus}, primitive root {root}: certifies {caught} new; "
            f"remaining={remaining}"
        )
    print("independent two-prime certificate covers all 1638 minors: PASS")
    require(centred_nonzero_count(WI, J) == 91,
            "a centred host entry vanished")
    require(colour_saturation(WI) == [12] * QMOD,
            "an owner cell missed a nonzero root colour")
    Wc = [[base.bsub(base.bscale(WI[k][ell], QMOD), J[k])
           for ell in range(QMOD)] for k in range(P)]
    for k in range(P):
        total = [0] * P
        for ell in range(QMOD):
            total = base.badd(total, Wc[k][ell])
        require(base.is_zero_b(total), "centred F_7 row sum failed")
    for ell in range(QMOD):
        total = [0] * P
        for k in range(P):
            total = base.badd(total, Wc[k][ell])
        require(base.is_zero_b(total), "centred F_13 column sum failed")
    k_dependent = any(
        base.canon(Wc[k][ell]) != base.canon(Wc[0][ell])
        for ell in range(QMOD) for k in range(1, P)
    )
    require(k_dependent, "centred host unexpectedly k-independent")
    raw_mod79 = [
        [eval_bucket(WI[k][ell], 79, 18) for ell in range(QMOD)]
        for k in range(P)
    ]
    centred_mod79 = [
        [eval_bucket(Wc[k][ell], 79, 18) for ell in range(QMOD)]
        for k in range(P)
    ]
    require(rank_mod(raw_mod79, 79) == 7,
            "raw host lost full column rank modulo 79")
    require(rank_mod(centred_mod79, 79) == 6,
            "centred host lost maximal row-zero rank modulo 79")
    universal, minimum = centred_denominator_data(WI, J, DENW)
    require(
        universal == 184240653122424862557594
        and minimum == Fraction(76, 11823941286254964867),
        "centred denominator/floor anchor changed",
    )
    print("Wc(k,ell)=W(k,ell)-J(k)/7: 91/91 entries nonzero")
    print("each owner cell fires all 12 nonzero root colours")
    print("rank_Q(zeta13)(W)=7; rank(Wc)=6 (maximal under row-zero): PASS")
    print("T1 row/column zero; T2 nonvacuous k-dependence: PASS")
    print(f"T3 universal denominator={universal}={base.factor_str(universal)}")
    print(f"T3 minimum nonzero coordinate={minimum} >= 1/D")

    print("\n== hostile comparison: separate a/r=3 fibre ==")
    QA = base.build_set(base.PAT_QA, base.ZELL)
    _, CA, _, WIA, _ = exact_host(E, QA, RPACK, DEPTH_A)
    require(all(CA[s][3] == CA[s][4] for s in range(P)),
            "stored a/r=3 central-column symmetry changed")
    records_a = minor_data(WIA)
    zero_a = [record[:4] for record in records_a
              if base.is_zero_b(record[4])]
    require(len(zero_a) == 90, "a/r=3 zero-minor anchor changed")
    print("a/r=3 has C_3(s)=C_4(s) and exactly 90 zero minors: PASS")
    print("b/r=5 replaces that even central profile by the nonzero odd Delta")

    print("\n== finite clock robustness; not part of the theorem quantifier ==")
    expected_rho = {
        3: Fraction(269990565109, 139817877040842),
        4: Fraction(4377228237904, 2726448602296419),
    }
    for K in (3, 4):
        rK, _, _, WIK, _ = exact_host(E, QB, 13**K, DEPTH_B)
        require(Fraction(rK["int_f_num"], 13**K * T_DEN)
                == expected_rho[K], "finite-clock return anchor changed")
        nzK = sum(not base.is_zero_b(record[4])
                  for record in minor_data(WIK))
        require(nzK == 1638, "finite-clock hostile lost mixedness")
        print(f"K={K}: stored return and 1638/1638 minors: PASS")

    print("\nSCOPE: canonical typed row, sigma={b}, K=2, r=5 host only")
    print("No cross-word multiplication, slaved contraction, physical current,")
    print("row exclusion, or LRC(14) conclusion.  All exact checks passed.")


if __name__ == "__main__":
    main()
