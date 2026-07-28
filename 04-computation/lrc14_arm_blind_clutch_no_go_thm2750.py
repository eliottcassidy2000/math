#!/usr/bin/env python3
"""Exact controls for THM-2750.

This canonical exact companion checks, over Q and F_13, that an
arm-blind scalar (indeed any arm-blind auxiliary operator) cannot move the
C3-invariant line into a charged summand.  It also finds the two smallest
escapes: a signed three-arm diagonal and the positive four-point S4 carrier.
"""

from fractions import Fraction
from itertools import permutations, product
from collections import Counter


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def eye(n):
    return tuple(
        tuple(Fraction(int(i == j)) for j in range(n)) for i in range(n)
    )


def zero(m, n):
    return tuple(tuple(Fraction(0) for _ in range(n)) for _ in range(m))


def add(a, b):
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
        for i in range(len(a))
    )


def sub(a, b):
    return tuple(
        tuple(a[i][j] - b[i][j] for j in range(len(a[0])))
        for i in range(len(a))
    )


def scale(c, a):
    return tuple(tuple(c * x for x in row) for row in a)


def mul(a, b):
    return tuple(
        tuple(
            sum((a[i][k] * b[k][j] for k in range(len(b))), Fraction(0))
            for j in range(len(b[0]))
        )
        for i in range(len(a))
    )


def transpose(a):
    return tuple(tuple(a[i][j] for i in range(len(a))) for j in range(len(a[0])))


def apply(a, v):
    return tuple(
        sum((a[i][j] * v[j] for j in range(len(v))), Fraction(0))
        for i in range(len(a))
    )


def kron(a, b):
    return tuple(
        tuple(
            a[i // len(b)][j // len(b[0])] * b[i % len(b)][j % len(b[0])]
            for j in range(len(a[0]) * len(b[0]))
        )
        for i in range(len(a) * len(b))
    )


def rank(a):
    rows = [list(row) for row in a]
    m = len(rows)
    n = len(rows[0])
    out = 0
    for col in range(n):
        pivot = next((i for i in range(out, m) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[out], rows[pivot] = rows[pivot], rows[out]
        z = rows[out][col]
        rows[out] = [x / z for x in rows[out]]
        for i in range(m):
            if i != out and rows[i][col]:
                z = rows[i][col]
                rows[i] = [rows[i][j] - z * rows[out][j] for j in range(n)]
        out += 1
    return out


def hs2(a):
    return sum((x * x for row in a for x in row), Fraction(0))


def diag(values):
    n = len(values)
    return tuple(
        tuple(Fraction(values[i]) if i == j else Fraction(0) for j in range(n))
        for i in range(n)
    )


def perm_matrix(p):
    n = len(p)
    a = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for source, target in enumerate(p):
        a[target][source] = Fraction(1)
    return tuple(tuple(row) for row in a)


def generated_matrix_group(gens):
    ident = eye(len(gens[0]))
    group = {ident}
    frontier = list(gens)
    while frontier:
        x = frontier.pop()
        if x in group:
            continue
        group.add(x)
        current = tuple(group)
        frontier.extend(mul(x, y) for y in current)
        frontier.extend(mul(y, x) for y in current)
    return frozenset(group)


def fmt_matrix(a):
    def fmt(x):
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    return "[" + ";".join(",".join(fmt(x) for x in row) for row in a) + "]"


def mod_matrix(a, p):
    return tuple(
        tuple((x.numerator * pow(x.denominator, -1, p)) % p for x in row)
        for row in a
    )


def mod_rank(a, p):
    rows = [list(row) for row in a]
    m = len(rows)
    n = len(rows[0])
    out = 0
    for col in range(n):
        pivot = next((i for i in range(out, m) if rows[i][col] % p), None)
        if pivot is None:
            continue
        rows[out], rows[pivot] = rows[pivot], rows[out]
        inv = pow(rows[out][col] % p, -1, p)
        rows[out] = [(x * inv) % p for x in rows[out]]
        for i in range(m):
            if i != out and rows[i][col] % p:
                z = rows[i][col] % p
                rows[i] = [(rows[i][j] - z * rows[out][j]) % p for j in range(n)]
        out += 1
    return out


def main():
    # The THM-2721 rational amplitude has denominator exactly 13^11.  Thus
    # its positive Q-value has no naive reduction in F_13; clearing the common
    # denominator produces residue 5, but that is a lattice choice, not a
    # canonical identification with THM-2744's normalized coefficient line.
    corolla_num = 2089891528601250990
    corolla_den = 1792160394037
    require(corolla_den == 13 ** 11, "THM2721 amplitude denominator")
    require(corolla_num % 13 == 5, "cleared THM2721 amplitude residue")

    # Regular C3 arm carrier.
    r = perm_matrix((1, 2, 0))
    r2 = mul(r, r)
    p3 = scale(Fraction(1, 3), add(add(eye(3), r), r2))
    q3 = sub(eye(3), p3)
    u3 = (Fraction(1),) * 3
    require(apply(p3, u3) == u3, "equal arm is invariant")
    require(mul(p3, p3) == p3, "C3 Reynolds projector")

    # THM-2744 supplies -1 on its normalized coefficient line.  Alone, or
    # tensored with any arm-independent coefficient operator, it is central.
    minus = scale(Fraction(-1), eye(3))
    scalar_off = mul(mul(q3, minus), p3)
    require(scalar_off == zero(3, 3), "scalar -1 cannot create C3 charge")
    require(sub(mul(minus, p3), mul(p3, minus)) == zero(3, 3), "central clutch")

    # Exhaust a box in Q[C3]: Q f(r) P always vanishes.
    group_ring_checks = 0
    for a0, a1, a2 in product(range(-3, 4), repeat=3):
        f = add(add(scale(Fraction(a0), eye(3)), scale(Fraction(a1), r)),
                scale(Fraction(a2), r2))
        require(mul(mul(q3, f), p3) == zero(3, 3), "group-ring no-go")
        require(mul(f, p3) == scale(Fraction(a0 + a1 + a2), p3),
                "augmentation action on invariant line")
        group_ring_checks += 1

    # Any arm-blind auxiliary clutch I_3 tensor T has the same no-go.
    aux = ((Fraction(1), Fraction(2)), (Fraction(3), Fraction(-1)))
    p6 = kron(p3, eye(2))
    q6 = sub(eye(6), p6)
    arm_blind = kron(eye(3), aux)
    require(mul(mul(q6, arm_blind), p6) == zero(6, 6),
            "arbitrary arm-blind auxiliary clutch")

    # Every unsigned permutation of the three arms fixes u3, hence cannot
    # charge it.  This includes every physical three-point involution.
    unsigned_off_ranks = []
    for p in permutations(range(3)):
        s = perm_matrix(p)
        unsigned_off_ranks.append(rank(mul(mul(q3, s), p3)))
    require(unsigned_off_ranks == [0] * 6, "all S3 reindexings are flat")

    # The first linear escape is an unequal diagonal arm gain.  For signs,
    # precisely the six nonconstant patterns have rank one and HS^2=8/9.
    sign_rows = []
    for signs in product((-1, 1), repeat=3):
        d = diag(signs)
        off = mul(mul(q3, d), p3)
        order = len(generated_matrix_group((r, d)))
        sign_rows.append((signs, rank(off), hs2(off), order))
        if len(set(signs)) == 1:
            require(rank(off) == 0, "constant signs are arm-blind")
        else:
            require(rank(off) == 1 and hs2(off) == Fraction(8, 9),
                    "nonconstant sign gate")

    d_one = diag((-1, 1, 1))
    off_one = mul(mul(q3, d_one), p3)
    require(apply(off_one, u3) ==
            (Fraction(-4, 3), Fraction(2, 3), Fraction(2, 3)),
            "one-arm sign charged vector")

    # F_13 control.  Here omega=3 is a primitive cube root and 1/3=9.
    off_one_13 = mod_matrix(off_one, 13)
    scalar_off_13 = mod_matrix(scalar_off, 13)
    require(mod_rank(scalar_off_13, 13) == 0, "F13 scalar no-go")
    require(mod_rank(off_one_13, 13) == 1, "F13 unequal-sign rank")
    omega = 3
    inv3 = pow(3, -1, 13)
    signs = (-1, 1, 1)
    fourier = tuple(
        sum(signs[j] * pow(omega, (-k * j) % 3, 13) for j in range(3))
        * inv3 % 13
        for k in range(3)
    )
    require(fourier == (9, 8, 8), "F13 sign Fourier coefficients")

    # Minimal positive permutation escape: add one C3-fixed coordinate and
    # swap it with one arm.  This is the standard four-point S4 carrier.
    r4 = perm_matrix((0, 2, 3, 1))
    s4 = perm_matrix((1, 0, 2, 3))
    p4 = scale(Fraction(1, 3), add(add(eye(4), r4), mul(r4, r4)))
    q4 = sub(eye(4), p4)
    off4 = mul(mul(q4, s4), p4)
    require(rank(off4) == 1 and hs2(off4) == Fraction(8, 9),
            "S4 rank-one leakage")
    require(len(generated_matrix_group((r4, s4))) == 24, "S4 group order")
    require(apply(off4, (Fraction(1),) * 4) == (Fraction(0),) * 4,
            "constant four-vector remains flat")
    require(apply(off4, (Fraction(0), Fraction(1), Fraction(1), Fraction(1))) ==
            (Fraction(0), Fraction(-2, 3), Fraction(1, 3), Fraction(1, 3)),
            "zero-reference corolla charges")
    require(apply(off4, (Fraction(3), Fraction(-1), Fraction(-1), Fraction(-1))) ==
            (Fraction(0), Fraction(8, 3), Fraction(-4, 3), Fraction(-4, 3)),
            "THM2743 w control")

    # Complete four-point involution classification.  Rank one is equivalent
    # to moving the C3-fixed point.  Three such involutions generate S4 and
    # three (the double transpositions) generate A4, so off-diagonal leakage
    # alone is not an abstract S4 certificate outside THM-2743's marked lift.
    four_involution_rows = []
    for p in permutations(range(4)):
        s = perm_matrix(p)
        if mul(s, s) != eye(4):
            continue
        off = mul(mul(q4, s), p4)
        four_involution_rows.append(
            (p[0] != 0, rank(off), hs2(off), len(generated_matrix_group((r4, s))))
        )
    four_involution_census = Counter(four_involution_rows)
    require(four_involution_census == Counter({
        (False, 0, Fraction(0), 3): 1,
        (False, 0, Fraction(0), 6): 3,
        (True, 1, Fraction(8, 9), 12): 3,
        (True, 1, Fraction(8, 9), 24): 3,
    }), "four-point involution census")

    # A central -1 can reverse an already-created S4 charged vector, but does
    # not change rank, support, or norm.  Creation still belongs to s4.
    central4 = scale(Fraction(-1), eye(4))
    signed_s4_off = mul(mul(q4, mul(central4, s4)), p4)
    require(signed_s4_off == scale(Fraction(-1), off4),
            "central sign only reverses S4 leakage")

    print("THM2750 ARM-BLIND C3 CLUTCH AND MINIMAL MARKED LEAKAGE")
    print("THM2721_A_denominator=13^11 cleared_numerator_mod13=5")
    print("naive_THM2721_to_F13=undefined_without_lattice_choice")
    print(f"group_ring_no_go_checks={group_ring_checks}")
    print(f"scalar_minus_offdiag_rank={rank(scalar_off)}")
    print("arbitrary_arm_blind_aux_offdiag_rank=0")
    print(f"unsigned_three_arm_offdiag_ranks={unsigned_off_ranks}")
    print("sign_pattern_rows=(signs,rank,HS2,generated_group_order)")
    for row in sign_rows:
        print(row)
    print(f"one_arm_sign_offdiag={fmt_matrix(off_one)}")
    print(f"one_arm_sign_F13_fourier={fourier}")
    print(f"S4_offdiag={fmt_matrix(off4)}")
    print("S4_general_input=(x,A,A,A)->(0,2(x-A)/3,(A-x)/3,(A-x)/3)")
    print(f"S4_offdiag_rank={rank(off4)} HS2={hs2(off4)} group_order=24")
    print("four_point_involutions=(moves_fixed,rank,HS2,group_order):count")
    for key in sorted(four_involution_census, key=str):
        print(f"{key}:{four_involution_census[key]}")
    print("central_minus_times_S4=-S4_offdiag")
    print("PASS")


if __name__ == "__main__":
    main()


