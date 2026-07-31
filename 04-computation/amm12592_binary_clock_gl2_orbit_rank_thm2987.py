#!/usr/bin/env python3
"""Exact referee for the candidate THM-2987 ledger-orbit trichotomy.

The homogeneous THM-2976 parity ledger is

    B_l = S_N l^d,
    S_N = p^N + q^N + l^N,       l=p+q,       N=M+1,

over F_2[p,q].  GL(2,F_2) permutes the three nonzero linear forms
{p,q,l}.  This companion checks the resulting orbit ranks, the corner-clock
specialization, the pointed-line stabilizer, and the positive integral
order-three hostile.  All decisions use exact integer/bit arithmetic.
"""


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_power_of_two(n):
    return n > 0 and (n & (n - 1)) == 0


def v2(n):
    require(n > 0, "v2 requires a positive integer")
    ans = 0
    while n % 2 == 0:
        n //= 2
        ans += 1
    return ans


def binom2(n, k):
    """Lucas: binomial(n,k) modulo two."""
    return int(0 <= k <= n and (k & ~n) == 0)


def linear_power(line, degree):
    """Homogeneous F2 polynomial as a bitset indexed by the q exponent.

    line=1,2,3 means p,q,p+q respectively.
    """
    if line == 1:
        return 1
    if line == 2:
        return 1 << degree
    require(line == 3, "unknown linear form")
    bits = 0
    for j in range(degree + 1):
        bits |= binom2(degree, j) << j
    return bits


def multiply(a, b):
    """Multiply two F2 polynomials represented by q-exponent bitsets."""
    out = 0
    aa = a
    i = 0
    while aa:
        if aa & 1:
            bb = b
            j = 0
            while bb:
                if bb & 1:
                    out ^= 1 << (i + j)
                bb >>= 1
                j += 1
        aa >>= 1
        i += 1
    return out


def s_polynomial(n):
    return linear_power(1, n) ^ linear_power(2, n) ^ linear_power(3, n)


def ledger_orbit(n, d):
    s_n = s_polynomial(n)
    return tuple(multiply(s_n, linear_power(line, d)) for line in (1, 2, 3))


def physical_beta(n, d):
    """THM-2976 homogeneous formula in p,q coordinates."""
    endpoint = linear_power(1, n) ^ linear_power(2, n)
    return linear_power(3, n + d) ^ multiply(endpoint, linear_power(3, d))


def rank_f2(vectors):
    basis = {}
    for value in vectors:
        x = value
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                break
    return len(basis)


def matmul(a, b, modulus=None):
    c = (
        (a[0][0] * b[0][0] + a[0][1] * b[1][0],
         a[0][0] * b[0][1] + a[0][1] * b[1][1]),
        (a[1][0] * b[0][0] + a[1][1] * b[1][0],
         a[1][0] * b[0][1] + a[1][1] * b[1][1]),
    )
    if modulus is None:
        return c
    return tuple(tuple(x % modulus for x in row) for row in c)


def matpow(a, exponent, modulus=None):
    out = ((1, 0), (0, 1))
    for _ in range(exponent):
        out = matmul(out, a, modulus)
    return out


def gl2_f2():
    mats = []
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    if (a * d - b * c) % 2 == 1:
                        mats.append(((a, b), (c, d)))
    return mats


def image_lines(matrix):
    """Images of p,q,p+q under row-wise variable substitution."""
    p_image = matrix[0]
    q_image = matrix[1]
    l_image = (p_image[0] ^ q_image[0], p_image[1] ^ q_image[1])
    return p_image, q_image, l_image


def matrix_order_mod2(matrix):
    identity = ((1, 0), (0, 1))
    power = identity
    for order in range(1, 7):
        power = matmul(power, matrix, 2)
        if power == identity:
            return order
    raise RuntimeError("GL2(F2) order exceeded six")


def main():
    rank_counts = {0: 0, 1: 0, 2: 0, 3: 0}
    cases = 0
    for n in range(1, 129):
        s_n = s_polynomial(n)
        require((s_n == 0) == is_power_of_two(n), f"S_N gate failed at N={n}")
        for d in range(65):
            orbit = ledger_orbit(n, d)
            require(physical_beta(n, d) == orbit[2], f"physical factorization failed {n},{d}")
            actual = rank_f2(orbit)
            if is_power_of_two(n):
                expected = 0
            elif d == 0:
                expected = 1
            elif is_power_of_two(d):
                expected = 2
            else:
                expected = 3
            require(actual == expected, f"rank failed N={n},d={d}: {actual}!={expected}")
            if expected == 2:
                require(orbit[0] ^ orbit[1] ^ orbit[2] == 0, "V4 sum failed")
                require(len({0, *orbit}) == 4, "V4 cardinality failed")
            cases += 1
            rank_counts[actual] += 1

    corner_counts = {"v0_rank1": 0, "v1_rank2": 0, "vge2_rank3": 0}
    for n in range(2, 257):
        if is_power_of_two(n):
            continue
        clock = v2(n)
        d = (1 << clock) - 1
        actual = rank_f2(ledger_orbit(n, d))
        expected = 1 if clock == 0 else 2 if clock == 1 else 3
        require(actual == expected, f"corner rank failed N={n},v={clock}")
        key = "v0_rank1" if clock == 0 else "v1_rank2" if clock == 1 else "vge2_rank3"
        corner_counts[key] += 1

    for odd_j in range(3, 32, 2):
        ladder = []
        for clock in range(9):
            n = odd_j << clock
            d = (1 << clock) - 1
            ladder.append(rank_f2(ledger_orbit(n, d)))
        require(ladder == [1, 2, 3, 3, 3, 3, 3, 3, 3], f"odd ladder failed J={odd_j}")

    mats = gl2_f2()
    require(len(mats) == 6, "GL2(F2) size changed")
    nonzero_lines = {(1, 0), (0, 1), (1, 1)}
    for matrix in mats:
        require(set(image_lines(matrix)) == nonzero_lines, "linear-form permutation failed")
    order_profile = {order: sum(matrix_order_mod2(g) == order for g in mats)
                     for order in (1, 2, 3)}
    require(order_profile == {1: 1, 2: 3, 3: 2}, "GL2(F2) order profile changed")
    stabilizer_l = sum(image_lines(g)[2] == (1, 1) for g in mats)
    require(stabilizer_l == 2, "pointed padding-line stabilizer is not C2")

    points = range(4)
    matchings = set()
    for direction in (1, 2, 3):
        matching = tuple(sorted({tuple(sorted((x, x ^ direction))) for x in points}))
        require(len(matching) == 2, "translation direction is not a perfect matching")
        matchings.add(matching)
    require(len(matchings) == 3, "V4 did not produce three K4 matchings")

    t_bar = ((0, 1), (1, 1))
    t_plus = t_bar
    t_signed = ((0, -1), (1, -1))
    require(tuple(tuple(x % 2 for x in row) for row in t_signed) == t_bar,
            "signed lift has wrong mod-two reduction")
    require(matpow(t_signed, 3) == ((1, 0), (0, 1)), "signed lift is not order three")
    require(matpow(t_plus, 2) == ((1, 1), (1, 2)), "positive-lift hostile changed")
    require(all(matpow(t_plus, k) not in (((1, 0), (0, 1)), ((-1, 0), (0, -1)))
                for k in range(1, 13)), "naive positive lift became projectively finite")

    # General exact gate used in the proof: A^3=I would force trace -1;
    # A^3=-I would force trace 1 and determinant 1.  For a,d>=0 with
    # a+d=1, ad=0 and det(A)=ad-bc<=0.  Both projective-order-three cases
    # are therefore impossible for an entrywise nonnegative integral A.
    trace_one_diagonals = [(a, 1 - a) for a in range(2)]
    require(all(a * d == 0 for a, d in trace_one_diagonals),
            "trace-one diagonal gate changed")

    print("THM-2987 binary-clock ledger GL2 orbit referee")
    print(f"box=N1..128,d0..64 cases={cases} rank_counts={rank_counts}")
    print("S_N_zero_iff_N_power_of_two=PASS; physical_beta=S_N*(p+q)^d=PASS")
    print("non_dyadic_rank=d0:1,d_positive_power2:2,otherwise:3=PASS")
    print(f"corner_scan_N2..256={corner_counts}")
    print("odd_unit_ladders_J3..31_clock0..8=ranks[1,2,3,3,3,3,3,3,3]=PASS")
    print(f"GL2_F2_size=6 order_profile={order_profile} pointed_line_stabilizer={stabilizer_l}")
    print("V4_translation_directions=three_K4_perfect_matchings=PASS")
    print(f"positive_lift_T2={matpow(t_plus, 2)} signed_lift_T3={matpow(t_signed, 3)}")
    print("entrywise_nonnegative_integral_projective_order3_lift=IMPOSSIBLE")
    print("ALL THM-2987 CANDIDATE REFEREE CHECKS PASSED")


if __name__ == "__main__":
    main()
