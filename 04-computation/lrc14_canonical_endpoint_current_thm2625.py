#!/usr/bin/env python3
"""Exact dual-field certificate for THM-2625.

This is a dependency-free reconstruction of one canonical THM-2334 marked
current before its endpoint-difference/Radon quotient.  For the typed control

    W=(1,14,27,40,53,66,13,13^3,2*13^5), sigma={a}, R=13^2,
    (X,m,Y)=(13,1,742599),

it computes the separately allocated left and right endpoint transforms

    L*(l)=sum_ell P_ell z13^(-<tau(ell),l>),
    R*(r)=sum_ell Q_ell z13^(+<tau(ell),r>),

and the determinant-sector pushforward

    S*(q,Delta)=sum_{det(q,r)=Delta} L*(r+q) R*(r).

The PLUS sign in R* is load-bearing: Q is the conjugate endpoint factor.  All
quantities are cyclotomic integers in Z[zeta_N], N=169*T_DEN.  Each exact
specialization zeta_N -> h in F_p is a ring homomorphism because p is proved
prime and h is checked to have exact order N.  A nonzero specialized value
therefore proves that the characteristic-zero cyclotomic integer is nonzero.

Two independent specializations are run.  No floating-point arithmetic enters
the construction or decision.  Every check uses require(), so ``python -O``
cannot erase the validity gate.

Script: 04-computation/lrc14_canonical_endpoint_current_thm2625.py
Output: 05-knowledge/results/lrc14_canonical_endpoint_current_thm2625.out
"""

from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# Canonical typed row and exact scales
# ---------------------------------------------------------------------------
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
UNIT_IDX = (1, 2, 3, 4, 5)
K_CLOCK = 2
RDIL = 13**K_CLOCK
X, M_DEEP, Y = 13, 1, 742599


def nu13(n):
    value = 0
    while n % 13 == 0:
        n //= 13
        value += 1
    return value


require(W[GUARD] % 2 == 1 and W[GUARD] % 13, "guard type")
require(all(W[i] % 13 for i in UNIT_IDX), "unit type")
require(len({W[i] for i in UNIT_IDX}) == 5, "distinct units")
PROFILE = (nu13(W[OWNER]), nu13(W[TA]), nu13(W[TB]))
require(PROFILE == (1, 3, 5), "valuation profile")
require(Y == X + M_DEEP * W[TB], "marked triangle")
require(nu13(X) == nu13(Y) == 1 and gcd(M_DEEP, 91) == 1, "triangle type")

LCM_W = 1
for value in W:
    LCM_W = LCM_W * value // gcd(LCM_W, value)
T_DEN = 182 * LCM_W
NN = RDIL * T_DEN
require(T_DEN == 297836897838480, "T_DEN drift")
require(NN == 50334435734703120, "cyclotomic order drift")
require(all(T_DEN % (182 * value) == 0 for value in W), "breakpoint scale")

NN_PRIMES = []
remaining = NN
for prime in (2, 3, 5, 7, 11, 13, 53):
    if remaining % prime == 0:
        NN_PRIMES.append(prime)
        while remaining % prime == 0:
            remaining //= prime
require(remaining == 1, "incomplete factorization of N")


# ---------------------------------------------------------------------------
# Lucas primality/order gate for the two exact specializations
# ---------------------------------------------------------------------------
def small_prime(n):
    if n < 2:
        return False
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            return False
        divisor += 1
    return True


def product_factorization(factors):
    value = 1
    for prime, exponent in factors.items():
        value *= prime**exponent
    return value


def verify_lucas_certificate(prime, factors, witnesses, label):
    """Lucas primality theorem with the complete factorization of p-1."""
    require(product_factorization(factors) == prime - 1, f"{label}: p-1 factorization")
    require(set(factors) == set(witnesses), f"{label}: Lucas witness set")
    for factor in factors:
        require(small_prime(factor), f"{label}: composite factor in p-1")
        witness = witnesses[factor]
        require(pow(witness, prime - 1, prime) == 1, f"{label}: Fermat condition")
        require(
            gcd(pow(witness, (prime - 1) // factor, prime) - 1, prime) == 1,
            f"{label}: Lucas gcd condition at {factor}",
        )


EMBEDDINGS = (
    (
        7,
        352341050142921841,
        435817657216,
        {2: 4, 3: 3, 5: 1, 7: 3, 11: 1, 13: 8, 53: 1},
        {2: 23, 3: 2, 5: 3, 7: 2, 11: 2, 13: 2, 53: 3},
    ),
    (
        19,
        956354278959359281,
        153943385426666320,
        {2: 4, 3: 3, 5: 1, 7: 2, 11: 1, 13: 8, 19: 1, 53: 1},
        {2: 23, 3: 2, 5: 2, 7: 2, 11: 2, 13: 2, 19: 2, 53: 2},
    ),
)
for embedding_index, (multiplier, prime, root, factors, witnesses) in enumerate(EMBEDDINGS):
    require(prime == multiplier * NN + 1, "p != c*N+1")
    verify_lucas_certificate(prime, factors, witnesses, f"p{embedding_index + 1}")
    require(pow(root, NN, prime) == 1, "root does not have N-power one")
    require(
        all(pow(root, NN // factor, prime) != 1 for factor in NN_PRIMES),
        "root does not have exact order N",
    )

MODS = tuple((prime, root) for _, prime, root, _, _ in EMBEDDINGS)


# ---------------------------------------------------------------------------
# Exact Boolean interval machinery on the T_DEN lattice
# ---------------------------------------------------------------------------
def in_comb(index, ell_i):
    speed = W[index]
    unit = T_DEN // (182 * speed)
    low = (-13 - 14 * ell_i) % 182
    intervals = []
    for number in range(speed):
        start = (low + 182 * number) * unit
        stop = start + 26 * unit
        if stop <= T_DEN:
            intervals.append((start, stop))
        else:
            intervals.append((start, T_DEN))
            intervals.append((0, stop - T_DEN))
    intervals.sort()
    return intervals


def subtract_comb(intervals, speed, period_den, low, high):
    unit = T_DEN // (period_den * speed)
    step = period_den * unit
    window_length = (high - low) * unit
    base = (low % period_den) * unit
    output = []
    for start, stop in intervals:
        first_possible = start - window_length + 1
        number = -((base - first_possible) // step)
        position = base + number * step
        cursor = start
        while position < stop:
            window_stop = position + window_length
            if window_stop > cursor:
                if position > cursor:
                    output.append((cursor, position))
                cursor = window_stop
                if cursor >= stop:
                    break
            position += step
        if cursor < stop:
            output.append((cursor, stop))
    return output


def build_set(pattern, ell):
    inside = [index for index, mode in pattern.items() if mode == "in"]
    start_index = min(inside, key=lambda index: W[index])
    intervals = in_comb(start_index, ell[start_index])
    for index, mode in pattern.items():
        if mode == "gout":
            intervals = subtract_comb(
                intervals, W[index], 91, -13 - 7 * ell[index], 13 - 7 * ell[index]
            )
    rest = sorted(
        (W[index], index)
        for index, mode in pattern.items()
        if index != start_index and mode in ("in", "out")
    )
    for _, index in rest:
        residue = ell[index]
        if pattern[index] == "out":
            intervals = subtract_comb(
                intervals, W[index], 182, -13 - 14 * residue, 13 - 14 * residue
            )
        else:
            intervals = subtract_comb(
                intervals, W[index], 182, 13 - 14 * residue, 169 - 14 * residue
            )
    return intervals


def interval_length(intervals):
    last = -1
    total = 0
    for start, stop in intervals:
        require(0 <= start < stop <= T_DEN and start >= last, "corrupt interval list")
        last = stop
        total += stop - start
    return total


PAT_E = {
    0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    6: "in", 7: "out", 8: "out",
}
PAT_QA = {
    0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    6: "out", 7: "in", 8: "out",
}
ZERO_ELL = (0,) * 9


# ---------------------------------------------------------------------------
# Exact endpoint Fourier sums
# ---------------------------------------------------------------------------
def make_tabs(q_intervals, frequency, mods):
    tables = []
    for prime, root in mods:
        lower = [pow(root, (-frequency * start) % NN, prime) for start, _ in q_intervals]
        upper = [pow(root, (-frequency * stop) % NN, prime) for _, stop in q_intervals]
        tables.append((lower, upper))
    return tables


def x_sweep(e_intervals, q_intervals, q_starts, frequency, mods, tables):
    """Return endpoint numerator of E^ell intersect T^-2 Q and overlap length."""
    count_q = len(q_intervals)
    sums = [0] * len(mods)
    overlap = 0
    wrap_forward = [
        pow(root, (-frequency * T_DEN) % NN, prime) for prime, root in mods
    ]
    wrap_backward = [
        pow(root, (frequency * T_DEN) % NN, prime) for prime, root in mods
    ]
    for a, b in e_intervals:
        lifted_a = RDIL * a
        source_a = lifted_a % T_DEN
        span = RDIL * (b - a)
        require(span < T_DEN, "E interval too long for one-wrap sweep")
        source_b = source_a + span
        index = bisect_right(q_starts, source_a) - 1
        offset = 0
        if index < 0:
            index = count_q - 1
            offset = -T_DEN
        base_exp = (-frequency * (lifted_a - source_a)) % NN
        base_values = [pow(root, base_exp, prime) for prime, root in mods]
        accumulators = [0] * len(mods)
        wrap = [1] * len(mods)
        if offset:
            wrap = list(wrap_backward)
        while True:
            q_a0, q_b0 = q_intervals[index]
            q_a, q_b = q_a0 + offset, q_b0 + offset
            if q_a >= source_b:
                break
            if q_b > source_a:
                low = max(source_a, q_a)
                high = min(source_b, q_b)
                if high > low:
                    overlap += high - low
                    for field_index, (prime, root) in enumerate(mods):
                        low_value = (
                            tables[field_index][0][index] * wrap[field_index] % prime
                            if low == q_a
                            else pow(root, (-frequency * low) % NN, prime)
                        )
                        high_value = (
                            tables[field_index][1][index] * wrap[field_index] % prime
                            if high == q_b
                            else pow(root, (-frequency * high) % NN, prime)
                        )
                        accumulators[field_index] = (
                            accumulators[field_index] + low_value - high_value
                        ) % prime
            index += 1
            if index == count_q:
                index = 0
                offset += T_DEN
                for field_index, (prime, _) in enumerate(mods):
                    wrap[field_index] = (
                        wrap[field_index] * wrap_forward[field_index] % prime
                    )
        for field_index, (prime, _) in enumerate(mods):
            sums[field_index] = (
                sums[field_index] + base_values[field_index] * accumulators[field_index]
            ) % prime
    return sums, overlap


def endpoint_sum(intervals, frequency, mods):
    sums = [0] * len(mods)
    for start, stop in intervals:
        start_exp = (-frequency * RDIL * start) % NN
        stop_exp = (-frequency * RDIL * stop) % NN
        for field_index, (prime, root) in enumerate(mods):
            sums[field_index] = (
                sums[field_index]
                + pow(root, start_exp, prime)
                - pow(root, stop_exp, prime)
            ) % prime
    return sums


# ---------------------------------------------------------------------------
# Relation quotient and endpoint coordinates
# ---------------------------------------------------------------------------
def gf13_rank(vectors):
    matrix = [list(vector) for vector in vectors]
    rank = 0
    for column in range(9):
        selected = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column] % 13),
            None,
        )
        if selected is None:
            continue
        matrix[rank], matrix[selected] = matrix[selected], matrix[rank]
        inverse = pow(matrix[rank][column], 11, 13)
        matrix[rank] = [(value * inverse) % 13 for value in matrix[rank]]
        for row in range(len(matrix)):
            if row != rank and matrix[row][column] % 13:
                factor = matrix[row][column]
                matrix[row] = [
                    (matrix[row][index] - factor * matrix[rank][index]) % 13
                    for index in range(9)
                ]
        rank += 1
    return rank


U0 = 5
PIVOTS = (OWNER, 0, 1, 2, 3, 4)
RELATION_ROWS = []
for pivot in PIVOTS:
    row = [0] * 9
    row[U0] = (row[U0] + W[pivot]) % 13
    row[pivot] = (row[pivot] - W[U0]) % 13
    RELATION_ROWS.append(row)
graft_a = RELATION_ROWS[PIVOTS.index(1)]
graft_a[U0] = (graft_a[U0] + W[TA]) % 13
graft_a[TA] = (graft_a[TA] - W[U0]) % 13
graft_b = RELATION_ROWS[PIVOTS.index(2)]
graft_b[U0] = (graft_b[U0] + W[TB]) % 13
graft_b[TB] = (graft_b[TB] - W[U0]) % 13

WMOD = tuple(value % 13 for value in W)
V1 = (0, 12, 0, 0, 0, 0, 0, 1, 0)
V2 = (0, 0, 12, 0, 0, 0, 0, 0, 1)
require(gf13_rank(RELATION_ROWS) == 6, "relation rank")
for vector in (WMOD, V1, V2):
    require(
        all(sum(a * b for a, b in zip(row, vector)) % 13 == 0 for row in RELATION_ROWS),
        "purported quotient vector is not in L-perp",
    )
require(gf13_rank((WMOD, V1, V2)) == 3, "quotient basis degeneracy")

REPS = {
    (alpha, beta): tuple(
        (alpha * V1[index] + beta * V2[index]) % 13 for index in range(9)
    )
    for alpha in range(13)
    for beta in range(13)
}
require(len(set(REPS.values())) == 169, "representative collision")
TAU = {address: (ell[TA], ell[TB]) for address, ell in REPS.items()}
require(len(set(TAU.values())) == 169, "endpoint pairing degeneracy")
require(all(TAU[(a, b)] == (a, b) for a in range(13) for b in range(13)), "tau order")


def digest_bank(bank, keys):
    payload = ",".join(str(bank[key]) for key in keys).encode()
    return sha256(payload).hexdigest()


EXPECTED_DIGESTS = (
    (
        "33d77c87c2970bf99bd1f0f139e74b8eef70677ba28844391f6d681ebf2e6c69",
        "755b046e39d71fbe83295dd7e43727eaa30a5220f359e2d389268d124e6a8f27",
        "d7e20ba35d1b25c1761b03c297f2449d74b29be0defbc7496b10eff8b9458b41",
    ),
    (
        "97471c735336996fab845da55ebbca7a19786297ebe5dd9f7f5026a0a816c9ee",
        "07bc100112b3cf83cfddec6a6354df5caf7d5a12edd7f9d10c921ff08153271b",
        "bb09105d839e5abdfadf7b4e1959bcc8e5a4af680a8dbe21ec485ab7bb7ec8f2",
    ),
)


def main():
    print("THM-2625 exact canonical endpoint-current sector certificate")
    print(f"W={W}; profile={PROFILE}; sigma={{a}}; R={RDIL}")
    print(f"(X,m,Y)=({X},{M_DEEP},{Y}); T_DEN={T_DEN}; N={NN}")
    print("validity gate: Lucas primality + exact N-order checks PASS")

    for index in range(9):
        require(interval_length(in_comb(index, 0)) == T_DEN // 7, "comb measure")
    e_zero = build_set(PAT_E, ZERO_ELL)
    q_intervals = build_set(PAT_QA, ZERO_ELL)
    e_measure = Fraction(interval_length(e_zero), T_DEN)
    q_measure = Fraction(interval_length(q_intervals), T_DEN)
    require(e_measure > 0 and q_measure > 0, "base Boolean sets must be positive")
    q_starts = [start for start, _ in q_intervals]
    empty_tabs = make_tabs(q_intervals, 1, ())
    _, overlap = x_sweep(e_zero, q_intervals, q_starts, 1, (), empty_tabs)
    stratum_measure = Fraction(overlap, NN)
    require(stratum_measure > 0, "chosen delayed word stratum is empty")
    print(f"base measures: E={e_measure}; Q_a={q_measure}; E meet T^-2 Q_a={stratum_measure}")

    tabs = make_tabs(q_intervals, X, MODS)
    p_banks = [dict(), dict()]
    q_banks = [dict(), dict()]
    z13 = [pow(root, NN // 13, prime) for prime, root in MODS]

    for address, ell in REPS.items():
        e_intervals = build_set(PAT_E, ell)
        a_values, _ = x_sweep(e_intervals, q_intervals, q_starts, X, MODS, tabs)
        # frequency -Y gives sum e(+Y*t), the conjugate endpoint numerator.
        b_values = endpoint_sum(e_intervals, -Y, MODS)
        for field_index, (prime, _) in enumerate(MODS):
            phase = pow(z13[field_index], M_DEEP * ell[TB] % 13, prime)
            p_banks[field_index][address] = phase * a_values[field_index] % prime
            q_banks[field_index][address] = b_values[field_index] % prime

    # Separate representative-independence hostile control.  This is stronger
    # than checking only the product P_ell Q_ell.
    for address in ((0, 0), (1, 0), (0, 1), (3, 7)):
        ell = REPS[address]
        shifted = tuple((ell[i] + WMOD[i]) % 13 for i in range(9))
        shifted_intervals = build_set(PAT_E, shifted)
        a_values, _ = x_sweep(
            shifted_intervals, q_intervals, q_starts, X, MODS, tabs
        )
        b_values = endpoint_sum(shifted_intervals, -Y, MODS)
        for field_index, (prime, _) in enumerate(MODS):
            phase = pow(z13[field_index], M_DEEP * shifted[TB] % 13, prime)
            require(
                phase * a_values[field_index] % prime == p_banks[field_index][address],
                "left allocation is not representative-independent",
            )
            require(
                b_values[field_index] % prime == q_banks[field_index][address],
                "right allocation is not representative-independent",
            )
    print("separate P/Q quotient-gauge checks: PASS")

    xy_keys = [(x, y) for x in range(13) for y in range(13)]
    sector_keys = [
        (q0, q1, delta)
        for q0 in range(13)
        for q1 in range(13)
        for delta in range(13)
    ]

    for field_index, (prime, _) in enumerate(MODS):
        powers = [pow(z13[field_index], exponent, prime) for exponent in range(13)]
        p_bank, q_bank = p_banks[field_index], q_banks[field_index]
        left = {}
        right = {}
        aggregate = {}

        for point in xy_keys:
            lx = 0
            rx = 0
            for address in xy_keys:
                tau0, tau1 = TAU[address]
                pairing = (tau0 * point[0] + tau1 * point[1]) % 13
                lx = (lx + p_bank[address] * powers[-pairing % 13]) % prime
                rx = (rx + q_bank[address] * powers[pairing]) % prime
            left[point] = lx
            right[point] = rx

        require(sum(value != 0 for value in left.values()) == 169, "left support hole")
        require(sum(value != 0 for value in right.values()) == 169, "right support hole")

        # All inverse-transform controls, with signs reversed on the right.
        for address in xy_keys:
            tau0, tau1 = TAU[address]
            recover_left = 0
            recover_right = 0
            for point in xy_keys:
                pairing = (tau0 * point[0] + tau1 * point[1]) % 13
                recover_left = (recover_left + left[point] * powers[pairing]) % prime
                recover_right = (recover_right + right[point] * powers[-pairing % 13]) % prime
            require(recover_left == 169 * p_bank[address] % prime, "left DFT inversion")
            require(recover_right == 169 * q_bank[address] % prime, "right DFT inversion")

        anum = {}
        for q in xy_keys:
            total = 0
            for address in xy_keys:
                tau0, tau1 = TAU[address]
                pairing = (tau0 * q[0] + tau1 * q[1]) % 13
                total = (
                    total
                    + p_bank[address] * q_bank[address] * powers[-pairing % 13]
                ) % prime
            anum[q] = total
        require(sum(value != 0 for value in anum.values()) == 169, "target support hole")

        sectors = {key: 0 for key in sector_keys}
        joint_support = 0
        nondegenerate_edge_support = 0
        for q in xy_keys:
            for r in xy_keys:
                l = ((r[0] + q[0]) % 13, (r[1] + q[1]) % 13)
                coefficient = left[l] * right[r] % prime
                if coefficient:
                    joint_support += 1
                delta = (l[0] * r[1] - l[1] * r[0]) % 13
                if q != (0, 0) and delta != 0 and coefficient:
                    nondegenerate_edge_support += 1
                key = (q[0], q[1], delta)
                sectors[key] = (sectors[key] + coefficient) % prime

        require(joint_support == 13**4, "joint endpoint coefficient hole")
        require(nondegenerate_edge_support == 26208, "nondegenerate edge hole")
        for q in xy_keys:
            line_sum = sum(sectors[(q[0], q[1], delta)] for delta in range(13)) % prime
            require(line_sum == 169 * anum[q] % prime, "Radon line-sum identity")
        for delta in range(13):
            expected_nonzero = delta == 0
            require(
                (sectors[(0, 0, delta)] != 0) == expected_nonzero,
                "q=0 structural sector classification",
            )
        for q in xy_keys:
            if q == (0, 0):
                continue
            require(
                all(sectors[(q[0], q[1], delta)] != 0 for delta in range(13)),
                "admissible determinant-sector cancellation",
            )

        sector_support = sum(value != 0 for value in sectors.values())
        nondegenerate_sectors = sum(
            sectors[(q[0], q[1], delta)] != 0
            for q in xy_keys if q != (0, 0)
            for delta in range(1, 13)
        )
        require(sector_support == 2185, "sector support census")
        require(nondegenerate_sectors == 2016, "transvection-sector census")

        digests = (
            digest_bank(left, xy_keys),
            digest_bank(right, xy_keys),
            digest_bank(sectors, sector_keys),
        )
        require(digests == EXPECTED_DIGESTS[field_index], "bank digest mismatch")
        print(
            f"field {field_index + 1}: p={prime}; "
            "support(L,R,A,J,S,S_nondeg,J_nondeg)="
            f"(169,169,169,{joint_support},{sector_support},"
            f"{nondegenerate_sectors},{nondegenerate_edge_support})"
        )
        print(f"field {field_index + 1}: sha256(L,R,S)={digests}")

    print("structural zeros: exactly q=(0,0), Delta=1,...,12")
    print("all admissible sectors, all 2,016 nondegenerate sectors, and all their 13 edges survive")
    print("Radon control sum_Delta S*(q,Delta)=169*Anum(q): PASS in both fields")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
