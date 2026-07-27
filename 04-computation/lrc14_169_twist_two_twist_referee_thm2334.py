#!/usr/bin/env python3
"""Independent deterministic referee for the THM-2334 169-twist certificate.

This file deliberately does not import the original 169-twist computation.
It checks only the trivial character and

    v2 = -e_(q2) + e_(c3),

under two exact embeddings of Z[zeta_NN] into finite fields.  Its interval
engine is also independent: it refines half-open intervals at exact comb
boundaries, tests each resulting cell at an exact midpoint, and intersects
E^ell directly with all 169 pullbacks of the fixed word interval set.

The certified consequence is deliberately narrow.  On the explicit typed
but non-cover word THM-2309 (25), owner c1, word {a}, clock k=2, and triangle
(X,m,Y)=(13,1,742599), H(v2) differs from H(0).  Hence the 169-bank is
nonconstant and THM-2334 (42)--(43) gives positive unrestricted mod-13
nonzero-target energy.  No all-91-unit, address/gain, terminal-phase, scalar
row, or LRC(14) conclusion is checked here.
"""

from fractions import Fraction
from math import gcd


SCRIPT = "04-computation/lrc14_169_twist_two_twist_referee_thm2334.py"
OUTPUT = "05-knowledge/results/lrc14_169_twist_two_twist_referee_thm2334.out"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def product_factorization(factors):
    value = 1
    for prime, exponent in factors.items():
        value *= prime**exponent
    return value


def small_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def verify_lucas_certificate(prime, factors, witnesses, label):
    """Lucas primality certificate using the complete factorization of p-1."""
    require(product_factorization(factors) == prime - 1,
            f"{label}: factorization of p-1 is incomplete")
    require(set(factors) == set(witnesses),
            f"{label}: witness set does not match prime divisors of p-1")
    for q in factors:
        require(small_prime(q), f"{label}: factor {q} is not prime")
        a = witnesses[q]
        require(pow(a, prime - 1, prime) == 1,
                f"{label}: Fermat condition failed at q={q}")
        require(gcd(pow(a, (prime - 1) // q, prime) - 1, prime) == 1,
                f"{label}: Lucas gcd condition failed at q={q}")


def rank_mod_13(vectors):
    matrix = [[entry % 13 for entry in row] for row in vectors]
    rank = 0
    for column in range(9):
        pivot = next((r for r in range(rank, len(matrix))
                      if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], 11, 13)
        matrix[rank] = [(inverse * x) % 13 for x in matrix[rank]]
        for r in range(len(matrix)):
            if r == rank or matrix[r][column] == 0:
                continue
            multiplier = matrix[r][column]
            matrix[r] = [(matrix[r][c] - multiplier * matrix[rank][c]) % 13
                         for c in range(9)]
        rank += 1
    return rank


def dot_mod_13(left, right):
    return sum(a * b for a, b in zip(left, right)) % 13


# Typed word and marked triangle.
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD = 0
OWNER = 6
TARGET_A = 7
TARGET_B = 8
UNITS = (1, 2, 3, 4, 5)
OMITTED_UNIT = 5
R = 13**2
X = 13
M = 1
Y = X + M * W[TARGET_B]


def valuation_13(n):
    value = 0
    while n % 13 == 0:
        value += 1
        n //= 13
    return value


LCM_W = 1
for speed in W:
    LCM_W = LCM_W * speed // gcd(LCM_W, speed)
T_DEN = 182 * LCM_W
NN = R * T_DEN
NN_FACTORS = {2: 4, 3: 3, 5: 1, 7: 2, 11: 1, 13: 8, 53: 1}

# The two exact embeddings used by the original computation.
P1 = 352341050142921841
H1 = 435817657216
P2 = 956354278959359281
H2 = 153943385426666320
EMBEDDINGS = ((P1, H1), (P2, H2))

P1_FACTORS = {2: 4, 3: 3, 5: 1, 7: 3, 11: 1, 13: 8, 53: 1}
P1_WITNESSES = {2: 23, 3: 2, 5: 3, 7: 2, 11: 2, 13: 2, 53: 3}
P2_FACTORS = {2: 4, 3: 3, 5: 1, 7: 2, 11: 1, 13: 8, 19: 1, 53: 1}
P2_WITNESSES = {2: 23, 3: 2, 5: 2, 7: 2, 11: 2, 13: 2, 19: 2, 53: 2}


def verify_embedding(prime, root, label):
    require(pow(root, NN, prime) == 1,
            f"{label}: proposed root does not have order dividing NN")
    for q in NN_FACTORS:
        require(pow(root, NN // q, prime) != 1,
                f"{label}: proposed root order loses prime factor {q}")


def push_interval(output, left, right):
    if left >= right:
        return
    if output and output[-1][1] == left:
        output[-1] = (output[-1][0], right)
    else:
        output.append((left, right))


def boundary_points(left, right, speed, shift, width_denominator):
    """Comb boundaries strictly inside [left,right] on the T_DEN grid."""
    period_denominator = 13 * width_denominator
    require(T_DEN % (period_denominator * speed) == 0,
            "comb boundary is not integral on T_DEN grid")
    unit = T_DEN // (period_denominator * speed)
    step = period_denominator * unit
    points = []
    for sign in (-1, 1):
        base = (sign * 13 - width_denominator * shift) * unit
        first_n = (left - base) // step + 1
        point = base + first_n * step
        while point < right:
            points.append(point)
            point += step
    return sorted(set(points))


def midpoint_inside_comb(left, right, speed, shift, width_denominator):
    """Exact test of ||speed*t+shift/13|| < 1/width_denominator."""
    midpoint_sum = left + right
    total_denominator = 26 * T_DEN
    numerator = (13 * speed * midpoint_sum + 2 * T_DEN * shift) % total_denominator
    distance = min(numerator, total_denominator - numerator)
    return width_denominator * distance < total_denominator


def refine(intervals, speed, shift, width_denominator, want_inside):
    output = []
    for left, right in intervals:
        cuts = [left]
        cuts.extend(boundary_points(left, right, speed, shift, width_denominator))
        cuts.append(right)
        for a, b in zip(cuts, cuts[1:]):
            inside = midpoint_inside_comb(a, b, speed, shift, width_denominator)
            if inside == want_inside:
                push_interval(output, a, b)
    return output


def validate_intervals(intervals, label):
    previous_right = -1
    for left, right in intervals:
        require(0 <= left < right <= T_DEN, f"{label}: invalid endpoint")
        require(left > previous_right, f"{label}: overlap or unmerged adjacency")
        previous_right = right
    return sum(right - left for left, right in intervals)


PATTERN_E = {
    GUARD: "guard_safe",
    1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    OWNER: "in", TARGET_A: "out", TARGET_B: "out",
}
PATTERN_QA = {
    GUARD: "guard_safe",
    1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    OWNER: "out", TARGET_A: "in", TARGET_B: "out",
}


def build_boolean_set(pattern, shift):
    # Apply positive danger constraints first, then the guard, then exclusions.
    positive = sorted((W[i], i) for i, mode in pattern.items() if mode == "in")
    require(positive, "Boolean set must have an anchoring positive danger comb")
    intervals = [(0, T_DEN)]
    for _, i in positive:
        intervals = refine(intervals, W[i], shift[i], 14, True)
    if pattern.get(GUARD) == "guard_safe":
        intervals = refine(intervals, W[GUARD], shift[GUARD], 7, False)
    exclusions = sorted((W[i], i) for i, mode in pattern.items() if mode == "out")
    for _, i in exclusions:
        intervals = refine(intervals, W[i], shift[i], 14, False)
    validate_intervals(intervals, "built Boolean set")
    return intervals


def endpoint_sum_on_t_den(intervals, frequency):
    """Return (2*pi*i*frequency)*hat(1_intervals)(frequency) in both fields."""
    totals = [0, 0]
    for left, right in intervals:
        exponent_left = (-frequency * R * left) % NN
        exponent_right = (-frequency * R * right) % NN
        for j, (prime, root) in enumerate(EMBEDDINGS):
            totals[j] = (totals[j] + pow(root, exponent_left, prime)
                         - pow(root, exponent_right, prime)) % prime
    return tuple(totals)


def word_preimages(word_intervals):
    """Intervals of {t: R*t mod 1 lies in Q}, on the NN=R*T_DEN grid."""
    for branch in range(R):
        offset = branch * T_DEN
        for left, right in word_intervals:
            yield offset + left, offset + right


def marked_endpoint_sum(present_intervals, word_intervals, frequency):
    """Direct two-pointer intersection on the NN grid, independently of x_sweep."""
    scaled_present = [(R * left, R * right) for left, right in present_intervals]
    present_index = 0
    totals = [0, 0]
    total_length = 0
    component_count = 0
    for q_left, q_right in word_preimages(word_intervals):
        while (present_index < len(scaled_present)
               and scaled_present[present_index][1] <= q_left):
            present_index += 1
        scan = present_index
        while scan < len(scaled_present) and scaled_present[scan][0] < q_right:
            e_left, e_right = scaled_present[scan]
            left = max(e_left, q_left)
            right = min(e_right, q_right)
            if left < right:
                total_length += right - left
                component_count += 1
                exponent_left = (-frequency * left) % NN
                exponent_right = (-frequency * right) % NN
                for j, (prime, root) in enumerate(EMBEDDINGS):
                    totals[j] = (totals[j] + pow(root, exponent_left, prime)
                                 - pow(root, exponent_right, prime)) % prime
            scan += 1
    return tuple(totals), total_length, component_count


def make_owner_packet_rows():
    pivot_labels = (OWNER, GUARD, 1, 2, 3, 4)
    rows = []
    for label in pivot_labels:
        row = [0] * 9
        row[OMITTED_UNIT] += W[label]
        row[label] -= W[OMITTED_UNIT]
        if label == 1:
            row[OMITTED_UNIT] += W[TARGET_A]
            row[TARGET_A] -= W[OMITTED_UNIT]
        if label == 2:
            row[OMITTED_UNIT] += W[TARGET_B]
            row[TARGET_B] -= W[OMITTED_UNIT]
        rows.append(tuple(x % 13 for x in row))
    return rows


def main():
    print("THM-2334 independent deterministic two-twist referee")
    print(f"script={SCRIPT}")
    print(f"stored_output={OUTPUT}")
    print("")

    # Row and triangle typing.
    require(gcd(*W) == 1, "typed word is not primitive")
    require(W[GUARD] % 2 == 1 and W[GUARD] % 13 != 0,
            "guard is not an odd 13-unit")
    require(len({W[i] for i in UNITS}) == 5, "unit speeds are not distinct")
    require(all(W[i] % 13 != 0 for i in UNITS), "a unit speed is divisible by 13")
    profile = tuple(valuation_13(W[i]) for i in (OWNER, TARGET_A, TARGET_B))
    require(profile == (1, 3, 5), "wrong strict valuation profile")
    require(R == 13 ** (profile[0] + 1), "word clock is not k=lambda_owner+1")
    require(Y == X + M * W[TARGET_B], "frequency triangle does not close")
    require(valuation_13(X) == valuation_13(Y) == 1,
            "triangle endpoints do not have the same shallow root colour")
    require(gcd(M, 91) == 1, "deep multiplier is not a 91-unit")
    require(product_factorization(NN_FACTORS) == NN, "wrong NN factorization")
    require(all(T_DEN % (182 * speed) == 0 for speed in W),
            "T_DEN misses a comb endpoint denominator")
    print("[1] typed row and triangle controls: PASS")
    print(f"    w={W}")
    print(f"    profile={profile}  owner=c1  word={{a}}  k=2  R={R}")
    print(f"    (X,m,Y)=({X},{M},{Y})")
    print(f"    T_DEN={T_DEN}  NN={NN}")
    print("")

    # Deterministic primality and exact-order certificates.
    verify_lucas_certificate(P1, P1_FACTORS, P1_WITNESSES, "p1")
    verify_lucas_certificate(P2, P2_FACTORS, P2_WITNESSES, "p2")
    verify_embedding(P1, H1, "embedding 1")
    verify_embedding(P2, H2, "embedding 2")
    print("[2] Lucas primality and exact-order embedding controls: PASS")
    print(f"    p1={P1}  h1={H1}")
    print(f"    p2={P2}  h2={H2}")
    print("")

    # Exact quotient typing.
    rows = make_owner_packet_rows()
    w_mod = tuple(speed % 13 for speed in W)
    require(all(dot_mod_13(row, w_mod) == 0 for row in rows),
            "owner-packet row leaves K_13")
    require(rank_mod_13(rows) == 6, "owner packet does not have rank six")
    v1 = (0, 12, 0, 0, 0, 0, 0, 1, 0)
    v2 = (0, 0, 12, 0, 0, 0, 0, 0, 1)
    require(all(dot_mod_13(row, vector) == 0 for row in rows
                for vector in (w_mod, v1, v2)),
            "claimed character generator is not in L-perp")
    require(rank_mod_13((w_mod, v1, v2)) == 3,
            "w,v1,v2 do not span the three-dimensional L-perp")
    quotient_size = 13 ** ((9 - rank_mod_13(rows)) - 1)
    require(quotient_size == 169, "target character quotient does not have size 169")
    zero = (0,) * 9
    require(v2 == tuple(((-1 if i == 2 else 1 if i == TARGET_B else 0) % 13)
                        for i in range(9)),
            "v2 is not -e_q2+e_c3")
    for ell in (zero, v2):
        require(all((R * coordinate) % 13 == 0 for coordinate in ell),
                "word translation is not target-neutral")
    print("[3] target quotient controls: PASS")
    print(f"    rank(L)=6  |L_perp/<w>|={quotient_size}")
    print(f"    v1={v1}")
    print(f"    v2={v2}=-e_q2+e_c3")
    print("")

    # Independent exact Boolean geometry.
    q_a = build_boolean_set(PATTERN_QA, zero)
    e_zero = build_boolean_set(PATTERN_E, zero)
    e_v2 = build_boolean_set(PATTERN_E, v2)
    q_measure = Fraction(validate_intervals(q_a, "Q_a"), T_DEN)
    e_measure = Fraction(validate_intervals(e_zero, "E_0"), T_DEN)
    require(e_measure == Fraction(1882176, 28589561), "E_1 measure control failed")
    require(q_measure == Fraction(143103830843, 5727632650740),
            "Q_{1,{a}} measure control failed")
    print("[4] independent exact interval refinement: PASS")
    print(f"    E(0): intervals={len(e_zero)}  measure={e_measure}")
    print(f"    E(v2): intervals={len(e_v2)}")
    print(f"    Q_{{1,{{a}}}}: intervals={len(q_a)}  measure={q_measure}")
    print("")

    # Direct marked/bare coefficients at the two twists.
    ax_zero, stratum_length, stratum_components = marked_endpoint_sum(e_zero, q_a, X)
    ax_v2, _, v2_components = marked_endpoint_sum(e_v2, q_a, X)
    by_zero = endpoint_sum_on_t_den(e_zero, -Y)
    by_v2 = endpoint_sum_on_t_den(e_v2, -Y)
    stratum_measure = Fraction(stratum_length, NN)
    require(stratum_measure == Fraction(21376087, 17907461390),
            "delayed {a}-word stratum measure control failed")
    require(ax_zero == (117079706151067649, 1326261470990946),
            "marked X coefficient disagrees with the independent control")
    require(all(value != 0 for value in ax_zero + ax_v2 + by_zero + by_v2),
            "one required endpoint coefficient vanishes in a certified embedding")
    require(M % 7 != 0, "deep danger-comb coefficient vanishes")
    print("[5] direct 169-fold word-preimage intersection: PASS")
    print(f"    measure(E(0) intersect T^-2 Q_{{1,{{a}}}})={stratum_measure}")
    print(f"    components: ell=0 -> {stratum_components}, ell=v2 -> {v2_components}")
    print(f"    AX(0) images={ax_zero}")
    print(f"    AX(v2) images={ax_v2}")
    print(f"    conjugate-BY(0) images={by_zero}")
    print(f"    conjugate-BY(v2) images={by_v2}")
    print("")

    gamma_zero = []
    gamma_v2 = []
    differences = []
    for j, (prime, root) in enumerate(EMBEDDINGS):
        phase_v2 = pow(root, (M * v2[TARGET_B] * (NN // 13)) % NN, prime)
        g0 = ax_zero[j] * by_zero[j] % prime
        gv2 = phase_v2 * ax_v2[j] % prime * by_v2[j] % prime
        difference = (gv2 - g0) % prime
        require(difference != 0, f"embedding {j + 1} does not distinguish v2 from zero")
        gamma_zero.append(g0)
        gamma_v2.append(gv2)
        differences.append(difference)
    require(gamma_zero[0] == 310354333794505177,
            "p1 trivial-twist gamma control failed")
    require(gamma_v2[0] == 174331800739176126,
            "p1 v2 gamma control failed")
    print("[6] exact two-embedding nonconstancy certificate: PASS")
    print(f"    p1: gamma(0)={gamma_zero[0]}  gamma(v2)={gamma_v2[0]}  diff={differences[0]}")
    print(f"    p2: gamma(0)={gamma_zero[1]}  gamma(v2)={gamma_v2[1]}  diff={differences[1]}")
    print("")
    print("EXACT CLAIM: H(v2) != H(0) for the stated typed row/word/clock/triangle.")
    print("THM-2334 (42)-(43) therefore gives sum_(q!=0)|A(q)|^2 > 0.")
    print("SCOPE: unrestricted mod-13 target aggregate on a typed non-cover word only;")
    print("no all-91-unit, address/gain, terminal-phase, scalar-row, or LRC closure.")
    print("ALL DETERMINISTIC TWO-TWIST REFEREE CHECKS PASSED")


if __name__ == "__main__":
    main()
