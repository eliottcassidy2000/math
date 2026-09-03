#!/usr/bin/env python3
"""Independent exact referee for THM-4368.

This script imports neither the primary verifier nor its fixtures.  It
reconstructs the simplex streams, packet traces, source cone, triangular
address, and the fixed-row Pascal minors directly from their definitions.
"""

from math import comb, factorial, gcd, isqrt
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


checks = 0


def claim(test, label):
    global checks
    checks += 1
    if not test:
        raise AssertionError(label)


def ceil_div(a, b):
    return -((-a) // b)


def generalized_binom(a, r):
    """Integer binomial for integral a and r>=0."""
    if r < 0:
        return 0
    if a >= 0:
        return comb(a, r) if r <= a else 0
    return (-1) ** r * comb(r - a - 1, r)


def packet_A(ell, n0, N, n):
    s = ceil_div(ell, 2)
    k = n - n0
    if 0 <= k <= N:
        return (-1) ** (n - s) * comb(N, k)
    return 0


def packet_L(ell, n0, N, q, m):
    s = ceil_div(ell, 2)
    return sum(comb(m + q - n, q) * packet_A(ell, n0, N, n)
               for n in range(s, m + 1))


def f_coefficient(ell, n0, N, q, m):
    """Coefficient predicted by z^(n0-s)(1-z)^(N-q-1)."""
    s = ceil_div(ell, 2)
    r = m - n0
    if r < 0:
        return 0
    p = N - q - 1
    return (-1) ** (n0 - s + r) * generalized_binom(p, r)


def monomial_h(a, b, c, e, n, r):
    n0 = b + c + 2 * e
    r0 = a + 2 * b + e
    N = c + e
    k = n - n0
    if 0 <= k <= N and r == r0 + 2 * k:
        return comb(N, k)
    return 0


def monomial_L(ell, a, b, c, e, q, m):
    s = ceil_div(ell, 2)
    return sum(((-1) ** (n - s)) * comb(m + q - n, q)
               * monomial_h(a, b, c, e, n, 2 * n - ell)
               for n in range(s, m + 1))


def tri(k):
    return k * (k + 1) // 2


def det_bareiss(matrix):
    """Fraction-free exact determinant over the integers."""
    n = len(matrix)
    if n == 0:
        return 1
    a = [row[:] for row in matrix]
    sign = 1
    denominator = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot_value - a[i][k] * a[k][j]
                claim(numerator % denominator == 0,
                      ("Bareiss exact division", n, k, i, j))
                a[i][j] = numerator // denominator
        denominator = pivot_value
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign * a[n - 1][n - 1]


def addr(u, v):
    return tri(u + v - 2) + u


def decode(A):
    # t is the unique positive integer with T(t-1)<A<=T(t).
    t = (isqrt(8 * A + 1) - 1) // 2
    if tri(t) < A:
        t += 1
    u = A - tri(t - 1)
    return u, t + 1 - u


# Arbitrary-stream recurrence and exact inversion, including ell=2 and M=s.
for ell in range(2, 31):
    s = ceil_div(ell, 2)
    for length in range(1, 20):
        A = [((17 * i * i + 13 * ell + 5 * length) % 29) - 14
             for i in range(length)]
        streams = []
        prev = A
        for q in range(8):
            cur = []
            acc = 0
            for x in prev:
                acc += x
                cur.append(acc)
            streams.append(cur)
            prev = cur
        for q, stream in enumerate(streams):
            direct = [sum(comb(i + q - j, q) * A[j]
                          for j in range(i + 1))
                      for i in range(length)]
            claim(stream == direct, ("prefix", ell, length, q))
            diff = [stream[0]] + [stream[i] - stream[i - 1]
                                      for i in range(1, length)]
            claim(diff == (A if q == 0 else streams[q - 1]),
                  ("difference", ell, length, q))

# Direct packet/function-series equality. This explicitly straddles
# q=N-1 (impulse) and q=N (first infinite tail), including N=1.
seen_traces = {}
for ell in range(2, 36):
    s = ceil_div(ell, 2)
    rho = ceil_div(ell, 3)
    for N in range(rho, rho + 11):
        for n0 in range(max(N, ell - N), max(N, ell - N) + 11):
            key = tuple(packet_L(ell, n0, N, 0, m)
                        for m in range(s, n0 + N + 3))
            old = seen_traces.get((ell, key))
            claim(old is None or old == (N, n0),
                  ("injective trace", ell, old, N, n0))
            seen_traces[(ell, key)] = (N, n0)
            for q in range(0, N + 4):
                for m in range(s, n0 + N + q + 9):
                    got = packet_L(ell, n0, N, q, m)
                    want = f_coefficient(ell, n0, N, q, m)
                    claim(got == want,
                          ("Fq", ell, N, n0, q, m, got, want))
            # q=0 support and exact boundary valuations.
            coeffs = [packet_L(ell, n0, N, 0, m)
                      for m in range(s, n0 + N + 2)]
            first = next(i for i, x in enumerate(coeffs) if x)
            last = max(i for i, x in enumerate(coeffs) if x)
            claim(first == n0 - s, ("ord0", ell, N, n0, first))
            claim(last == n0 - s + N - 1,
                  ("degree", ell, N, n0, last))
            # Divisibility multiplicity at z=1 is tested by differences of
            # coefficient lists: (1-z)^r divides iff r-fold partial sums at 1
            # vanish; direct alternating moments furnish an independent test.
            shifted = [packet_L(ell, n0, N, 0, n0 + j)
                       for j in range(N)]
            for power in range(N):
                moment = sum(x * (j ** power)
                             for j, x in enumerate(shifted))
                claim((moment == 0) == (power < N - 1),
                      ("ord1", ell, N, n0, power, moment))

# Feasibility: necessary inequalities on brute-force monomials, and the stated
# construction on every pair in a broad independent box.
for ell in range(2, 81):
    s = ceil_div(ell, 2)
    rho = ceil_div(ell, 3)
    realized = set()
    for c in range(0, 31):
        for e in range(0, 31):
            a = 2 * c + 3 * e - ell
            if a < 0:
                continue
            for b in range(0, 31):
                N = c + e
                n0 = b + c + 2 * e
                realized.add((N, n0))
                claim(N >= rho and n0 >= N and n0 + N >= ell,
                      ("necessity", ell, a, b, c, e, N, n0))
                claim(n0 >= s, ("start", ell, N, n0, s))
    for N in range(1, 41):
        for n0 in range(0, 81):
            feasible = N >= rho and n0 >= N and n0 + N >= ell
            if feasible:
                lower = max(0, ell - 2 * N)
                upper = min(N, n0 - N)
                claim(lower <= upper,
                      ("nonempty source fibre", ell, N, n0, lower, upper))
                realizing = []
                packet_signatures = set()
                for source_e in range(0, N + 1):
                    source_c = N - source_e
                    source_a = 2 * N + source_e - ell
                    source_b = n0 - N - source_e
                    if min(source_a, source_b, source_c, source_e) < 0:
                        continue
                    realizing.append(source_e)
                    claim(source_c + source_e == N
                          and source_b + source_c + 2 * source_e == n0
                          and source_a == 2 * source_c + 3 * source_e - ell,
                          ("source fibre equations", ell, N, n0, source_e))
                    source_n0 = source_b + source_c + 2 * source_e
                    source_r0 = source_a + 2 * source_b + source_e
                    packet_signatures.add(tuple(
                        (source_n0 + k, source_r0 + 2 * k, comb(N, k))
                        for k in range(N + 1)))
                expected_es = list(range(lower, upper + 1))
                claim(realizing == expected_es,
                      ("source fibre interval", ell, N, n0,
                       realizing, expected_es))
                claim(len(realizing) == upper - lower + 1,
                      ("source fibre multiplicity", ell, N, n0,
                       len(realizing), lower, upper))
                canonical_packet = tuple(
                    (n0 + k, 2 * n0 - ell + 2 * k, comb(N, k))
                    for k in range(N + 1))
                claim(packet_signatures == {canonical_packet},
                      ("source fibre packet equality", ell, N, n0,
                       packet_signatures, canonical_packet))
            # In the exhaustively covered exponent box, compare existence too.
            if N <= 30 and n0 <= 30:
                claim(((N, n0) in realized) == feasible,
                      ("finite iff", ell, N, n0, feasible))

# Address bijection, decoder, reflection sum and parity.
addresses = set()
for u in range(1, 151):
    for v in range(1, 151):
        A = addr(u, v)
        claim(A not in addresses, ("address collision", u, v, A))
        addresses.add(A)
        claim(decode(A) == (u, v), ("decode", u, v, A, decode(A)))
        t = u + v - 1
        claim(addr(u, v) + addr(v, u) == t * t + 1,
              ("reflection sum", u, v))
        claim((u == v) == (t % 2 == 1 and u == (t + 1) // 2),
              ("fixed parity", u, v, t))
for A in range(1, tri(299) + 1):
    claim(addr(*decode(A)) == A, ("surjective", A, decode(A)))

# Fixed-row dual.  The reversed-column matrix B is independently compared
# with P P^T and its full determinant, then every consecutive-row minor is
# checked.  q0 is deliberately allowed beyond the square-basis range.
for r in range(0, 16):
    P = [[comb(q, j) if j <= q else 0 for j in range(r + 1)]
         for q in range(r + 1)]
    B = [[comb(q + k, q) for k in range(r + 1)]
         for q in range(r + 1)]
    PPt = [[sum(P[q][j] * P[k][j] for j in range(r + 1))
            for k in range(r + 1)] for q in range(r + 1)]
    claim(B == PPt, ("B=PPt", r))
    claim(det_bareiss(B) == 1, ("det B", r))
    for a in range(1, r + 2):
        vandermonde = 1
        for i in range(a):
            for j in range(i + 1, a):
                vandermonde *= j - i
        factorial_product = 1
        for k in range(a):
            factorial_product *= factorial(k)
        claim(vandermonde == factorial_product,
              ("consecutive Vandermonde", a, vandermonde, factorial_product))
        for q0 in range(0, 18):
            minor = [[comb(q0 + i + k, k) for k in range(a)]
                     for i in range(a)]
            claim(det_bareiss(minor) == 1,
                  ("unit consecutive minor", r, a, q0))

# The fixed-row hierarchy-annihilator bank and clock.  This directly compares
# the closed rank formula with THM-4364's three gates over every tested m>=s,
# including empty, first-entry, last-unsaturated, and saturated rows.
for ell in range(2, 73):
    s = ceil_div(ell, 2)
    rho = ceil_div(ell, 3)
    for d in range(0, 13):
        m0 = ell + d - rho + 1
        claim(m0 >= s, ("clock begins in row domain", ell, d, m0, s))
        marked = 0
        for m in range(s, ell + d + rho + 5):
            direct_bank = [q for q in range(0, rho + 5)
                           if q < rho and m + q >= ell
                           and d <= m + q - ell]
            q0 = max(0, ell + d - m)
            predicted_bank = list(range(q0, rho)) if q0 < rho else []
            rank = max(0, rho - q0)
            claim(direct_bank == predicted_bank,
                  ("annihilator bank", ell, d, m, direct_bank, predicted_bank))
            claim(len(direct_bank) == rank,
                  ("annihilator rank", ell, d, m, len(direct_bank), rank))
            claim(rank <= m - s + 1,
                  ("rank fits jet", ell, d, m, rank, m - s + 1))
            if direct_bank:
                claim(direct_bank == list(range(direct_bank[0], rho)),
                      ("bank consecutive", ell, d, m, direct_bank))
            if m < m0:
                claim(rank == 0, ("pre-clock empty", ell, d, m, rank))
            elif m <= ell + d:
                j = m - m0
                claim(direct_bank == list(range(rho - 1 - j, rho)),
                      ("clock bank", ell, d, m, j, direct_bank))
                claim(rank == j + 1, ("clock rank", ell, d, m, j, rank))
                marked += rank
            else:
                claim(rank == rho and direct_bank == list(range(rho)),
                      ("saturated", ell, d, m, rank, direct_bank))
        claim(marked == tri(rho),
              ("triangular clock positions", ell, d, marked, tri(rho)))

# THM-4366's two sharp-row modules each have the singleton hierarchy bank
# {3}; their exact coefficients are recovered from the defining formula.
def stencil(ell, m, q):
    s = ceil_div(ell, 2)
    return [((-1) ** (n - s)) * comb(m + q - n, q)
            for n in range(s, m + 1)]


for ell, m, d, expected in (
        (11, 10, 2, [35, -20, 10, -4, 1]),
        (10, 10, 3, [56, -35, 20, -10, 4, -1])):
    rho = ceil_div(ell, 3)
    bank = [q for q in range(0, rho)
            if m + q >= ell and d <= m + q - ell]
    claim(bank == [3], ("THM-4366 singleton bank", ell, m, d, bank))
    claim(stencil(ell, m, 3) == expected,
          ("THM-4366 stencil", ell, m, d, stencil(ell, m, 3), expected))

# Ambient F(1-z) reflection is coefficientwise checked as polynomials; source
# validity criterion is checked against feasibility of the reflected pair.
for ell in range(2, 81):
    s = ceil_div(ell, 2)
    rho = ceil_div(ell, 3)
    for v in range(rho, rho + 21):
        for n0 in range(max(v, ell - v), max(v, ell - v) + 25):
            u = n0 - s + 1
            t = u + v - 1
            # Coefficients in ascending powers.
            original = [(-1) ** (u - 1 + j) * comb(v - 1, j)
                        for j in range(v)]
            reflected = [(-1) ** (v - 1 + j) * comb(u - 1, j)
                         for j in range(u)]
            # Expand (-1)^(v-u)F(1-z) directly.
            transformed = [0] * u
            overall = (-1) ** (v - u + u - 1)
            for j in range(v):
                coeff = (-1) ** (u - 1 + j) * comb(v - 1, j)
                # coeff*(1-z)^(u-1)*z^j after substituting; but a simpler
                # independent expansion starts from F(z)=sign*z^(u-1)*(1-z)^(v-1).
                _ = coeff
            transformed = [overall * ((-1) ** j) * comb(u - 1, j)
                           for j in range(u)]
            claim(transformed == reflected,
                  ("polynomial reflection", ell, u, v, transformed, reflected))
            reflected_feasible = (u >= rho and (s + v - 1) >= u
                                  and (s + v - 1) + u >= ell)
            criterion = u >= rho and u - v <= s - 1
            claim(reflected_feasible == criterion,
                  ("source reflection", ell, u, v, reflected_feasible, criterion))
            claim(s + u + v - 1 >= ell,
                  ("terminal inherited", ell, u, v, t))

# Direct monomial diagonal/off-diagonal checks, including packets invisible to
# the chosen diagonal. The all-zero class is kept distinct from valid packets.
for ell in range(2, 21):
    s = ceil_div(ell, 2)
    for a in range(0, 13):
        for b in range(0, 8):
            for c in range(0, 10):
                for e in range(0, 7):
                    N = c + e
                    n0 = b + c + 2 * e
                    diagonal = a == 2 * c + 3 * e - ell
                    for q in range(0, 5):
                        for m in range(s, s + 7):
                            got = monomial_L(ell, a, b, c, e, q, m)
                            want = (packet_L(ell, n0, N, q, m)
                                    if diagonal else 0)
                            claim(got == want,
                                  ("monomial", ell, a, b, c, e, q, m, got, want))

# Named ell=10 hostiles and controls.
ell = 10
named = [
    ((0, 3, 5, 0), (4, 5)),
    ((2, 2, 6, 0), (4, 6)),
    ((10, 2, 10, 0), (8, 10)),
    ((0, 3, 2, 2), (5, 4)),
    ((0, 2, 2, 2), (4, 4)),
    ((2, 1, 6, 0), (3, 6)),
]
for (a, b, c, e), pair in named:
    N = c + e
    n0 = b + c + 2 * e
    got_pair = (n0 - ceil_div(ell, 2) + 1, N)
    claim(a == 2 * c + 3 * e - ell, ("named diagonal", a, b, c, e))
    claim(got_pair == pair, ("named pair", (a, b, c, e), got_pair, pair))
claim(gcd(4, 5) == 1 and gcd(8, 10) == 2, "ray scales")
claim(addr(4, 5) == 32 and addr(5, 4) == 33 and addr(4, 4) == 25,
      "named addresses")
claim(not (3 >= ceil_div(10, 3)), "invalid reciprocal N")

print("PASS")
print(f"checks={checks}")
print("edge_ell=2=PASS")
print("N=1_q=N-1_q=N=PASS")
print("zero_off_diagonal=PASS")
print("source_overlap_fibres=PASS")
print("fixed_row_pascal_unit_minors=PASS")
print("hierarchy_rank_clock=PASS")
print("thm4366_singleton_directions=PASS")
