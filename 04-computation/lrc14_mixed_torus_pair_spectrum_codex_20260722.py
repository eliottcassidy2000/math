#!/usr/bin/env python3
"""Exact pair spectrum and a hostile rank-eight equality family (THM-2120)."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm, pi, sin


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def circle_distance(value):
    residue = F(value) % 1
    return min(residue, 1 - residue)


def primitive_relation(g, f, fp):
    raw = (det(f, fp), det(fp, g), det(g, f))
    divisor = gcd(gcd(abs(raw[0]), abs(raw[1])), abs(raw[2]))
    require(divisor > 0, "zero relation")
    relation = tuple(value // divisor for value in raw)
    if relation[0] < 0:
        relation = tuple(-value for value in relation)
    return raw, divisor, relation


def positive_square(value):
    return max(value, F(0)) ** 2


def square_linear_cdf(t, P, Q):
    """Area in [-1,1]^2 on which P*x+Q*y <= t, for P,Q>0."""
    t = F(t)
    s = t + P + Q
    return (
        positive_square(s)
        - positive_square(s - 2 * P)
        - positive_square(s - 2 * Q)
        + positive_square(s - 2 * P - 2 * Q)
    ) / (2 * P * Q)


def independent_pair_weight(A, B, C):
    """Exact weight for primitive A*g+B*f+C*f'=0."""
    A, P, Q = abs(A), abs(B), abs(C)
    require(A > 0 and P > 0 and Q > 0, "full-rank transverse relation required")
    cutoff = (P + Q + 2 * A) // 14 + 2
    guard_area_sum = F(0)
    for n in range(-cutoff, cutoff + 1):
        guard_area_sum += square_linear_cdf(14 * n + 2 * A, P, Q)
        guard_area_sum -= square_linear_cdf(14 * n - 2 * A, P, Q)
    return F(1, 49) - guard_area_sum / (196 * A)


def fourier_weight(A, B, C, terms=20000):
    A, B, C = abs(A), abs(B), abs(C)
    correction = 0.0
    for m in range(1, terms + 1):
        guard = sin(2 * pi * m * A / 7) / (pi * m * A)
        first = sin(pi * m * B / 7) / (pi * m * B)
        second = sin(pi * m * C / 7) / (pi * m * C)
        correction += 2 * guard * first * second
    return 5 / 343 - correction


def line_weight(A):
    """Closed form for w(A,1,-1)."""
    A = abs(A)
    residue = A % 7
    if residue == 0:
        return F(5, 343)
    return F(5, 343) + F(2 * residue - 7, 343 * A)


def max_spanning_tree(weights, n):
    used = {0}
    total = F(0)
    edges = []
    while len(used) < n:
        best = None
        for i in used:
            for j in range(n):
                if j in used:
                    continue
                candidate = (weights[tuple(sorted((i, j)))], i, j)
                if best is None or candidate > best:
                    best = candidate
        require(best is not None, "disconnected")
        weight, i, j = best
        used.add(j)
        total += weight
        edges.append((i, j, weight))
    return total, edges


def main():
    require(independent_pair_weight(1, 1, 1) == 0, "pair-rigidity zero")
    require(independent_pair_weight(1, 1, 2) == F(1, 392), "first positive")
    for A in range(1, 401):
        require(
            line_weight(A) == independent_pair_weight(A, 1, -1),
            (A, line_weight(A), independent_pair_weight(A, 1, -1)),
        )

    saturated = ((0, 1), (1, 0), (1, 2))
    nonsaturated = ((0, 1), (2, 0), (2, 2))
    raw1, index1, relation1 = primitive_relation(*saturated)
    raw2, index2, relation2 = primitive_relation(*nonsaturated)
    require(relation1 == relation2 == (2, 1, -1), (raw1, raw2, relation1, relation2))
    require(index1 == 1 and index2 == 2, "finite-kernel indices")
    finite_kernel_weight = independent_pair_weight(*relation1)

    neutrality_checks = 0
    for A in range(1, 30):
        for B in range(1, 16):
            for C in range(1, 16):
                if gcd(gcd(A, B), C) > 1:
                    continue
                if A % 7 == 0 or B % 7 == 0 or C % 7 == 0:
                    require(independent_pair_weight(A, B, C) == F(5, 343), (A, B, C))
                    neutrality_checks += 1

    fourier_cases = ((1, 1, 2), (2, 1, 1), (3, 2, 5), (5, 4, 3), (8, 3, 2), (11, 5, 6))
    max_fourier_error = 0.0
    for relation in fourier_cases:
        exact = float(independent_pair_weight(*relation))
        approximate = fourier_weight(*relation)
        max_fourier_error = max(max_fourier_error, abs(exact - approximate))
    require(max_fourier_error < 2e-10, max_fourier_error)

    guard = (0, 1)
    characters = [(1, 7 * i) for i in range(8)]
    direction = (2, 1)
    h = guard[0] * direction[0] + guard[1] * direction[1]
    speeds = [c[0] * direction[0] + c[1] * direction[1] for c in characters]
    require(h == 1 and speeds == [2 + 7 * i for i in range(8)], (h, speeds))
    require(len(set(speeds)) == 8 and all(q > 0 for q in speeds), "positive specialization")
    weights = {}
    for i, j in combinations(range(8), 2):
        raw, index, relation = primitive_relation(guard, characters[i], characters[j])
        require(
            relation[0] == 7 * (j - i) and abs(relation[1]) == abs(relation[2]) == 1,
            (i, j, raw, index, relation),
        )
        weight = independent_pair_weight(*relation)
        require(weight == F(5, 343), (i, j, relation, weight))
        weights[(i, j)] = weight
    tau, tree = max_spanning_tree(weights, 8)
    require(tau == F(5, 49), (tau, tree))

    # THM-2114's two new carriers both detect that equality is non-live.
    # The mod-seven vector (1,1) misses every character kernel.  Independently,
    # this guard-safe point has active set {0,1}, which cannot be connected in
    # every spanning tree when every tree is maximum.
    needle = (1, 1)
    needle_residues = (
        (guard[0] * needle[0] + guard[1] * needle[1]) % 7,
        *tuple((c[0] * needle[0] + c[1] * needle[1]) % 7 for c in characters),
    )
    require(all(needle_residues), needle_residues)
    equality_point = (F(3, 50), F(47, 175))
    active = tuple(
        i
        for i, c in enumerate(characters)
        if circle_distance(c[0] * equality_point[0] + c[1] * equality_point[1])
        < F(1, 14)
    )
    require(circle_distance(equality_point[1]) > F(1, 7), equality_point)
    require(active == (0, 1), active)

    asymptotic = []
    for N in (8, 15, 29, 50, 99, 197, 400):
        family_weights = {}
        for i, j in combinations(range(8), 2):
            family_weights[(i, j)] = line_weight(N * (j - i))
        family_tau, _ = max_spanning_tree(family_weights, 8)
        asymptotic.append((N, family_tau, family_tau - F(5, 49)))
    require(abs(float(asymptotic[-1][2])) < F(1, 1000), asymptotic[-1])

    # The vanishing margin survives the actual terminal arithmetic sidecars.
    # With L=lcm(2,...,14), N=1+kL, and specialization d=(7,1), the speeds
    # Q={7,7+N,...,7+7N} are divisor-complete through 14.  Every one-deletion
    # retains an adjacent pair whose gcd is gcd(7,N)=1.
    clock_lcm = 1
    for modulus in range(2, 15):
        clock_lcm = lcm(clock_lcm, modulus)
    require(clock_lcm == 360360, "clock lcm")
    live_rows = []
    for k in (0, 1, 2, 10, 100):
        N = 1 + k * clock_lcm
        speeds_live = [7 + N * i for i in range(8)]
        require(N % 7 == 1, "live residue")
        require(
            all(any(q % modulus == 0 for q in speeds_live) for modulus in range(2, 15)),
            (k, N, speeds_live, "divisor completeness"),
        )
        require(
            all(
                gcd(*[q for j, q in enumerate(speeds_live) if j != i]) == 1
                for i in range(8)
            ),
            (k, N, speeds_live, "hereditary primitivity"),
        )
        live_weights = {
            (i, j): line_weight(N * (j - i)) for i, j in combinations(range(8), 2)
        }
        live_tau, _ = max_spanning_tree(live_weights, 8)
        live_margin = live_tau - F(5, 49)
        require(live_margin == F(17, 1470 * N), (k, N, live_tau, live_margin))
        live_rows.append((k, N, speeds_live[-1], live_margin))

    print(f"zero_layer={independent_pair_weight(1, 1, 1)}")
    print(f"first_positive_layer={independent_pair_weight(1, 1, 2)}")
    print(
        f"finite_kernel_relation={relation1} saturated_index={index1} "
        f"nonsaturated_index={index2}"
    )
    print(f"finite_kernel_weight={finite_kernel_weight}")
    print(f"neutrality_checks={neutrality_checks}")
    print(f"max_fourier_error={max_fourier_error:.3e}")
    print(f"equality_guard={h} equality_speeds={tuple(speeds)}")
    print(f"equality_pair_weight={F(5, 343)} equality_tree_weight={tau}")
    print(f"equality_mod7_needle_residues={needle_residues}")
    print(f"equality_maxbasis_control_point={equality_point} active={active}")
    print("asymptotic_tree_margins:")
    for N, family_tau, margin in asymptotic:
        print(f"  N={N} tau={family_tau} margin={margin} ({float(margin):+.6e})")
    print(f"live_clock_lcm={clock_lcm}")
    print("divisor_complete_hereditary_margins:")
    for k, N, maximum, margin in live_rows:
        print(f"  k={k} N={N} maxQ={maximum} margin={margin}")
    print("PASS")


if __name__ == "__main__":
    main()
