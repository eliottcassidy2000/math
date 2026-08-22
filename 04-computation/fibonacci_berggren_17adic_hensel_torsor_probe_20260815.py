#!/usr/bin/env python3
"""Unnumbered exact probe for the Fibonacci/Berggren 17-adic Hensel torsor.

Dependency-free, integer-only, deterministic, and safe under ``python -O``.
The companion verifies finite universes and emits the first-layer certificate
used by the coefficient-independent all-exponent proof in the reflection.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path


P = 17
TEST_EXPONENTS = (1, 2, 3, 4)
G = ((0, 1), (1, 1))
IDENTITY = ((1, 0), (0, 1))
BRANCHES = {
    "A": ((0, 1), (-1, 2)),
    "B": ((0, 1), (1, 2)),
    "C": ((1, 0), (2, 1)),
}
EXPECTED_SEMANTIC_SHA256 = "bc6795aa3a4647b44cfa0b15f5112503cace4c9228dbdf18afb9e5aa8217fff0"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def lf_sha256(path):
    return sha256(lf_bytes(path)).hexdigest()


def mmul(left, right):
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def mvec(matrix, vector):
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def madd(left, right):
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(2)) for i in range(2)
    )


def mscale(scalar, matrix):
    return tuple(tuple(scalar * matrix[i][j] for j in range(2)) for i in range(2))


def mmod(matrix, modulus):
    return tuple(tuple(entry % modulus for entry in row) for row in matrix)


def mpow(matrix, exponent):
    require(exponent >= 0, exponent)
    out = IDENTITY
    base = matrix
    while exponent:
        if exponent & 1:
            out = mmul(out, base)
        base = mmul(base, base)
        exponent >>= 1
    return out


def mpow_mod(matrix, exponent, modulus):
    require(exponent >= 0 and modulus > 0, (exponent, modulus))
    out = IDENTITY
    base = mmod(matrix, modulus)
    while exponent:
        if exponent & 1:
            out = mmod(mmul(out, base), modulus)
        base = mmod(mmul(base, base), modulus)
        exponent >>= 1
    return out


def determinant(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def discriminant(matrix, modulus):
    trace = (matrix[0][0] + matrix[1][1]) % modulus
    det = determinant(matrix) % modulus
    return (trace * trace - 4 * det) % modulus


def legendre(value, prime=P):
    value %= prime
    require(value != 0, ("Legendre zero", prime))
    out = pow(value, (prime - 1) // 2, prime)
    require(out in (1, prime - 1), (value, prime, out))
    return 1 if out == 1 else -1


def normalize_projective(vector, modulus):
    x, y = vector[0] % modulus, vector[1] % modulus
    if gcd(y, P) == 1:
        return (x * pow(y, -1, modulus) % modulus, 1)
    require(gcd(x, P) == 1, ("nonprimitive vector", vector, modulus))
    out = (1, y * pow(x, -1, modulus) % modulus)
    require(out[1] % P == 0, (vector, modulus, out))
    return out


def projective_points(exponent):
    modulus = P**exponent
    return tuple(
        [(x, 1) for x in range(modulus)]
        + [(1, P * y) for y in range(P ** (exponent - 1))]
    )


def projective_action(matrix, point, modulus):
    return normalize_projective(mvec(matrix, point), modulus)


def reduce_to_base(point):
    return normalize_projective((point[0] % P, point[1] % P), P)


def cycle_decomposition(matrix, points, modulus):
    unseen = set(points)
    cycles = []
    for start in points:
        if start not in unseen:
            continue
        cycle = []
        current = start
        while current in unseen:
            unseen.remove(current)
            cycle.append(current)
            current = projective_action(matrix, current, modulus)
        require(current == start, ("non-cycle return", start, current))
        cycles.append(tuple(cycle))
    require(not unseen, len(unseen))
    return tuple(cycles)


def norm_form(point):
    m, n = point
    return n * n - m * n - m * m


def fib_pair(index, modulus=None):
    require(index >= 0, index)
    a, b = 0, 1
    for _ in range(index):
        a, b = b, a + b
        if modulus is not None:
            a %= modulus
            b %= modulus
    return a, b


def apply_word(word, vector):
    out = vector
    for letter in word:
        out = mvec(BRANCHES[letter], out)
    return out


def normalized_fibonacci_parameter(index):
    m, n = fib_pair(index)
    if m % 2 == 1 and n % 2 == 1:
        return ((n - m) // 2, (n + m) // 2)
    return m, n


def ray_word(index):
    require(index >= 2, index)
    if index % 3 == 2:
        r = (index - 2) // 3
        return "BA" * r, "R2"
    if index % 3 == 0:
        r = (index - 3) // 3
        return "A" + "BA" * r, "R0"
    r = (index - 4) // 3
    return "C" + "BC" * r, "R1"


def euclid_triple(parameter):
    m, n = parameter
    return n * n - m * m, 2 * m * n, n * n + m * m


def local_t4(parameter):
    m, n = parameter
    a, b, _ = euclid_triple(parameter)
    delta = b - a
    require(delta != 0, ("T4 tie", parameter))
    epsilon = int(delta < 0)
    vertices = {"P": parameter}
    vertices.update({letter: mvec(matrix, parameter) for letter, matrix in BRANCHES.items()})
    weights = {label: u[0] * u[0] + u[1] * u[1] for label, u in vertices.items()}
    require(len(set(weights.values())) == 4, (parameter, weights))
    order = tuple(sorted(weights, key=weights.get))
    expected = ("P", "A", "C", "B") if epsilon == 0 else ("P", "C", "A", "B")
    require(order == expected, (parameter, delta, weights, order, expected))
    return epsilon, order


EDGE_LABELS = ("01", "02", "03", "12", "13", "23")


def fibonacci_t6(index):
    require(index >= 2, index)
    f0, _ = fib_pair(index - 1)
    f1, _ = fib_pair(index)
    f2, _ = fib_pair(index + 1)
    f3, _ = fib_pair(index + 2)
    window = (f0, f1, f2, f3)
    weights = {
        "01": window[0] * window[1],
        "02": window[0] * window[2],
        "03": window[0] * window[3],
        "12": window[1] * window[2],
        "13": window[1] * window[3],
        "23": window[2] * window[3],
    }
    order = tuple(sorted(EDGE_LABELS, key=lambda edge: (weights[edge], edge)))
    return window, weights, order


def oriented_farey_state(index):
    left = fib_pair(index)
    right = fib_pair(index + 1)
    det = left[0] * right[1] - left[1] * right[0]
    require(abs(det) == 1, (index, left, right, det))
    x, y = (left, right) if det == 1 else (right, left)
    xbar = (x[0] % 2, x[1] % 2)
    ybar = (y[0] % 2, y[1] % 2)
    labels = {(1, 0): "a", (0, 1): "b", (1, 1): "c"}
    middle = ((xbar[0] + ybar[0]) % 2, (xbar[1] + ybar[1]) % 2)
    require(xbar in labels and ybar in labels and middle in labels, (index, xbar, ybar))
    return labels[xbar] + labels[middle] + labels[ybar]


def base_orbit_certificate():
    points = projective_points(1)
    cycles = cycle_decomposition(G, points, P)
    require(tuple(len(cycle) for cycle in cycles) == (9, 9), cycles)
    square_set = {x * x % P for x in range(1, P)}
    require(5 not in square_set and 3 not in square_set and -1 % P in square_set, square_set)

    for point in points:
        image = mvec(G, point)
        require((norm_form(image) + norm_form(point)) % P == 0, (point, image))
        require(norm_form(point) % P != 0, ("isotropic line", point))

    signed_cycles = []
    for cycle in cycles:
        signs = {legendre(norm_form(point)) for point in cycle}
        require(len(signs) == 1, (cycle, signs))
        signed_cycles.append((next(iter(signs)), cycle))
    require(tuple(sign for sign, _ in signed_cycles) == (1, -1), signed_cycles)
    require(tuple(len(cycle) for _, cycle in signed_cycles) == (9, 9), signed_cycles)

    fibonacci_lines = {
        normalize_projective(fib_pair(index, P), P) for index in range(2, 11)
    }
    require(fibonacci_lines == set(signed_cycles[0][1]), fibonacci_lines)
    for index in range(2, 101):
        u = fib_pair(index)
        require(norm_form(u) == (-1) ** index, (index, u, norm_form(u)))
        require(legendre(norm_form(u)) == 1, (index, u))
    return tuple((sign, cycle) for sign, cycle in signed_cycles)


def first_layer_certificate(g9, k, c, b):
    require(g9 == ((21, 34), (34, 55)), g9)
    require(k == ((1597, 2584), (2584, 4181)), k)
    require(g9 == madd(mscale(4, IDENTITY), mscale(P, c)), (g9, c))
    require(k == madd(mscale(-1, IDENTITY), mscale(P, b)), (k, b))
    require(mmod(k, P) == ((16, 0), (0, 16)), mmod(k, P))
    require(mmod(k, 2) == IDENTITY, mmod(k, 2))
    require(discriminant(c, P) == 3 and legendre(3) == -1, c)
    require(discriminant(b, P) == 5 and legendre(5) == -1, b)

    u = mscale(-1, k)
    require(u == madd(IDENTITY, mscale(-P, b)), u)
    congruence_rows = []
    for s in range(4):
        modulus = P ** (s + 2)
        for q in (1, 2, 16, 18):
            require(q % P != 0, q)
            exponent = (P**s) * q
            observed = mpow_mod(u, exponent, modulus)
            predicted = mmod(madd(IDENTITY, mscale(-q * P ** (s + 1), b)), modulus)
            require(observed == predicted, (s, q, observed, predicted))
        congruence_rows.append((s, modulus, "4/4"))

    split_b = ((1, 0), (0, 2))
    split_u = madd(IDENTITY, mscale(P, split_b))
    modulus = P**2
    split_fibre = {point for point in projective_points(2) if reduce_to_base(point) == (1, 0)}
    split_seed = (1, 0)
    require(len(split_fibre) == P, len(split_fibre))
    require(projective_action(split_u, split_seed, modulus) == split_seed, split_u)
    require(discriminant(split_b, P) == 1, split_b)
    return tuple(congruence_rows), (len(split_fibre), 1, discriminant(split_b, P))


def exhaustive_hensel_audit(k):
    rows = []
    for exponent in TEST_EXPONENTS:
        modulus = P**exponent
        fibre_size = P ** (exponent - 1)
        points = projective_points(exponent)
        require(len(points) == (P + 1) * fibre_size, (exponent, len(points)))
        require(len(set(points)) == len(points), exponent)

        g_cycles = cycle_decomposition(G, points, modulus)
        k_cycles = cycle_decomposition(k, points, modulus)
        require(len(g_cycles) == 2, (exponent, tuple(map(len, g_cycles))))
        require({len(cycle) for cycle in g_cycles} == {9 * fibre_size}, exponent)
        require(len(k_cycles) == P + 1, (exponent, len(k_cycles)))
        require({len(cycle) for cycle in k_cycles} == {fibre_size}, exponent)

        fibres = {}
        for point in points:
            fibres.setdefault(reduce_to_base(point), set()).add(point)
        require(len(fibres) == P + 1, (exponent, len(fibres)))
        require({len(fibre) for fibre in fibres.values()} == {fibre_size}, exponent)

        seen_cycle_fibres = set()
        for cycle in k_cycles:
            reductions = {reduce_to_base(point) for point in cycle}
            require(len(reductions) == 1, (exponent, reductions))
            base = next(iter(reductions))
            require(set(cycle) == fibres[base], (exponent, base))
            seen_cycle_fibres.add(base)
        require(seen_cycle_fibres == set(fibres), exponent)

        rows.append(
            (
                exponent,
                len(points),
                len(g_cycles),
                next(iter({len(cycle) for cycle in g_cycles})),
                len(k_cycles),
                next(iter({len(cycle) for cycle in k_cycles})),
                len(fibres),
                fibre_size,
            )
        )
    return tuple(rows)


def fibonacci_berggren_and_tournament_audit(k):
    for index in range(2, 201):
        word, ray = ray_word(index)
        parameter = normalized_fibonacci_parameter(index)
        require(apply_word(word, (1, 2)) == parameter, (index, word, ray, parameter))
        shifted_word, shifted_ray = ray_word(index + 18)
        suffix = "BA" * 6 if ray in ("R0", "R2") else "BC" * 6
        require(shifted_ray == ray and shifted_word == word + suffix, (index, word, shifted_word))

        epsilon, t4_order = local_t4(parameter)
        shifted_epsilon, shifted_t4_order = local_t4(normalized_fibonacci_parameter(index + 18))
        require(epsilon == int(index % 3 == 1), (index, parameter, epsilon))
        require((epsilon, t4_order) == (shifted_epsilon, shifted_t4_order), index)

        if index >= 3:
            _, weights, order = fibonacci_t6(index)
            _, shifted_weights, shifted_order = fibonacci_t6(index + 18)
            require(len(set(weights.values())) == 6, (index, weights))
            expected = (
                ("01", "02", "03", "12", "13", "23")
                if index % 2 == 1
                else ("01", "02", "12", "03", "13", "23")
            )
            require(order == expected and shifted_order == expected, (index, order, shifted_order))
            require(weights["03"] - weights["12"] == (-1) ** index, (index, weights))
            require(shifted_weights["03"] - shifted_weights["12"] == (-1) ** index, index)

        require(oriented_farey_state(index + 18) == oriented_farey_state(index), index)
        raw = fib_pair(index)
        shifted_raw = fib_pair(index + 18)
        require(mvec(k, raw) == shifted_raw, (index, raw, shifted_raw))

    farey_cycle = tuple(oriented_farey_state(index) for index in range(2, 8))
    require(farey_cycle == ("bca", "bac", "abc", "acb", "cab", "cba"), farey_cycle)

    root_window, root_weights, _ = fibonacci_t6(2)
    tied_values = tuple(sorted(value for value in set(root_weights.values()) if list(root_weights.values()).count(value) > 1))
    require(root_window == (1, 1, 2, 3), root_window)
    require(tied_values == (2, 3), (root_weights, tied_values))

    modulus = P**2
    state_3 = normalize_projective(fib_pair(3, modulus), modulus)
    state_21 = normalize_projective(fib_pair(21, modulus), modulus)
    require(state_3 == (97, 1) and state_21 == (63, 1), (state_3, state_21))
    require(reduce_to_base(state_3) == reduce_to_base(state_21), (state_3, state_21))
    require(projective_action(k, state_3, modulus) == state_21, (state_3, state_21))
    require(local_t4(normalized_fibonacci_parameter(3)) == local_t4(normalized_fibonacci_parameter(21)), "T4 hostile")
    require(fibonacci_t6(3)[2] == fibonacci_t6(21)[2], "T6 hostile")

    return (
        (2, 200, "all_three_exact_rays", "shift_18_appends_six_branch_pairs"),
        farey_cycle,
        (root_window, tied_values),
        ((3, state_3), (21, state_21), reduce_to_base(state_3)),
    )


def periodic_support(base_index, modulus_index, subset):
    period = 18 * modulus_index
    return {(base_index + 18 * residue) % period for residue in subset}


def xor_harmonic_audit():
    rows = []
    base_index = 3
    for exponent in TEST_EXPONENTS:
        torsor_order = P ** (exponent - 1)
        period = 18 * torsor_order
        left = {j for j in range(torsor_order) if j % 2 == 0}
        right = {j for j in range(torsor_order) if j % 3 == 0}
        carrier = periodic_support(base_index, torsor_order, set(range(torsor_order)))
        s_left = periodic_support(base_index, torsor_order, left)
        s_right = periodic_support(base_index, torsor_order, right)
        s_xor = periodic_support(base_index, torsor_order, left ^ right)
        require(len(carrier) == torsor_order, (exponent, len(carrier)))
        require(s_left ^ s_right == s_xor, exponent)
        require(s_left & s_right == periodic_support(base_index, torsor_order, left & right), exponent)
        require(s_left | s_right == periodic_support(base_index, torsor_order, left | right), exponent)
        require(carrier - s_left == periodic_support(base_index, torsor_order, set(range(torsor_order)) - left), exponent)
        pole = Fraction(len(left), period)
        xor_pole = Fraction(len(left ^ right), period)
        require(xor_pole == Fraction(len(left) + len(right) - 2 * len(left & right), period), exponent)
        rows.append((exponent, torsor_order, period, len(left), len(right), len(left ^ right), pole, xor_pole))
    return tuple(rows)


def affine_lift_transducer_audit(k):
    rows = []
    previous_states = None
    seed_vector = fib_pair(3)
    for exponent in TEST_EXPONENTS:
        modulus = P**exponent
        order = P ** (exponent - 1)
        seed = normalize_projective(seed_vector, modulus)
        states = tuple(
            projective_action(mpow_mod(k, lift, modulus), seed, modulus)
            for lift in range(order)
        )
        require(len(set(states)) == order, (exponent, len(set(states)), order))
        for lift in range(order):
            require(projective_action(k, states[lift], modulus) == states[(lift + 1) % order], (exponent, lift))
        if previous_states is not None:
            for lift, state in enumerate(states):
                reduced = normalize_projective((state[0] % (P ** (exponent - 1)), state[1] % (P ** (exponent - 1))), P ** (exponent - 1))
                require(reduced == previous_states[lift % len(previous_states)], (exponent, lift, reduced))
        rows.append((exponent, order, len(set(states)), states[0], states[1 % order]))
        previous_states = states
    return tuple(rows)


def half_twist_fibre_section(modulus, residue, prime, base_sheet):
    base_modulus = modulus // prime
    hits = []
    for t in range(prime):
        sheet = base_sheet + base_modulus * t
        phase_numerator = (residue * (2 * sheet + 1)) % (2 * modulus)
        distance_numerator = min(phase_numerator, 2 * modulus - phase_numerator)
        if 7 * distance_numerator < modulus:
            hits.append(t)
    return tuple(hits)


def affine_lift_hostile():
    modulus = P**2
    residues = (1, 35)
    require(residues[0] % (2 * P) == residues[1] % (2 * P) == 1, residues)
    fixed = (modulus // P - 1) // 2
    fixed_sections = tuple(half_twist_fibre_section(modulus, residue, P, fixed) for residue in residues)
    moving_sections = tuple(half_twist_fibre_section(modulus, residue, P, 1) for residue in residues)
    require(fixed_sections == ((0, 16), (0, 16)), fixed_sections)
    require(moving_sections == ((0, 1, 16), (13, 14, 15)), moving_sections)
    require(set(moving_sections[0]).isdisjoint(moving_sections[1]), moving_sections)
    return (modulus, residues, (0, 1), fixed, fixed_sections, 1, moving_sections)


def main():
    g9 = mpow(G, 9)
    k = mpow(G, 18)
    c = tuple(tuple((g9[i][j] - 4 * IDENTITY[i][j]) // P for j in range(2)) for i in range(2))
    b = tuple(tuple((k[i][j] + IDENTITY[i][j]) // P for j in range(2)) for i in range(2))

    base_cycles = base_orbit_certificate()
    first_layer = first_layer_certificate(g9, k, c, b)
    hensel_rows = exhaustive_hensel_audit(k)
    tournament_rows = fibonacci_berggren_and_tournament_audit(k)
    xor_rows = xor_harmonic_audit()
    transducer_rows = affine_lift_transducer_audit(k)
    lift_hostile = affine_lift_hostile()

    semantic_surface = (
        P,
        TEST_EXPONENTS,
        G,
        g9,
        k,
        c,
        b,
        base_cycles,
        first_layer,
        hensel_rows,
        tournament_rows,
        xor_rows,
        transducer_rows,
        lift_hostile,
        (
            "T4:labelled parent/A/B/C hypotenuse order",
            "T6:labelled Fibonacci-window edge-product order",
            "XOR:symmetric_difference_on_calibrated_K_torsor",
            "harmonic_pole=periodic_support_density",
            "no_LRC14_decrement",
        ),
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("UNNUMBERED fibonacci/Berggren 17-adic Hensel torsor exact companion")
    print(f"script_sha256_lf={lf_sha256(__file__)}")
    print(f"matrices=G:{G};G9:{g9};K=G18:{k}")
    print(f"first_layers=G9=4I+17C:C={c}:disc(C)={discriminant(c, P)};K=-I+17B:B={b}:disc(B)={discriminant(b, P)}")
    print(f"P1_F17_cycles_by_Legendre_D={base_cycles}")
    print(f"hensel_exact_rows=(a,|P1|,#G,Glen,#K,Klen,#fibres,fibre_size)={hensel_rows}")
    print(f"all_a_proof_certificate=U=-K=I-17B;disc(B)=5_nonsquare;binomial_layers_and_split_hostile={first_layer}")
    print(f"berggren_farey_tournaments={tournament_rows}")
    print(f"affine_lift_transducer=(a,order,distinct,state0,state1)={transducer_rows}")
    print(f"xor_harmonic=(a,N,period,|J2|,|J3|,|xor|,pole_J2,pole_xor)={xor_rows}")
    print(f"Q289_same_mod34_affine_lift_hostile={lift_hostile}")
    print("tournament_gate=T4/T6_are_labelled_transitive_comparison_tournaments_on_the_stated_tie-free_domains;no_missing_or_bidirected_edges")
    print("preserved=Legendre_D_orbit,three_ray_class,labelled_T4_bit,labelled_T6_order,Cassini_sign,oriented_mod2_Farey_frame,XOR_after_calibration,harmonic_pole")
    print("lost=exact_fraction,weights,triple,word_depth,17adic_lift,subset_phase,affine_owner,current;sidecar=base_lift_plus_K_exponent_and_exact_weights/current")
    print(f"semantic_sha256={semantic_digest}")
    print("status=FINITE-EXACT_a=1..4_PLUS_COEFFICIENT-INDEPENDENT_ALL-a_PROOF_PACKET;unnumbered;indexed_as_result_only;no_theorem_or_LRC14_decrement")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
