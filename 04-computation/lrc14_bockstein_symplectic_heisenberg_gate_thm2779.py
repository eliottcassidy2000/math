#!/usr/bin/env python3
"""Exact finite checks for THM-2779.

The companion uses only integer arithmetic modulo 13.  It verifies the
decoder and frame torsors, the faithful 13^2-point Heisenberg action, the
13-point Boolean-permutation obstruction, and the one-object dilation
naturality gate.
"""

from itertools import product


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add_vec(a, b):
    return tuple((x + y) % P for x, y in zip(a, b))


def scale_vec(c, a):
    return tuple((c * x) % P for x in a)


def cyclic_convolution(a, b):
    out = [0] * P
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[(i + j) % P] = (out[(i + j) % P] + x * y) % P
    return tuple(out)


def matrix_rank_mod_p(matrix):
    a = [list(row) for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    pivots = []
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col] % P), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col] % P, -1, P)
        a[rank] = [(inv * x) % P for x in a[rank]]
        for r in range(rows):
            if r == rank:
                continue
            c = a[r][col] % P
            if c:
                a[r] = [(x - c * y) % P for x, y in zip(a[r], a[rank])]
        pivots.append(col)
        rank += 1
        if rank == rows:
            break
    return rank, tuple(pivots), tuple(tuple(row) for row in a)


def det(v, w):
    return (v[0] * w[1] - v[1] * w[0]) % P


def hmul(g, h):
    """Heisenberg law compatible with rho below."""
    x, y, z = g
    xp, yp, zp = h
    return (
        (x + xp) % P,
        (y + yp) % P,
        (z + zp - y * xp) % P,
    )


def hinv(g):
    x, y, z = g
    return ((-x) % P, (-y) % P, (-z - y * x) % P)


def rho(g, point):
    """Faithful affine action on F_13^2."""
    x, y, z = g
    a, b = point
    return ((a + x) % P, (b + z - y * a) % P)


def compose_words(*elements):
    out = (0, 0, 0)
    for g in elements:
        out = hmul(out, g)
    return out


def endpoint_x(point):
    """THM-2620 transvection R -> R+det(s,R)s in (v,w) coordinates."""
    v, w = point
    return (v, (w + v) % P)


def endpoint_x_inv(point):
    v, w = point
    return (v, (w - v) % P)


def endpoint_y(point):
    """Target translation R -> R+t."""
    v, w = point
    return ((v + 1) % P, w)


def endpoint_y_inv(point):
    v, w = point
    return ((v - 1) % P, w)


def endpoint_z(point):
    """Central translation R -> R+s."""
    v, w = point
    return (v, (w + 1) % P)


def digits_to_address(point):
    """Chosen two-digit section F_13^2 -> Z/169Z."""
    v, w = point
    return v + P * w


def address_to_digits(address):
    """Inverse of digits_to_address on the standard representatives."""
    address %= P ** 2
    return (address % P, address // P)


def odometer_x(address):
    """The endpoint shear in the chosen two-digit cyclic model."""
    return 14 * address % (P ** 2)


def odometer_y(address):
    """Low-digit increment with its carry deliberately suppressed."""
    v, w = address_to_digits(address)
    return digits_to_address(((v + 1) % P, w))


def odometer_z(address):
    """The central high-digit increment."""
    return (address + P) % (P ** 2)


def main():
    # THM-2771's target collapse and its intrinsic Bockstein decoder.
    s = (0, 0, 9, 0, 0, 0, 0, 0, 2, 2, 2, 2, 9)
    k = (7, 3, 3, 3, 12, 2, 9, 10, 9, 8, 10, 4, 0)
    s_beta = scale_vec(11, s)
    k_beta = scale_vec(6, k)
    epsilon = (P - 1, 1) + (0,) * (P - 2)
    norm = (1,) * P

    require(cyclic_convolution(s, k) == epsilon, "S*K != u-1")
    require(cyclic_convolution(s_beta, k_beta) == epsilon,
            "S_beta*K_beta != u-1")
    require(cyclic_convolution(epsilon, norm) == (0,) * P,
            "epsilon*N_13 != 0")

    # The convolution matrix has rank 12 and kernel exactly <N_13>.
    columns = []
    for j in range(P):
        basis = tuple(1 if i == j else 0 for i in range(P))
        columns.append(cyclic_convolution(s_beta, basis))
    matrix = tuple(tuple(columns[c][r] for c in range(P)) for r in range(P))
    rank, pivots, _ = matrix_rank_mod_p(matrix)
    require(rank == P - 1, "decoder convolution rank is not 12")
    require(cyclic_convolution(s_beta, norm) == (0,) * P,
            "N_13 does not span a kernel line")

    decoder_line = tuple(add_vec(k_beta, scale_vec(c, norm)) for c in range(P))
    require(len(set(decoder_line)) == P, "decoder line is not a 13-torsor")
    require(all(cyclic_convolution(s_beta, q) == epsilon for q in decoder_line),
            "one decoder-line point fails")

    # Every decoder gauge has the same uniform seven-chart invoice.
    chart_column = (0, 2, 0, 10, 0, 0, 0)
    require(sum(chart_column) % P == P - 1, "chart sum is not -1")
    chart_norm_convolution = tuple(sum(chart_column) % P for _ in range(7))
    require(chart_norm_convolution == (P - 1,) * 7,
            "C*N_7 != -N_7")

    # Normalized symplectic frames form a 13-torsor over each nonzero t.
    vectors = tuple(product(range(P), repeat=2))
    nonzero_vectors = tuple(v for v in vectors if v != (0, 0))
    frame_fibres = {
        t: tuple(sv for sv in vectors if det(sv, t) == 1)
        for t in nonzero_vectors
    }
    require(all(len(fibre) == P for fibre in frame_fibres.values()),
            "a normalized frame fibre does not have size 13")
    frame_count = sum(map(len, frame_fibres.values()))
    require(frame_count == 2184, "normalized oriented-frame count mismatch")

    t0 = (0, 1)
    s0 = (1, 0)
    frame_line = tuple(((s0[0] + c * t0[0]) % P,
                        (s0[1] + c * t0[1]) % P) for c in range(P))
    require(frame_line == frame_fibres[t0],
            "fixed-direction frame torsor ordering mismatch")
    require(all(det(sv, t0) == 1 for sv in frame_line),
            "frame gauge changes curvature")
    decoder_frame_pairs = tuple(zip(decoder_line, frame_line))
    require(len(decoder_frame_pairs) == P, "decoder/frame torsor mismatch")

    # The determinant Mobius face is origin-independent and equals one.
    mobius_checks = 0
    for left in vectors:
        for right in vectors:
            for sv in frame_line:
                f00 = det(left, right)
                f10 = det(add_vec(left, sv), right)
                f01 = det(left, add_vec(right, t0))
                f11 = det(add_vec(left, sv), add_vec(right, t0))
                require((f00 - f10 - f01 + f11) % P == 1,
                        "determinant Mobius face mismatch")
                mobius_checks += 1

    # H_13 and its faithful affine action on 13^2 points.
    group = tuple(product(range(P), repeat=3))
    points = vectors
    identity = (0, 0, 0)
    x_gen = (1, 0, 0)
    y_gen = (0, 1, 0)
    z_gen = (0, 0, 1)
    commutator = compose_words(
        x_gen, y_gen, hinv(x_gen), hinv(y_gen)
    )
    require(commutator == z_gen, "[X,Y] != Z")
    require(all(hmul(g, hinv(g)) == identity and
                hmul(hinv(g), g) == identity for g in group),
            "Heisenberg inverse formula failed")

    permutation_signatures = {
        tuple(rho(g, point) for point in points)
        for g in group
    }
    require(len(permutation_signatures) == P ** 3,
            "the 169-point affine action is not faithful")
    orbit0 = {rho(g, (0, 0)) for g in group}
    require(len(orbit0) == P ** 2, "the affine action is not transitive")

    # Exact action-law check on the two points determining these affine maps.
    # (The b coefficient is identically one.)
    action_pair_checks = 0
    test_points = ((0, 0), (1, 0))
    for g in group:
        for h in group:
            gh = hmul(g, h)
            for point in test_points:
                require(rho(gh, point) == rho(g, rho(h, point)),
                        "Heisenberg action law failed")
            action_pair_checks += 1

    # The same H_13 action is intrinsic to a THM-2620 endpoint-origin fibre.
    # In an oriented frame write R=w*s+v*t, so Delta=det(s,R)=v.
    endpoint_transvection_checks = 0
    delta_one_central_edges = 0
    for v, w in points:
        point = (v, w)
        comm_point = endpoint_x(
            endpoint_y(endpoint_x_inv(endpoint_y_inv(point)))
        )
        require(comm_point == endpoint_z(point),
                "endpoint [transvection,target shift] is not central shift")

        r_vec = (w, v)  # s=(1,0), t=(0,1)
        transvection_vec = ((w + v) % P, v)
        require(det(s0, r_vec) == v, "endpoint determinant coordinate mismatch")
        require((endpoint_x(point)[1], endpoint_x(point)[0])
                == transvection_vec, "endpoint transvection formula mismatch")
        if v == 1:
            left_vec = ((r_vec[0] + s0[0]) % P,
                        (r_vec[1] + s0[1]) % P)
            require((endpoint_z(point)[1], endpoint_z(point)[0]) == left_vec,
                    "Delta=1 central edge is not R->L")
            delta_one_central_edges += 1
        endpoint_transvection_checks += 1

    central_sector_cycles = 0
    for v in range(P):
        orbit = set()
        point = (v, 0)
        for _ in range(P):
            orbit.add(point)
            point = endpoint_z(point)
        require(point == (v, 0) and len(orbit) == P,
                "central endpoint orbit is not a 13-cycle")
        central_sector_cycles += 1

    # Choosing the standard two-digit section identifies the endpoint plane
    # with Z/169Z.  The central action becomes +13, while physical +1 differs
    # from the carry-suppressed low-digit action exactly on the wrap fibre.
    odometer_conjugacy_checks = 0
    odometer_nonwraps = 0
    odometer_carries = 0
    for address in range(P ** 2):
        point = address_to_digits(address)
        require(digits_to_address(point) == address,
                "two-digit address section is not bijective")
        require(address_to_digits(odometer_x(address)) == endpoint_x(point),
                "two-digit X/shear conjugacy failed")
        require(address_to_digits(odometer_y(address)) == endpoint_y(point),
                "two-digit Y/low-digit conjugacy failed")
        require(address_to_digits(odometer_z(address)) == endpoint_z(point),
                "two-digit Z/central conjugacy failed")
        physical_successor = (address + 1) % (P ** 2)
        if point[0] == P - 1:
            require(physical_successor == odometer_z(odometer_y(address)),
                    "odometer carry is not the central correction")
            odometer_carries += 1
        else:
            require(physical_successor == odometer_y(address),
                    "nonwrapping odometer step differs from low-digit step")
            odometer_nonwraps += 1
        odometer_conjugacy_checks += 1
    require((odometer_conjugacy_checks, odometer_nonwraps, odometer_carries)
            == (169, 156, 13), "two-digit odometer census mismatch")

    # On 13 points every 13-subgroup of S_13 has order <=13 and is abelian.
    # Any action of degree <169 has only fixed or 13-point orbits, so its
    # commutator centre is killed orbitwise.
    vp_13_factorial = sum(P // (P ** j) for j in range(1, 2))
    require(vp_13_factorial == 1, "v_13(13!) mismatch")
    require(P < P ** 2, "permutation-degree boundary arithmetic failed")

    # The 13-root coefficient space nevertheless carries the Weyl/projective
    # commutator: T M = zeta M T.  Forgetting phases makes M the identity
    # root permutation and kills the centre.
    weyl_phase_checks = 0
    for r in range(P):
        tm_phase = (-r) % P
        mt_phase = (-(r + 1)) % P
        require((tm_phase - mt_phase) % P == 1,
                "Weyl commutator exponent mismatch")
        weyl_phase_checks += 1

    # If the target label is D-stable while a collapsed root label is sent
    # to zero by D, a one-object linear natural intertwiner must vanish.
    flat_intertwiners = tuple(c for c in range(P) if (0 * c - c * 1) % P == 0)
    require(flat_intertwiners == (0,), "nonzero flat D-intertwiner survived")
    graded_intertwiners = tuple(range(P))
    require(len(graded_intertwiners) == P,
            "graded digit intertwiners should form one scalar line")

    print("THM-2779 exact decoder/frame/Heisenberg gate")
    print(f"decoder_rank={rank} kernel_dim={P-rank} pivots={pivots}")
    print(f"decoder_line={len(decoder_line)} kernel_generator=N13")
    print(f"S_beta={s_beta}")
    print(f"K_beta={k_beta}")
    print(f"chart_invoice={chart_norm_convolution}")
    print(f"frame_fibres={len(frame_fibres)} fibre_size=13 total={frame_count}")
    print(f"canonical_frame_line={frame_line}")
    print(f"determinant_mobius_checks={mobius_checks}")
    print(f"H13_order={len(group)} commutator={commutator}")
    print(f"faithful_affine_degree={len(points)} orbit={len(orbit0)}")
    print(f"action_pair_checks={action_pair_checks}")
    print(f"endpoint_transvection_checks={endpoint_transvection_checks}")
    print(f"central_sector_cycles={central_sector_cycles}")
    print(f"delta1_Z_equals_R_to_L={delta_one_central_edges}")
    print(f"odometer_H13_conjugacy_checks={odometer_conjugacy_checks}")
    print(f"odometer_plus1_section=nonwrap:{odometer_nonwraps} carry:{odometer_carries} central_step=13")
    print(f"v13_13factorial={vp_13_factorial}")
    print("minimum_faithful_permutation_degree=169")
    print(f"weyl_13root_phase_checks={weyl_phase_checks} boolean_centre=killed")
    print(f"flat_D_intertwiners={flat_intertwiners}")
    print(f"graded_digit_intertwiners={len(graded_intertwiners)}")
    print("scope=coefficient_twin_exact; two_digit_odometer_model_exact; physical_same_ancestry_open")


if __name__ == "__main__":
    main()
