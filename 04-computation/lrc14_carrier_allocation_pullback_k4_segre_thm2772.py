#!/usr/bin/env python3
"""Exact companion for reserved proof-complete THM-2772.

The companion verifies the finite parts of the carrier-allocation pullback
theorem over F_13:

* the Boolean bare/source/target/both allocation-lift fibre census;
* the 2,185-sector affine-square dictionary and marked K4 transform;
* every endpoint-origin determinant fibre;
* the endpoint-parallelogram determinant curvature and ordered-frame census;
* all 13^4 Segre-Hadamard endpoint quadruples and the sharp hostiles; and
* the 13^4/13^4/13^6 cospan cardinalities.

The universal property of the pullback and the Cech boundary invoice are
proved algebraically in the theorem.  This script uses explicit exceptions,
not Python assertions, so optimized mode retains every truth-bearing check.
"""

from collections import Counter


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def det(x, y):
    return (x[0] * y[1] - x[1] * y[0]) % P


def marked_difference(values, modulus=P):
    """THM-2757 difference for mark 00 and arms 10,01,11."""
    v00, v10, v01, v11 = values
    return (
        (v00 + v10 - v01 - v11) % modulus,
        (v00 - v10 + v01 - v11) % modulus,
        (v00 - v10 - v01 + v11) % modulus,
    )


def mixed_rectangle(values):
    v00, v10, v01, v11 = values
    return (v00 - v10 - v01 + v11) % P


def compatible_bits(k_source, k_target):
    """Allocation states whose absent endpoint has harmonic zero."""
    return tuple(
        (source_present, target_present)
        for source_present in (0, 1)
        for target_present in (0, 1)
        if (source_present or k_source == 0)
        and (target_present or k_target == 0)
    )


def dft2(table, zeta, prime):
    """Unnormalized two-dimensional DFT on F_13^2."""
    return tuple(tuple(
        sum(
            table[a][b]
            * pow(zeta, (-k * a - ell * b) % P, prime)
            for a in range(P) for b in range(P)
        ) % prime
        for ell in range(P)
    ) for k in range(P))


def main():
    fibre_histogram = Counter()
    lifted_harmonic_points = 0
    for k_source in range(P):
        for k_target in range(P):
            size = len(compatible_bits(k_source, k_target))
            fibre_histogram[size] += 1
            lifted_harmonic_points += size
    require(fibre_histogram == Counter({1: 144, 2: 24, 4: 1}),
            "allocation-lift fibre histogram")
    require(lifted_harmonic_points == 196,
            "allocation-lift cardinality over the harmonic plane")
    require(set(compatible_bits(0, 0)) == {
        (0, 0), (1, 0), (0, 1), (1, 1)
    }, "central allocation K4 fibre")
    allocation_lift_with_q = P**2 * lifted_harmonic_points
    require(allocation_lift_with_q == 33124,
            "allocation lift with target-difference coordinate")

    # Use the ordered endpoint frame s=(1,0), t=(0,1).  Any other ordered
    # basis gives the same theorem after an invertible linear reparametrization.
    source_axis = (1, 0)
    target_axis = (0, 1)
    sector_profiles = {}
    charged_profiles = 0
    for q0 in range(P):
        for q1 in range(P):
            q = (q0, q1)
            deltas = (0,) if q == (0, 0) else range(P)
            a = det(q, source_axis)
            b = det(q, target_axis)
            for delta in deltas:
                profile = (
                    delta,
                    (delta + a) % P,
                    (delta + b) % P,
                    (delta + a + b) % P,
                )
                require(mixed_rectangle(profile) == 0,
                        "affine determinant profile has nonzero curvature")
                difference = marked_difference(profile)
                require(difference == ((-2 * b) % P,
                                       (-2 * a) % P, 0),
                        "determinant-profile Hadamard formula")
                is_charged = len(set(difference)) != 1
                require(is_charged == (q != (0, 0)),
                        "marked charge does not detect nonzero q")
                charged_profiles += is_charged
                require(profile not in sector_profiles,
                        "sector-to-profile map is not injective")
                sector_profiles[profile] = (q, delta)
    require(len(sector_profiles) == 2185,
            "admissible determinant-sector count")
    require(charged_profiles == 2184,
            "charged nonzero-target sector count")

    # Explicit inverse for the chosen frame.
    for profile, (q, delta) in sector_profiles.items():
        a = (profile[1] - profile[0]) % P
        b = (profile[2] - profile[0]) % P
        recovered_q = (b, (-a) % P)
        require(recovered_q == q and profile[0] == delta,
                "sector-profile inverse")

    # The endpoint origin is a genuine two-dimensional sidecar.  For q!=0,
    # every determinant occurs on exactly thirteen of the 169 origins.
    for q0 in range(P):
        for q1 in range(P):
            q = (q0, q1)
            counts = Counter(
                det(q, (r0, r1))
                for r0 in range(P) for r1 in range(P)
            )
            if q == (0, 0):
                require(counts == Counter({0: P**2}),
                        "zero-target determinant fibre")
            else:
                require(counts == Counter({delta: P for delta in range(P)}),
                        "nonzero-target determinant fibre")

    # A genuinely two-sided endpoint-allocation square has a bilinear term
    # that the fixed-q affine dictionary cannot see.  For
    #
    #   F(epsilon,eta)=det(L+epsilon*s,R+eta*t),
    #
    # its Boolean mixed rectangle is det(s,t), independently of L and R.
    # Check all increment pairs at one base and every base at one normalized
    # frame; the theorem gives the division-free bilinear proof in general.
    frame_area_counts = Counter()
    fixed_left = (4, 7)
    fixed_right = (9, 3)
    for s0 in range(P):
        for s1 in range(P):
            source_increment = (s0, s1)
            for t0 in range(P):
                for t1 in range(P):
                    target_increment = (t0, t1)
                    area = det(source_increment, target_increment)
                    values = (
                        det(fixed_left, fixed_right),
                        det(((fixed_left[0] + s0) % P,
                             (fixed_left[1] + s1) % P), fixed_right),
                        det(fixed_left,
                            ((fixed_right[0] + t0) % P,
                             (fixed_right[1] + t1) % P)),
                        det(((fixed_left[0] + s0) % P,
                             (fixed_left[1] + s1) % P),
                            ((fixed_right[0] + t0) % P,
                             (fixed_right[1] + t1) % P)),
                    )
                    require(mixed_rectangle(values) == area,
                            "endpoint-parallelogram determinant curvature")
                    if area:
                        frame_area_counts[area] += 1
    require(frame_area_counts
            == Counter({area: 2184 for area in range(1, P)}),
            "ordered endpoint-frame area census")

    for l0 in range(P):
        for l1 in range(P):
            left = (l0, l1)
            for r0 in range(P):
                for r1 in range(P):
                    right = (r0, r1)
                    values = (
                        det(left, right),
                        det(((l0 + 1) % P, l1), right),
                        det(left, (r0, (r1 + 1) % P)),
                        det(((l0 + 1) % P, l1),
                            (r0, (r1 + 1) % P)),
                    )
                    require(mixed_rectangle(values) == 1,
                            "normalized frame curvature depends on endpoint")

    marker_root = 5
    normalized_corrections = tuple((-marker_root) % P for _ in range(7))
    require(all((marker_root + correction) % P == 0
                for correction in normalized_corrections),
            "normalized endpoint-square pointwise root correction")
    require(sum(normalized_corrections) % P
            == (-7 * marker_root) % P,
            "normalized endpoint-square Cech invoice")

    # A common endpoint atom gives a rank-one 2x2 allocation square.  Check
    # the Segre and marked-Hadamard identities on every F_13 quadruple.
    census = Counter()
    nonzero_endpoint_census = Counter()
    pluecker_gate_count = 0
    pluecker_gate_all_endpoint_factors_nonzero = 0
    for p0 in range(P):
        for p1 in range(P):
            for q0 in range(P):
                for q1 in range(P):
                    values = (
                        p0 * q0 % P,
                        p1 * q0 % P,
                        p0 * q1 % P,
                        p1 * q1 % P,
                    )
                    require((values[0] * values[3]
                             - values[1] * values[2]) % P == 0,
                            "Segre rank-one identity")
                    difference = marked_difference(values)
                    invariant = sum(values) % P
                    expected = (
                        (p0 + p1) * (q0 - q1) % P,
                        (p0 - p1) * (q0 + q1) % P,
                        (p0 - p1) * (q0 - q1) % P,
                    )
                    require(difference == expected,
                            "Segre-Hadamard factorization")
                    require(invariant == (p0 + p1) * (q0 + q1) % P,
                            "Segre-Hadamard invariant coordinate")
                    require(invariant * difference[2] % P
                            == difference[0] * difference[1] % P,
                            "Segre-Hadamard Pluecker identity")
                    require(difference[2] == mixed_rectangle(values),
                            "third marked difference is not the mixed face")
                    pluecker_gate = (
                        invariant != 0
                        and difference[0] != 0
                        and difference[1] != 0
                    )
                    require(not pluecker_gate or difference[2] != 0,
                            "three-coordinate Pluecker gate lost mixed face")
                    pluecker_gate_count += pluecker_gate
                    if (pluecker_gate
                            and all(value != 0
                                    for value in (p0, p1, q0, q1))):
                        pluecker_gate_all_endpoint_factors_nonzero += 1
                    charge = len(set(difference)) != 1
                    mixed = difference[2] != 0
                    key = ("charged" if charge else "flat",
                           "mixed" if mixed else "mixed_zero")
                    census[key] += 1
                    if all(value != 0 for value in (p0, p1, q0, q1)):
                        nonzero_endpoint_census[key] += 1
                    flat_formula = (
                        p1 * (q0 - q1) % P == 0
                        and q1 * (p0 - p1) % P == 0
                    )
                    require((not charge) == flat_formula,
                            "sharp marked-flatness criterion")

    require(census == Counter({
        ("charged", "mixed"): 24192,
        ("charged", "mixed_zero"): 3744,
        ("flat", "mixed"): 144,
        ("flat", "mixed_zero"): 481,
    }), "complete Segre-Hadamard census")
    require(nonzero_endpoint_census == Counter({
        ("charged", "mixed"): 17424,
        ("charged", "mixed_zero"): 3168,
        ("flat", "mixed_zero"): 144,
    }), "all-endpoint-factors-nonzero census")
    require(pluecker_gate_count == 20736,
            "Pluecker mixed-face gate census")
    require(pluecker_gate_all_endpoint_factors_nonzero == 14400,
            "all-endpoint-factors-nonzero Pluecker gate census")

    # Exact finite control for the virtual-face/common-cosupport distinction.
    # F_53 contains a primitive thirteenth root.  Taking p=q=delta_0 and
    # nonzero bare constants realizes full one-dimensional DFT support while
    # keeping the four corner arrays elementary.
    control_prime = 53
    zeta = pow(2, (control_prime - 1) // P, control_prime)
    require(pow(zeta, P, control_prime) == 1
            and all(pow(zeta, exponent, control_prime) != 1
                    for exponent in range(1, P)),
            "F53 thirteenth root order")
    bare = tuple(tuple(1 for _b in range(P)) for _a in range(P))
    source = tuple(tuple(1 if a == 0 else 0 for _b in range(P))
                   for a in range(P))
    target = tuple(tuple(1 if b == 0 else 0 for b in range(P))
                   for _a in range(P))
    both = tuple(tuple(1 if a == b == 0 else 0 for b in range(P))
                 for a in range(P))
    virtual = tuple(tuple(
        (bare[a][b] - source[a][b] - target[a][b] + both[a][b])
        % control_prime
        for b in range(P)
    ) for a in range(P))
    transforms = tuple(
        dft2(table, zeta, control_prime)
        for table in (bare, source, target, both, virtual)
    )
    supports = tuple(tuple(
        (k, ell) for k in range(P) for ell in range(P)
        if transform[k][ell]
    ) for transform in transforms)
    require(supports[0] == ((0, 0),),
            "bare-corner transform support")
    require(supports[1] == tuple((k, 0) for k in range(P)),
            "source-corner transform support")
    require(supports[2] == tuple((0, ell) for ell in range(P)),
            "target-corner transform support")
    require(len(supports[3]) == P**2,
            "both-corner transform support")
    require(set(supports[0]).intersection(
        supports[1], supports[2], supports[3]
    ) == {(0, 0)}, "four-corner common transform support")
    mixed_addresses = tuple(
        (k, ell) for k in range(1, P) for ell in range(1, P)
    )
    require(all(
        transforms[0][k][ell] == transforms[1][k][ell]
        == transforms[2][k][ell] == 0
        and transforms[4][k][ell] == transforms[3][k][ell] != 0
        for k, ell in mixed_addresses
    ), "virtual mixed face is not the both-corner transform")
    require(all(
        (lambda h: (
            h != 0
            and marked_difference((0, 0, 0, h), control_prime)
                == ((-h) % control_prime, (-h) % control_prime, h)
            and h * h % control_prime
                == ((-h) % control_prime) * ((-h) % control_prime)
                   % control_prime
        ))(transforms[3][k][ell])
        for k, ell in mixed_addresses
    ), "one-corner spectral Pluecker hostile")

    charged_mixed_zero = (1, 1, 2, 2)
    require(marked_difference(charged_mixed_zero) == (11, 0, 0)
            and mixed_rectangle(charged_mixed_zero) == 0,
            "charged-but-mixed-zero hostile")
    mixed_positive = (1, 2, 3, 6)
    require(marked_difference(mixed_positive) == (7, 9, 2)
            and mixed_rectangle(mixed_positive) == 2,
            "two-sided mixed positive control")

    full_carrier_address = P**4
    endpoint_plane = P**4
    common_pullback = P**6
    pullback_origins_per_carrier_address = P**2
    require((full_carrier_address, endpoint_plane, common_pullback,
             pullback_origins_per_carrier_address)
            == (28561, 28561, 4826809, 169),
            "carrier/endpoint pullback cardinalities")

    print("THM-2772 carrier-allocation pullback and mixed-face certificate")
    print(f"allocation_lift_fibre_histogram={dict(sorted(fibre_histogram.items()))}")
    print(f"allocation_lift_harmonic_points={lifted_harmonic_points}")
    print(f"allocation_lift_with_q={allocation_lift_with_q}")
    print("central_zero_harmonic_fibre=K4=(bare,source,target,both)")
    print(f"admissible_sector_profiles={len(sector_profiles)}")
    print(f"charged_affine_sector_profiles={charged_profiles}")
    print("affine_profile=(Delta,Delta+a,Delta+b,Delta+a+b)")
    print("marked_difference_affine=(-2b,-2a,0); mixed_rectangle=0")
    print("endpoint_parallelogram_mixed_face=det(source_increment,target_increment)")
    print("ordered_frame_area_census="
          f"{dict(sorted(frame_area_counts.items()))}")
    print("normalized_det1_frames=2184; c_j=-a; sum_j_c_j=-7a")
    print(f"Segre_Hadamard_census={dict(sorted(census.items()))}")
    print("Segre_Hadamard_all_endpoint_factors_nonzero="
          f"{dict(sorted(nonzero_endpoint_census.items()))}")
    print(f"Pluecker_gate_D0_D1_D2_nonzero={pluecker_gate_count}")
    print("Pluecker_gate_all_endpoint_factors_nonzero="
          f"{pluecker_gate_all_endpoint_factors_nonzero}")
    print("virtual_face_control_supports_(bare,source,target,both,Omega)="
          f"{tuple(len(support) for support in supports)}")
    print(f"virtual_face_mixed_addresses={len(mixed_addresses)}; "
          "at_each: bare=source=target=0 and Omega=both!=0")
    print("spectral_one_corner_Pluecker_gate=144/144; "
          "common_atom_required_before_DFT")
    print("charged_mixed_zero_hostile="
          f"{charged_mixed_zero}; D={marked_difference(charged_mixed_zero)}")
    print("two_sided_mixed_positive="
          f"{mixed_positive}; D={marked_difference(mixed_positive)}")
    print(f"carrier_address_size={full_carrier_address}")
    print(f"endpoint_plane_size={endpoint_plane}")
    print(f"canonical_common_pullback_size={common_pullback}")
    print("pullback_endpoint_origins_per_carrier_address="
          f"{pullback_origins_per_carrier_address}")
    print("POSITIVE CONDITIONAL: a common normalized endpoint-translation "
          "square has determinant mixed face 1 and gives c_j=-a, so seven "
          "coherent clocks pay sum_j c_j=-7a exactly.")
    print("CONCLUSION: the determinant-sector affine K4 is charged on every q!=0 sector but has zero mixed curvature.  A physical THM-2591 correction must retain the Boolean allocation bits and endpoint origin on one common atom and make the third Segre-Hadamard component survive before paying the seven-clock Cech invoice.")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
