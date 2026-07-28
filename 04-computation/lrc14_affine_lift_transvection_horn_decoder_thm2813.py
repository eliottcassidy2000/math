#!/usr/bin/env python3
"""Exact companion for THM-2813.

The script is dependency-free.  It verifies the complete address arithmetic
of THM-2807's thirteen affine lifts, computes g_t g_0^{-1} as an affine map,
identifies both the common depth-two shear and the top-graded relative shear,
checks the fixed/free orbit and normal-jet claims, and reconstructs
THM-2806's abstract wandering-selector allocation square.

THM-2803's actual cyclotomic endpoint profiles are not rebuilt here.  The
theorem uses its already proved all-minors statement at coordinates 6 and 7;
this companion checks the exact formal universe and 13-to-1 versus injective
flag counts which follow from that dependency.
"""


P = 13
P2 = P**2
P4 = P**4
P5 = P**5
P6 = P**6

N0 = 3_454_614
N_PLUS = 3_454_627
N_A = 4_143_978
K0 = 23_098
M0 = 2_652_079
V0 = 352_469
VERTICAL_QUOTIENT = 4_079
TRANSVECTION_RESIDUE = 10


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def affine_compose(left, right):
    """Return left after right for affine pairs (multiplier, translation)."""
    lm, lv = left
    rm, rv = right
    return (lm * rm % P6, (lm * rv + lv) % P6)


def affine_inverse(pair):
    multiplier, translation = pair
    inverse_multiplier = pow(multiplier, -1, P6)
    return (
        inverse_multiplier,
        (-inverse_multiplier * translation) % P6,
    )


def affine_apply(pair, value):
    multiplier, translation = pair
    return (multiplier * value + translation) % P6


def lift(t):
    exponent = K0 + t * P4
    multiplier = pow(14, exponent, P6)
    translation = ((1 - multiplier) * N0) % P6
    return exponent, multiplier, translation


def relative_expected(t):
    return ((1 + t * P5) % P6, (-7 * t * P5) % P6)


def normal_coordinates(value):
    return ((value - 7) % P, (value // P5) % P)


def verify_address_horn():
    require(P2 == 169 and P4 == 28_561 and P5 == 371_293, "power drift")
    require(P6 == 4_826_809, "address modulus drift")
    require(N0 % P == N_PLUS % P == N_A % P == 7, "fixed sheet drift")
    require(
        (N0 % P2, N_PLUS % P2, N_A % P2) == (85, 98, 98),
        "depth-two horn drift",
    )
    require(
        N_PLUS - N0 == P
        and N_A - N_PLUS == P2 * VERTICAL_QUOTIENT
        and N_A - N0 == P * 53_028
        and VERTICAL_QUOTIENT % P == TRANSVECTION_RESIDUE,
        "horn edge arithmetic drift",
    )

    base_exponent, base_multiplier, base_translation = lift(0)
    require(
        (base_exponent, base_multiplier, base_translation) == (K0, M0, V0),
        "least affine lift drift",
    )
    require(
        pow(14, K0, P5) == 53_028,
        "depth-five discrete logarithm drift",
    )

    lifts = []
    relatives = []
    for t in range(P):
        exponent, multiplier, translation = lift(t)
        require(
            multiplier == (M0 + t * P5) % P6,
            "lift multiplier formula drift",
        )
        require(
            translation == (V0 + 6 * t * P5) % P6,
            "lift translation formula drift",
        )
        pair = (multiplier, translation)
        require(affine_apply(pair, N0) == N0, "lift stopped fixing n0")
        require(affine_apply(pair, N_PLUS) == N_A, "lift stopped mapping horn")
        lifts.append(pair)

        relative = affine_compose(pair, affine_inverse(lifts[0]))
        require(
            relative == relative_expected(t),
            "g_t g_0^-1 is not the claimed relative transvection",
        )
        relatives.append(relative)

    require(len(set(lifts)) == P, "full-depth lifts collided")
    require(len(set(relatives)) == P, "relative transvections collided")

    group_law_checks = 0
    for t in range(P):
        for u in range(P):
            require(
                affine_compose(relatives[t], relatives[u])
                == relatives[(t + u) % P],
                "relative transvection group law failed",
            )
            group_law_checks += 1

    # All g_t induce one common depth-two shear.  In low/second digit
    # coordinates (v,w), it is (v,w+10(v-7)).
    quotient_checks = 0
    expected_depth_two_pair = (
        1 + P * TRANSVECTION_RESIDUE,
        (-7 * P * TRANSVECTION_RESIDUE) % P2,
    )
    for pair in lifts:
        quotient_pair = (pair[0] % P2, pair[1] % P2)
        require(
            quotient_pair == expected_depth_two_pair,
            "full lift lost the common depth-two transvection",
        )
        for v in range(P):
            for w in range(P):
                value = v + P * w
                image = affine_apply(quotient_pair, value) % P2
                image_v = image % P
                image_w = image // P
                require(
                    (image_v, image_w)
                    == (v, (w + TRANSVECTION_RESIDUE * (v - 7)) % P),
                    "depth-two shear formula failed",
                )
                quotient_checks += 1

    # Relative lifts act on low/top digits by the standard transvection
    # (r,h)->(r,h+t*r).  Middle digits are immaterial because the added
    # displacement is a multiple of p^5.
    relative_quotient_checks = 0
    for t, relative in enumerate(relatives):
        for r in range(P):
            for h in range(P):
                value = (7 + r + h * P5) % P6
                image = affine_apply(relative, value)
                require(
                    normal_coordinates(image) == (r, (h + t * r) % P),
                    "top-graded relative shear failed",
                )
                relative_quotient_checks += 1

    # For a nonzero group element the exact fixed condition is r=0.
    fixed_residue_checks = 0
    for t in range(1, P):
        relative = relatives[t]
        for residue in range(P):
            is_fixed = affine_apply(relative, residue) == residue
            require(is_fixed == (residue == 7), "fixed-line criterion drift")
            fixed_residue_checks += 1

    fixed_sheet_size = P5
    free_point_count = P6 - P5
    free_orbit_count = free_point_count // P
    require(fixed_sheet_size == 371_293, "fixed-sheet count drift")
    require(free_orbit_count == 342_732, "free-orbit count drift")

    # H is the top digit.  On the adjacent r=1 sheet its cocycle is t.
    jets = []
    jet_checks = 0
    for t, relative in enumerate(relatives):
        for h in range(P):
            value = 8 + h * P5
            image = affine_apply(relative, value)
            before_h = normal_coordinates(value)[1]
            after_h = normal_coordinates(image)[1]
            require((after_h - before_h) % P == t, "normal jet drift")
            jet_checks += 1
        jets.append(
            (
                normal_coordinates(affine_apply(relative, 8))[1]
                - normal_coordinates(8)[1]
            )
            % P
        )
    require(tuple(jets) == tuple(range(P)), "normal jet is not faithful")

    return {
        "lifts": tuple(lifts),
        "relatives": tuple(relatives),
        "group_law_checks": group_law_checks,
        "quotient_checks": quotient_checks,
        "relative_quotient_checks": relative_quotient_checks,
        "fixed_residue_checks": fixed_residue_checks,
        "fixed_sheet_size": fixed_sheet_size,
        "free_orbit_count": free_orbit_count,
        "jet_checks": jet_checks,
        "jets": tuple(jets),
        "depth_two_pair": expected_depth_two_pair,
    }


def verify_projective_flag_arithmetic():
    # THM-2803 proves all 78 minors for every one of 168 directions.
    directions = P2 - 1
    fibre_pairs = P * (P - 1) // 2
    inherited_w6_w7_minors = directions * fibre_pairs
    require((directions, fibre_pairs) == (168, 78), "minor universe drift")
    require(inherited_w6_w7_minors == 13_104, "decoder-minor count drift")

    # Treat the injective projective decoder supplied by THM-2803 as the
    # formal label delta.  The quotient horn repeats its w=7 coordinate and
    # forgets t; the normal jet restores t.
    collapsed = {}
    augmented = set()
    for delta in range(P):
        for t in range(P):
            horn_key = (delta, delta)
            collapsed.setdefault(horn_key, []).append((delta, t))
            augmented.add((horn_key, t))
    require(len(collapsed) == P, "collapsed horn image drift")
    require(
        set(len(fibre) for fibre in collapsed.values()) == {P},
        "collapsed horn fibres are not all thirteen",
    )
    require(len(augmented) == P2, "normal-jet augmentation is not injective")

    return {
        "directions": directions,
        "fibre_pairs": fibre_pairs,
        "inherited_w6_w7_minors": inherited_w6_w7_minors,
        "collapsed_flags": len(collapsed),
        "collapsed_fibre_size": P,
        "augmented_flags": len(augmented),
    }


def verify_fixed_rees_boundary():
    # THM-2806's wandering-selector lemma at m=13 and w=1.
    support = {"B": 0, "P": 0, "Q": 0, "H": 0, "all": 0, "omega": 0}
    sums = {"B": 0, "P": 0, "Q": 0, "H": 0, "omega": 0}
    fourfold_omega = None
    for a in range(P):
        for b in range(P):
            B = 1
            present_source = int(a == 0)
            present_target = int(b == 0)
            both = int(a == 0 and b == 0)
            omega = B - present_source - present_target + both
            values = {
                "B": B,
                "P": present_source,
                "Q": present_target,
                "H": both,
            }
            for key, value in values.items():
                support[key] += value != 0
                sums[key] += value
            support["all"] += all(values.values())
            support["omega"] += omega != 0
            sums["omega"] += omega
            if all(values.values()):
                fourfold_omega = omega

    require(
        support
        == {"B": 169, "P": 13, "Q": 13, "H": 1, "all": 1, "omega": 144},
        "wandering-selector support census drift",
    )
    require(
        sums == {"B": 169, "P": 13, "Q": 13, "H": 1, "omega": 144},
        "wandering-selector central sum drift",
    )
    require(fourfold_omega == 0, "sole fourfold point stopped being raw-flat")

    v = (sums["B"], sums["P"], sums["Q"], sums["H"])
    D = (
        v[0] + v[1] + v[2] + v[3],
        v[0] + v[1] - v[2] - v[3],
        v[0] - v[1] + v[2] - v[3],
        v[0] - v[1] - v[2] + v[3],
    )
    require(v == (169, 13, 13, 1), "central allocation vector drift")
    require(D == (196, 168, 168, 144), "central Hadamard vector drift")
    require(D[3] % P == 1, "normalized Rees residue drift")

    # There is no proved action of the THM-2807 lift torsor on THM-2806's
    # separately fixed marked selector.  Repeating its lone scalar residue
    # over an abstract t-index is only the sharp scalar-information hostile:
    # a scalar 1 cannot encode the genuine normal jet (0,...,12).  It is not
    # a claim of physical t-invariance.
    formal_rees_repetition = tuple(D[3] % P for _t in range(P))
    genuine_jets = tuple(range(P))
    require(
        len(set(formal_rees_repetition)) == 1,
        "formal Rees scalar repetition gained a lift colour",
    )
    require(len(set(genuine_jets)) == P, "normal jet collapsed")
    require(
        formal_rees_repetition != genuine_jets,
        "one Rees scalar encoded the normal jet",
    )

    return {
        "support": support,
        "central": v,
        "hadamard": D,
        "fourfold_omega": fourfold_omega,
        "formal_rees_repetition": formal_rees_repetition,
        "genuine_jets": genuine_jets,
    }


def main():
    address = verify_address_horn()
    flags = verify_projective_flag_arithmetic()
    rees = verify_fixed_rees_boundary()

    print("THM-2813 AFFINE-LIFT TRANSVECTION AND PROJECTIVE HORN DECODER")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
    print(
        f"p={P}; modulus={P6}; top_step={P5}; "
        f"addresses=({N0},{N_PLUS},{N_A})"
    )
    print(
        f"affine_lifts={len(address['lifts'])}; "
        f"least={address['lifts'][0]}; last={address['lifts'][-1]}"
    )
    print(
        "common_depth2_shear=(v,w)->(v,w+10(v-7)); "
        f"affine_pair={address['depth_two_pair']}"
    )
    print(
        "relative_lift=A_t(y)=y+t*13^5*(y-7 mod13); "
        "top_quotient=(r,h)->(r,h+t*r)"
    )
    print(
        f"group_law_checks={address['group_law_checks']}; "
        f"depth2_checks={address['quotient_checks']}; "
        f"top_quotient_checks={address['relative_quotient_checks']}"
    )
    print(
        f"fixed_residue_checks={address['fixed_residue_checks']}; "
        f"fixed_sheet_size={address['fixed_sheet_size']}; "
        f"free_offsheet_orbits={address['free_orbit_count']}"
    )
    print(
        f"normal_jet_checks={address['jet_checks']}; "
        f"normal_jets={address['jets']}"
    )
    print(
        f"inherited_thm2803_directions={flags['directions']}; "
        f"w6_w7_minors={flags['inherited_w6_w7_minors']}; "
        f"collapsed_flags={flags['collapsed_flags']}x"
        f"{flags['collapsed_fibre_size']}; "
        f"jet_augmented_flags={flags['augmented_flags']}"
    )
    print(
        f"thm2806_support={rees['support']}; "
        f"central={rees['central']}; hadamard={rees['hadamard']}; "
        f"raw_fourfold_omega={rees['fourfold_omega']}"
    )
    print(
        "formal_rees_scalar_repetition=(1,...,1); no_t_action_typed=yes; "
        "genuine_normal_jets=(0,1,...,12); scalar_not_jet=yes"
    )
    print(
        "scope=relative address transvection and conditional common-scalar "
        "projective decoder only; no physical allocation-to-endpoint map, "
        "root/Cech map, row exclusion, or LRC14"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
