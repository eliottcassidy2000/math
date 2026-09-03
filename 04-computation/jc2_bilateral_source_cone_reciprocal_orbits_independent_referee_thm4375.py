"""Import-free clean-room verifier for the proposed THM-4375.

This script deliberately imports nothing from the repository or the Python
standard library.  It reconstructs source fibres from exponent tuples and
checks the proposed cone, orbit, block, generating-function, primitive-ray,
and fibre statements against those tuples on an explicit finite universe.

The proof of the unbounded statements is recorded separately in REPORT.md;
the finite checks here are hostile controls and independent formula audits.
"""


CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise RuntimeError("FAIL: " + label)
    CHECKS += 1


def ceil_div(a, b):
    return -((-a) // b)


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def tri(k):
    return k * (k + 1) // 2


def s_rho(ell):
    return (ell + 1) // 2, (ell + 2) // 3


def address(u, v):
    tau = u + v - 1
    return tri(tau - 1) + u


def decode_address(value):
    tau = 1
    while tri(tau) < value:
        tau += 1
    u = value - tri(tau - 1)
    v = tau + 1 - u
    return u, v, tau


def source_tuples(ell, u, v):
    """Reconstruct all (a,b,c,e) from the defining packet equations."""
    s, unused_rho = s_rho(ell)
    del unused_rho
    n_value = v
    n_zero = s + u - 1
    tuples = []
    if n_value <= 0:
        return tuples
    for e in range(n_value + 1):
        c = n_value - e
        a = 2 * c + 3 * e - ell
        b = n_zero - c - 2 * e
        if a >= 0 and b >= 0:
            tuples.append((a, b, c, e))
    return tuples


def source_wall_formula(ell, u, v):
    s, rho = s_rho(ell)
    return v >= rho and u - v >= 1 - s and u + v >= ell - s + 1


def bilateral_direct(ell, u, v):
    return bool(source_tuples(ell, u, v)) and bool(source_tuples(ell, v, u))


def bilateral_band(ell, u, v):
    s, rho = s_rho(ell)
    return u >= rho and v >= rho and abs(u - v) <= s - 1


def fibre_closed(ell, u, v):
    s, unused_rho = s_rho(ell)
    del unused_rho
    n_value = v
    n_zero = s + u - 1
    return min(n_value, n_zero - n_value) - max(0, ell - 2 * n_value) + 1


def coefficient_by_orbits(ell, tau):
    s, rho = s_rho(ell)
    tau_zero = 2 * rho - 1
    count = 0
    for d in range(s):
        remainder = tau - tau_zero - d
        if remainder >= 0 and remainder % 2 == 0:
            count += 1
    return count


def block_types(ell, tau):
    total = tau + 1
    return [(u, total - u) for u in range(1, total) if bilateral_band(ell, u, total - u)]


def coeff(values, index):
    if index < 0 or index >= len(values):
        return 0
    return values[index]


def check_inversion_and_source_cone():
    # A direct forward scan of exponent triples is independent of the inverse
    # e-parametrization used by source_tuples.
    cap = 26
    for ell in range(2, 21):
        s, unused_rho = s_rho(ell)
        del unused_rho
        seen = set()
        max_n_zero = s + cap - 1
        for c in range(cap + 1):
            for e in range(cap + 1 - c):
                if c + e == 0:
                    continue
                a = 2 * c + 3 * e - ell
                if a < 0:
                    continue
                for b in range(max_n_zero + 1):
                    n_value = c + e
                    n_zero = b + c + 2 * e
                    u = n_zero - s + 1
                    v = n_value
                    if 1 <= u <= cap and 1 <= v <= cap:
                        seen.add((u, v, a, b, c, e))
        for u in range(1, cap + 1):
            for v in range(1, cap + 1):
                forward = sorted(item[2:] for item in seen if item[0] == u and item[1] == v)
                inverse = source_tuples(ell, u, v)
                check(forward == inverse, "forward/inverse exponent reconstruction")

    # Larger wall and bilateral audit, including points on and just outside
    # each of the three source walls.
    for ell in range(2, 61):
        s, rho = s_rho(ell)
        check(2 * rho >= ell - s + 1, "bilateral sum wall is redundant")
        for u in range(1, 81):
            for v in range(1, 81):
                direct_source = bool(source_tuples(ell, u, v))
                check(direct_source == source_wall_formula(ell, u, v), "source cone from tuples")
                direct_bilateral = direct_source and bool(source_tuples(ell, v, u))
                check(direct_bilateral == bilateral_band(ell, u, v), "bilateral cone equivalence")
        # Native relation has loops and a missing long pair, and hence is not
        # a tournament.  When s >= 2 it also has both orientations of a
        # neighboring pair rather than exactly one.
        check(bilateral_band(ell, rho, rho), "looped path-power loop")
        check(not bilateral_band(ell, rho, rho + s), "path-power missing long edge")
        if s >= 2:
            check(bilateral_band(ell, rho, rho + 1), "path-power forward neighbor")
            check(bilateral_band(ell, rho + 1, rho), "path-power reverse neighbor")


def check_orbits_blocks_and_series():
    max_tau = 240
    for ell in range(2, 71):
        s, rho = s_rho(ell)
        tau_zero = 2 * rho - 1
        orbit_counts = [0] * (max_tau + 1)
        oriented_counts = [0] * (max_tau + 1)
        fixed_counts = [0] * (max_tau + 1)

        for tau in range(1, max_tau + 1):
            total = tau + 1
            types = block_types(ell, tau)
            lower = max(rho, ceil_div(total - (s - 1), 2))
            upper = min(total - rho, (total + (s - 1)) // 2)
            interval = [] if lower > upper else [(u, total - u) for u in range(lower, upper + 1)]
            check(types == interval, "centered block interval")
            check(lower + upper == total, "formal reciprocal block bounds")
            if types:
                addresses = [address(u, v) for u, v in types]
                check(addresses == list(range(tri(tau - 1) + lower, tri(tau - 1) + upper + 1)),
                      "consecutive block addresses")
                check(types[0] == (lower, upper) and types[-1] == (upper, lower),
                      "reciprocal nonempty endpoints")
                check(addresses[0] + addresses[-1] == tau * tau + 1,
                      "reciprocal endpoint address sum")
            for u, v in types:
                decoded_u, decoded_v, decoded_tau = decode_address(address(u, v))
                check((decoded_u, decoded_v, decoded_tau) == (u, v, tau), "address decoder")
                check(address(u, v) + address(v, u) == tau * tau + 1, "address reflection")

            oriented = len(types)
            fixed = sum(1 for u, v in types if u == v)
            orbits = sum(1 for u, v in types if u <= v)
            orbit_counts[tau] = orbits
            oriented_counts[tau] = oriented
            fixed_counts[tau] = fixed

            parameter_orbits = coefficient_by_orbits(ell, tau)
            k = tau - tau_zero
            closed_orbits = 0
            if k >= 0:
                closed_orbits = sum(1 for d in range(min(s - 1, k) + 1) if d % 2 == k % 2)
            expected_fixed = 1 if k >= 0 and k % 2 == 0 else 0
            check(orbits == parameter_orbits == closed_orbits, "orbit coefficient extraction")
            check(fixed == expected_fixed, "fixed coefficient extraction")
            check(oriented == 2 * orbits - expected_fixed, "oriented coefficient extraction")

        # Verify all three displayed rational generating functions as exact
        # coefficient recurrences, not through symbolic-algebra imports.
        for tau in range(max_tau + 1):
            lhs_o = (coeff(orbit_counts, tau) - coeff(orbit_counts, tau - 1)
                     - coeff(orbit_counts, tau - 2) + coeff(orbit_counts, tau - 3))
            rhs_o = (1 if tau == tau_zero else 0) - (1 if tau == tau_zero + s else 0)
            check(lhs_o == rhs_o, "orbit generating function")

            lhs_m = (coeff(oriented_counts, tau) - coeff(oriented_counts, tau - 1)
                     - coeff(oriented_counts, tau - 2) + coeff(oriented_counts, tau - 3))
            rhs_m = ((1 if tau == tau_zero else 0)
                     + (1 if tau == tau_zero + 1 else 0)
                     - 2 * (1 if tau == tau_zero + s else 0))
            check(lhs_m == rhs_m, "oriented generating function")

            lhs_f = coeff(fixed_counts, tau) - coeff(fixed_counts, tau - 2)
            rhs_f = 1 if tau == tau_zero else 0
            check(lhs_f == rhs_f, "fixed generating function")

        # After every d-stream has begun, two more blocks add exactly one
        # point on each d-stream.  This is the exact periodic-error mechanism
        # behind both cumulative asymptotics.
        cumulative_o = [0] * (max_tau + 1)
        cumulative_m = [0] * (max_tau + 1)
        for tau in range(1, max_tau + 1):
            cumulative_o[tau] = cumulative_o[tau - 1] + orbit_counts[tau]
            cumulative_m[tau] = cumulative_m[tau - 1] + oriented_counts[tau]
        threshold = tau_zero + s - 1
        for tau in range(threshold, max_tau - 1):
            check(cumulative_o[tau + 2] - cumulative_o[tau] == s,
                  "orbit cumulative two-block slope")
            check(cumulative_m[tau + 2] - cumulative_m[tau] == 2 * s - 1,
                  "oriented cumulative two-block slope")

        # Canonical orbit parametrization, exact address separation/midpoint,
        # and the fixed-address formula.
        for j in range(81):
            for d in range(s):
                w = rho + j
                u, v = w + d, w
                tau = 2 * rho - 1 + 2 * j + d
                check(bilateral_band(ell, u, v), "canonical orbit is bilateral")
                check(min(u, v) == w and abs(u - v) == d, "unique orbit parameters")
                a_low = tri(tau - 1) + rho + j
                a_high = a_low + d
                check({address(u, v), address(v, u)} == {a_low, a_high}, "orbit address pair")
                check(a_low + a_high == tau * tau + 1, "orbit address midpoint")
                check(2 * a_low == tau * tau + 1 - d, "lower midpoint formula")
                check(2 * a_high == tau * tau + 1 + d, "upper midpoint formula")

        for w in range(rho, rho + 121):
            tau = 2 * w - 1
            fixed_address = address(w, w)
            check(fixed_address == w * w + (w - 1) * (w - 1), "fixed address formula")
            check(tau % 2 == 1 and bilateral_band(ell, w, w), "fixed odd block")
            check(block_types(ell, tau).count((w, w)) == 1, "unique fixed type in odd block")

        # Exact finite path-power orbit count, including loops.
        for n in range(1, 61):
            vertices = list(range(rho, rho + n))
            direct = sum(1 for i, u in enumerate(vertices) for v in vertices[i:]
                         if bilateral_band(ell, u, v))
            diameter = min(s - 1, n - 1)
            closed = (diameter + 1) * n - diameter * (diameter + 1) // 2
            check(direct == closed, "finite path-power orbit count")


def check_primitive_rays():
    for ell in range(2, 61):
        s, rho = s_rho(ell)
        for p in range(1, 31):
            for q in range(1, 31):
                if gcd(p, q) != 1:
                    continue
                if p == q:
                    check((p, q) == (1, 1), "only primitive diagonal ray")
                    for scale in range(1, 91):
                        check(bilateral_band(ell, scale * p, scale * q) == (scale >= rho),
                              "diagonal ray scale tail")
                else:
                    lower = ceil_div(rho, min(p, q))
                    upper = (s - 1) // abs(p - q)
                    direct_scales = [scale for scale in range(1, 91)
                                     if bilateral_band(ell, scale * p, scale * q)]
                    closed_scales = [] if lower > upper else list(range(lower, upper + 1))
                    check(direct_scales == closed_scales, "primitive nonfixed scale interval")
                    check(len(closed_scales) == max(0, upper - lower + 1), "scale interval size")


def check_fibres_and_imbalance():
    for ell in range(2, 101):
        s, rho = s_rho(ell)
        sharp_blocks = []
        upper_w = max(rho + 12, 2 * s + 8)
        for w in range(rho, upper_w + 1):
            for d in range(s):
                plus_tuples = source_tuples(ell, w + d, w)
                minus_tuples = source_tuples(ell, w, w + d)
                mu_plus = len(plus_tuples)
                mu_minus = len(minus_tuples)

                closed_plus = min(w, s - 1 + d) - max(0, ell - 2 * w) + 1
                closed_minus = (min(w + d, s - 1 - d)
                                - max(0, ell - 2 * w - 2 * d) + 1)
                defect_plus = s + d - max(0, s - 1 + d - w) - max(0, ell - 2 * w)
                defect_minus = (s - d - max(0, s - 1 - w - 2 * d)
                                - max(0, ell - 2 * w - 2 * d))

                check(bilateral_band(ell, w + d, w), "fibre orbit lies in bilateral cone")
                check(mu_plus == fibre_closed(ell, w + d, w) == closed_plus == defect_plus,
                      "plus fibre formulas")
                check(mu_minus == fibre_closed(ell, w, w + d) == closed_minus == defect_minus,
                      "minus fibre formulas")
                check(mu_plus > 0 and mu_minus > 0, "bilateral fibres are positive")
                check(mu_plus + mu_minus <= 2 * s, "fibre conservation ceiling")

                sum_equality = mu_plus + mu_minus == 2 * s
                closed_sum_equality = w >= max(s, s - 1 + d)
                check(sum_equality == closed_sum_equality, "fibre-sum equality classification")
                if closed_sum_equality:
                    check((mu_plus, mu_minus) == (s + d, s - d), "stable endpoint fibres")

                difference = abs(mu_plus - mu_minus)
                check(difference <= 2 * s - 2, "sharp asymmetry bound")
                if s == 1:
                    check(d == 0 and difference == 0, "ell=2 fixed-only asymmetry")
                else:
                    sharp = difference == 2 * s - 2
                    closed_sharp = d == s - 1 and w >= 2 * s - 2
                    check(sharp == closed_sharp, "sharp asymmetry equality classification")
                    if sharp:
                        check((mu_plus, mu_minus) == (2 * s - 1, 1), "sharp ordered fibres")
                        sharp_blocks.append(2 * w + d - 1)

        if s >= 2:
            check(min(sharp_blocks) == 5 * s - 6, "first sharp-asymmetry block")


def check_low_order_and_named_hostiles():
    # ell=2: only fixed points, identical generating functions, singleton fibres.
    ell = 2
    for w in range(1, 81):
        check(bilateral_band(ell, w, w), "ell=2 fixed ray")
        check(len(source_tuples(ell, w, w)) == 1, "ell=2 singleton fibre")
        check(not bilateral_band(ell, w, w + 1), "ell=2 no nonfixed arc")
    for tau in range(1, 101):
        expected = 1 if tau % 2 == 1 else 0
        types = block_types(ell, tau)
        check(len(types) == expected, "ell=2 common series coefficient")

    # ell=3: one orbit per block; oriented counts alternate 1,2; named fibres.
    ell = 3
    for tau in range(1, 101):
        types = block_types(ell, tau)
        check(sum(1 for u, v in types if u <= v) == 1, "ell=3 one orbit per block")
        check(len(types) == (1 if tau % 2 == 1 else 2), "ell=3 oriented alternation")
    check((len(source_tuples(3, 1, 2)), len(source_tuples(3, 2, 1))) == (1, 1),
          "ell=3 first nonfixed balanced fibres")
    check((len(source_tuples(3, 2, 3)), len(source_tuples(3, 3, 2))) == (1, 3),
          "ell=3 first sharp imbalance")

    # ell=10 fixed control and first reciprocal-fibre hostile.
    ell = 10
    fixed = source_tuples(ell, 4, 4)
    left = source_tuples(ell, 4, 5)
    right = source_tuples(ell, 5, 4)
    check(address(4, 4) == 25, "ell=10 first fixed address")
    check(fixed == [(0, 2, 2, 2), (1, 1, 1, 3), (2, 0, 0, 4)],
          "ell=10 fixed fibre tuples")
    check((address(4, 5), address(5, 4)) == (32, 33), "ell=10 hostile addresses")
    check(left == [(0, 3, 5, 0), (1, 2, 4, 1), (2, 1, 3, 2), (3, 0, 2, 3)],
          "ell=10 address-32 fibre")
    check(right == [(0, 3, 2, 2), (1, 2, 1, 3), (2, 1, 0, 4)],
          "ell=10 address-33 fibre")

    expected_ray = [
        (1, 4, 5, 32, 33, 4, 3),
        (2, 8, 10, 144, 146, 3, 7),
        (3, 12, 15, 337, 340, 2, 8),
        (4, 16, 20, 611, 615, 1, 9),
    ]
    actual_ray = []
    for scale in range(1, 10):
        u, v = 4 * scale, 5 * scale
        if bilateral_band(ell, u, v):
            actual_ray.append((scale, u, v, address(u, v), address(v, u),
                               len(source_tuples(ell, u, v)),
                               len(source_tuples(ell, v, u))))
    check(actual_ray == expected_ray, "ell=10 primitive 4/5 scale ladder")


check_inversion_and_source_cone()
check_orbits_blocks_and_series()
check_primitive_rays()
check_fibres_and_imbalance()
check_low_order_and_named_hostiles()

OUTPUT_LINES = [
    "THM-4375 CLEAN-ROOM INDEPENDENT REFEREE",
    "universe=ell_cone_2..60,u_v_1..80;ell_blocks_2..70,tau_1..240;",
    "         primitive_p_q_1..30,scale_1..90;ell_fibres_2..100",
    "method=import-free exponent tuples + direct path/block/orbit counts + coefficient recurrences",
    "ell=2 control=fixed-only,singleton-fibres,common-series z/(1-z^2)",
    "ell=3 hostile=(2,3)<->(3,2),fibres=1<->3,first-sharp-block=4",
    "ell=10 hostile=Addr32<->Addr33,fibres=4<->3",
    "ell=10 ray=4/5,scales=1..4,fibres=(4,3),(3,7),(2,8),(1,9)",
    "carrier=looped path-power with absent long pairs and honest fixed ties; not a tournament",
    "scope=no JC(2) or DC(2) consequence",
    "checks=" + str(CHECKS),
    "PASS",
]
open(1, "wb", buffering=0, closefd=False).write(("\n".join(OUTPUT_LINES) + "\n").encode("ascii"))
