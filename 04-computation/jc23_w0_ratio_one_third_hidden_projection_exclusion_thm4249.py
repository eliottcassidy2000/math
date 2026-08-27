#!/usr/bin/env python3
"""Exact good-reduction exclusion of THM-4249's R=1/3 lane.

THM-4249 leaves 132 degree-42 map orbits over the unique common admissible
torsion-envelope ratio R=U/Z=1/3.  If a candidate collapses the
attachments, its hidden projector H=A*f+B*g does too and hence H(Q0)=O.

The load-bearing compression is uniform.  If delta=A^2+omega*B^2 were
divisible by pi=omega^2-1, reduction modulo pi=F_3 would force pi|A,B.
Then H=pi*H' would have degree q(H')=4K with 3 not dividing K, impossible
because every hidden degree is divisible by six.  Thus 3 does not divide
N(delta).  At one exact good
specialization to F_397, this script proves

    Ann_O(f(Q0)) = (6+12*omega),       g(Q0)=[15]f(Q0),

and `f(Q0)` has exact additive order 18.  Collapse would give
delta*f(Q0)=O, hence 18|N(delta), contradicting 3 not dividing N(delta).
The explicit 132-row annihilator check is retained as a hostile control.
Thus R=1/3 is not in S_42.  The other 34 envelope ratios,
the remaining 1,512 map-ratio incidences, W=0, M=12, and JC(2) stay open.

The shell/orbit universe is imported from THM-4249's primary exact
certificate.  The companion independent audit reconstructs it without that
import and uses a separate finite-group implementation.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


HERE = Path(__file__).resolve().parent
PARENT_PATH = HERE / "jc23_w0_cyclic_projector_squeeze_thm4249.py"
SPEC = importlib.util.spec_from_file_location("thm4249_parent", PARENT_PATH)
require(SPEC is not None and SPEC.loader is not None, "parent import failed")
parent = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(parent)


Point = tuple[int, int] | None

Q = 397
ZETA12 = 157
RHO = 161
SCALE = 27
SOURCE_X = 15
SOURCE_Y = 28
ORIGIN: Point = None


def inv(value: int) -> int:
    value %= Q
    require(value != 0, "attempted inversion of zero")
    return pow(value, Q - 2, Q)


def point_neg(point: Point) -> Point:
    if point is None:
        return None
    return point[0], (-point[1]) % Q


def point_add(left: Point, right: Point) -> Point:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 + y2) % Q == 0:
        return None
    if left == right:
        require(y1 != 0, "doubling a two-torsion point missed O")
        slope = 3 * x1 * x1 * inv(2 * y1) % Q
    else:
        slope = (y2 - y1) * inv(x2 - x1) % Q
    x3 = (slope * slope - x1 - x2) % Q
    y3 = (slope * (x1 - x3) - y1) % Q
    answer = x3, y3
    require((y3 * y3 - x3**3 - 1) % Q == 0,
            "elliptic addition left E0")
    return answer


def point_mul(multiplier: int, point: Point) -> Point:
    if multiplier < 0:
        return point_mul(-multiplier, point_neg(point))
    answer = ORIGIN
    addend = point
    remaining = multiplier
    while remaining:
        if remaining & 1:
            answer = point_add(answer, addend)
        addend = point_add(addend, addend)
        remaining //= 2
    return answer


OMEGA = pow(ZETA12, 4, Q)


def omega_action(point: Point) -> Point:
    if point is None:
        return None
    return OMEGA * point[0] % Q, point[1]


def cm_action(coefficient: parent.E, point: Point) -> Point:
    m, n = coefficient
    return point_add(point_mul(m, point), point_mul(n, omega_action(point)))


def explicit_f(x_value: int, y_value: int, t_value: int) -> Point:
    half = inv(2)
    x_target = (
        SCALE**2 * half * x_value
        * (t_value**2 - RHO**2) * inv(t_value)
    ) % Q
    y_target = (
        SCALE**3 * half * y_value
        * (t_value**2 + RHO**3) * inv(t_value)
    ) % Q
    require((y_target**2 - x_target**3 - 1) % Q == 0,
            "explicit hidden point missed E0")
    return x_target, y_target


def eisenstein_divides(divisor: parent.E, dividend: parent.E) -> bool:
    norm = parent.e_norm(divisor)
    numerator = parent.e_mul(dividend, parent.e_conjugate(divisor))
    return numerator[0] % norm == 0 and numerator[1] % norm == 0


def point_order(point: Point, bound: int) -> int:
    for order in range(1, bound + 1):
        if point_mul(order, point) is None:
            return order
    raise AssertionError("point order exceeded certified bound")


def verify_good_specialization() -> tuple[Point, Point, int]:
    require(sp.isprime(Q), "specialization modulus is not prime")
    require(Q not in (2, 3), "source or target has bad characteristic")

    require((ZETA12**4 - ZETA12**2 + 1) % Q == 0,
            "zeta_12 cyclotomic equation failed")
    require(pow(ZETA12, 12, Q) == 1, "zeta_12 lost order dividing twelve")
    require(all(pow(ZETA12, divisor, Q) != 1
                for divisor in (1, 2, 3, 4, 6)),
            "zeta_12 lost exact order twelve")
    require((OMEGA**2 + OMEGA + 1) % Q == 0 and OMEGA != 1,
            "target CM root changed")

    sqrt_three = (2 * ZETA12 - ZETA12**3) % Q
    require(sqrt_three**2 % Q == 3, "sqrt(3) embedding changed")
    quartic = RHO**4 - 2 * RHO**3 - 2 * RHO + 1
    quadratic = RHO**2 - (1 + 2 * ZETA12 - ZETA12**3) * RHO + 1
    require(quartic % Q == 0 and quadratic % Q == 0,
            "quartic root is incompatible with zeta_12")
    require((4 * RHO**3 - 6 * RHO**2 - 2) % Q != 0,
            "quartic root is inseparable")
    scale_denominator = (2 * RHO**3 + 3 * RHO**2 - 1) % Q
    require(scale_denominator != 0, "hidden-map scale denominator vanished")
    require(pow(SCALE, 6, Q) == 4 * inv(scale_denominator) % Q,
            "hidden-map sixth-root scale changed")
    require(6 * pow(SCALE, 5, Q) % Q != 0,
            "sixth-root scale is inseparable")

    half = inv(2)
    require(pow(SOURCE_X, 3, Q) == half,
            "chosen attachment x^3 is not 1/2")
    require(pow(SOURCE_Y, 2, Q) == sqrt_three * half % Q,
            "chosen attachment y^2 is not sqrt(3)/2")
    require((pow(SOURCE_X, 6, Q) + pow(SOURCE_Y, 4, Q)) % Q == 1,
            "chosen point missed C0")
    require(3 * SOURCE_X**2 % Q != 0 and 2 * SOURCE_Y % Q != 0,
            "attachment radicals became inseparable")

    ratio = Fraction(1, 3)
    require(pow(SOURCE_X, 6, Q) * inv(pow(SOURCE_Y, 4, Q)) % Q
            == ratio.numerator * inv(ratio.denominator) % Q,
            "marked ratio is not 1/3")
    t_value = (1 + SOURCE_Y**2) * inv(SOURCE_X**3) % Q
    require(t_value == 379 and t_value != 0,
            "attachment t-coordinate changed or became a pole")
    require((t_value**2 - 2 * sqrt_three * t_value - 1) % Q == 0,
            "R=1/3 t-equation changed")

    f_point = explicit_f(SOURCE_X, SOURCE_Y, t_value)
    tau_x = pow(ZETA12, 2, Q) * SOURCE_X % Q
    tau_y = pow(ZETA12, 3, Q) * SOURCE_Y % Q
    tau_t = -inv(t_value) % Q
    require(tau_t == (1 + tau_y**2) * inv(tau_x**3) % Q,
            "tau no longer sends t to -1/t")
    g_point = explicit_f(tau_x, tau_y, tau_t)

    require(f_point == (340, 181), "f(Q0) changed")
    require(g_point == (327, 3), "g(Q0) changed")
    require(point_mul(6, f_point) == (0, 396), "[6]f(Q0) changed")
    require(point_mul(9, f_point) == (35, 0), "[9]f(Q0) changed")
    require(point_mul(18, f_point) is None, "[18]f(Q0) is nonzero")
    require(point_order(f_point, 108) == 18, "integer order of f(Q0) changed")
    require(point_order(g_point, 108) == 6, "integer order of g(Q0) changed")
    require(g_point == point_mul(15, f_point), "g(Q0) is no longer [15]f(Q0)")

    return f_point, g_point, t_value


def verify_uniform_pi_obstruction() -> parent.E:
    pi = parent.e_sub(parent.OMEGA2, parent.ONE)
    require(parent.e_norm(pi) == 3, "pi lost norm three")
    zero_pairs = {
        (left, right)
        for left in range(3)
        for right in range(3)
        if (left * left + right * right) % 3 == 0
    }
    require(zero_pairs == {(0, 0)},
            "A^2+B^2=0 acquired a nonzero solution over F3")

    # If pi divides delta=A^2+omega*B^2, then omega=1 modulo pi and the
    # preceding residue check forces pi|A,B.  For every residual profile,
    # q(H)=12K would become q(H')=4K, but 3 does not divide any listed K,
    # whereas the hidden Gram has every degree divisible by 6.
    hidden_indices = {5, 7, 10, 11, 13}
    require(all(index % 3 != 0 and (4 * index) % 6 != 0
                for index in hidden_indices),
            "uniform hidden-degree contradiction changed")
    return pi


def verify_annihilator(f_point: Point) -> parent.E:
    annihilator = (6, 12)
    require(parent.e_norm(annihilator) == 108,
            "annihilator norm changed")
    require(cm_action(annihilator, f_point) is None,
            "claimed annihilator does not kill f(Q0)")

    # The 6-by-18 box gives 108 distinct O-multiples.  Since the principal
    # ideal (6+12omega) already has index/norm 108 and kills the point, this
    # proves equality Ann_O(f(Q0))=(6+12omega), not merely containment.
    orbit = {
        cm_action((m, n), f_point)
        for m in range(6)
        for n in range(18)
    }
    require(len(orbit) == 108, "O-orbit of f(Q0) changed")
    require(None in orbit, "O-orbit lost the origin")
    return annihilator


def verify_complete_ratio_lane(
    f_point: Point, g_point: Point, annihilator: parent.E, pi: parent.E
) -> tuple[str, Counter[int], Counter[tuple[int, int]]]:
    _, _, residual, _, _, _, _ = parent.enumerate_full_residual()
    lane = residual[42]
    require(len(lane) == 3168, "degree-42 residual vector lane changed")
    representatives, sizes = parent.symmetry_orbits(lane)
    require(len(representatives) == 132 and sizes == Counter({24: 132}),
            "degree-42 residual orbit quotient changed")

    determinant_norms: Counter[int] = Counter()
    profiles: Counter[tuple[int, int]] = Counter()
    rows: list[str] = []
    for index, vector in enumerate(representatives):
        _, b, c, d = vector
        hidden_a = parent.e_add(
            parent.e_scale(2, b), parent.e_mul(parent.OMEGA2, d)
        )
        hidden_b = parent.e_add(parent.e_scale(2, c), d)
        hidden_degree = parent.hidden_degree(hidden_a, hidden_b)
        require(hidden_degree % 12 == 0,
                "residual hidden degree lost divisibility by twelve")
        hidden_k = hidden_degree // 12
        profiles[(parent.e_norm(d), hidden_k)] += 1
        require(hidden_k % 3 != 0 and (4 * hidden_k) % 6 != 0,
                "uniform pi obstruction acquired a divisible-by-six quotient")
        require(eisenstein_divides(pi, d),
                "R=1/3 is not common to a residual d-kernel")

        determinant = parent.e_add(
            parent.e_mul(hidden_a, hidden_a),
            parent.e_mul(parent.OMEGA,
                         parent.e_mul(hidden_b, hidden_b)),
        )
        determinant_norms[parent.e_norm(determinant)] += 1
        require(parent.e_norm(determinant) % 3 != 0,
                "a residual cyclic determinant became divisible by three")

        reduced_coefficient = parent.e_add(
            hidden_a, parent.e_scale(15, hidden_b)
        )
        require(not eisenstein_divides(annihilator, reduced_coefficient),
                "a residual hidden projection entered Ann(f(Q0))")

        hidden_value = point_add(
            cm_action(hidden_a, f_point),
            cm_action(hidden_b, g_point),
        )
        require(hidden_value == cm_action(reduced_coefficient, f_point),
                "g=[15]f reduction failed on a residual vector")
        require(hidden_value is not None,
                "a residual hidden projection vanished at Q0")
        rows.append(
            f"orbit{index:02d} A={hidden_a} B={hidden_b} "
            f"Aplus15B={reduced_coefficient} H397={hidden_value}"
        )

    require(profiles == Counter({
        (3, 13): 28, (9, 11): 24, (12, 10): 36,
        (21, 7): 32, (27, 5): 12,
    }), "degree-42 ratio-one-third profile changed")
    require(determinant_norms == Counter({
        4: 2, 52: 12, 100: 4, 148: 8, 196: 8, 208: 12,
        244: 8, 292: 12, 388: 4, 400: 12, 436: 4, 592: 12,
        628: 4, 676: 14, 724: 4, 772: 4, 916: 4, 964: 4,
    }), "cyclic determinant hostile distribution changed")

    # For fixed sixth/fourth powers the 24 radical choices split into two
    # C12 orbits.  rho:(x,y)->(-x,-y), exponent shift (3,2), exchanges them.
    # Precomposition by rho permutes the full Hom shell, so the tested node
    # orbit is complete for all 132 map-orbit representatives.
    root_choices = {(left, right) for left in range(6) for right in range(4)}
    canonical = {(j % 6, j % 4) for j in range(12)}
    rho_image = {((left + 3) % 6, (right + 2) % 4)
                 for left, right in canonical}
    require(len(canonical) == len(rho_image) == 12,
            "attachment root orbit changed size")
    require(canonical.isdisjoint(rho_image)
            and canonical | rho_image == root_choices,
            "rho no longer exchanges the two root orbits")

    row_digest = hashlib.sha256("\n".join(rows).encode("ascii")).hexdigest()
    return row_digest, determinant_norms, profiles


def main() -> None:
    f_point, g_point, t_value = verify_good_specialization()
    pi = verify_uniform_pi_obstruction()
    annihilator = verify_annihilator(f_point)
    row_digest, determinant_norms, profiles = verify_complete_ratio_lane(
        f_point, g_point, annihilator, pi
    )

    print("THM-4249 R=1/3 hidden-projection exclusion exact certificate")
    print("field=F_397 good_characteristic=PASS")
    print(
        "embedding=zeta12:157 quartic_root:161 scale:27 "
        "source_x:15 source_y:28 t:379"
    )
    print("relations=cyclotomic+quartic_pair+sixth_root+C0+R_one_third PASS")
    print(f"fQ={f_point} order18 gQ={g_point} order6 relation_g=[15]f")
    print(
        f"uniform_pi={pi} norm3 delta_mod_pi=A2+B2 "
        "zero_only_at_A_B_zero all_K_nonzero_mod3"
    )
    print("order_witnesses=[6]f:(0,396) [9]f:(35,0) [18]f:O")
    print(f"Ann_O(fQ)=({annihilator[0]}+{annihilator[1]}omega) norm108 orbit108")
    print("degree42_common_ratio_lane=vectors3168 symmetry_orbits132 all_size24")
    print(f"Nd_K_orbit_profiles={dict(sorted(profiles.items()))}")
    print(f"cyclic_determinant_norms={dict(sorted(determinant_norms.items()))}")
    print(f"orbit_row_sha256={row_digest}")
    print("annihilator_divisibility=0_of_132 direct_hidden_zero=0_of_132")
    print("root_choices=24 two_C12_orbits rho_exchanges one_canonical_orbit_complete")
    print("consequence=R_42_common_ratio_1/3_EXCLUDED S42_envelope_at_most34")
    print("incidence_frontier=degree34:864 degree42:648 total:1512")
    print("scope=remaining_ratios_mixed_tests_W0_M12_entry_JC2_OPEN")
    print("result=PASS")


if __name__ == "__main__":
    main()
