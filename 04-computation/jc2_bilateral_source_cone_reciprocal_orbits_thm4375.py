"""Exact finite audit for the reciprocal source-cone theorem candidate.

This is a scratch verifier.  It imports no repository code.  The checked
identities are symbolic formulas whose finite tests include all residue and
boundary regimes used in the accompanying proof note.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
import sys


CHECKS = 0


def check(condition: bool, message: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def params(ell: int) -> tuple[int, int]:
    return (ell + 1) // 2, (ell + 2) // 3


def triangular(k: int) -> int:
    return k * (k + 1) // 2


def address(u: int, v: int) -> int:
    tau = u + v - 1
    return triangular(tau - 1) + u


def source_valid(ell: int, u: int, v: int) -> bool:
    s, rho = params(ell)
    n0 = s + u - 1
    n = v
    return u >= 1 and v >= 1 and n >= rho and n0 >= n and n0 + n >= ell


def bilateral(ell: int, u: int, v: int) -> bool:
    return source_valid(ell, u, v) and source_valid(ell, v, u)


def bilateral_band(ell: int, u: int, v: int) -> bool:
    s, rho = params(ell)
    return u >= rho and v >= rho and abs(u - v) <= s - 1


def mu_formula(ell: int, u: int, v: int) -> int:
    s, _ = params(ell)
    n = v
    n0 = s + u - 1
    return min(n, n0 - n) - max(0, ell - 2 * n) + 1


def source_fibre(ell: int, u: int, v: int) -> tuple[tuple[int, int, int, int], ...]:
    s, _ = params(ell)
    n = v
    n0 = s + u - 1
    ans = []
    for e in range(n + 1):
        c = n - e
        a = 2 * n + e - ell
        b = n0 - n - e
        if min(a, b, c, e) >= 0:
            ans.append((a, b, c, e))
    return tuple(ans)


def block_points(ell: int, tau: int) -> tuple[tuple[int, int], ...]:
    total = tau + 1
    return tuple((u, total - u) for u in range(1, total) if bilateral(ell, u, total - u))


def block_interval(ell: int, tau: int) -> tuple[int, int] | None:
    s, rho = params(ell)
    a = s - 1
    total = tau + 1
    low = max(rho, (total - a + 1) // 2)
    high = min(total - rho, (total + a) // 2)
    return None if low > high else (low, high)


def orbit_count_coefficient(ell: int, tau: int) -> int:
    s, rho = params(ell)
    k = tau - (2 * rho - 1)
    if k < 0:
        return 0
    return sum(1 for d in range(s) if d <= k and (d - k) % 2 == 0)


def oriented_count_coefficient(ell: int, tau: int) -> int:
    s, rho = params(ell)
    k = tau - (2 * rho - 1)
    if k < 0:
        return 0
    return sum((1 if d == 0 else 2) for d in range(s) if d <= k and (d - k) % 2 == 0)


def primitive_scale_interval(ell: int, p: int, q: int) -> tuple[int, int | None]:
    s, rho = params(ell)
    low = (rho + min(p, q) - 1) // min(p, q)
    if p == q:
        return low, None
    return low, (s - 1) // abs(p - q)


def run() -> None:
    # The one elementary inequality which makes the two-sided cone a band.
    for ell in range(2, 301):
        s, rho = params(ell)
        check(2 * rho >= ell - s + 1, f"redundant sum wall, ell={ell}")

    # Directly compare THM-4368 feasibility on both ends with the claimed band.
    for ell in range(2, 91):
        s, rho = params(ell)
        for u in range(1, 5 * s + 18):
            for v in range(1, 5 * s + 18):
                check(
                    bilateral(ell, u, v) == bilateral_band(ell, u, v),
                    f"band mismatch ell={ell}, pair={(u, v)}",
                )

    # Every orbit has the unique (w+d,w), w=rho+j, 0<=d<s representative.
    for ell in range(2, 91):
        s, rho = params(ell)
        seen = set()
        for j in range(80):
            for d in range(s):
                u, v = rho + j + d, rho + j
                check(bilateral(ell, u, v), f"bad orbit parameter ell={ell}, {(j, d)}")
                key = (min(u, v), max(u, v))
                check(key not in seen, f"duplicate orbit parameter ell={ell}, {key}")
                seen.add(key)
                check((min(u, v) - rho, abs(u - v)) == (j, d))

    # Triangular-block address intervals, reflection, orbit counts, and the
    # coefficient interpretation of the two rational generating functions.
    for ell in range(2, 91):
        s, rho = params(ell)
        tau0 = 2 * rho - 1
        for tau in range(max(1, tau0 - 4), tau0 + 5 * s + 45):
            pts = block_points(ell, tau)
            interval = block_interval(ell, tau)
            if interval is None:
                check(not pts)
            else:
                low, high = interval
                check(pts == tuple((u, tau + 1 - u) for u in range(low, high + 1)))
                check(low + high == tau + 1)
                addrs = tuple(address(u, v) for u, v in pts)
                check(addrs == tuple(range(triangular(tau - 1) + low,
                                               triangular(tau - 1) + high + 1)))
                for u, v in pts:
                    check(address(u, v) + address(v, u) == tau * tau + 1)
                    d = abs(u - v)
                    check(
                        {address(u, v), address(v, u)}
                        == {(tau * tau + 1 - d) // 2, (tau * tau + 1 + d) // 2}
                    )

            orbit_keys = {(min(u, v), max(u, v)) for u, v in pts}
            check(len(pts) == oriented_count_coefficient(ell, tau))
            check(len(orbit_keys) == orbit_count_coefficient(ell, tau))
            fixed = [(u, v) for u, v in pts if u == v]
            expected_fixed = int(tau >= tau0 and (tau - tau0) % 2 == 0)
            check(len(fixed) == expected_fixed)

    # Finite initial-segment path-power-with-loops count.
    for ell in range(2, 91):
        s, rho = params(ell)
        for n in range(1, 80):
            vertices = range(rho, rho + n)
            orbits = {
                (min(u, v), max(u, v))
                for u in vertices
                for v in vertices
                if bilateral(ell, u, v)
            }
            dmax = min(s - 1, n - 1)
            predicted = (dmax + 1) * n - dmax * (dmax + 1) // 2
            check(len(orbits) == predicted)

    # Primitive Stern--Brocot rays: exact scale interval.  The diagonal ray
    # is truncated only for the purpose of this finite comparison.
    for ell in range(2, 91):
        for p in range(1, 50):
            for q in range(1, 50):
                if gcd(p, q) != 1:
                    continue
                low, high = primitive_scale_interval(ell, p, q)
                brute = [g for g in range(1, 121) if bilateral(ell, g * p, g * q)]
                if high is None:
                    predicted = list(range(low, 121))
                elif low <= high:
                    predicted = list(range(low, min(high, 120) + 1))
                else:
                    predicted = []
                check(brute == predicted, f"ray mismatch ell={ell}, ray={(p, q)}")

    # Source fibres, endpoint formulas, defect decomposition, conservation
    # ceiling, equality boundary, and sharp maximum asymmetry.
    for ell in range(2, 91):
        s, rho = params(ell)
        for w in range(rho, 5 * s + 30):
            for d in range(s):
                up = mu_formula(ell, w + d, w)
                down = mu_formula(ell, w, w + d)
                up_fibre = source_fibre(ell, w + d, w)
                down_fibre = source_fibre(ell, w, w + d)
                check(up == len(up_fibre) and down == len(down_fibre))
                check(up >= 1 and down >= 1)
                up_defect = (
                    s + d
                    - max(0, s - 1 + d - w)
                    - max(0, ell - 2 * w)
                )
                down_defect = (
                    s - d
                    - max(0, s - 1 - w - 2 * d)
                    - max(0, ell - 2 * w - 2 * d)
                )
                check(up == up_defect)
                check(down == down_defect)
                check(up <= s + d and down <= s - d and up + down <= 2 * s)
                equality = w >= max(s, s - 1 + d)
                check((up + down == 2 * s) == equality)
                if equality:
                    check((up, down) == (s + d, s - d))
                check(abs(up - down) <= 2 * s - 2)
                if s >= 2:
                    check(
                        (abs(up - down) == 2 * s - 2)
                        == (d == s - 1 and w >= 2 * s - 2)
                    )
                else:
                    check(up == down == 1)

        if s == 1:
            check(mu_formula(ell, rho, rho) == 1)
        else:
            w, d = 2 * s - 2, s - 1
            check((mu_formula(ell, w + d, w), mu_formula(ell, w, w + d))
                  == (2 * s - 1, 1))

    # Named ell=10 sharp controls.  This is the earliest bilateral block and
    # the earliest nontrivial reciprocal block, respectively.
    ell = 10
    check(block_points(ell, 7) == ((4, 4),))
    check(address(4, 4) == 25)
    check(source_fibre(ell, 4, 4) == ((0, 2, 2, 2), (1, 1, 1, 3), (2, 0, 0, 4)))
    check(block_points(ell, 8) == ((4, 5), (5, 4)))
    check((address(4, 5), address(5, 4)) == (32, 33))
    check(source_fibre(ell, 5, 4) == ((0, 3, 2, 2), (1, 2, 1, 3), (2, 1, 0, 4)))
    check(source_fibre(ell, 4, 5) == ((0, 3, 5, 0), (1, 2, 4, 1),
                                              (2, 1, 3, 2), (3, 0, 2, 3)))
    check((mu_formula(ell, 5, 4), mu_formula(ell, 4, 5)) == (3, 4))

    ray45 = []
    for g in range(1, 6):
        u, v = 4 * g, 5 * g
        ray45.append(
            (
                g,
                bilateral(ell, u, v),
                address(u, v),
                address(v, u),
                mu_formula(ell, u, v) if bilateral(ell, u, v) else None,
                mu_formula(ell, v, u) if bilateral(ell, u, v) else None,
            )
        )
    check(
        ray45
        == [
            (1, True, 32, 33, 4, 3),
            (2, True, 144, 146, 3, 7),
            (3, True, 337, 340, 2, 8),
            (4, True, 611, 615, 1, 9),
            (5, False, 966, 971, None, None),
        ]
    )

    # Rational identity is only a coordinate check; it asserts no transfer to
    # Stern--Brocot dynamics or to any open problem.
    check(Fraction(5, 4) == 1 / Fraction(4, 5))

    # Exact low-order edge cases.
    check(params(2) == (1, 1))
    for w in range(1, 80):
        check(block_points(2, 2 * w - 1) == ((w, w),))
        check(not block_points(2, 2 * w))
        check(address(w, w) == w * w + (w - 1) * (w - 1))
        check(mu_formula(2, w, w) == 1)
    check(params(3) == (2, 1))
    for tau in range(1, 100):
        check(orbit_count_coefficient(3, tau) == 1)
        check(oriented_count_coefficient(3, tau) == (1 if tau % 2 else 2))
    check((mu_formula(3, 2, 1), mu_formula(3, 1, 2)) == (1, 1))
    check((mu_formula(3, 3, 2), mu_formula(3, 2, 3)) == (3, 1))

    # Eventual two-block coefficient sums prove the asymptotic constants in
    # the note, independently of any floating-point limit calculation.
    for ell in range(2, 91):
        s, rho = params(ell)
        tau0 = 2 * rho - 1
        for tau in range(tau0 + 2 * s + 10, tau0 + 2 * s + 80, 2):
            check(
                orbit_count_coefficient(ell, tau)
                + orbit_count_coefficient(ell, tau + 1)
                == s
            )
            check(
                oriented_count_coefficient(ell, tau)
                + oriented_count_coefficient(ell, tau + 1)
                == 2 * s - 1
            )

    lines = [
        "RECIPROCAL SOURCE-CONE AUDIT: PASS",
        f"assertions={CHECKS}",
        "ell_range=2..90; pair boxes through 5*s+17",
        "orbit_parameter_j=0..79; block window=tau0-4..tau0+5*s+44",
        "primitive_coprime_rays=1..49 squared; scale_window=1..120",
        "fibre_window=w=rho..5*s+29, d=0..s-1",
        "edge_ell2=loop-only, fixed addresses w^2+(w-1)^2, mu=1",
        "edge_ell3=one orbit per block; oriented counts alternate 1,2; sharp fibre starts (3,2)",
        "ell10_fixed_control=(u,v,Addr,mu)=(4,4,25,3)",
        "ell10_first_nontrivial_orbit=(4,5)<->(5,4), Addr=32<->33, mu=4<->3",
        "ell10_ray_4_over_5_scales=1..4 with mu ladder (4,3),(3,7),(2,8),(1,9)",
    ]
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("ascii"))


if __name__ == "__main__":
    run()
