#!/usr/bin/env python3
"""Independent recurrence referee for THM-4378.

Unlike the primary coefficient-matrix audit, this script computes in the
quadratic free module Z[q] + z Z[q], using only z^2=z-q.  It reconstructs
band packets by two multiplication recurrences and checks the quotient,
reciprocal action, eigensectors, parity gluing, norm/difference images, the
ell=2 exception, and the unsigned-reflection hostile.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("independent check failed: " + label)


def clean(poly: tuple[int, ...]) -> tuple[int, ...]:
    work = list(poly)
    while work and work[-1] == 0:
        work.pop()
    return tuple(work)


def plus(left: tuple[int, ...], right: tuple[int, ...], sign: int = 1) -> tuple[int, ...]:
    size = max(len(left), len(right))
    return clean(tuple(
        (left[k] if k < len(left) else 0)
        + sign * (right[k] if k < len(right) else 0)
        for k in range(size)
    ))


def q_times(poly: tuple[int, ...]) -> tuple[int, ...]:
    return (0,) + poly if poly else ()


def scalar(value: int, poly: tuple[int, ...]) -> tuple[int, ...]:
    return clean(tuple(value * coefficient for coefficient in poly))


Pair = tuple[tuple[int, ...], tuple[int, ...]]


def pair_plus(left: Pair, right: Pair, sign: int = 1) -> Pair:
    return plus(left[0], right[0], sign), plus(left[1], right[1], sign)


def times_z(pair: Pair) -> Pair:
    """z(A+zB)=-qB+z(A+B)."""
    a, b = pair
    return scalar(-1, q_times(b)), plus(a, b)


def times_one_minus_z(pair: Pair) -> Pair:
    """(1-z)(A+zB)=A+qB-zA."""
    a, b = pair
    return plus(a, q_times(b)), scalar(-1, a)


@lru_cache(maxsize=None)
def packet_pair(i: int, j: int) -> Pair:
    check(i >= 0 and j >= 0, "packet quadrant")
    pair: Pair = ((1,), ())
    for _ in range(i):
        pair = times_z(pair)
    for _ in range(j):
        pair = times_one_minus_z(pair)
    return pair


def reciprocal(pair: Pair) -> Pair:
    """A(q)+zB(q) -> (A+B)(q)-zB(q)."""
    a, b = pair
    return plus(a, b), scalar(-1, b)


def params(ell: int) -> tuple[int, int]:
    return (ell + 1) // 2, (ell + 2) // 3


def direct_source(ell: int, u: int, v: int) -> bool:
    s, rho = params(ell)
    n0 = s + u - 1
    return u >= 1 and v >= rho and n0 >= v and n0 + v >= ell


def band(ell: int, u: int, v: int) -> bool:
    s, rho = params(ell)
    return u >= rho and v >= rho and abs(u - v) <= s - 1


def finite_band(ell: int, depth: int) -> tuple[tuple[int, int], ...]:
    s, _ = params(ell)
    return tuple(
        (i, total - i)
        for total in range(depth + 1)
        for i in range(total + 1)
        if abs(2 * i - total) <= s - 1
    )


def fibre_size(ell: int, u: int, v: int) -> int:
    s, _ = params(ell)
    return sum(
        1
        for e in range(v + 1)
        if min(2 * v + e - ell, s + u - v - 1 - e, v - e, e) >= 0
    )


def unit(power: int) -> tuple[int, ...]:
    return (0,) * power + (1,)


def quotient_checks() -> tuple[int, int]:
    blocks = 0
    packet_count = 0
    for ell in range(2, 61):
        s, rho = params(ell)
        for u in range(1, 4 * s + 12):
            for v in range(1, 4 * s + 12):
                check(
                    (direct_source(ell, u, v) and direct_source(ell, v, u))
                    == band(ell, u, v),
                    "bilateral band from two source tests",
                )

        for depth in range(25):
            points = finite_band(ell, depth)
            if ell == 2:
                check(points == tuple((k, k) for k in range(depth // 2 + 1)),
                      "ell two diagonal only")
                for k, point in enumerate(points):
                    check(packet_pair(*point) == (unit(k), ()), "ell two q power")
                blocks += 1
                packet_count += len(points)
                continue

            # The central spine is literally the standard pair-module basis.
            for k in range(depth // 2 + 1):
                check(packet_pair(k, k) == (unit(k), ()), "diagonal q spine")
                if 2 * k + 1 <= depth:
                    check(packet_pair(k + 1, k) == ((), unit(k)), "oriented zq spine")

            # Every band packet has legal truncated A and B coordinates.
            for i, j in points:
                a, b = packet_pair(i, j)
                check(len(a) <= depth // 2 + 1, "invariant coordinate bound")
                check(len(b) <= (depth + 1) // 2, "oriented coordinate bound")
                check(reciprocal(packet_pair(i, j)) == packet_pair(j, i),
                      "plain residual polynomial reciprocity")
                packet_count += 1

            # Direct recurrence verification of every internal Pascal cell.
            for i, j in points:
                if i + j == depth or abs(i - j) > s - 2:
                    continue
                cell = pair_plus(
                    packet_pair(i, j),
                    pair_plus(packet_pair(i + 1, j), packet_pair(i, j + 1)),
                    -1,
                )
                check(cell == ((), ()), "quadratic-module Pascal cell")
            blocks += 1
    return blocks, packet_count


def eigen_and_smith_checks() -> tuple[int, int]:
    blocks = 0
    coefficient_boxes = 0
    for depth in range(33):
        na = depth // 2 + 1
        nb = (depth + 1) // 2

        # Exhaust small boxes in one q-coordinate at a time. This validates
        # the equations Jf=f iff B=0 and Jf=-f iff B=-2A without importing
        # a linear-algebra package.
        for k in range(na):
            for a_value, b_value in product(range(-3, 4), repeat=2):
                if k >= nb and b_value:
                    continue
                a = scalar(a_value, unit(k))
                b = scalar(b_value, unit(k)) if k < nb else ()
                pair = (a, b)
                image = reciprocal(pair)
                is_plus = image == pair
                is_minus = image == (scalar(-1, a), scalar(-1, b))
                check(is_plus == (not b), "exact plus equation")
                check(is_minus == (b == scalar(-2, a)), "exact minus equation")
                coefficient_boxes += 1

        # Parity of every B coefficient is the complete eigen-sum quotient.
        for k in range(nb):
            qk: Pair = (unit(k), ())
            zqk: Pair = ((), unit(k))
            xqk: Pair = (scalar(-1, unit(k)), scalar(2, unit(k)))
            check(reciprocal(qk) == qk, "q invariant")
            check(reciprocal(xqk) == (scalar(-1, xqk[0]), scalar(-1, xqk[1])),
                  "xq anti-invariant")
            check(pair_plus(qk, xqk) == ((), scalar(2, unit(k))),
                  "two times orientation eigensplits")

            difference = pair_plus(zqk, reciprocal(zqk), -1)
            norm = pair_plus(zqk, reciprocal(zqk))
            check(difference == xqk, "difference reaches anti generator")
            check(norm == qk, "norm reaches paired invariant")

        if depth % 2 == 0:
            top: Pair = (unit(depth // 2), ())
            check(pair_plus(top, reciprocal(top))
                  == (scalar(2, top[0]), ()), "unpaired top norm doubles")
        blocks += 1
    return blocks, coefficient_boxes


def hostile_checks() -> None:
    # Use the signed F convention at ell=3. In pair coordinates:
    # F11=1, F12=1-z, F21=-z.
    f11: Pair = ((1,), ())
    f12: Pair = times_one_minus_z(f11)
    f21: Pair = scalar(-1, times_z(f11)[0]), scalar(-1, times_z(f11)[1])
    circuit = pair_plus(f11, pair_plus(f12, f21, -1), -1)
    reflected = pair_plus(f11, pair_plus(f21, f12, -1), -1)
    check(circuit == ((), ()), "smallest circuit zero")
    check(reflected == ((2,), ()), "bare swap has trace two")


def fibre_bridge_checks() -> int:
    cases = 0
    for ell in range(3, 51):
        s, rho = params(ell)
        for d in range(1, s):
            # This obtains B_d from a wholly separate recurrence rather than
            # from the primary's monomial coefficient conversion.
            _, b_d = packet_pair(d, 0)
            check(b_d and b_d[0] == 1, "B_d constant term")
            for w in range(rho, s + d + 3):
                k = w - rho
                plus_pair = packet_pair(k + d, k)
                minus_pair = packet_pair(k, k + d)
                check(reciprocal(plus_pair) == minus_pair,
                      "bridge reciprocal packet pair")
                plus_gamma = tuple(value % 2 for value in plus_pair[1])
                minus_gamma = tuple(value % 2 for value in minus_pair[1])
                check(plus_gamma == minus_gamma and any(plus_gamma),
                      "same nonzero reciprocal gamma class")

                mu_plus = fibre_size(ell, w + d, w)
                mu_minus = fibre_size(ell, w, w + d)
                aggregate = tuple(
                    (mu_plus * a + mu_minus * b) % 2
                    for a, b in zip(plus_gamma, minus_gamma)
                )
                check(any(aggregate) == ((mu_plus - mu_minus) % 2 == 1),
                      "aggregate gamma iff odd imbalance")
                if w >= s:
                    beta = min(w - s + d + 1, 2 * d)
                    check(mu_plus - mu_minus == beta, "jet rank equals imbalance")
                    check(any(aggregate) == (beta % 2 == 1), "jet parity toggle")
                    if w >= s + d - 1:
                        check(beta == 2 * d and not any(aggregate),
                              "full jet gamma vanishes")
                cases += 1
    return cases


def main() -> None:
    blocks, packets = quotient_checks()
    eigen_blocks, boxes = eigen_and_smith_checks()
    hostile_checks()
    bridges = fibre_bridge_checks()
    print("THM-4378 independent quadratic-module referee")
    print("design=Z[q] + z Z[q] recurrence with z^2=z-q; no primary import")
    print("scope=finite bilateral packet quotient; JC(2), DC(2), bracket and seam OPEN")
    print()
    print(f"band quotient blocks={blocks}, reconstructed packets={packets}")
    print(f"eigen/Smith depths={eigen_blocks}, local coefficient boxes={boxes}")
    print(f"source-fibre parity bridge cases={bridges}")
    print("verified quotient spine=1,z,q,zq,... and reciprocal (A,B)->(A+B,-B)")
    print("verified eigen sum parity cokernel=(Z/2)^ceil(D/2), no 4-torsion")
    print("verified difference surjects Q_minus; norm has one even-top Z/2")
    print("verified weighted gamma nonzero iff odd imbalance; full 2d jet kills it")
    print("verified ell=2 diagonal J=id boundary")
    print("smallest hostile ell=3,R=2: trace(C)=0, trace(bare-swap C)=2")
    print()
    print(f"PASS checks={CHECKS}")


if __name__ == "__main__":
    main()
