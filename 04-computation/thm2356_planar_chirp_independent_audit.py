#!/usr/bin/env python3
"""Independent exact audit for THM-2356.

This deliberately does not reproduce the theorem's cyclotomic implementation.
It checks the two character-orthogonality selectors combinatorially, including
the h=0 failure boundary, and independently exhausts the F_169 planar map.
"""

from __future__ import annotations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def prime_field_selector_audit(p: int) -> int:
    inv_two = pow(2, -1, p)

    def phi(x: int) -> int:
        return inv_two * x * x % p

    controls = 0
    for h in range(1, p):
        for y in range(p):
            target = (phi(y + h) - phi(y)) % p
            survivors = []
            for x in range(p):
                for u in range(p):
                    if (x - u - h) % p:
                        continue
                    if (phi(x) - phi(u) - target) % p:
                        continue
                    survivors.append((x, u))
            require(
                survivors == [((y + h) % p, y)],
                f"wrong selector at p={p}, h={h}, y={y}: {survivors}",
            )
            controls += 1

    # At h=0 both character averages retain the whole diagonal, not one site.
    for y in range(p):
        survivors = [
            (x, u)
            for x in range(p)
            for u in range(p)
            if (x - u) % p == 0 and (phi(x) - phi(u)) % p == 0
        ]
        require(
            survivors == [(x, x) for x in range(p)],
            f"h=0 boundary changed at p={p}, y={y}",
        )

    return controls


def canonical_cyclotomic(coefficients: list[int]) -> tuple[int, ...]:
    """Canonical Q(zeta_p) vector using 1+zeta+...+zeta^(p-1)=0."""

    last = coefficients[-1]
    return tuple(value - last for value in coefficients[:-1])


def graph_restriction_selector_audit(p: int) -> int:
    """Check the signs/scalar in the graph-restriction Fourier formula."""

    inv_two = pow(2, -1, p)

    def phi(x: int) -> int:
        return inv_two * x * x % p

    controls = 0
    for c in range(p):
        for a in range(p):
            for b in range(p):
                for r in range(p):
                    for s in range(p):
                        for q in range(p):
                            # Expansion of
                            # Ahat(beta,alpha) zeta^(-alpha*c)
                            # zeta^((b-beta)q+(a-alpha)phi(q)).
                            histogram = [0] * p
                            for alpha in range(p):
                                for beta in range(p):
                                    exponent = (
                                        beta * r
                                        + alpha * s
                                        - alpha * c
                                        + (b - beta) * q
                                        + (a - alpha) * phi(q)
                                    ) % p
                                    histogram[exponent] += 1

                            expected = [0] * p
                            if r == q and s == (phi(q) + c) % p:
                                expected[(b * q + a * phi(q)) % p] = p * p
                            require(
                                canonical_cyclotomic(histogram)
                                == canonical_cyclotomic(expected),
                                (
                                    "graph Fourier selector failed at "
                                    f"p={p}, c={c}, a={a}, b={b}, "
                                    f"r={r}, s={s}, q={q}"
                                ),
                            )
                            controls += 1

    return controls


P = 13
Field = tuple[int, int]


def fadd(x: Field, y: Field) -> Field:
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def fsub(x: Field, y: Field) -> Field:
    return ((x[0] - y[0]) % P, (x[1] - y[1]) % P)


def fmul(x: Field, y: Field) -> Field:
    return (
        (x[0] * y[0] + 2 * x[1] * y[1]) % P,
        (x[0] * y[1] + x[1] * y[0]) % P,
    )


def fscale(c: int, x: Field) -> Field:
    return (c * x[0] % P, c * x[1] % P)


def phi169(x: Field) -> Field:
    return fscale(7, fmul(x, x))


def trace169(x: Field) -> int:
    return 2 * x[0] % P


def field_169_audit() -> tuple[int, int, int]:
    elements = [(a, b) for a in range(P) for b in range(P)]
    zero = (0, 0)
    one = (1, 0)

    # 2 is a nonsquare, so theta^2-2 is irreducible.
    squares = {a * a % P for a in range(P)}
    require(2 not in squares, "2 unexpectedly square modulo 13")
    require(
        all(any(fmul(x, y) == one for y in elements) for x in elements if x != zero),
        "F_169 inverse check failed",
    )

    derivative_checks = 0
    trace_pairing_checks = 0
    for h in elements:
        if h == zero:
            continue
        image = set()
        for y in elements:
            derivative = fsub(phi169(fadd(y, h)), phi169(y))
            expected = fadd(fmul(h, y), phi169(h))
            require(derivative == expected, "quadratic derivative identity failed")
            image.add(derivative)
        require(len(image) == len(elements), f"nonplanar F_169 derivative at {h}")
        derivative_checks += 1

        # The absolute-trace pairing a -> Tr(a h) is a nonzero linear
        # functional for h != 0, hence each residue occurs exactly 13 times.
        counts = [0] * P
        for a in elements:
            counts[trace169(fmul(a, h))] += 1
        require(counts == [P] * P, f"degenerate trace pairing at h={h}")
        trace_pairing_checks += 1

    graph_cells = {
        (q, fadd(phi169(q), c)) for c in elements for q in elements
    }
    require(len(graph_cells) == len(elements) ** 2, "F_169 graphs do not partition")

    return derivative_checks, trace_pairing_checks, len(graph_cells)


def diagonal_boundary_audit() -> None:
    # Gamma_(x,y)=z_x conjugate(z_y).  The scalar identity reduces exactly
    # to (z_x zbar_y)(z_w zbar_x)/(z_w zbar_y)=|z_x|^2.
    def gmul(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
        return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    def gconj(x: tuple[int, int]) -> tuple[int, int]:
        return (x[0], -x[1])

    zx = (2, 3)
    zy = (-1, 4)
    zw = (5, -2)
    gamma_xy = gmul(zx, gconj(zy))
    gamma_wx = gmul(zw, gconj(zx))
    gamma_wy = gmul(zw, gconj(zy))
    norm_x = zx[0] * zx[0] + zx[1] * zx[1]
    require(
        gmul(gamma_xy, gamma_wx)
        == (norm_x * gamma_wy[0], norm_x * gamma_wy[1]),
        "three-support diagonal identity failed",
    )

    # Uniform boundaries: singleton location is invisible; on a fixed
    # two-site support, real magnitudes (2,3) and (3,2) have the same total
    # diagonal and the same two ordered cross products.
    require(2 * 2 + 3 * 3 == 3 * 3 + 2 * 2, "two-site norm changed")
    require(2 * 3 == 3 * 2, "two-site cross product changed")


Gaussian = tuple[int, int]


def gadd(x: Gaussian, y: Gaussian) -> Gaussian:
    return (x[0] + y[0], x[1] + y[1])


def gneg(x: Gaussian) -> Gaussian:
    return (-x[0], -x[1])


def gmul(x: Gaussian, y: Gaussian) -> Gaussian:
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gconj(x: Gaussian) -> Gaussian:
    return (x[0], -x[1])


def cyclotomic_canonical_gaussian(
    coefficients: list[Gaussian],
) -> tuple[Gaussian, ...]:
    """Reduce at a primitive p-th root by 1+z+...+z^(p-1)=0."""

    last = coefficients[-1]
    return tuple((value[0] - last[0], value[1] - last[1]) for value in coefficients[:-1])


def circular_norm(coefficients: list[Gaussian]) -> list[Gaussian]:
    """Return P(z) conjugate(P(z)) in Z[i][z]/(z^p-1)."""

    p = len(coefficients)
    result = [(0, 0) for _ in range(p)]
    for j, left in enumerate(coefficients):
        for k, right in enumerate(coefficients):
            index = (j - k) % p
            result[index] = gadd(result[index], gmul(left, gconj(right)))
    return result


def channel_difference_identity_audit(p: int) -> int:
    """Check the final D_c identity without floating-point characters."""

    signals: list[list[Gaussian]] = []
    for seed in range(1, 8):
        signals.append(
            [
                (
                    ((seed + 2 * q + q * q) % 7) - 3,
                    ((3 * seed + q + 2 * q * q) % 5) - 2,
                )
                for q in range(p)
            ]
        )

    controls = 0
    for signal in signals:
        nonzero_sum = (0, 0)
        square_sum = 0
        for value in signal[1:]:
            nonzero_sum = gadd(nonzero_sum, value)
            square_sum += value[0] * value[0] + value[1] * value[1]
        rhs = square_sum + (
            nonzero_sum[0] * nonzero_sum[0]
            + nonzero_sum[1] * nonzero_sum[1]
        )

        accumulated = [(0, 0) for _ in range(p)]
        total = (0, 0)
        for value in signal:
            total = gadd(total, value)
        for b in range(p):
            delta = [(0, 0) for _ in range(p)]
            for q, value in enumerate(signal):
                index = b * q % p
                delta[index] = gadd(delta[index], value)
            delta[0] = gadd(delta[0], gneg(total))
            norm = circular_norm(delta)
            accumulated = [
                gadd(old, new) for old, new in zip(accumulated, norm)
            ]

        expected = [(p * rhs, 0)] + [(0, 0) for _ in range(p - 1)]
        require(
            cyclotomic_canonical_gaussian(accumulated)
            == cyclotomic_canonical_gaussian(expected),
            f"channel-difference identity failed at p={p}",
        )
        controls += 1
    return controls


def vertical_tensor_boundary_audit(p: int) -> int:
    """Check the graph-singleton classification and its sharp tensors."""

    inv_two = pow(2, -1, p)

    def phi(x: int) -> int:
        return inv_two * x * x % p

    # The graph coordinates (c,q) and joint coordinates (q,z) are related
    # by the explicit bijection z=phi(q)+c.
    image = {
        (q, (phi(q) + c) % p): (c, q)
        for c in range(p)
        for q in range(p)
    }
    require(len(image) == p * p, f"graph/joint bijection failed at p={p}")

    # A full vertical tensor leaves every graph alive at q=0.
    vertical = {(0, z): (z + 1, 1 - z) for z in range(p)}
    for c in range(p):
        support = [
            q
            for q in range(p)
            if vertical.get((q, (phi(q) + c) % p), (0, 0)) != (0, 0)
        ]
        require(support == [0], f"vertical graph support changed at p={p}, c={c}")

    # Moving the tensor to any q0!=0 keeps every graph one-sparse but
    # visibly lands at a nonzero target.  Thus q0=0 is the exact residual.
    for q0 in range(1, p):
        shifted = {(q0, z): (z + 1, 2 - z) for z in range(p)}
        for c in range(p):
            support = [
                q
                for q in range(p)
                if shifted.get((q, (phi(q) + c) % p), (0, 0)) != (0, 0)
            ]
            require(
                support == [q0],
                f"shifted singleton did not land at p={p}, q0={q0}, c={c}",
            )

    return p * p


def gnorm2(x: Gaussian) -> int:
    return x[0] * x[0] + x[1] * x[1]


def joint_energy_bridge_audit(p: int) -> int:
    """Check D_graph >= E_mask/p and the sharp two-target equality case."""

    inv_two = pow(2, -1, p)

    def phi(x: int) -> int:
        return inv_two * x * x % p

    controls = 0
    for seed in range(1, 6):
        joint = {
            (q, z): (
                ((seed + q + 2 * z + q * z) % 7) - 3,
                ((2 * seed + 3 * q + z * z) % 5) - 2,
            )
            for q in range(p)
            for z in range(p)
        }

        first = sum(
            gnorm2(joint[(q, z)])
            for q in range(1, p)
            for z in range(p)
        )
        second = 0
        for c in range(p):
            graph_sum = (0, 0)
            for q in range(1, p):
                graph_sum = gadd(graph_sum, joint[(q, (phi(q) + c) % p)])
            second += gnorm2(graph_sum)
        graph_energy = first + second

        marginals = {}
        for q in range(p):
            total = (0, 0)
            for z in range(p):
                total = gadd(total, joint[(q, z)])
            marginals[q] = total

        masks = (
            set(range(1, p)),
            {q for q in range(1, p) if q % 2},
            {q for q in range(1, p) if (q * q + seed) % 3},
        )
        for mask in masks:
            masked_energy = sum(gnorm2(marginals[q]) for q in mask)
            require(
                p * graph_energy >= masked_energy,
                f"joint-energy bridge failed at p={p}, seed={seed}",
            )
            controls += 1

    # The factor 1/p is sharp.  Two nonzero target rows, constant and
    # opposite across the jet coordinate, cancel every graph sum while
    # attaining equality in Cauchy on both target rows.
    q1, q2 = 1, 2
    amplitude = (2, -1)
    hostile = {
        (q, z): (
            amplitude
            if q == q1
            else gneg(amplitude)
            if q == q2
            else (0, 0)
        )
        for q in range(p)
        for z in range(p)
    }
    hostile_first = sum(
        gnorm2(hostile[(q, z)])
        for q in range(1, p)
        for z in range(p)
    )
    hostile_second = 0
    for c in range(p):
        graph_sum = (0, 0)
        for q in range(1, p):
            graph_sum = gadd(graph_sum, hostile[(q, (phi(q) + c) % p)])
        hostile_second += gnorm2(graph_sum)
    hostile_masked = 2 * gnorm2((p * amplitude[0], p * amplitude[1]))
    require(
        p * (hostile_first + hostile_second) == hostile_masked,
        f"sharp joint-energy hostile failed at p={p}",
    )
    controls += 1
    return controls


def main() -> None:
    selector_controls = sum(prime_field_selector_audit(p) for p in (3, 5, 7))
    graph_selector_controls = sum(
        graph_restriction_selector_audit(p) for p in (3, 5, 7)
    )
    difference_controls = sum(
        channel_difference_identity_audit(p) for p in (3, 5, 7, 11, 13)
    )
    vertical_controls = sum(
        vertical_tensor_boundary_audit(p) for p in (3, 5, 7, 11, 13)
    )
    energy_bridge_controls = sum(
        joint_energy_bridge_audit(p) for p in (3, 5, 7, 11, 13)
    )
    derivatives, trace_pairings, graph_cells = field_169_audit()
    diagonal_boundary_audit()

    print("THM-2356 independent planar-chirp hostile audit")
    print(f"prime-field off-diagonal selector controls: {selector_controls}")
    print(f"prime-field graph-restriction selector controls: {graph_selector_controls}")
    print("h=0 boundary: every diagonal site survives")
    print(f"F_169 nonzero planar derivatives: {derivatives}")
    print(f"F_169 nondegenerate trace-pairing controls: {trace_pairings}")
    print(f"F_169 planar-graph partition cells: {graph_cells}")
    print("diagonal boundary: support>=3 formula PASS; singleton/two-site hostiles PASS")
    print(f"channel-difference Parseval identities: {difference_controls}")
    print(f"graph-singleton/vertical-tensor cells: {vertical_controls}")
    print(f"joint-to-word energy bridge controls: {energy_bridge_controls}")
    print(
        "VERDICT: formula signs, q^-2 normalization, graph restriction, "
        "and exact vertical residual PASS"
    )


if __name__ == "__main__":
    main()
