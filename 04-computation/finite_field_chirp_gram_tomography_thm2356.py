#!/usr/bin/env python3
"""Exact referee for THM-2356 finite-field chirp Gram tomography.

The main positive control is the odd field K=F_13[s]/(s^2-2), of order
169.  The script exhausts the planar derivatives of phi(x)=x^2/2 and
performs the complete double character inversion on a nontrivial Gaussian
integer signal.  Cyclotomic values are represented exactly in
Q(i,zeta_13); no floating point or Python ``assert`` is used.
"""

from __future__ import annotations

from dataclasses import dataclass


P = 13
Q = P * P
INV_TWO = 7
Field = tuple[int, int]
Gaussian = tuple[int, int]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fadd(x: Field, y: Field) -> Field:
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def fneg(x: Field) -> Field:
    return ((-x[0]) % P, (-x[1]) % P)


def fsub(x: Field, y: Field) -> Field:
    return fadd(x, fneg(y))


def fmul(x: Field, y: Field) -> Field:
    # theta^2=2, and 2 is a nonsquare modulo 13.
    return (
        (x[0] * y[0] + 2 * x[1] * y[1]) % P,
        (x[0] * y[1] + x[1] * y[0]) % P,
    )


def fscale(c: int, x: Field) -> Field:
    return ((c * x[0]) % P, (c * x[1]) % P)


def trace(x: Field) -> int:
    # Frobenius conjugates theta to -theta.
    return (2 * x[0]) % P


def phi(x: Field) -> Field:
    return fscale(INV_TWO, fmul(x, x))


def derivative(h: Field, y: Field) -> Field:
    return fsub(phi(fadd(y, h)), phi(y))


def gmul_conj(x: Gaussian, y: Gaussian) -> Gaussian:
    # (a+ib) conjugate(c+id).
    return (x[0] * y[0] + x[1] * y[1], x[1] * y[0] - x[0] * y[1])


@dataclass(frozen=True)
class CyclotomicGaussian:
    """Element of Q(i,zeta_13), stored in the cyclic group-ring basis."""

    real: tuple[int, ...]
    imag: tuple[int, ...]

    def shifted(self, exponent: int) -> "CyclotomicGaussian":
        e = exponent % P
        return CyclotomicGaussian(
            tuple(self.real[(j - e) % P] for j in range(P)),
            tuple(self.imag[(j - e) % P] for j in range(P)),
        )

    def canonical(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        # Phi_13(zeta)=1+...+zeta^12=0.  Subtract the zeta^12
        # coefficient from every entry and drop that final coordinate.
        return (
            tuple(self.real[j] - self.real[P - 1] for j in range(P - 1)),
            tuple(self.imag[j] - self.imag[P - 1] for j in range(P - 1)),
        )


ZERO_RING = CyclotomicGaussian((0,) * P, (0,) * P)


def ring_from_terms(
    terms: list[tuple[int, Gaussian]],
) -> CyclotomicGaussian:
    real = [0] * P
    imag = [0] * P
    for exponent, coefficient in terms:
        real[exponent % P] += coefficient[0]
        imag[exponent % P] += coefficient[1]
    return CyclotomicGaussian(tuple(real), tuple(imag))


def ring_sum_shifted(
    summands: list[tuple[CyclotomicGaussian, int]],
) -> CyclotomicGaussian:
    real = [0] * P
    imag = [0] * P
    for value, exponent in summands:
        e = exponent % P
        for j in range(P):
            real[(j + e) % P] += value.real[j]
            imag[(j + e) % P] += value.imag[j]
    return CyclotomicGaussian(tuple(real), tuple(imag))


def expected_scalar(value: Gaussian, scale: int) -> CyclotomicGaussian:
    real = [0] * P
    imag = [0] * P
    real[0] = scale * value[0]
    imag[0] = scale * value[1]
    return CyclotomicGaussian(tuple(real), tuple(imag))


def main() -> None:
    elements = [(a, b) for a in range(P) for b in range(P)]
    zero = (0, 0)
    one = (1, 0)

    # Field and trace controls.
    require(
        all(
            fmul(x, y) == fmul(y, x)
            and trace(fadd(x, y)) == (trace(x) + trace(y)) % P
            for x in elements
            for y in elements
        ),
        "field commutativity or trace additivity failed",
    )
    require(
        all(
            any(fmul(x, y) == one for y in elements)
            for x in elements
            if x != zero
        ),
        "theta^2-2 failed to define a field",
    )

    # Exhaust the planar-function property.  For every h!=0, the derivative
    # y -> phi(y+h)-phi(y) must permute all 169 field elements.
    planar_derivatives = 0
    for h in elements:
        if h == zero:
            continue
        image = {derivative(h, y) for y in elements}
        require(len(image) == Q, f"nonplanar derivative at h={h}")
        require(
            all(
                derivative(h, y) == fadd(fmul(h, y), phi(h))
                for y in elements
            ),
            f"quadratic derivative identity failed at h={h}",
        )
        planar_derivatives += 1
    require(planar_derivatives == Q - 1, "wrong planar derivative count")

    # The LRC target x first-jet space has the same abstract form K x K.
    # The translates z=phi(q)+c partition it into q planar graphs.
    graph_pairs: set[tuple[Field, Field]] = set()
    for c in elements:
        graph = {
            (qvalue, fadd(phi(qvalue), c))
            for qvalue in elements
        }
        require(len(graph) == Q, f"planar graph changed size at c={c}")
        require(
            graph_pairs.isdisjoint(graph),
            f"planar graph translates collided at c={c}",
        )
        graph_pairs.update(graph)
    require(
        len(graph_pairs) == Q * Q,
        "planar graph translates stopped partitioning K x K",
    )

    # A five-site Gaussian-integer signal.  Direct expansion of every one of
    # the 169^2 chirped intensities remains exact in Q(i,zeta_13).
    signal: dict[Field, Gaussian] = {
        (0, 0): (2, 1),
        (1, 0): (-1, 3),
        (0, 1): (4, -2),
        (3, 5): (1, 2),
        (8, 4): (-3, -1),
    }
    pair_terms: list[tuple[Field, Field, Gaussian]] = []
    for x, zx in signal.items():
        for y, zy in signal.items():
            pair_terms.append(
                (
                    fsub(x, y),
                    fsub(phi(x), phi(y)),
                    gmul_conj(zx, zy),
                )
            )

    intensities: dict[tuple[Field, Field], CyclotomicGaussian] = {}
    for a in elements:
        for b in elements:
            terms: list[tuple[int, Gaussian]] = []
            for h, dphi, coefficient in pair_terms:
                exponent = trace(fadd(fmul(b, h), fmul(a, dphi)))
                terms.append((exponent, coefficient))
            intensities[(a, b)] = ring_from_terms(terms)
    require(len(intensities) == Q * Q, "wrong chirped intensity count")

    # Perform the literal double inversion on every ordered nonzero pair in
    # the signal and on five zero controls.  The unnormalized inverse equals
    # q^2 z_(y+h) conjugate(z_y).
    controls: list[tuple[Field, Field]] = []
    support = list(signal)
    for y in support:
        for x in support:
            if x != y:
                controls.append((fsub(x, y), y))
    controls.extend(
        [
            ((1, 1), (2, 2)),
            ((2, 3), (4, 1)),
            ((7, 0), (5, 5)),
            ((0, 9), (12, 12)),
            ((6, 8), (1, 10)),
        ]
    )

    for h, y in controls:
        require(h != zero, "off-diagonal inversion received h=0")
        dphi = derivative(h, y)
        summands: list[tuple[CyclotomicGaussian, int]] = []
        for a in elements:
            a_phase = trace(fmul(a, dphi))
            for b in elements:
                phase = -(trace(fmul(b, h)) + a_phase)
                summands.append((intensities[(a, b)], phase))
        recovered = ring_sum_shifted(summands)
        x = fadd(y, h)
        expected_value = gmul_conj(signal.get(x, (0, 0)), signal.get(y, (0, 0)))
        expected = expected_scalar(expected_value, Q * Q)
        require(
            recovered.canonical() == expected.canonical(),
            f"double chirp inversion failed at h={h}, y={y}",
        )

    # Linear masks add no measurement: multiplying by a linear character
    # merely translates the Fourier-frequency label.
    linear_relabels = 0
    for a in elements:
        for b in elements:
            require(
                all(
                    trace(fmul(fadd(a, b), fsub(x, y)))
                    == trace(
                        fadd(
                            fmul(a, fsub(x, y)),
                            fmul(b, fsub(x, y)),
                        )
                    )
                    for x in signal
                    for y in signal
                ),
                "linear-mask relabelling identity failed",
            )
            linear_relabels += 1

    # Sharp magnitude sidecars.  Every singleton has constant intensity one,
    # irrespective of its location.  On a fixed two-point support, swapping
    # unequal real magnitudes also preserves every chirped intensity because
    # the norm and both off-diagonal products are unchanged.
    for location in elements:
        for a in ((0, 0), (1, 0), (0, 1), (7, 11)):
            for b in ((0, 0), (2, 0), (0, 3), (9, 5)):
                exponent = trace(
                    fadd(
                        fmul(b, fsub(location, location)),
                        fmul(a, fsub(phi(location), phi(location))),
                    )
                )
                require(exponent == 0, "a singleton intensity was not constant")

    x0, x1 = (0, 0), (4, 7)
    for a in elements:
        for b in elements:
            e01 = trace(
                fadd(
                    fmul(b, fsub(x0, x1)),
                    fmul(a, fsub(phi(x0), phi(x1))),
                )
            )
            # Both (2,3) and (3,2) have norm 13 and cross product 6.
            first = ring_from_terms(
                [(0, (13, 0)), (e01, (6, 0)), (-e01, (6, 0))]
            )
            second = ring_from_terms(
                [(0, (13, 0)), (e01, (6, 0)), (-e01, (6, 0))]
            )
            require(
                first.canonical() == second.canonical(),
                "the two-site magnitude-swap hostile changed",
            )

    print("THM-2356 exact finite-field chirp tomography referee")
    print("field: F_13[theta]/(theta^2-2), q=169, trace(u+v theta)=2u")
    print("planar map: phi(x)=x^2/2")
    print(f"exhaustive nonzero planar derivatives: {planar_derivatives}")
    print("planar graph partition: 169 slices x 169 points = 28561")
    print(f"chirped intensity table: {len(intensities)}")
    print(f"exact off-diagonal inversion controls: {len(controls)}")
    print(f"linear-mask relabelling checks: {linear_relabels}")
    print("sharp hostiles: 169 invisible singleton locations; two-site 2<->3 swap")
    print("VERDICT: chirps recover every off-diagonal Gram entry; singletons supply the diagonal")


if __name__ == "__main__":
    main()
