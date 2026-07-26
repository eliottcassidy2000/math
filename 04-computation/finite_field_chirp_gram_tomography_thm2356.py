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
from fractions import Fraction


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


def gadd(x: Gaussian, y: Gaussian) -> Gaussian:
    return (x[0] + y[0], x[1] + y[1])


def gneg(x: Gaussian) -> Gaussian:
    return (-x[0], -x[1])


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


def signal_intensity(
    signal: dict[Field, Gaussian],
    a: Field,
    b: Field,
) -> CyclotomicGaussian:
    terms: list[tuple[int, Gaussian]] = []
    for x, zx in signal.items():
        for y, zy in signal.items():
            exponent = trace(
                fadd(
                    fmul(b, fsub(x, y)),
                    fmul(a, fsub(phi(x), phi(y))),
                )
            )
            terms.append((exponent, gmul_conj(zx, zy)))
    return ring_from_terms(terms)


def signal_transform(
    signal: dict[Field, Gaussian],
    b: Field,
) -> CyclotomicGaussian:
    return ring_from_terms(
        [(trace(fmul(b, x)), value) for x, value in signal.items()]
    )


def joint_transform(
    joint: dict[tuple[Field, Field], Gaussian],
    eta: Field,
    b: Field,
) -> CyclotomicGaussian:
    return ring_from_terms(
        [
            (
                trace(fadd(fmul(eta, zvalue), fmul(b, qvalue))),
                value,
            )
            for (qvalue, zvalue), value in joint.items()
        ]
    )


def coarse_sum(
    graph_signals: dict[Field, dict[Field, Gaussian]],
) -> dict[Field, Gaussian]:
    result: dict[Field, Gaussian] = {}
    for signal in graph_signals.values():
        for qvalue, value in signal.items():
            result[qvalue] = gadd(result.get(qvalue, (0, 0)), value)
    return {
        qvalue: value
        for qvalue, value in result.items()
        if value != (0, 0)
    }


def signal_total(signal: dict[Field, Gaussian]) -> Gaussian:
    total = (0, 0)
    for value in signal.values():
        total = gadd(total, value)
    return total


def singleton_ledger(
    signal: dict[Field, Gaussian],
) -> dict[Field, int]:
    return {
        location: gmul_conj(value, value)[0]
        for location, value in signal.items()
    }


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

    # The exact refined-zero locus is the vertical tensor
    # delta_0(q) B(z).  Tensor the THM-2333 full-support convolution
    # inverse with a deterministic nowhere-zero jet profile.  Every joint
    # fibre then contains 169 nonzero atomic terms, every planar graph
    # survives, and every surviving target is nevertheless q=0.
    p_shift = (3, 4)

    def endpoint_u(x: Field) -> Fraction:
        return Fraction(2 if x == zero else 1)

    def endpoint_v(x: Field) -> Fraction:
        return Fraction(Q if x == p_shift else -1, Q + 1)

    target_hostile: dict[Field, Fraction] = {}
    for qvalue in elements:
        terms = [
            endpoint_u(u)
            * endpoint_v(fadd(fsub(u, qvalue), p_shift))
            for u in elements
        ]
        require(
            all(term != 0 for term in terms),
            f"vertical-tensor hostile lost an atomic term at q={qvalue}",
        )
        target_hostile[qvalue] = sum(terms, Fraction(0))
        require(
            target_hostile[qvalue]
            == (Fraction(1) if qvalue == zero else Fraction(0)),
            f"target convolution inverse failed at q={qvalue}",
        )

    def jet_profile(zvalue: Field) -> Fraction:
        return Fraction(1 + zvalue[0] + P * zvalue[1])

    require(
        all(jet_profile(zvalue) != 0 for zvalue in elements),
        "deterministic jet profile stopped being full support",
    )
    graph_singletons = 0
    for c in elements:
        for qvalue in elements:
            joint_value = (
                target_hostile[qvalue]
                * jet_profile(fadd(phi(qvalue), c))
            )
            expected_joint = (
                jet_profile(c) if qvalue == zero else Fraction(0)
            )
            require(
                joint_value == expected_joint,
                f"vertical graph restriction failed at c={c}, q={qvalue}",
            )
        graph_singletons += 1
    require(graph_singletons == Q, "wrong vertical singleton graph count")

    # B=1 saturates support uncertainty on K x K: delta_0(q) has all
    # 169 target-character frequencies, while 1_K(z) has only the
    # trivial jet-character frequency.  Check the latter character sum
    # exactly in Q(zeta_13).
    jet_fourier_support = 0
    for d in elements:
        character_sum = ring_from_terms(
            [(trace(fmul(d, zvalue)), (1, 0)) for zvalue in elements]
        )
        expected_character_sum = expected_scalar(
            (1, 0),
            Q if d == zero else 0,
        )
        require(
            character_sum.canonical()
            == expected_character_sum.canonical(),
            f"jet character orthogonality failed at d={d}",
        )
        if d == zero:
            jet_fourier_support += 1
    joint_spatial_support = Q
    joint_fourier_support = Q * jet_fourier_support
    require(
        joint_spatial_support * joint_fourier_support == Q * Q,
        "vertical tensor stopped saturating support uncertainty",
    )

    # The same target delta is the sharp separate-degree footprint
    # extremizer on the 13 x 13 target grid.
    footprint_support = 0
    for x0_value in range(P):
        for x1_value in range(P):
            polynomial_delta = (
                (1 - pow(x0_value, P - 1, P))
                * (1 - pow(x1_value, P - 1, P))
            ) % P
            indicator_delta = int(x0_value == 0 and x1_value == 0)
            require(
                polynomial_delta == indicator_delta,
                "finite-field delta footprint identity failed",
            )
            footprint_support += indicator_delta
    require(footprint_support == 1, "wrong target footprint support")
    # Refined-to-coarse hostile.  One planar graph has two target points,
    # including a nonzero target, but a singleton on a second graph cancels
    # that target in the coarser sum over the jet coordinate.
    target_a = (1, 0)
    graph_c = (0, 1)
    minimal_graphs = {
        zero: {zero: (1, 0), target_a: (1, 0)},
        graph_c: {target_a: (-1, 0)},
    }
    minimal_joint = {
        (zero, zero): (1, 0),
        (target_a, phi(target_a)): (1, 0),
        (target_a, fadd(phi(target_a), graph_c)): (-1, 0),
    }
    minimal_coarse = coarse_sum(minimal_graphs)
    require(
        minimal_coarse == {zero: (1, 0)},
        "the three-atom refined/coarse hostile stopped cancelling",
    )
    require(
        gmul_conj(
            minimal_graphs[zero][target_a],
            minimal_graphs[zero][zero],
        )
        == (1, 0),
        "the two-supported graph lost its off-diagonal Gram entry",
    )
    require(
        all(
            signal_transform(minimal_coarse, b).canonical()
            == expected_scalar((1, 0), 1).canonical()
            for b in elements
        ),
        "the coarse target twists of the minimal hostile were not constant",
    )
    separating_eta = (0, 1)
    require(
        trace(fmul(separating_eta, graph_c)) != 0,
        "the chosen jet character did not separate the hostile graph labels",
    )
    require(
        all(
            joint_transform(minimal_joint, zero, b).canonical()
            == expected_scalar((1, 0), 1).canonical()
            for b in elements
        ),
        "the trivial jet polarizer did not reproduce the coarse hostile",
    )
    require(
        joint_transform(minimal_joint, separating_eta, zero).canonical()
        != joint_transform(minimal_joint, separating_eta, one).canonical(),
        "the separating jet polarizer had no target-edge defect",
    )

    # Independent graph phases are not glued by graphwise intensities or
    # singleton energies.  These two families differ only by the global sign
    # on graph c.  They have the same graph data and the same total current,
    # but one coarse target is delta_0 and the other is 3 delta_0-2 delta_a.
    graph_d = (1, 1)
    base_graph = {zero: (1, 0), target_a: (-1, 0)}
    cancel_graphs = {
        zero: base_graph,
        graph_c: {
            location: gneg(value)
            for location, value in base_graph.items()
        },
        graph_d: {zero: (1, 0)},
    }
    escape_graphs = {
        zero: base_graph,
        graph_c: dict(base_graph),
        graph_d: {zero: (1, 0)},
    }
    phase_gauge_checks = 0
    for a in elements:
        for b in elements:
            require(
                signal_intensity(cancel_graphs[graph_c], a, b).canonical()
                == signal_intensity(escape_graphs[graph_c], a, b).canonical(),
                "a global graph sign changed its chirp intensity",
            )
            phase_gauge_checks += 1
    require(
        all(
            singleton_ledger(cancel_graphs[label])
            == singleton_ledger(escape_graphs[label])
            for label in (zero, graph_c, graph_d)
        ),
        "the graph-phase pair changed a singleton ledger",
    )
    cancel_coarse = coarse_sum(cancel_graphs)
    escape_coarse = coarse_sum(escape_graphs)
    require(
        cancel_coarse == {zero: (1, 0)},
        "the same-data cancel family lost its zero-only coarse target",
    )
    require(
        escape_coarse == {zero: (3, 0), target_a: (-2, 0)},
        "the same-data escape family changed its coarse target",
    )
    require(
        signal_total(cancel_coarse) == (1, 0)
        and signal_total(escape_coarse) == (1, 0),
        "the same-data pair stopped having the same total current",
    )
    require(
        signal_transform(escape_coarse, zero).canonical()
        != signal_transform(escape_coarse, one).canonical(),
        "the escape family lost its nonconstant coarse target response",
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
    print(
        "vertical tensor hostile: 28561 termwise-full joint fibres; "
        "169 one-sparse zero-target graphs"
    )
    print(
        "uncertainty/footprint extremizer: supports 169 x 169; "
        "target grid support 1"
    )
    print("refined/coarse hostile: 3 atoms, one support-2 graph, coarse C=delta_0")
    print("jet polarizer: trivial target edge zero; separating-jet edge nonzero")
    print(
        "independent graph-phase checks: "
        f"{phase_gauge_checks}; same graph data and total, coarse energy 0 vs >0"
    )
    print(
        "VERDICT: chirps recover within-graph Gram data; "
        "coarse target needs cross-graph phase gluing"
    )


if __name__ == "__main__":
    main()
