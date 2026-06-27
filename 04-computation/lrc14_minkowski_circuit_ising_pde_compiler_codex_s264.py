#!/usr/bin/env python3
"""S264 scout: Minkowski/circuit/Ising/PDE/De Moivre carrier compiler.

This is a dependency-free synthesis script.  It checks the exact De Moivre
quintic cancellation

  x = z - a/z
  x^5 + 5 a x^3 + 5 a^2 x = z^5 - a^5/z^5

as a Laurent-polynomial identity, then builds a small tournament over proof
carriers rather than runners.  The goal is to make the merged analogy auditable:
which sidecar retains which LRC14 proof payload, and which sidecar is only a
stress test unless another coordinate is restored.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations

Monomial = tuple[int, int]  # (z exponent, a exponent)
Poly = dict[Monomial, int]


def clean(poly: Poly) -> Poly:
    return {m: c for m, c in poly.items() if c}


def add(p: Poly, q: Poly) -> Poly:
    out = dict(p)
    for m, c in q.items():
        out[m] = out.get(m, 0) + c
    return clean(out)


def scale(c: int, p: Poly) -> Poly:
    return clean({m: c * v for m, v in p.items()})


def mul(p: Poly, q: Poly) -> Poly:
    out: Poly = {}
    for (z1, a1), c1 in p.items():
        for (z2, a2), c2 in q.items():
            m = (z1 + z2, a1 + a2)
            out[m] = out.get(m, 0) + c1 * c2
    return clean(out)


def pow_poly(p: Poly, n: int) -> Poly:
    out: Poly = {(0, 0): 1}
    for _ in range(n):
        out = mul(out, p)
    return out


def fmt_poly(poly: Poly) -> str:
    if not poly:
        return "0"
    bits: list[str] = []
    for (z, a), coeff in sorted(poly.items(), reverse=True):
        term_bits = []
        if z:
            term_bits.append("z" if z == 1 else f"z^{z}")
        if a:
            term_bits.append("a" if a == 1 else f"a^{a}")
        term = "*".join(term_bits) or "1"
        bits.append(f"{coeff:+d}*{term}")
    text = " ".join(bits)
    return text[1:] if text.startswith("+") else text


def de_moivre_identity() -> tuple[bool, Poly, Poly]:
    x = {(1, 0): 1, (-1, 1): -1}
    a = {(0, 1): 1}
    lhs = add(add(pow_poly(x, 5), scale(5, mul(a, pow_poly(x, 3)))), scale(5, mul(pow_poly(a, 2), x)))
    rhs = {(5, 0): 1, (-5, 5): -1}
    residual = add(lhs, scale(-1, rhs))
    return residual == {}, lhs, residual


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves_lrc: int
    names_lost_coordinate: int
    finite_checkable: int
    tail_control: int
    zero_mode_control: int
    transfer_or_partition: int
    algebraic_cancellation: int
    circuit_uniformity: int
    lean_facing: int
    packet_exit: int

    @property
    def score(self) -> int:
        return sum(
            [
                self.preserves_lrc,
                self.names_lost_coordinate,
                self.finite_checkable,
                self.tail_control,
                self.zero_mode_control,
                self.transfer_or_partition,
                self.algebraic_cancellation,
                self.circuit_uniformity,
                self.lean_facing,
                self.packet_exit,
            ]
        )


CARRIERS = [
    Carrier("finite_address_packet", 1, 1, 1, 0, 1, 0, 0, 1, 1, 1),
    Carrier("endpoint_phi_gap_operator", 1, 1, 1, 0, 1, 0, 0, 1, 1, 1),
    Carrier("pde_weak_form_stiffness", 1, 1, 1, 0, 1, 0, 0, 0, 1, 1),
    Carrier("lee_yang_ising_transfer", 1, 1, 1, 0, 1, 1, 0, 0, 0, 1),
    Carrier("minkowski_successive_minima", 1, 1, 1, 1, 0, 0, 0, 0, 0, 1),
    Carrier("circuit_complexity_audit", 0, 1, 1, 0, 0, 0, 0, 1, 1, 1),
    Carrier("de_moivre_quintic_fold", 0, 1, 1, 0, 0, 0, 1, 1, 0, 0),
    Carrier("jacobi_theta_signed_tail", 1, 1, 0, 1, 0, 1, 0, 0, 0, 0),
    Carrier("raw_scalar_p0", 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
]


def orient(a: Carrier, b: Carrier) -> str:
    if a.score != b.score:
        return a.name if a.score > b.score else b.name
    # Tie path: preserve the explicit order in CARRIERS as a Hamiltonian
    # retention path when payload scores match.
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    return a.name if order[a.name] < order[b.name] else b.name


def tournament_fingerprint() -> None:
    wins = Counter({carrier.name: 0 for carrier in CARRIERS})
    margins: list[tuple[int, str, str]] = []
    adjacency = {carrier.name: set() for carrier in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        winner = orient(a, b)
        loser = b.name if winner == a.name else a.name
        wins[winner] += 1
        adjacency[winner].add(loser)
        margins.append((abs(a.score - b.score), a.name, b.name))

    directed_3cycles = 0
    for a, b, c in combinations([x.name for x in CARRIERS], 3):
        ab = b in adjacency[a]
        bc = c in adjacency[b]
        ca = a in adjacency[c]
        ba = a in adjacency[b]
        cb = b in adjacency[c]
        ac = c in adjacency[a]
        if (ab and bc and ca) or (ba and cb and ac):
            directed_3cycles += 1

    ordered = sorted(CARRIERS, key=lambda x: (-x.score, [c.name for c in CARRIERS].index(x.name)))
    print("TOURNAMENT ANALYSIS: proof carriers as vertices, not runners")
    print("  axes: preserves_lrc, names_lost_coordinate, finite_checkable, tail_control,")
    print("        zero_mode_control, transfer_or_partition, algebraic_cancellation,")
    print("        circuit_uniformity, lean_facing, packet_exit")
    print("  scores:")
    for carrier in ordered:
        print(f"    {carrier.name:31s} score={carrier.score}")
    print(f"  score_hist={dict(sorted(Counter(c.score for c in CARRIERS).items()))}")
    print(f"  directed_3cycles={directed_3cycles}")
    print("  Hamilton retention path:")
    print("    " + " -> ".join(c.name for c in ordered))
    print("  low-margin edge-flip risks:")
    for margin, a, b in sorted(margins)[:10]:
        print(f"    margin={margin}: {a} vs {b}")


def compiler_routes() -> None:
    print("\nMERGED ROUTES")
    routes = [
        (
            "support-six tail",
            [
                "low-height wall ledger",
                "relation lattice",
                "minkowski_successive_minima",
                "jacobi_theta_signed_tail",
                "finite_address_packet",
            ],
        ),
        (
            "coverage/free-energy packet",
            [
                "lee_yang_ising_transfer",
                "pde_weak_form_stiffness",
                "endpoint_phi_gap_operator",
                "observer_gluing_certificate",
                "finite_address_packet",
            ],
        ),
        (
            "finite-depth algebraic fold",
            [
                "de_moivre_quintic_fold",
                "circuit_complexity_audit",
                "lost-coordinate ledger",
                "finite_address_packet",
            ],
        ),
    ]
    for label, route in routes:
        print(f"  {label:28s}: " + " -> ".join(route))


def main() -> None:
    print("S264 Minkowski/circuit/Ising/PDE/De Moivre proof-carrier compiler")
    ok, lhs, residual = de_moivre_identity()
    print("\nDE MOIVRE QUINTIC FOLD")
    print("  identity: x=z-a/z => x^5+5*a*x^3+5*a^2*x = z^5-a^5/z^5")
    print(f"  verified={ok}")
    print(f"  folded polynomial={fmt_poly(lhs)}")
    print(f"  residual={fmt_poly(residual)}")
    print("  readout: this is a finite-depth cancellation detector, not a global")
    print("  quintic solver for arbitrary LRC packets.")
    compiler_routes()
    print("\nPDE LOOKBACK")
    print("  repo anchors: tournament_matrix_expansion_atlas entries 103-110,")
    print("  HYP-2108 endpoint-cover circuit positivity, HYP-2112 Phi gap functional,")
    print("  HYP-3101 normal-fan/Cech component bound, HYP-3108/HYP-3109 phi4/Lee-Yang.")
    print("  PDE translation: weak forms keep mass/stiffness/zero-mode data;")
    print("  any quotient that drops them is like solving a PDE after deleting")
    print("  boundary conditions.")
    print()
    tournament_fingerprint()


if __name__ == "__main__":
    main()
