#!/usr/bin/env python3
"""HYP-3229 modular/arithmetic sidecar audit for the LRC(14) Fejer route.

The point of this scout is not to prove LRC(14).  It separates certificate
payload from attractive arithmetic coincidences around the HYP-3214 Fejer /
de-Moivre magic function.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Sequence, Tuple


Poly = List[Fraction]  # low-to-high coefficients


def trim(p: Poly) -> Poly:
    p = p[:]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def poly_add(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    out = [Fraction(0) for _ in range(n)]
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
    return trim(out)


def poly_sub(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    out = [Fraction(0) for _ in range(n)]
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
    return trim(out)


def poly_mul(a: Poly, b: Poly) -> Poly:
    out = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return trim(out)


def poly_divmod(a: Poly, b: Poly) -> Tuple[Poly, Poly]:
    a = trim(a)
    b = trim(b)
    if b == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(a) < len(b):
        return [Fraction(0)], a
    q = [Fraction(0) for _ in range(len(a) - len(b) + 1)]
    r = a[:]
    while len(r) >= len(b) and r != [0]:
        coeff = r[-1] / b[-1]
        shift = len(r) - len(b)
        q[shift] = coeff
        subtractor = [Fraction(0)] * shift + [coeff * c for c in b]
        r = poly_sub(r, subtractor)
    return trim(q), trim(r)


def poly_compose(p: Poly, q: Poly) -> Poly:
    out = [Fraction(0)]
    power = [Fraction(1)]
    for coeff in p:
        out = poly_add(out, [coeff * c for c in power])
        power = poly_mul(power, q)
    return trim(out)


def v_n_poly(n: int) -> Poly:
    """V_n where V_n(2 cos t) = 2 cos(nt)."""
    if n == 0:
        return [Fraction(2)]
    if n == 1:
        return [Fraction(0), Fraction(1)]
    v0 = [Fraction(2)]
    v1 = [Fraction(0), Fraction(1)]
    x = [Fraction(0), Fraction(1)]
    for _ in range(2, n + 1):
        v0, v1 = v1, poly_sub(poly_mul(x, v1), v0)
    return v1


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def poly_str(p: Poly, var: str = "u") -> str:
    p = trim(p)
    pieces: List[str] = []
    for degree in range(len(p) - 1, -1, -1):
        coeff = p[degree]
        if coeff == 0:
            continue
        sign = "-" if coeff < 0 else "+"
        abs_coeff = abs(coeff)
        if degree == 0:
            body = fmt_frac(abs_coeff)
        elif degree == 1:
            body = var if abs_coeff == 1 else f"{fmt_frac(abs_coeff)}*{var}"
        else:
            body = f"{var}^{degree}" if abs_coeff == 1 else f"{fmt_frac(abs_coeff)}*{var}^{degree}"
        pieces.append((sign, body))
    if not pieces:
        return "0"
    first_sign, first_body = pieces[0]
    out = first_body if first_sign == "+" else f"-{first_body}"
    for sign, body in pieces[1:]:
        out += f" {sign} {body}"
    return out


def cubic_discriminant(a: int, b: int, c: int, d: int) -> int:
    return b * b * c * c - 4 * a * c**3 - 4 * b**3 * d - 27 * a * a * d * d + 18 * a * b * c * d


def sigma1(n: int) -> int:
    return sum(d for d in range(1, n + 1) if n % d == 0)


def gamma0_7_eisenstein_coeff(n: int) -> int:
    """Coefficient of q^n in (7 E_2(7 tau) - E_2(tau))/6, constant term omitted."""
    base = 4 * sigma1(n)
    if n % 7 == 0:
        base -= 28 * sigma1(n // 7)
    return base


@dataclass(frozen=True)
class Qomega:
    """Element a + b*omega with omega^2 + omega + 1 = 0."""

    a: Fraction
    b: Fraction = Fraction(0)

    def __add__(self, other: "Qomega") -> "Qomega":
        return Qomega(self.a + other.a, self.b + other.b)

    def __neg__(self) -> "Qomega":
        return Qomega(-self.a, -self.b)

    def __sub__(self, other: "Qomega") -> "Qomega":
        return self + (-other)

    def scale(self, c: Fraction) -> "Qomega":
        return Qomega(c * self.a, c * self.b)

    def __mul__(self, other: "Qomega") -> "Qomega":
        # (a+bw)(c+dw)=ac+(ad+bc)w+bd w^2 = (ac-bd)+(ad+bc-bd)w.
        return Qomega(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    def denominator_lcm(self) -> int:
        return math.lcm(self.a.denominator, self.b.denominator)

    def __str__(self) -> str:
        a = fmt_frac(self.a)
        b = fmt_frac(abs(self.b))
        if self.b == 0:
            return a
        if self.a == 0:
            return f"{'-' if self.b < 0 else ''}{b}*omega"
        sign = "-" if self.b < 0 else "+"
        return f"{a} {sign} {b}*omega"


OMEGA_POWERS = [
    Qomega(Fraction(1), Fraction(0)),
    Qomega(Fraction(0), Fraction(1)),
    Qomega(Fraction(-1), Fraction(-1)),
]


def bernoulli2(x: Fraction) -> Fraction:
    return x * x - x + Fraction(1, 6)


def chi_mod7_even_value(k: int, a: int) -> Qomega:
    """Character chi_k on (Z/7)^*, restricted to even k=0,2,4 in Q(omega)."""
    generator_logs = {1: 0, 3: 1, 2: 2, 6: 3, 4: 4, 5: 5}
    if k == 0:
        return Qomega(Fraction(1), Fraction(0))
    if k not in (2, 4):
        raise ValueError("this exact Q(omega) helper only handles k=0,2,4")
    exponent_mod_3 = ((k // 2) * generator_logs[a]) % 3
    return OMEGA_POWERS[exponent_mod_3]


def generalized_bernoulli_B2_mod7(k: int) -> Qomega:
    total = Qomega(Fraction(0), Fraction(0))
    for a in range(1, 7):
        total += chi_mod7_even_value(k, a).scale(bernoulli2(Fraction(a, 7)))
    return total.scale(Fraction(7))


def tournament_fingerprint(vertices: Sequence[str], score: Dict[str, Tuple[int, int, int]]) -> Dict[str, object]:
    """Orient a proof-carrier tournament by lexicographic score."""
    edges: Dict[Tuple[str, str], str] = {}
    for a, b in itertools.combinations(vertices, 2):
        winner = a if score[a] > score[b] else b
        edges[(a, b)] = winner
    outdegrees = {v: 0 for v in vertices}
    for (a, b), winner in edges.items():
        outdegrees[winner] += 1
    cycles = []
    for triple in itertools.combinations(vertices, 3):
        wins = {v: 0 for v in triple}
        for a, b in itertools.combinations(triple, 2):
            winner = edges[(a, b)] if (a, b) in edges else edges[(b, a)]
            wins[winner] += 1
        if sorted(wins.values()) == [1, 1, 1]:
            cycles.append(triple)
    path = sorted(vertices, key=lambda v: score[v], reverse=True)
    return {
        "score_hist": {d: list(outdegrees.values()).count(d) for d in sorted(set(outdegrees.values()))},
        "directed_3_cycles": len(cycles),
        "hamiltonian_path": path,
        "outdegrees": outdegrees,
    }


def main() -> None:
    print("=== HYP-3229 modular magic sidecar audit for LRC(14) ===")

    de_moivre = [Fraction(-1), Fraction(-2), Fraction(1), Fraction(1)]
    v7 = v_n_poly(7)
    quotient, remainder = poly_divmod(poly_sub(v7, [Fraction(2)]), [Fraction(-2), Fraction(1)])
    square = poly_mul(de_moivre, de_moivre)
    print("\n[1] Exact Fejer / de-Moivre / Chebyshev identity")
    print(f"V_7(u) = {poly_str(v7)}")
    print(f"m(u) = {poly_str(de_moivre)}")
    print(f"(V_7(u)-2)/(u-2) = {poly_str(quotient)}")
    print(f"m(u)^2 = {poly_str(square)}")
    print(f"identity_exact = {quotient == square and remainder == [0]}")
    print(f"disc(m) = {cubic_discriminant(1, 1, -2, -1)} = 7^2")
    print("Fejer coefficients: " + ", ".join(f"{n}:{7 - abs(n)}" for n in range(-6, 7)))

    b = [Fraction(-2), Fraction(1)]  # u = B - 2
    beraha_poly = poly_compose(de_moivre, b)
    roots = [2 * math.cos(2 * math.pi * j / 7) for j in (1, 2, 3)]
    mahler_measure = math.prod(max(1.0, abs(r)) for r in roots)
    b7 = 2 + 2 * math.cos(2 * math.pi / 7)
    print("\n[2] Beraha/Tutte and Mahler sidecars")
    print(f"minimal polynomial for B_7 = 2 + 2cos(2pi/7): {poly_str(beraha_poly, 'B')}")
    print(f"B_7 numeric = {b7:.15f}")
    print(f"Mahler(m) numeric = {mahler_measure:.15f}")
    print(f"Mahler(m) - (B_7 - 1) = {mahler_measure - (b7 - 1):.3e}")
    print("status = arithmetic height sidecar; not yet a certificate inequality")

    print("\n[3] Explicit Gamma0(7) Eisenstein candidate")
    print("E_7(tau) := (7 E_2(7 tau) - E_2(tau))/6")
    coeffs = [gamma0_7_eisenstein_coeff(n) for n in range(1, 29)]
    print("q^1..q^28 coefficients:")
    print(coeffs)
    print("multiples of 7 retain the local conductor correction a_n=4*sigma1(n)-28*sigma1(n/7)")
    print("certificate use = q-expansion generator for a level-7 LP dual sidecar")

    print("\n[4] Dirichlet-L / Stark sidecar at conductor 7")
    raw_grid_denoms = [bernoulli2(Fraction(a, 7)).denominator for a in range(1, 7)]
    print(f"raw B_2(a/7) denominator lcm = {math.lcm(*raw_grid_denoms)}")
    for k in (0, 2, 4):
        b2 = generalized_bernoulli_B2_mod7(k)
        lminus1 = b2.scale(Fraction(-1, 2))
        label = "principal" if k == 0 else f"even primitive chi_{k}"
        print(
            f"{label}: B_2,chi = {b2}; L(-1,chi) = {lminus1}; "
            f"den_lcm={lminus1.denominator_lcm()}"
        )
    print("guardrail = 7^2 appears in the sampling grid/discriminant, but L(-1,chi) denominators contract")
    print("status = inconclusive for a closed cap formula; retain as conductor-7 sidecar, not proof core")

    print("\n[5] Comb-overlap Gram kernel bridge (mac-mini S75 integration)")
    print("K(p,q) := meas(D_p cap D_q) = <1_Dp, 1_Dq> is PSD by construction")
    print("S75 proved the single-arc row K(1,q)=1/(7q), q<=13:")
    print(", ".join(f"K(1,{q})={fmt_frac(Fraction(1, 7 * q))}" for q in range(1, 14)))
    print(f"K(7,q) constant sidecar = {fmt_frac(Fraction(1, 49))}")
    print("single-arc peeling: cap(P)=cap(P\\{1})-(1/7)*(1-1/min(P\\{1})) for speeds <=13")
    print("certificate use = finite spatial Gram basis dual to the Fejer Fourier kernel")
    print("remaining debt = order-3 triple-overlap constants, not the order-2 PSD kernel")

    print("\n[6] Subshift / transfer-operator compression")
    autocorr = {lag: 7 - abs(lag) for lag in range(-6, 7)}
    print(f"word autocorrelation of P(z)=1+...+z^6: {autocorr}")
    print("P(z)P(z^-1) has the same triangular Fejer coefficients and double 7th-root zeros")
    print("rank-one block transfer Perron eigenvalue = 7; square mass at theta=0 = 49")
    print("use = model AP as the unique rank-one/equioscillating transfer state before perturbation")

    print("\n[7] Tournament Analysis over proof carriers")
    carriers = [
        "fejer_demoivre_kernel",
        "comb_overlap_gram_kernel",
        "johnson_14_cap_kernel",
        "gamma0_7_eisenstein_sidecar",
        "toeplitz_green_conductance_bridge",
        "subshift_transfer_operator",
        "dirichlet_l_stark_sidecar",
        "beraha_tutte_sidecar",
        "mahler_lehmer_height_sidecar",
        "raw_arithmetic_numerology_guardrail",
    ]
    # Tuple = (certificate payload retained, formal checkability, negative risk).
    carrier_score = {
        "fejer_demoivre_kernel": (10, 9, 0),
        "comb_overlap_gram_kernel": (9, 9, 0),
        "johnson_14_cap_kernel": (9, 8, 0),
        "gamma0_7_eisenstein_sidecar": (7, 6, -2),
        "toeplitz_green_conductance_bridge": (7, 7, -1),
        "subshift_transfer_operator": (5, 5, -2),
        "dirichlet_l_stark_sidecar": (4, 4, -5),
        "beraha_tutte_sidecar": (4, 4, -4),
        "mahler_lehmer_height_sidecar": (3, 4, -4),
        "raw_arithmetic_numerology_guardrail": (1, 1, -8),
    }
    fp = tournament_fingerprint(carriers, carrier_score)
    print(f"vertices = {carriers}")
    print("pairwise observable = certificate payload, then formal checkability, then risk")
    print(f"score_hist = {fp['score_hist']}")
    print(f"directed_3_cycles = {fp['directed_3_cycles']}")
    print("hamiltonian_path = " + " -> ".join(fp["hamiltonian_path"]))
    print("\nAssumption challenge:")
    print("- Rejected vertices: runners, arcs, roots, q-coefficients, and famous constants alone.")
    print("- Chosen vertices: proof carriers, because the LRC predicate is certificate payload retention.")
    print("- Preserves: Fejer positivity/PSD, comb-Gram positivity, double-zero sharpness, cap-vs-sector split, conductor-7 labels.")
    print("- Destroys: literal runner identity, raw height/numerology, and any unlabelled scalar extremality.")

    print("\nNext proof tasks emitted by this audit:")
    print("1. Match Gamma0(7) E_7 coefficients to the S75 comb-overlap Gram basis for finite LP rows.")
    print("2. Prove Fejer/Gamma0/Gram slack dominates green-current trap weights and isolates order-3 debt.")
    print("3. Feed successful rows into HYP-3215 Gap A while keeping the induction-base flag separate.")
    print("4. Treat Dirichlet-L, Beraha, Mahler, and subshift data as sidecars until they preserve an LRC predicate.")


if __name__ == "__main__":
    main()
