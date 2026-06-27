#!/usr/bin/env python3
"""S162: packet-anchored Fejer interval certificate scaffold for LRC14.

The S157 Fourier-Toeplitz audit found floating Fejer-vector PSD violations for
all positive rows in the HYP-2963 packet bank, with AP/GW as the only zero-safe
atoms.  The remaining proof obligation is to convert those floating values into
rigorous interval certificates anchored to labelled packet fibers P(S).

This script is a finite proof-object scaffold.  For selected hard rows it:

1. imports the existing exact safe-component and labelled-packet classifiers;
2. chooses the S157 Fejer center and degree;
3. recomputes the Fejer quadratic form with rational interval arithmetic;
4. prints the labelled packet fiber key that the certificate should attach to.

The interval backend uses exact Fraction endpoints and a fixed rational
enclosure of pi.  Its output is intended as a checkable certificate format and
as a prototype for the Lean/arb-style backend, not as a final formal proof of
the hard-coded pi enclosure.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from math import factorial
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
N = 14


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s157 = load_module(
    "s162_fejer_float_source",
    REPO / "04-computation" / "lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.py",
)
packets_mod = load_module(
    "s162_packets",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)


def dec_fraction(text: str) -> Fraction:
    if "." not in text:
        return Fraction(int(text), 1)
    whole, frac = text.split(".")
    sign = -1 if whole.startswith("-") else 1
    whole_digits = whole[1:] if sign == -1 else whole
    numerator = int(whole_digits + frac)
    denominator = 10 ** len(frac)
    return sign * Fraction(numerator, denominator)


# pi is between these adjacent 78-decimal rationals.
PI = None
PI_LOW = dec_fraction(
    "3.141592653589793238462643383279502884197169399375105820974944592307816406286"
)
PI_HIGH = dec_fraction(
    "3.141592653589793238462643383279502884197169399375105820974944592307816406287"
)


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError((self.lo, self.hi))

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other: "Interval") -> "Interval":
        return Interval(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __mul__(self, other: "Interval") -> "Interval":
        vals = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(vals), max(vals))

    def scale(self, q: Fraction) -> "Interval":
        if q >= 0:
            return Interval(q * self.lo, q * self.hi)
        return Interval(q * self.hi, q * self.lo)

    def div_pos_int(self, q: int) -> "Interval":
        if q <= 0:
            raise ValueError(q)
        return Interval(self.lo / q, self.hi / q)

    def div_pos_interval(self, other: "Interval") -> "Interval":
        if other.lo <= 0:
            raise ValueError(other)
        return self * Interval(1 / other.hi, 1 / other.lo)

    def pow_nonneg(self, n: int) -> "Interval":
        if n < 0 or self.lo < 0:
            raise ValueError((self, n))
        return Interval(self.lo**n, self.hi**n)

    @property
    def width(self) -> Fraction:
        return self.hi - self.lo

    def contains_negative(self) -> bool:
        return self.hi < 0


ZERO = Interval(Fraction(0), Fraction(0))
ONE = Interval(Fraction(1), Fraction(1))
PI = Interval(PI_LOW, PI_HIGH)


@lru_cache(maxsize=None)
def sin_small_pi(beta: Fraction, terms: int) -> Interval:
    """Enclose sin(pi*beta) for 0 <= beta <= 1/2."""
    if beta == 0:
        return ZERO
    if not (0 <= beta <= Fraction(1, 2)):
        raise ValueError(beta)
    y = PI.scale(beta)
    total = ZERO
    for n in range(terms):
        term = y.pow_nonneg(2 * n + 1).div_pos_int(factorial(2 * n + 1))
        total = total + term if n % 2 == 0 else total - term
    tail = y.pow_nonneg(2 * terms + 1).div_pos_int(factorial(2 * terms + 1))
    return total + Interval(-tail.hi, tail.hi)


@lru_cache(maxsize=None)
def cos_small_pi(beta: Fraction, terms: int) -> Interval:
    """Enclose cos(pi*beta) for 0 <= beta <= 1/2."""
    if beta == 0:
        return ONE
    if not (0 <= beta <= Fraction(1, 2)):
        raise ValueError(beta)
    y = PI.scale(beta)
    total = ZERO
    for n in range(terms):
        term = y.pow_nonneg(2 * n).div_pos_int(factorial(2 * n))
        total = total + term if n % 2 == 0 else total - term
    tail = y.pow_nonneg(2 * terms).div_pos_int(factorial(2 * terms))
    return total + Interval(-tail.hi, tail.hi)


def mod_two(alpha: Fraction) -> Fraction:
    return alpha - 2 * (alpha // 2)


@lru_cache(maxsize=None)
def sin_pi(alpha: Fraction, terms: int) -> Interval:
    r = mod_two(alpha)
    if r == 0 or r == 1:
        return ZERO
    if r <= Fraction(1, 2):
        return sin_small_pi(r, terms)
    if r <= 1:
        return sin_small_pi(1 - r, terms)
    if r <= Fraction(3, 2):
        return -sin_small_pi(r - 1, terms)
    return -sin_small_pi(2 - r, terms)


@lru_cache(maxsize=None)
def cos_pi(alpha: Fraction, terms: int) -> Interval:
    r = mod_two(alpha)
    if r <= Fraction(1, 2):
        return cos_small_pi(r, terms)
    if r <= 1:
        return -cos_small_pi(1 - r, terms)
    if r <= Fraction(3, 2):
        return -cos_small_pi(r - 1, terms)
    return cos_small_pi(2 - r, terms)


@lru_cache(maxsize=None)
def sinc_seventh(m: int, terms: int) -> Interval:
    numerator = sin_pi(Fraction(m, 7), terms)
    denominator = PI.scale(Fraction(m, 1))
    return numerator.div_pos_interval(denominator)


def fejer_interval(
    speeds: tuple[int, ...],
    degree: int,
    center: Fraction,
    terms: int,
) -> Interval:
    total = Interval(Fraction(len(speeds), 7) - 1, Fraction(len(speeds), 7) - 1)
    for k in range(1, degree + 1):
        coeff = ZERO
        for v in speeds:
            if k % v == 0:
                coeff = coeff + sinc_seventh(k // v, terms)
        weight = Fraction(2 * (degree + 1 - k), degree + 1)
        cos_term = cos_pi(2 * k * center, terms)
        total = total + (coeff * cos_term).scale(weight)
    return total


@dataclass(frozen=True)
class RowSpec:
    name: str
    source_family: str
    speeds: tuple[int, ...]
    degree: int


def row_replace(holes: tuple[int, ...], adds: tuple[int, ...], name: str, family: str, degree: int) -> RowSpec:
    return RowSpec(name, family, tuple(sorted((set(AP) - set(holes)) | set(adds))), degree)


def default_rows() -> list[RowSpec]:
    return [
        row_replace((12,), (36,), "near/K33 12->36", "K33 state-lift", 159),
        row_replace((10, 12), (20, 24), "P10+GW", "two-block splice", 280),
        row_replace((12,), (168,), "covering 12->168", "covering comb", 63),
        row_replace((12, 13), (14, 29), "two drop(12, 13)->add(14, 29)", "two-swap", 41),
        row_replace((6,), (63,), "single swap 6->63", "one-swap AP bank", 266),
    ]


def packet_for(row: RowSpec):
    packets = packets_mod.compute_packets([(row.name, row.source_family, row.speeds)], workers=1)
    return packets[0]


def safe_center(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction, Fraction, int]:
    safe_measure, comps = s157.safe_components(speeds)
    if not comps:
        return safe_measure, Fraction(0), Fraction(0), 0
    largest = max(comps, key=lambda comp: comp.width)
    return safe_measure, largest.width, largest.midpoint, len(comps)


def as_decimal(q: Fraction, digits: int = 36) -> str:
    getcontext().prec = digits + 10
    value = Decimal(q.numerator) / Decimal(q.denominator)
    return format(+value, "f")


def interval_line(iv: Interval) -> str:
    return f"[{as_decimal(iv.lo)}, {as_decimal(iv.hi)}]"


def print_header(terms: int) -> None:
    print("S162 LRC14 PACKET-ANCHORED FEJER INTERVAL SCAFFOLD")
    print("=" * 78)
    print("[0] Certificate object")
    print("  Fejer form:")
    print("    Q_d(x)=6/7+2*sum_{1<=k<=d}(1-k/(d+1))*c_k*cos(2*pi*k*x)")
    print("    c_k=sum_{v|k, v in S} sin(pi*(k/v)/7)/(pi*(k/v))")
    print("  interval backend:")
    print(f"    rational pi enclosure width={PI.width}")
    print(f"    Taylor terms per sine/cosine={terms}")
    print("  packet anchor:")
    print("    P(S)=(route, family, q_class, packet_route, state_lift, q_threshold)")
    print()


def print_assumption_challenge() -> None:
    print("[1] Assumption challenge / quotient guardrail")
    print("  considered vertices:")
    print("    runners, divisor fibers v|k, Fourier modes, Fejer centers, endpoint")
    print("    owners, primitive q-phases, labelled packet routes, and proof obligations.")
    print("  chosen vertices:")
    print("    labelled packet fibers P(S) carrying explicit Fejer certificates.")
    print("  preserved LRC predicate:")
    print("    if danger arcs cover then F_S=C_S-1>=0, so every Fejer quadratic")
    print("    form is nonnegative; an interval with upper endpoint <0 certifies")
    print("    a strict safe interval for that packet fiber.")
    print("  destroyed information:")
    print("    raw endpoint order inside the packet; the certificate must record the")
    print("    center, degree, and packet key so this loss is explicit.")
    print()


def print_rows(rows: list[RowSpec], terms: int) -> None:
    print("[2] Packet-row certificates")
    for row in rows:
        safe_mu, width, center, comp_count = safe_center(row.speeds)
        packet = packet_for(row)
        iv = fejer_interval(row.speeds, row.degree, center, terms)
        packet_key = (
            packet.route,
            packet.family,
            packet.q_class,
            packet.packet_route,
            packet.state_lift,
            packet.q_threshold,
        )
        print(f"  row={row.name}")
        print(f"    speeds={row.speeds}")
        print(f"    packet_key={packet_key}")
        print(
            f"    safe_mu={safe_mu} components={comp_count} "
            f"largest_width={width} center={center}"
        )
        print(f"    Fejer degree={row.degree}")
        print(f"    interval={interval_line(iv)}")
        print(f"    interval_width={as_decimal(iv.width, 18)}")
        print(f"    certified_negative={iv.contains_negative()}")
    print()


def print_robin_robbins_synthesis() -> None:
    print("[3] Robin/Robbins synthesis for the proof target")
    print("  Robin's number-theory criterion is a scalar extremal inequality for")
    print("  sigma(n)/n; it is useful here as a warning that divisor averages do not")
    print("  remember phase labels by themselves.")
    print("  Robbins' graph theorem is the better proof analogy: strong orientation")
    print("  is possible exactly when the graph has no bridge.  For LRC14, the")
    print("  certificate graph should be orientable only after packet labels retain")
    print("  the load-bearing bridges: endpoint owner, q-phase, Toeplitz degree,")
    print("  and K33/state-lift debt.")
    print()


def print_next() -> None:
    print("[4] Next proof obligation")
    print("  Promote this scaffold to a production certificate generator by:")
    print("    (a) replacing the hard-coded pi enclosure with a formally cited backend;")
    print("    (b) grouping rows by packet_key and proving shared interval templates;")
    print("    (c) emitting Lean/arb-compatible rational interval certificates;")
    print("    (d) extending from the selected hard rows to all HYP-2963 packet fibers.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--terms", type=int, default=36)
    args = parser.parse_args()
    print_header(args.terms)
    print_assumption_challenge()
    print_rows(default_rows(), args.terms)
    print_robin_robbins_synthesis()
    print_next()


if __name__ == "__main__":
    main()
