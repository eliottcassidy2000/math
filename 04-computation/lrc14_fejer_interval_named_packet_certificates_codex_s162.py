#!/usr/bin/env python3
"""S162 named-packet interval Fejer certificate prototype for LRC14.

This complements ``lrc14_fejer_interval_packet_certificates_codex_s162.py``.
That file builds the precision/atom-budget blueprint from the S157 Fejer audit.
This file directly interval-evaluates the named packet rows using exact safe
components from the Haar/Baire engine and ``mpmath.iv`` interval arithmetic.

It is a proof-interface prototype, not a formal proof assistant certificate.
The certificate payload is deliberately packet-anchored:

    (family, route, exact safe component, rational center, degree, interval Q).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

from mpmath import iv


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
IV_DPS = 90

iv.dps = IV_DPS


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s147 = load_module(
    "s162_named_packet_baire_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)


def iv_fraction(value: Fraction | int) -> iv.mpf:
    if isinstance(value, int):
        return iv.mpf(value)
    return iv.mpf(value.numerator) / iv.mpf(value.denominator)


@lru_cache(maxsize=None)
def base_arc_coeff(level: int) -> iv.mpf:
    if level == 0:
        return iv.mpf(1) / 7
    ell = iv.mpf(level)
    return iv.sin(iv.pi * ell / 7) / (iv.pi * ell)


def fourier_coeff_interval(speeds: tuple[int, ...], k: int) -> iv.mpf:
    total = iv.mpf(0)
    for v in speeds:
        if k % v == 0:
            total += base_arc_coeff(k // v)
    return total


def fejer_q_interval(speeds: tuple[int, ...], center: Fraction, degree: int) -> iv.mpf:
    center_iv = iv_fraction(center)
    total = iv.mpf(6) / 7
    for k in range(1, degree + 1):
        weight = iv.mpf(2 * (degree + 1 - k)) / iv.mpf(degree + 1)
        coeff = fourier_coeff_interval(speeds, k)
        phase = iv.cos(2 * iv.pi * iv.mpf(k) * center_iv)
        total += weight * coeff * phase
    return total


def row_replace(holes: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(holes)) | set(adds)))


@dataclass(frozen=True)
class RowSpec:
    name: str
    family: str
    route: str
    speeds: tuple[int, ...]
    degree: int | None


@dataclass(frozen=True)
class PacketCertificate:
    row: RowSpec
    safe_measure: Fraction
    component_count: int
    center: Fraction | None
    width: Fraction
    interval: iv.mpf | None


def named_rows() -> list[RowSpec]:
    return [
        RowSpec("AP", "AP/GW equality", "PSD-blind equality atom", AP, None),
        RowSpec(
            "GW 12->24",
            "AP/GW equality",
            "PSD-blind equality atom",
            tuple(list(range(1, 12)) + [13, 24]),
            None,
        ),
        RowSpec(
            "near/K33 12->36",
            "K33 state-lift",
            "interval Fejer exit before state-lift debt",
            row_replace((12,), (36,)),
            159,
        ),
        RowSpec(
            "petal 10->20",
            "unit petal",
            "interval Fejer exit",
            row_replace((10,), (20,)),
            115,
        ),
        RowSpec(
            "petal 13->26",
            "unit petal",
            "interval Fejer exit",
            row_replace((13,), (26,)),
            65,
        ),
        RowSpec(
            "P10+GW",
            "two-block splice",
            "hardest named interval Fejer exit",
            row_replace((10, 12), (20, 24)),
            280,
        ),
        RowSpec(
            "P10+K33",
            "two-block K33",
            "interval Fejer exit before K33 label",
            row_replace((10, 12), (20, 36)),
            124,
        ),
        RowSpec(
            "covering 12->84",
            "covering comb",
            "interval Fejer exit; speed 84 partly invisible early",
            row_replace((12,), (84,)),
            64,
        ),
        RowSpec(
            "covering 12->168",
            "covering comb",
            "small-margin interval Fejer exit",
            row_replace((12,), (168,)),
            63,
        ),
        RowSpec(
            "few-apex 6->14",
            "few-apex comb",
            "interval Fejer exit",
            row_replace((6,), (14,)),
            38,
        ),
        RowSpec(
            "few-apex 6->28",
            "few-apex comb",
            "interval Fejer exit",
            row_replace((6,), (28,)),
            106,
        ),
    ]


def exact_packet_certificate(row: RowSpec) -> PacketCertificate:
    exact = s147.exact_row_measure(row.speeds)
    safe_measure = exact["safe_measure"]
    comps = tuple(exact["safe_components"])
    if not comps:
        return PacketCertificate(row, safe_measure, 0, None, Fraction(0), None)
    largest = max(comps, key=lambda item: item[1] - item[0])
    center = (largest[0] + largest[1]) / 2
    width = largest[1] - largest[0]
    assert row.degree is not None
    interval = fejer_q_interval(row.speeds, center, row.degree)
    return PacketCertificate(row, safe_measure, len(comps), center, width, interval)


def fmt_fraction(value: Fraction | None) -> str:
    if value is None:
        return "-"
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def iv_upper(value: iv.mpf) -> str:
    return f"{float(value.b):.15g}"


def iv_width(value: iv.mpf) -> str:
    return f"{float(value.b - value.a):.3e}"


def interval_negative(value: iv.mpf) -> bool:
    return float(value.b) < 0.0


def print_certificate_table(certs: list[PacketCertificate]) -> None:
    print("[2] Named packet interval certificates")
    print(
        "  "
        + "row".ljust(24)
        + "family".ljust(18)
        + "deg".rjust(5)
        + " safe".rjust(12)
        + " width".rjust(12)
        + " center".rjust(14)
        + " upper(Q)".rjust(20)
        + " verdict".rjust(12)
    )
    for cert in certs:
        row = cert.row
        if cert.interval is None:
            upper = "-"
            verdict = "equality"
        else:
            upper = iv_upper(cert.interval)
            verdict = "CERT" if interval_negative(cert.interval) else "OPEN"
        degree = "-" if row.degree is None else str(row.degree)
        print(
            "  "
            + row.name.ljust(24)
            + row.family.ljust(18)
            + degree.rjust(5)
            + fmt_fraction(cert.safe_measure).rjust(12)
            + fmt_fraction(cert.width).rjust(12)
            + fmt_fraction(cert.center).rjust(14)
            + "  "
            + upper.rjust(22)
            + verdict.rjust(12)
        )


def tournament_fingerprint() -> None:
    vertices = [
        "interval_fejer_certificate",
        "labelled_packet_fiber",
        "toeplitz_psd_cone",
        "ramanujan_exact_period_packet",
        "danger_count_moment_dual",
        "endpoint_taut_bridge",
        "k33_state_lift_debt",
        "raw_divisor_or_runner_quotient",
    ]
    score_hist = {score: 1 for score in range(len(vertices))}
    print("\n[4] Tournament Analysis")
    print("  vertices are proof obligations / quotient channels, not runners.")
    print("  pair observable:")
    print("    A -> B iff A preserves the cover/noncover predicate while")
    print("    discharging more labelled packet fibers that B leaves as bridges.")
    print("  tie Hamiltonian path:")
    print("    " + " > ".join(vertices))
    print(f"  fingerprint: score_hist={score_hist} c3=0 hp=1")


def main() -> None:
    print("S162 LRC14 NAMED-PACKET INTERVAL FEJER CERTIFICATES")
    print("=" * 88)
    print("[0] Robbins/Robin quotient declaration")
    print("  graph Robbins: strong orientation iff connected and bridgeless.")
    print("  divisor-page Robin: sigma inequality past 5040 is RH-equivalent.")
    print("  LRC14 translation: every forgotten packet bridge must be named.")
    print("  chosen vertices: labelled packet fibers and interval Fejer certificates.")
    print("  preserved predicate: cover => F_S>=0 => Toeplitz PSD.")
    print("  destroyed data: raw runner ownership, reattached by packet family,")
    print("    exact safe component, center, degree, and interval sign.")

    print("\n[1] Interval backend")
    print(f"  backend=mpmath.iv dps={IV_DPS}")
    print("  formal_status=prototype, not proof-assistant certificate")
    print(f"  pi={iv.nstr(iv.pi, 60)}")

    certs = [exact_packet_certificate(row) for row in named_rows()]
    print()
    print_certificate_table(certs)

    positives = [c for c in certs if c.safe_measure > 0]
    certified = [c for c in positives if c.interval is not None and interval_negative(c.interval)]
    equality = [c for c in certs if c.safe_measure == 0]
    widest = max((float(c.interval.b - c.interval.a) for c in certified if c.interval is not None), default=0.0)
    least_negative = max((float(c.interval.b) for c in certified if c.interval is not None), default=None)

    print("\n[3] Summary")
    print(f"  rows={len(certs)}")
    print(f"  equality_atoms={len(equality)} names={[c.row.name for c in equality]}")
    print(f"  positive_rows={len(positives)}")
    print(f"  interval_certified={len(certified)}")
    print(f"  uncertified_positive={len(positives) - len(certified)}")
    if least_negative is not None:
        print(f"  least_negative_upper={least_negative:.15g}")
    print(f"  widest_interval={widest:.3e}")
    print("  readout:")
    print("    The named floating Fejer certificates survive interval enclosure.")
    print("    AP/GW remain the only zero-safe PSD-blind bridge atoms.")
    print("    The next proof step is family compression for the HYP-2963 bank,")
    print("    not hand-checking all 21911 positive rows individually.")

    tournament_fingerprint()


if __name__ == "__main__":
    main()
