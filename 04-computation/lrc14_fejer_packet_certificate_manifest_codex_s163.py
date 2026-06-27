#!/usr/bin/env python3
"""S163: manifest for packet-anchored Fejer interval certificates.

S162 showed selected hard LRC14 rows have rational interval Fejer certificates
with upper endpoint below zero.  This follow-up turns those printed rows into a
stable JSONL-style manifest.  Each line records the labelled packet key P(S),
the exact rational center, the Fejer degree, the interval bound, and the proof
fields that a quotient is not allowed to forget.

The manifest is not a new proof of LRC14.  It is a bridge from exploratory
interval arithmetic to a future arb/Lean certificate importer.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import argparse
import json
import sys


REPO = Path(__file__).resolve().parents[1]
S162 = REPO / "04-computation" / "lrc14_packet_fejer_interval_scaffold_codex_s162.py"

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s162 = load_module("s163_s162_scaffold", S162)


def qstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def clipped(text: str, max_chars: int = 180) -> str:
    if len(text) <= max_chars:
        return text
    head = max_chars // 2
    tail = max_chars - head - 18
    return f"{text[:head]}...{text[-tail:]} (chars={len(text)})"


def fraction_digits(q: Fraction) -> tuple[int, int]:
    return len(str(abs(q.numerator))), len(str(q.denominator))


@dataclass(frozen=True)
class CertificateRecord:
    row: str
    speeds: tuple[int, ...]
    source_family: str
    packet_key: tuple[str, str, str, str, str, int]
    safe_mu: str
    component_count: int
    largest_component_width: str
    center: str
    degree: int
    interval_lo: str
    interval_hi: str
    interval_hi_sign: str
    interval_hi_digits: tuple[int, int]
    interval_width: str
    certified_negative: bool
    robbins_bridges: tuple[str, ...]
    quotient_may_forget: tuple[str, ...]
    quotient_must_retain: tuple[str, ...]


def bridge_labels(packet_key: tuple[str, str, str, str, str, int]) -> tuple[str, ...]:
    return (
        "exact_rational_center",
        "Fejer_degree",
        "divisor_curried_atom_formula",
        "trig_interval_enclosure",
        "signed_negative_margin",
        "packet_route",
        "packet_family",
        "q_class",
        "state_lift_debt",
        f"q_threshold={packet_key[-1]}",
    )


def certificate_for(row, terms: int) -> CertificateRecord:
    safe_mu, width, center, comp_count = s162.safe_center(row.speeds)
    packet = s162.packet_for(row)
    packet_key = (
        str(packet.route),
        str(packet.family),
        str(packet.q_class),
        str(packet.packet_route),
        str(packet.state_lift),
        int(packet.q_threshold),
    )
    iv = s162.fejer_interval(row.speeds, row.degree, center, terms)
    return CertificateRecord(
        row=row.name,
        speeds=tuple(row.speeds),
        source_family=row.source_family,
        packet_key=packet_key,
        safe_mu=qstr(safe_mu),
        component_count=comp_count,
        largest_component_width=qstr(width),
        center=qstr(center),
        degree=row.degree,
        interval_lo=qstr(iv.lo),
        interval_hi=qstr(iv.hi),
        interval_hi_sign="negative" if iv.hi < 0 else "nonnegative",
        interval_hi_digits=fraction_digits(iv.hi),
        interval_width=qstr(iv.width),
        certified_negative=iv.hi < 0,
        robbins_bridges=bridge_labels(packet_key),
        quotient_may_forget=(
            "raw_runner_names_after_packet_key_is_fixed",
            "floating_decimal_shadow_after_interval_bound_is_recorded",
        ),
        quotient_must_retain=(
            "packet_key",
            "center",
            "degree",
            "interval_hi",
            "divisor_curried_atom_formula",
            "route_handoff",
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--terms", type=int, default=36)
    parser.add_argument("--format", choices=("text", "jsonl"), default="text")
    args = parser.parse_args()

    records = [certificate_for(row, args.terms) for row in s162.default_rows()]
    failures = [rec.row for rec in records if not rec.certified_negative]

    if args.format == "jsonl":
        for rec in records:
            print(json.dumps(asdict(rec), sort_keys=True))
        return

    print("S163 LRC14 FEJER PACKET CERTIFICATE MANIFEST")
    print("=" * 72)
    print(f"source_scaffold={S162.relative_to(REPO)}")
    print(f"taylor_terms={args.terms}")
    print(f"records={len(records)} certified_negative={len(records)-len(failures)} failures={failures or '-'}")
    print()
    print("Robbins bridge contract:")
    print("  exact center -> Fejer degree -> divisor atom formula -> trig interval")
    print("  -> negative upper bound -> packet fiber -> route handoff")
    print("  A quotient may pass through this manifest only if those bridges are")
    print("  retained or formally reconstructed.")
    print()
    for rec in records:
        print(f"row={rec.row}")
        print(f"  packet_key={rec.packet_key}")
        print(f"  safe_mu={rec.safe_mu} components={rec.component_count} largest_width={rec.largest_component_width}")
        print(f"  center={rec.center} degree={rec.degree}")
        print(f"  interval_hi_sign={rec.interval_hi_sign} digits={rec.interval_hi_digits}")
        print(f"  interval_hi={clipped(rec.interval_hi)}")
        print(f"  certified_negative={rec.certified_negative}")
        print(f"  must_retain={', '.join(rec.quotient_must_retain)}")
    print()
    print("Theorem-facing readout:")
    print("  These are packet-fiber certificates, not scalar row scores.  Robin-style")
    print("  sigma/tau shadows may summarize divisor pressure, but Robbins-style")
    print("  no-bridge assembly says the proof must retain the interval and packet")
    print("  bridges until a formal checker verifies upper(Q)<0.")


if __name__ == "__main__":
    main()
