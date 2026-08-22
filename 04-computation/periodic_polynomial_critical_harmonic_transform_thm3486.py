#!/usr/bin/env python3
"""Exact structural companion for proved and independently audited THM-3486.

The analytic theorem is proved in the theorem file.  This companion audits
its finite algebraic inputs: the trivial top Fourier coefficient, the
mean-zero top discrepancy, the absolute/conditional boundary, the THM-3484
word, the two THM-3455 Boolean periods, and a declared K4/XOR address period.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from json import dumps
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/periodic_polynomial_critical_harmonic_transform_thm3486.py"
OUTPUT = "05-knowledge/results/periodic_polynomial_critical_harmonic_transform_thm3486.out"
EXPECTED_SEMANTIC_SHA256 = (
    "d30150749708eee29dc4f1773fc1e757e6c8f0df0b7b707b509023e8e06db748"
)

PINS = (
    (
        "THM-3485",
        ROOT / "01-canon/theorems/THM-3485-periodic-polynomial-fourier-jordan-recurrence-classification.md",
        "814ad59925079c9c69c3b90ea1ee7947dbf01ec1c8855bb760fdd6dc5098fafc",
        "PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED",
    ),
    (
        "THM-3455",
        ROOT / "01-canon/theorems/THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum.md",
        "e357aaf8036f485c4abbfa5969a0d1f9761372dcd7337b658daaadf4cd049c09",
        "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED",
    ),
)

THM3484_LANES = (
    (0, 0, 0, 0, -2048, 4096, 8192, -16384),
    (0, 0, 0, 0, 2048, 4096, -24576, -16384),
    (0, -256, 3072, -15360, 40960, -61440, 49152, -16384),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def trim(lane: tuple[int | Fraction, ...]) -> tuple[Fraction, ...]:
    values = [Fraction(value) for value in lane]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def critical_packet(
    name: str, lanes: tuple[tuple[int | Fraction, ...], ...]
) -> tuple[object, ...]:
    normalized = tuple(trim(lane) for lane in lanes)
    degree = max(len(lane) - 1 for lane in normalized)
    leading = tuple(
        lane[degree] if degree < len(lane) else Fraction(0)
        for lane in normalized
    )
    trivial_top = sum(leading, Fraction(0)) / len(leading)
    residual_top = tuple(value - trivial_top for value in leading)
    require(sum(residual_top, Fraction(0)) == 0,
            (name, "top residual is not mean zero"))

    # Sequence indices begin at one, so traverse residues 1,...,p-1,0.
    prefix = Fraction(0)
    prefixes = []
    for index in range(1, len(residual_top) + 1):
        prefix += residual_top[index % len(residual_top)]
        prefixes.append(prefix)
    require(prefix == 0, (name, "one-period discrepancy did not close"))
    discrepancy = max((abs(value) for value in prefixes), default=Fraction(0))
    convergence = "ABSOLUTE_AFTER_LOG" if all(value == 0 for value in residual_top) else "CONDITIONAL_AFTER_LOG"
    if trivial_top == 0:
        convergence = "CONDITIONAL_NO_LOG"
        require(any(value != 0 for value in residual_top),
                (name, "zero word entered nonzero packet census"))
    return (
        name,
        len(normalized),
        degree,
        leading,
        trivial_top,
        residual_top,
        discrepancy,
        convergence,
    )


def abel_identity(values: tuple[Fraction, ...], horizon: int) -> bool:
    period = len(values)
    require(sum(values, Fraction(0)) == 0, "Abel input is not mean zero")
    direct = sum(
        (values[index % period] / index for index in range(1, horizon + 1)),
        Fraction(0),
    )
    prefixes = []
    running = Fraction(0)
    for index in range(1, horizon + 1):
        running += values[index % period]
        prefixes.append(running)
    transformed = prefixes[-1] / horizon
    transformed += sum(
        (
            prefixes[index - 1]
            * (Fraction(1, index) - Fraction(1, index + 1))
            for index in range(1, horizon)
        ),
        Fraction(0),
    )
    return direct == transformed


def boolean_certificate(
    name: str, values: tuple[int, ...]
) -> tuple[object, ...]:
    period = len(values)
    require(set(values) <= {0, 1}, (name, "non-Boolean word"))
    count = sum(values)
    density = Fraction(count, period)
    centered = tuple(Fraction(value) - density for value in values)
    require(sum(centered, Fraction(0)) == 0, (name, "density mismatch"))

    running = Fraction(0)
    discrepancy = Fraction(0)
    for index in range(1, period + 1):
        running += centered[index % period]
        discrepancy = max(discrepancy, abs(running))
    require(running == 0, (name, "period did not close"))
    horizons = (1, period - 1 if period > 1 else 1, period, 2 * period + 7)
    require(all(abel_identity(centered, horizon) for horizon in horizons),
            (name, "Abel identity"))
    return name, period, count, density, discrepancy, horizons


def berggren_spine_word() -> tuple[int, ...]:
    period = 1683
    values = [0] * period
    for residue in range(period):
        t = period if residue == 0 else residue
        values[residue] = int(
            t % 9 in (2, 6)
            or t % 11 in (1, 9)
            or t % 51 in (3, 20, 30, 47)
        )
    require(sum(values) == 684, "THM-3455 full-spine count drift")
    return tuple(values)


def fibonacci_sample_word() -> tuple[int, ...]:
    period = 360
    values = [0] * period
    previous, current = 0, 1
    for index in range(1, period + 1):
        previous, current = current, previous + current
        fibonacci = previous
        values[index % period] = int(
            fibonacci % 9 in (2, 6)
            or fibonacci % 11 in (1, 9)
            or fibonacci % 51 in (3, 20, 30, 47)
        )
    require(sum(values) == 172, "THM-3455 Fibonacci count drift")
    return tuple(values)


def harmonic(number: int) -> Fraction:
    return sum((Fraction(1, index) for index in range(1, number + 1)),
               Fraction(0))


def alternating_harmonic_controls() -> tuple[object, ...]:
    samples = []
    for half_horizon in range(1, 65):
        direct = sum(
            (Fraction((-1) ** index, index)
             for index in range(1, 2 * half_horizon + 1)),
            Fraction(0),
        )
        closed = harmonic(half_horizon) - harmonic(2 * half_horizon)
        require(direct == closed, (half_horizon, direct, closed))
        if half_horizon in (1, 2, 3, 8):
            samples.append((half_horizon, direct))
    return tuple(samples) + (("all_half_horizons_through", 64),)


def security_gate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    forbidden = {
        "assert", "eval", "exec", "compile", "open", "system", "popen",
        "run", "Popen", "unlink", "remove", "rename", "write_text",
        "write_bytes",
    }
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            bad.append("assert")
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id in forbidden:
                bad.append(node.func.id)
            if isinstance(node.func, ast.Attribute) and node.func.attr in forbidden:
                bad.append(node.func.attr)
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), tuple(sorted(set(bad)))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    pins = []
    for label, path, expected_hash, status in PINS:
        actual_hash = lf_sha256(path)
        require(actual_hash == expected_hash, (label, "hash drift", actual_hash))
        require(status in path.read_text(encoding="utf-8"),
                (label, "status drift"))
        pins.append((label, actual_hash))

    packets = (
        critical_packet("THM-3484-ternary-word", THM3484_LANES),
        critical_packet(
            "shared-leading-p3",
            ((1, 2, 5), (-3, 7, 5), (4, -1, 5)),
        ),
        critical_packet(
            "mean-zero-top-p2",
            ((0, 0, 0, 0, 1), (0, 0, 0, 0, -1)),
        ),
        critical_packet(
            "mixed-leading-p3",
            ((2, 1, 4), (0, -3, 7), (5, 2, 1)),
        ),
    )
    require(packets[0][4:8] == (
        Fraction(-16384),
        (Fraction(0), Fraction(0), Fraction(0)),
        Fraction(0),
        "ABSOLUTE_AFTER_LOG",
    ), packets[0])
    require(packets[2][4] == 0 and packets[2][7] == "CONDITIONAL_NO_LOG",
            packets[2])

    k4_star_addresses = tuple(int(index.bit_count() % 2 == 1)
                              for index in range(8))
    booleans = (
        boolean_certificate("berggren-cap7-spine-index", berggren_spine_word()),
        boolean_certificate("berggren-cap7-fibonacci-index", fibonacci_sample_word()),
        boolean_certificate("declared-k4-xor-star-address", k4_star_addresses),
        boolean_certificate("one-residue-mod3", (1, 0, 0)),
    )
    require(booleans[0][3] == Fraction(76, 187), booleans[0])
    require(booleans[1][3] == Fraction(43, 90), booleans[1])
    require(booleans[2][2:4] == (4, Fraction(1, 2)), booleans[2])

    alternating = alternating_harmonic_controls()
    security = security_gate(ROOT / SCRIPT)
    semantic_payload = {
        "pins": tuple(pins),
        "packets": packets,
        "booleans": booleans,
        "alternating": alternating,
        "security": security,
        "theorem_boundary": (
            "trivial-top=log coefficient",
            "nontrivial-top=conditional renormalized channel",
            "shared-top=absolute renormalized channel",
        ),
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"),
              default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("THM-3486 PERIODIC-POLYNOMIAL CRITICAL HARMONIC TRANSFORM COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"PINS: {tuple(pins)}")
    print(f"CRITICAL_PACKETS: {packets}")
    print(f"BOOLEAN_PERIODS: {booleans}")
    print(f"ALTERNATING_HARMONIC_EXACT_CONTROLS: {alternating}")
    print("FOURIER_BOUNDARY: trivial top colour is the unique logarithmic channel; nonconstant mean-zero top data give conditional convergence after log subtraction; shared leading data give absolute convergence after subtraction")
    print("K4_SCOPE: the 4/4 star/triangle row is one declared Boolean-address period, not a tournament orientation or invariant temporal order")
    print("SUBSET_SCOPE: periodic Boolean subsets have a harmonic coefficient; arbitrary subsets of N need not have density and their reciprocal subseries may converge or diverge")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print(f"SECURITY_AST_NODES_AND_FORBIDDEN: {security}")
    print("VERDICT: finite algebraic controls agree with the proved theorem's top-Fourier-layer trichotomy and sharpen the two THM-3455 O(1) harmonic statements to convergent constants")


if __name__ == "__main__":
    main()
