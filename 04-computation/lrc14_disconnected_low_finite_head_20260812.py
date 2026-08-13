#!/usr/bin/env python3
"""Build and referee the exact disconnected-low finite-head C++ compiler.

The Python layer freezes the canonical context universe, compiles the disjoint
__int128 engine, cross-checks deterministic probes against THM-3352, and hashes
the ordered per-channel semantic summary.  It does not implement an affine
tail and makes no claim outside raw p<264 at primitive common scale g<=3.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc14_disconnected_low_finite_head_20260812.cpp"
TAIL = ROOT / "04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py"
MASS = ROOT / "04-computation/lrc_general_reflected_pair_mass_thm3352.py"
EXPECTED_TAIL = "78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9"
EXPECTED_MASS = "afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b"
EXPECTED_CONTEXT = "efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4"
EXPECTED_CHANNEL = "04b3b49cc072e73aa7472e9dd8573e8240c715d4bc5badb1cbd90fa92b233c28"
EXPECTED_SOURCE = "d7d6e97a760a2c03e124aa24c1f7f3adaec52f19491dd4ce5a5fa54e3e2af151"
EXPECTED_PROBE = "779775924578a481977044abcfd2f859190712585d48adfdc583e22add3c4392"
EXPECTED_SUMMARY = "85ff8655e647585c9113bdf204807d62e1c9e0243e50ec93971b592f9b46949e"
EXPECTED_RULERS = (
    168, 336, 420, 504, 560, 588, 784, 840, 980, 1008, 1176, 1260,
    1680, 1764, 1848, 1960, 2184, 2352, 2520, 2772, 2940, 3080,
    3276, 3528, 3640, 3696, 3920, 4312, 4368,
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def normalized_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module spec", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def contexts(tail):
    answer = set()
    for body, ruler in tail.SEL.MS.body_universe():
        cell, *_ = tail.SEL.body_geometry(body, ruler)
        if ruler < 4592:
            for e in body:
                for f in body:
                    if e != f:
                        answer.add((ruler, cell, e, f))
    answer = tuple(sorted(answer))
    payload = "".join(f"{L} {cell} {e} {f}\n" for L, cell, e, f in answer).encode()
    require(len(answer) == 2530, ("context count", len(answer)))
    require(tuple(sorted({row[0] for row in answer})) == EXPECTED_RULERS, "ruler atlas")
    require(sha256(payload).hexdigest() == EXPECTED_CONTEXT, ("context digest", sha256(payload).hexdigest()))
    return answer, payload


def channels():
    rows = tuple(
        (g, P, Q)
        for g in (1, 2, 3)
        for P in range(1, 264)
        if g * P < 264
        for Q in range(P + 1, 8 * P)
        if gcd(P, Q) == 1 and P + Q >= 8 and (P, Q) != (3, 5)
    )
    counts = tuple((g, sum(row[0] == g for row in rows)) for g in (1, 2, 3))
    require(counts == ((1, 148110), (2, 36978), (3, 16286)), counts)
    payload = "".join(f"{g} {P} {Q}\n" for g, P, Q in rows).encode()
    require(sha256(payload).hexdigest() == EXPECTED_CHANNEL, ("channel digest", sha256(payload).hexdigest()))
    return rows


def compiler() -> str:
    configured = os.environ.get("CXX")
    candidates = tuple(x for x in (configured, "g++", "clang++") if x)
    for candidate in candidates:
        found = shutil.which(candidate)
        if found:
            return found
    raise RuntimeError(("no C++ compiler", candidates))


def probe_rows(context_rows, channel_rows):
    probes = [
        (336, 174, 12, 3, 3, 5),  # 3:5 hostile positive control, 158/46397
        (168, 90, 1, 10, 3, 2),   # swapped low-channel zero control
        (168, 90, 12, 4, 6, 5),   # included candidate-band seed, 92/7645
    ]
    # Deterministic, spread-out finite-head controls.  The two engines share
    # no arithmetic implementation: one is Python Fraction, one __int128.
    for index in range(96):
        L, cell, e, f = context_rows[(index * 977 + 31) % len(context_rows)]
        g, P, Q = channel_rows[(index * 65537 + 17) % len(channel_rows)]
        probes.append((L, cell, e, g * P, f, g * Q))
    return tuple(probes)


def literal_tooth_mass(L, cell, e, p, f, q):
    """Definition-level clipped-interval intersection, independent of floor sums."""

    def intervals(endpoint, level):
        modulus = L * level - endpoint
        residue = endpoint * cell % L
        radius = Fraction(L, 14 * modulus)
        answer = []
        for index in range(-1, level + 1):
            centre = Fraction(residue + index * L, modulus)
            left, right = max(Fraction(0), centre - radius), min(Fraction(1), centre + radius)
            if left < right:
                answer.append((left, right))
        return answer

    first, second = intervals(e, p), intervals(f, q)
    i = j = 0
    answer = Fraction(0)
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            answer += right - left
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return answer


def run_probe(executable: Path, path: Path, probes, mass):
    path.write_text("".join(" ".join(map(str, row)) + "\n" for row in probes), encoding="ascii", newline="\n")
    completed = subprocess.run(
        (str(executable), "--probe", str(path)), check=True, capture_output=True, text=True
    )
    output = tuple(tuple(map(int, line.split())) for line in completed.stdout.splitlines())
    require(len(output) == len(probes), ("probe output count", len(output)))
    semantic = sha256()
    for expected, row in zip(probes, output):
        require(row[:6] == expected, ("probe address", expected, row))
        reference = mass.mass(*expected)
        literal = literal_tooth_mass(*expected)
        observed = Fraction(row[6], row[7])
        require(observed == reference, ("probe mismatch", expected, observed, reference))
        require(observed == literal, ("literal-tooth mismatch", expected, observed, literal))
        semantic.update((repr((expected, observed.numerator, observed.denominator)) + "\n").encode())
    require(Fraction(output[0][6], output[0][7]) == Fraction(158, 46397), output[0])
    require(Fraction(output[1][6], output[1][7]) == 0, output[1])
    require(Fraction(output[2][6], output[2][7]) == Fraction(92, 7645), output[2])
    return semantic.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=max(1, min(os.cpu_count() or 1, 12)))
    parser.add_argument("--keep-summary", type=Path)
    args = parser.parse_args()
    require(args.threads >= 1, ("threads", args.threads))

    require(normalized_hash(TAIL) == EXPECTED_TAIL, ("tail hash", normalized_hash(TAIL)))
    require(normalized_hash(MASS) == EXPECTED_MASS, ("mass hash", normalized_hash(MASS)))
    require(normalized_hash(SOURCE) == EXPECTED_SOURCE, ("source hash", normalized_hash(SOURCE)))
    tail = load("finite_head_tail_20260812", TAIL)
    mass = load("finite_head_mass_20260812", MASS)
    context_rows, context_payload = contexts(tail)
    channel_rows = channels()

    with tempfile.TemporaryDirectory(prefix="lrc14_finite_head_20260812_") as raw_directory:
        directory = Path(raw_directory)
        executable = directory / ("finite_head.exe" if os.name == "nt" else "finite_head")
        context_path = directory / "contexts.txt"
        probe_path = directory / "probes.txt"
        summary_path = directory / "summary.txt"
        context_path.write_bytes(context_payload)
        command = (
            compiler(), "-O3", "-DNDEBUG", "-std=c++20", "-pthread",
            str(SOURCE), "-o", str(executable),
        )
        subprocess.run(command, check=True)
        probe_digest = run_probe(executable, probe_path, probe_rows(context_rows, channel_rows), mass)
        require(probe_digest == EXPECTED_PROBE, ("probe digest", probe_digest))
        completed = subprocess.run(
            (str(executable), str(context_path), str(summary_path), str(args.threads)),
            check=True, capture_output=True, text=True,
        )
        summary = summary_path.read_bytes()
        require(summary.count(b"\n") == len(channel_rows), ("summary rows", summary.count(b"\n")))
        summary_digest = sha256(summary).hexdigest()
        require(summary_digest == EXPECTED_SUMMARY, ("summary digest", summary_digest))
        if args.keep_summary:
            args.keep_summary.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(summary_path, args.keep_summary)

    print(completed.stdout, end="")
    print(f"context_sha256={EXPECTED_CONTEXT}")
    print(f"channel_sha256={EXPECTED_CHANNEL}")
    print(f"probe_count=99;probe_semantic_sha256={probe_digest}")
    print(f"channel_summary_sha256={summary_digest}")
    print(f"compiler_source_sha256={normalized_hash(SOURCE)}")


if __name__ == "__main__":
    main()
