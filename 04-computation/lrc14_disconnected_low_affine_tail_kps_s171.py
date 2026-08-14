#!/usr/bin/env python3
"""Build and referee the exact disconnected-low affine-tail compiler.

This wrapper freezes the 79 context/scale lanes not closed by the monotone
midpoint envelope, compiles the independent integer C++ engine, and checks
deterministic physical-mass and homogenized-limit probes against exact Python
implementations.  The analytic input is the averaging lemma documented in
the companion reflection; this script certifies its complete finite residue.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import ceil, floor, gcd
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc14_disconnected_low_affine_tail_kps_s171.cpp"
TAIL = ROOT / "04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py"
MASS = ROOT / "04-computation/lrc_general_reflected_pair_mass_thm3352.py"
EXPECTED_TAIL = "78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9"
EXPECTED_MASS = "afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b"
EXPECTED_SOURCE = "0804e5067a0bbaa3c96f6f3d69915b7baf55a4a8a6b4a35792ec0b4b9ed75e64"
EXPECTED_CONTEXT = "0a36a95b906bf3fb9911d0d234fe9be8032a1d49c090bd21062c4fa9ec74a31f"
EXPECTED_PROBE = "701992b74f33b263545136c87f896f239e7fb7787e18dddec33c2a3ad188a530"
EXPECTED_SUMMARY = "86c21305cc7bf3f86dff2b99d29756855b6cb65bb2c20b1f7148a37d6ac45c58"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def normalized_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gamma(k: int) -> F:
    # Conservative parity-free version, decreasing in real k.
    return F(k * k + 1, 2 * k * k)


def envelope_margin(P: int, g: int, L: int, e: int, f: int) -> F:
    determinant = max(F((e + f) ** 2), F((8 * e + f) ** 2, 8))
    error = (
        gamma(P) * F(e * P, g * L * P - e)
        + gamma(P) * F(f * P, g * L * P - f)
        + determinant / (2 * g * g * L * L)
        + F(e + f, 2 * g * g * L * P)
    )
    return F(1, 49) - F(12, 49 * P * (P + 1)) - error - F(1, 294)


def hostile_contexts(tail):
    all_contexts = set()
    for body, L in tail.SEL.MS.body_universe():
        cell, *_ = tail.SEL.body_geometry(body, L)
        if L < 4592:
            for e in body:
                for f in body:
                    if e != f:
                        all_contexts.add((L, cell, e, f))
    require(len(all_contexts) == 2530, len(all_contexts))
    rows = []
    counts = []
    for L, cell, e, f in sorted(all_contexts):
        mask = 0
        for g in (1, 2, 3):
            P = (264 + g - 1) // g
            if envelope_margin(P, g, L, e, f) <= 0:
                mask |= 1 << (g - 1)
        if mask:
            rows.append((L, cell, e, f, mask))
    for g in (1, 2, 3):
        counts.append(sum(bool(mask & (1 << (g - 1))) for *_, mask in rows))
    require(len(rows) == 79 and tuple(counts) == (79, 10, 4), (len(rows), counts))
    require(tuple(sorted({x[0] for x in rows})) == (168, 336, 420, 504, 560, 588), "hostile rulers")
    payload = "".join(" ".join(map(str, row)) + "\n" for row in rows).encode("ascii")
    digest = sha256(payload).hexdigest()
    if EXPECTED_CONTEXT:
        require(digest == EXPECTED_CONTEXT, ("context digest", digest))
    return tuple(rows), payload, digest


def rays():
    answer = []
    for d in range(1, 9):
        for a in range(0, 7 * d + 1):
            for c in tuple(range(-46, 0)) + tuple(range(1, 47)):
                if (a == 0 and c < 0) or (a == 7 * d and c > 0):
                    continue
                for p0 in range(1, d + 1):
                    if (a * p0 + c) % d == 0:
                        answer.append((d, a, c, p0, p0 + (a * p0 + c) // d))
    require(len(answer) == 22890, len(answer))
    require(len({row[:3] for row in answer}) == 17206, "unique triples")
    return tuple(answer)


def triangle_primitive(R: F, u: F) -> F:
    if u <= -R:
        return F(0)
    if u <= 0:
        return (u + R) ** 2 / 2
    if u < R:
        return R * R / 2 + R * u - u * u / 2
    return R * R


def integrated_tent(R: F, left: F, right: F) -> F:
    require(left <= right, "tent interval")
    lo = floor(-right - R) - 2
    hi = floor(-left + R) + 2
    return sum(
        (triangle_primitive(R, right + n) - triangle_primitive(R, left + n)
         for n in range(lo, hi + 1)),
        F(0),
    )


def tent_value(R: F, x: F) -> F:
    lo = floor(-x - R) - 2
    hi = floor(-x + R) + 2
    return sum((max(F(0), R - abs(x + n)) for n in range(lo, hi + 1)), F(0))


def ray_limit(context, d: int, a: int, c: int) -> F:
    L, cell, e, f = context[:4]
    D = d + a
    k = gcd(d, D)
    P, Q = d // k, D // k
    r, s = e * cell % L, f * cell % L
    A = L * c + D * e - d * f
    origin = F(-D * r + d * s, k * L)
    RA, RB = F(P + Q, 14), F(abs(P - Q), 14)
    if A == 0:
        return (tent_value(RA, origin) - tent_value(RB, origin)) / (P * Q)
    other = origin - F(A, k * L)
    left, right = sorted((origin, other))
    return (integrated_tent(RA, left, right) - integrated_tent(RB, left, right)) / (
        P * Q * abs(right - left)
    )


def cone_start(ray) -> int:
    d, a, _, p0, q0 = ray
    D = d + a
    n = max(0, ceil(F(264 - p0, d)), ceil(F(1 - q0, D)))
    while p0 + d * n >= q0 + D * n:
        n += 1
    while q0 + D * n >= 8 * (p0 + d * n):
        n += 1
    return n


def probe_rows(contexts, ray_rows):
    selected = [((168, 90, 12, 1, 1), (3, 8, -1, 2, 7), 1_000_000)]
    for index in range(98):
        context = contexts[(977 * index + 31) % len(contexts)]
        ray = ray_rows[(65537 * index + 17) % len(ray_rows)]
        selected.append((context, ray, cone_start(ray) + (7919 * index + 13) % 211))
    return tuple(selected)


def compiler() -> str:
    for candidate in tuple(x for x in (os.environ.get("CXX"), "g++", "clang++") if x):
        found = shutil.which(candidate)
        if found:
            return found
    raise RuntimeError("no C++ compiler")


def run_probes(executable: Path, path: Path, probes, mass):
    path.write_text(
        "".join(
            f"{L} {cell} {e} {f} {d} {a} {c} {p0} {q0} {n}\n"
            for (L, cell, e, f, _), (d, a, c, p0, q0), n in probes
        ),
        encoding="ascii", newline="\n",
    )
    completed = subprocess.run((str(executable), "--probe", str(path)), check=True, capture_output=True, text=True)
    output = tuple(tuple(map(int, line.split())) for line in completed.stdout.splitlines())
    require(len(output) == len(probes), len(output))
    semantic = sha256()
    for (context, ray, n), row in zip(probes, output):
        L, cell, e, f, _ = context
        d, a, c, p0, q0 = ray
        D = d + a
        address = (L, cell, e, f, d, a, c, p0, q0, n)
        require(row[:10] == address, (row[:10], address))
        p, q = p0 + d * n, q0 + D * n
        expected_mass = mass.mass(L, cell, e, p, f, q)
        expected_limit = ray_limit(context, d, a, c)
        observed_mass, observed_limit = F(row[10], row[11]), F(row[12], row[13])
        require(observed_mass == expected_mass, ("mass", address, observed_mass, expected_mass))
        require(observed_limit == expected_limit, ("limit", address, observed_limit, expected_limit))
        semantic.update((repr((address, observed_mass, observed_limit)) + "\n").encode())
    require(F(output[0][12], output[0][13]) == F(709, 48048), output[0])
    digest = semantic.hexdigest()
    if EXPECTED_PROBE:
        require(digest == EXPECTED_PROBE, ("probe digest", digest))
    return digest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=max(1, min(os.cpu_count() or 1, 12)))
    parser.add_argument("--keep-summary", type=Path)
    args = parser.parse_args()
    require(args.threads >= 1, args.threads)
    require(normalized_hash(TAIL) == EXPECTED_TAIL, "tail hash")
    require(normalized_hash(MASS) == EXPECTED_MASS, "mass hash")
    source_hash = normalized_hash(SOURCE)
    if EXPECTED_SOURCE:
        require(source_hash == EXPECTED_SOURCE, ("source hash", source_hash))
    tail, mass = load("affine_tail_tail", TAIL), load("affine_tail_mass", MASS)
    contexts, context_payload, context_digest = hostile_contexts(tail)
    ray_rows = rays()

    with tempfile.TemporaryDirectory(prefix="lrc14_affine_tail_s171_") as raw:
        directory = Path(raw)
        executable = directory / ("affine-tail.exe" if os.name == "nt" else "affine-tail")
        context_path, probe_path, summary_path = (directory / name for name in ("contexts.txt", "probes.txt", "summary.txt"))
        context_path.write_bytes(context_payload)
        subprocess.run((compiler(), "-O3", "-DNDEBUG", "-std=c++20", "-pthread", str(SOURCE), "-o", str(executable)), check=True)
        probe_digest = run_probes(executable, probe_path, probe_rows(contexts, ray_rows), mass)
        completed = subprocess.run((str(executable), str(context_path), str(summary_path), str(args.threads)), check=True, capture_output=True, text=True)
        summary = summary_path.read_bytes()
        require(summary.count(b"\n") == 22890, ("summary rows", summary.count(b"\n")))
        summary_digest = sha256(summary).hexdigest()
        if EXPECTED_SUMMARY:
            require(summary_digest == EXPECTED_SUMMARY, ("summary digest", summary_digest))
        if args.keep_summary:
            args.keep_summary.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(summary_path, args.keep_summary)

    print(completed.stdout, end="")
    print(f"context_mask_sha256={context_digest};mask_counts=(79,10,4)")
    print(f"probe_count=99;probe_semantic_sha256={probe_digest}")
    print(f"ray_summary_sha256={summary_digest}")
    print(f"compiler_source_sha256={source_hash}")


if __name__ == "__main__":
    main()
