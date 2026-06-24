#!/usr/bin/env python3
"""HYP-2977: spectral-shadow dual probe for the LRC14 proof.

This is deliberately a different route from the endpoint-owner / aperture
packet work.  Given an exact strict-safe open set

    U(S) = {t in R/Z : ||v t|| > 1/14 for every v in S},

the script treats 1_U as the primary object and asks which Fourier bands see
its mass.  The proof target is not a local pinch.  It is a finite
trigonometric/dual certificate:

    positive low-band shadow or controlled high-frequency tail
        -> U(S) has positive Haar mass.

Tournament Analysis is over frequency bands.  The pairwise observable is which
band captures more Parseval energy of 1_U across the audited positive rows.
This quotient preserves Haar mass and spectral detectability while destroying
endpoint-owner labels, which is the point of this session's alternative angle.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import cmath
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
TAU = 2.0 * 3.141592653589793238462643383279502884


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s146 = load_module(
    "spectral_shadow_s146",
    REPO / "04-computation" / "lrc14_haar_baire_taut_boundary_s146.py",
)
s152 = load_module(
    "spectral_shadow_s152",
    REPO / "04-computation" / "lrc14_few_apex_lift_packet_probe_codex_s152.py",
)
s124 = load_module(
    "spectral_shadow_s124",
    REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py",
)


@dataclass(frozen=True)
class Row:
    name: str
    source: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class SpectralAudit:
    row: Row
    components: tuple[tuple[Fraction, Fraction], ...]
    mass: Fraction
    max_component: Fraction
    coeffs: tuple[complex, ...]
    captured_by_h: dict[int, float]
    first_capture: dict[str, int | None]
    band_energy: tuple[float, ...]
    dominant: tuple[tuple[int, float], ...]
    fejer_mid: dict[int, float]


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def qdiv(speeds: tuple[int, ...], cap: int = 240) -> int:
    for d in range(2, cap + 1):
        if all(v % d for v in speeds):
            return d
    return cap + 1


def replace_one(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - {drop}) | {add}))


def replace_many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(drops)) | set(adds)))


def named_rows() -> list[Row]:
    return [
        Row("AP", "boundary", AP, "strict-Haar-zero equality atom"),
        Row("GW 12->24", "boundary", replace_one(12, 24), "strict-Haar-zero equality atom"),
        Row("near 12->36", "named", replace_one(12, 36), "K33 near-miss"),
        Row("petal 10->20", "named", replace_one(10, 20), "unit-petal row"),
        Row("petal 13->26", "named", replace_one(13, 26), "unit-petal row"),
        Row(
            "two-swap 10,12->20,24",
            "named",
            replace_many((10, 12), (20, 24)),
            "P10 plus GW splice",
        ),
        Row(
            "two-swap 10,12->20,36",
            "named",
            replace_many((10, 12), (20, 36)),
            "P10 plus K33 splice",
        ),
        Row("covering 6->98", "covering", replace_one(6, 98), "HYP-2965 smallest covering mass"),
        Row("covering 12->84", "covering", replace_one(12, 84), "first-order apex blocked row"),
        Row("covering 12->168", "covering", replace_one(12, 168), "same lcm-tail packet as 12->84"),
    ]


def few_apex_rows(max_k14: int, max_multiplier: int, count: int) -> list[Row]:
    if count <= 0:
        return []
    bank = s152.build_bank(max_k14, max_multiplier, include_stretch=True)
    audits = [s152.lift_packet(row) for row in bank]
    audits.sort(key=lambda a: (a.safe_t_measure, a.open_lifts, max(a.row.speeds), a.row.speeds))
    out: list[Row] = []
    seen = {tuple(r.speeds) for r in named_rows()}
    for audit in audits:
        speeds = tuple(audit.row.speeds)
        if speeds in seen:
            continue
        seen.add(speeds)
        out.append(
            Row(
                audit.row.name,
                "few-apex-worst",
                speeds,
                f"k14={audit.k14}, lift_mass={fmt(audit.safe_t_measure)}, open_lifts={audit.open_lifts}",
            )
        )
        if len(out) >= count:
            break
    return out


def interval_coeff(intervals: tuple[tuple[Fraction, Fraction], ...], n: int) -> complex:
    total = 0j
    denom = 1j * TAU * n
    for a, b in intervals:
        af = float(a)
        bf = float(b)
        total += (cmath.exp(-1j * TAU * n * af) - cmath.exp(-1j * TAU * n * bf)) / denom
    return total


def fourier_coeffs(intervals: tuple[tuple[Fraction, Fraction], ...], max_h: int) -> tuple[complex, ...]:
    # coeffs[n] = integral_U exp(-2*pi*i*n*t) dt.  coeffs[0] is stored too.
    mass = float(sum((b - a for a, b in intervals), Fraction(0)))
    coeffs = [complex(mass, 0.0)]
    for n in range(1, max_h + 1):
        coeffs.append(interval_coeff(intervals, n))
    return tuple(coeffs)


def captured_ratio(coeffs: tuple[complex, ...], mass: Fraction, h: int) -> float:
    mu = float(mass)
    if mu <= 0.0:
        return 0.0
    energy = abs(coeffs[0]) ** 2
    energy += 2.0 * sum(abs(coeffs[n]) ** 2 for n in range(1, h + 1))
    return min(energy / mu, 1.0)


def first_capture_levels(coeffs: tuple[complex, ...], mass: Fraction, max_h: int) -> dict[str, int | None]:
    levels = (("50%", 0.50), ("75%", 0.75), ("90%", 0.90), ("99%", 0.99))
    out: dict[str, int | None] = {name: None for name, _ in levels}
    for h in range(1, max_h + 1):
        r = captured_ratio(coeffs, mass, h)
        for name, target in levels:
            if out[name] is None and r >= target:
                out[name] = h
    return out


def band_ranges(max_h: int) -> tuple[tuple[int, int], ...]:
    raw = ((1, 7), (8, 14), (15, 28), (29, 56), (57, 112), (113, 224), (225, 448))
    return tuple((a, min(b, max_h)) for a, b in raw if a <= max_h)


def band_energy(coeffs: tuple[complex, ...], bands: tuple[tuple[int, int], ...]) -> tuple[float, ...]:
    out = []
    for lo, hi in bands:
        out.append(2.0 * sum(abs(coeffs[n]) ** 2 for n in range(lo, hi + 1)))
    return tuple(out)


def dominant_modes(coeffs: tuple[complex, ...], count: int) -> tuple[tuple[int, float], ...]:
    vals = [(n, abs(coeffs[n])) for n in range(1, len(coeffs))]
    vals.sort(key=lambda item: (-item[1], item[0]))
    return tuple(vals[:count])


def largest_component_midpoint(components: tuple[tuple[Fraction, Fraction], ...]) -> Fraction | None:
    if not components:
        return None
    a, b = max(components, key=lambda iv: (iv[1] - iv[0], -iv[0]))
    return (a + b) / 2


def fejer_value(coeffs: tuple[complex, ...], h: int, x: Fraction) -> float:
    xf = float(x)
    value = coeffs[0].real
    for n in range(1, h + 1):
        weight = 1.0 - n / (h + 1.0)
        value += 2.0 * weight * (coeffs[n] * cmath.exp(1j * TAU * n * xf)).real
    return value


def audit_row(row: Row, max_h: int, capture_h: tuple[int, ...], bands: tuple[tuple[int, int], ...]) -> SpectralAudit:
    comps = tuple(s146.safe_open_components(row.speeds))
    mass = s146.interval_measure(list(comps))
    max_component = max((b - a for a, b in comps), default=Fraction(0))
    coeffs = fourier_coeffs(comps, max_h)
    captured = {h: captured_ratio(coeffs, mass, min(h, max_h)) for h in capture_h if h <= max_h}
    first = first_capture_levels(coeffs, mass, max_h)
    energies = band_energy(coeffs, bands)
    dominant = dominant_modes(coeffs, 6)
    mid = largest_component_midpoint(comps)
    fejer = {}
    for h in capture_h:
        if h <= max_h:
            fejer[h] = 0.0 if mid is None else fejer_value(coeffs, h, mid)
    return SpectralAudit(
        row=row,
        components=comps,
        mass=mass,
        max_component=max_component,
        coeffs=coeffs,
        captured_by_h=captured,
        first_capture=first,
        band_energy=energies,
        dominant=dominant,
        fejer_mid=fejer,
    )


def edge_from_adj(adj: list[list[bool]], i: int, j: int) -> bool:
    return adj[i][j]


def tournament_fingerprint(adj: list[list[bool]]) -> dict[str, object]:
    n = len(adj)
    outdeg = [sum(1 for j in range(n) if i != j and adj[i][j]) for i in range(n)]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            c3 += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            c3 += 1

    graph = [[j for j in range(n) if i != j and adj[i][j]] for i in range(n)]
    rgraph = [[j for j in range(n) if i != j and adj[j][i]] for i in range(n)]
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    seen.clear()
    comps = []
    for v in reversed(order):
        if v in seen:
            continue
        comp = set()
        stack = [v]
        seen.add(v)
        while stack:
            x = stack.pop()
            comp.add(x)
            for w in rgraph[x]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(comp)

    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for _ in range(1, n):
        nxt: dict[tuple[int, int], int] = defaultdict(int)
        for (used, last), val in dp.items():
            for v in range(n):
                if used & (1 << v):
                    continue
                if adj[last][v]:
                    nxt[(used | (1 << v), v)] += val
        dp = nxt
    full = (1 << n) - 1
    hp = sum(val for (used, _), val in dp.items() if used == full)
    return {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "c3": c3,
        "scc": tuple(sorted((len(c) for c in comps), reverse=True)),
        "hp": hp,
    }


def band_tournament(
    audits: list[SpectralAudit], bands: tuple[tuple[int, int], ...]
) -> tuple[list[list[bool]], dict[tuple[int, int], tuple[int, int]]]:
    positive = [a for a in audits if a.mass > 0]
    n = len(bands)
    adj = [[False] * n for _ in range(n)]
    ledger: dict[tuple[int, int], tuple[int, int]] = {}
    for i, j in combinations(range(n), 2):
        iwins = 0
        jwins = 0
        for audit in positive:
            ei = audit.band_energy[i]
            ej = audit.band_energy[j]
            if ei > ej:
                iwins += 1
            elif ej > ei:
                jwins += 1
        ledger[(i, j)] = (iwins, jwins)
        if iwins > jwins or (iwins == jwins and i < j):
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj, ledger


def print_assumption_challenge(bands: tuple[tuple[int, int], ...]) -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, gaps, fixed sections, section boundaries, wall-crossing events,")
    print("    residues, cover arcs, lift packets, Fourier modes, and proof obligations.")
    print("  chosen vertices:")
    print("    Fourier bands " + ", ".join(f"{lo}-{hi}" for lo, hi in bands) + ".")
    print("  LRC predicate preserved:")
    print("    positive strict Haar mass is the L2 mass of 1_U, and a positive")
    print("    Fejer convolution against 1_U is a spectral witness of nonempty U.")
    print("  information destroyed:")
    print("    endpoint-owner labels, C27 shell addresses, and exact packet family.")
    print("  challenged assumption:")
    print("    that the proof must locate a local boundary pinch.  This route asks")
    print("    whether every live packet has a finite trigonometric shadow, or else")
    print("    a high-frequency tail controlled by relation-lattice structure.")


def print_audit_table(audits: list[SpectralAudit], capture_h: tuple[int, ...]) -> None:
    print()
    print("[1] Exact safe set and Fourier visibility")
    header = (
        f"  {'row':34s} {'src':15s} {'qdiv':>5s} {'mass':>14s} {'max_comp':>12s} "
        + " ".join(f"E<={h:>3d}" for h in capture_h)
        + "  H90"
    )
    print(header)
    for audit in audits:
        caps = " ".join(f"{audit.captured_by_h.get(h, 0.0):6.3f}" for h in capture_h)
        h90 = audit.first_capture["90%"]
        h90_text = "-" if h90 is None else str(h90)
        print(
            f"  {audit.row.name[:34]:34s} {audit.row.source[:15]:15s} "
            f"{qdiv(audit.row.speeds):5d} {fmt(audit.mass):>14s} {fmt(audit.max_component):>12s} "
            f"{caps} {h90_text:>4s}"
        )


def print_dominant_modes(audits: list[SpectralAudit]) -> None:
    print()
    print("[2] Dominant Fourier modes")
    for audit in audits:
        if audit.mass == 0:
            print(f"  {audit.row.name}: strict-safe mass zero; spectral shadow is zero.")
            continue
        modes = ", ".join(f"{n}:{amp:.5g}" for n, amp in audit.dominant[:5])
        print(f"  {audit.row.name}: {modes}")


def print_fejer_table(audits: list[SpectralAudit], capture_h: tuple[int, ...]) -> None:
    print()
    print("[3] Fejer midpoint shadow")
    print("  Values are Fejer_H * 1_U at the midpoint of the largest exact safe component.")
    print("  Positive values are spectral witnesses for the already exact open set;")
    print("  the proof target is to lower-bound one of these without enumerating endpoints.")
    header = f"  {'row':34s} " + " ".join(f"F{h:>3d}" for h in capture_h)
    print(header)
    for audit in audits:
        vals = " ".join(f"{audit.fejer_mid.get(h, 0.0):8.5f}" for h in capture_h)
        print(f"  {audit.row.name[:34]:34s} {vals}")


def print_band_tournament(audits: list[SpectralAudit], bands: tuple[tuple[int, int], ...]) -> None:
    print()
    print("[4] Tournament Analysis: frequency-band majority relation")
    adj, ledger = band_tournament(audits, bands)
    fp = tournament_fingerprint(adj)
    names = [f"{lo}-{hi}" for lo, hi in bands]
    outdeg = [sum(1 for j in range(len(bands)) if i != j and adj[i][j]) for i in range(len(bands))]
    order = sorted(range(len(bands)), key=lambda i: (-outdeg[i], i))
    print("  vertices:")
    print("    " + ", ".join(names))
    print("  pair observable:")
    print("    band A -> band B iff A has larger Parseval band energy than B on")
    print("    more audited positive rows; ties follow low-frequency Hamiltonian path.")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  Hamiltonian path by outdegree:")
    print("    " + " > ".join(names[i] for i in order))
    print("  close pair ledgers:")
    for (i, j), (iwins, jwins) in ledger.items():
        if abs(iwins - jwins) <= 2:
            print(f"    {names[i]} vs {names[j]}: {iwins}-{jwins}")


def print_theorem_readout(audits: list[SpectralAudit], max_h: int) -> None:
    print()
    print("[5] Proof readout")
    zero = [a.row.name for a in audits if a.mass == 0]
    positive = [a for a in audits if a.mass > 0]
    unresolved_90 = [a.row.name for a in positive if a.first_capture["90%"] is None]
    min_mass = min((a.mass for a in positive), default=Fraction(0))
    min_fejer14 = min((a.fejer_mid.get(14, 0.0) for a in positive), default=0.0)
    print(f"  zero strict-mass rows in this audit: {zero}")
    print(f"  positive rows audited: {len(positive)}")
    print(f"  smallest exact positive mass: {fmt(min_mass)}")
    print(f"  smallest Fejer_14 midpoint value among positive rows: {min_fejer14:.8f}")
    if unresolved_90:
        print(f"  rows not 90%-captured by H={max_h}: {unresolved_90[:8]}{' ...' if len(unresolved_90) > 8 else ''}")
    else:
        print(f"  every positive row reached 90% Parseval capture by H={max_h}.")
    print("  candidate theorem form:")
    print("    In the post-THM-571 Moon core, either a bounded frequency packet")
    print("    gives a positive Fejer/Beurling-Selberg minorant, or the energy is")
    print("    forced into high bands whose relation-lattice tail is one of the")
    print("    already labelled lift/K33/covering packets.  AP/GW are exactly the")
    print("    zero-shadow atoms.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-frequency", type=int, default=224)
    parser.add_argument("--few-apex-worst", type=int, default=18)
    parser.add_argument("--max-k14", type=int, default=6)
    parser.add_argument("--max-multiplier", type=int, default=180)
    args = parser.parse_args()

    max_h = args.max_frequency
    captures = tuple(h for h in (14, 28, 56, 112, 224, 448) if h <= max_h)
    bands = band_ranges(max_h)

    rows = named_rows() + few_apex_rows(args.max_k14, args.max_multiplier, args.few_apex_worst)
    seen: set[tuple[int, ...]] = set()
    unique_rows: list[Row] = []
    for row in rows:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        unique_rows.append(row)

    print("HYP-2977 / T1061 spectral-shadow dual probe for LRC14")
    print(
        f"parameters: max_frequency={max_h} few_apex_worst={args.few_apex_worst} "
        f"max_k14={args.max_k14} max_multiplier={args.max_multiplier}"
    )
    print_assumption_challenge(bands)

    audits = [audit_row(row, max_h, captures, bands) for row in unique_rows]
    print_audit_table(audits, captures)
    print_dominant_modes(audits)
    print_fejer_table(audits, captures)
    print_band_tournament(audits, bands)
    print_theorem_readout(audits, max_h)


if __name__ == "__main__":
    main()
