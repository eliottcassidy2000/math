#!/usr/bin/env python3
"""Fourier-Toeplitz dual scout for LRC14.

This is deliberately different from the labelled-packet / boundary-moment
route.  It starts from the covering formulation itself:

    D_S(t) = sum_{v in S} 1_{||v t|| < 1/14}.

If S were a strict LRC14 counterexample, then D_S(t) >= 1 almost everywhere,
so F_S(t)=D_S(t)-1 would be a nonnegative integrable function.  Therefore every
Toeplitz moment matrix

    T_K(S) = (hat F_S(i-j))_{0 <= i,j <= K}

must be positive semidefinite.  A negative eigenvalue is a Farkas/Fourier dual
certificate that the danger arcs do not cover the circle, hence S has a safe
open interval.

The output is a scouting ledger: which hard rows are discharged by low-degree
Toeplitz negativity, which rows remain invisible to this dual lens, and which
harmonic bands dominate the negative eigenvectors.

Tournament Analysis uses harmonic proof lenses as vertices, not runners,
packets, or cover arcs.  The pairwise observable is preservation of the
nonnegative-multiplicity obstruction before scalarizing to exact safe mass.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, pi, sin
from pathlib import Path
import argparse
import sys

import numpy as np


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
U14 = frozenset({1, 3, 5, 9, 11, 13})


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s150 = load_module(
    "ft_s150_packet_migration",
    REPO / "04-computation" / "lrc14_packet_migration_gauntlet_codex_s150.py",
)


@dataclass(frozen=True)
class Row:
    name: str
    family: str
    speeds: tuple[int, ...]
    qdiv: int


@dataclass(frozen=True)
class ExactHint:
    status: str
    mass: Fraction
    atom_keys: tuple[str, ...]
    transfer: str


@dataclass(frozen=True)
class DualResult:
    row: Row
    exact: ExactHint
    first_negative_degree: int | None
    min_eig: float
    min_eig_degree: int
    dominant_band: str
    band_energy: tuple[tuple[str, float], ...]


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def row_tuple(items: list[int] | tuple[int, ...]) -> tuple[int, ...]:
    out = tuple(sorted(set(items)))
    assert len(out) == 13 and primitive(out)
    return out


def curated_rows() -> list[Row]:
    rows = [
        Row("AP", "tight-boundary", AP, s150.qdiv(AP)),
        Row("GW 12->24", "tight-boundary", row_tuple(list(range(1, 12)) + [13, 24]), 14),
        Row("K33 near 12->36", "K33-state-lift", row_tuple(list(range(1, 12)) + [13, 36]), 14),
        Row("unit petal 10->20", "unit-petal", row_tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 20]), 14),
        Row("unit petal 13->26", "unit-petal", row_tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 26]), 14),
        Row("two-hole P10+GW", "two-hole-splice", row_tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 24]), 14),
        Row("two-hole P10+K33", "two-hole-splice", row_tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 36]), 14),
        Row("covering comb 12->84", "covering-tail", row_tuple(list(range(1, 12)) + [13, 84]), 15),
        Row("covering comb 12->168", "covering-tail", row_tuple(list(range(1, 12)) + [13, 168]), 15),
        Row("magnitude liar 12->96", "magnitude-liar", row_tuple(list(range(1, 12)) + [13, 96]), 15),
        Row("boundary-gap small 6->98", "boundary-gap", row_tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 98]), 15),
        Row("few-apex lift 6->210", "few-apex-lift", row_tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 210]), 15),
        Row("few-apex lift 5->196", "few-apex-lift", row_tuple([1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 196]), 15),
    ]
    # Recompute qdiv for every curated row to avoid trusting the hand labels.
    return [Row(r.name, r.family, r.speeds, s150.qdiv(r.speeds)) for r in rows]


def generated_rows(single_limit: int, two_swap_limit: int, max_generated: int) -> list[Row]:
    selected: list[Row] = []
    for k, add_max in ((1, single_limit), (2, two_swap_limit)):
        specs, _ = s150.generate_bank(k, add_max)
        hard = [s for s in specs if s.qdiv >= 14]
        # Prefer rows in the live covering zone, then deterministic spread over names.
        hard.sort(key=lambda s: (0 if s.qdiv > 14 else 1, s.qdiv, s.name))
        stride = max(1, len(hard) // max(1, max_generated // 2))
        for spec in hard[::stride]:
            selected.append(Row(f"{k}-swap {spec.name}", f"{k}-swap-qdiv>=14", spec.speeds, spec.qdiv))
            if len(selected) >= max_generated:
                return selected
    return selected


def unique_rows(rows: list[Row]) -> list[Row]:
    seen: set[tuple[int, ...]] = set()
    out: list[Row] = []
    for row in rows:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        out.append(row)
    return out


def exact_hint(row: Row) -> ExactHint:
    holes = tuple(v for v in AP if v not in row.speeds)
    adds = tuple(v for v in row.speeds if v not in AP)
    try:
        exact = s150.classify_exact(
            s150.RowSpec("dual-scout", row.name, row.speeds, holes, adds, row.qdiv, "qdiv>=14 exact")
        )
        return ExactHint(exact.status, exact.mass, exact.atom_keys, exact.transfer)
    except Exception:
        return ExactHint("not-classified", Fraction(-1), (), "-")


def arc_coeff(k: int) -> float:
    """Fourier coefficient of 1_{||t||<1/14} at nonzero frequency k."""
    return sin(pi * k / 7.0) / (pi * k)


def fhat_cache(speeds: tuple[int, ...], max_m: int) -> np.ndarray:
    vals = np.zeros(max_m + 1, dtype=float)
    vals[0] = len(speeds) / 7.0 - 1.0
    for m in range(1, max_m + 1):
        total = 0.0
        for v in speeds:
            if m % v == 0:
                total += arc_coeff(m // v)
        vals[m] = total
    return vals


def toeplitz(vals: np.ndarray, degree: int) -> np.ndarray:
    idx = np.arange(degree + 1)
    return vals[np.abs(idx[:, None] - idx[None, :])]


def band_for_index(i: int) -> str:
    if i == 0:
        return "constant"
    if i % 14 == 0:
        return "apex-14-multiple"
    if i % 27 == 0:
        return "C27-multiple"
    if i % 14 in U14:
        return "unit-apex-residue"
    if i <= 13:
        return "low-AP-band"
    if i <= 54:
        return "mid-transfer-band"
    return "high-tail-band"


def band_energy(vec: np.ndarray) -> tuple[tuple[str, float], ...]:
    acc: Counter[str] = Counter()
    total = float(np.vdot(vec, vec).real)
    if total <= 0:
        return ()
    for i, z in enumerate(vec):
        acc[band_for_index(i)] += float(abs(z) ** 2) / total
    return tuple(sorted(acc.items(), key=lambda kv: (-kv[1], kv[0])))


def scan_row(row: Row, degrees: tuple[int, ...], tol: float) -> DualResult:
    vals = fhat_cache(row.speeds, max(degrees))
    first_negative: int | None = None
    best_eig = float("inf")
    best_degree = -1
    best_vec: np.ndarray | None = None
    for degree in degrees:
        mat = toeplitz(vals, degree)
        eigvals, eigvecs = np.linalg.eigh(mat)
        e0 = float(eigvals[0])
        if e0 < best_eig:
            best_eig = e0
            best_degree = degree
            best_vec = eigvecs[:, 0]
        if first_negative is None and e0 < -tol:
            first_negative = degree
    energies = band_energy(best_vec) if best_vec is not None else ()
    dominant = energies[0][0] if energies else "-"
    return DualResult(row, exact_hint(row), first_negative, best_eig, best_degree, dominant, energies)


def fmt_frac(frac: Fraction) -> str:
    if frac < 0:
        return "-"
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  Vertex-set choice:")
    print("    harmonic modes / Toeplitz proof lenses, not runners, arcs, or packets.")
    print("  LRC predicate preserved:")
    print("    danger-covering implies nonnegative multiplicity F_S=D_S-1, hence PSD")
    print("    Toeplitz moment matrices for every degree.")
    print("  Information destroyed:")
    print("    endpoint owners, exact interval components, C27/K33 addresses, and")
    print("    packet-family ancestry are invisible unless reattached afterwards.")
    print("  Challenged assumption:")
    print("    the hard residual may be easier to disprove as a positive-measure")
    print("    trigonometric moment problem than as a labelled covering packet.")


def print_dual_lemma() -> None:
    print("[1] Fourier-Toeplitz dual lemma")
    print("  Let D_S(t)=sum_v 1_{||v t||<1/14} and F_S=D_S-1.")
    print("  If S is a strict counterexample, then F_S(t)>=0 a.e.")
    print("  Therefore for every polynomial q(t)=sum_{j=0}^K c_j exp(2*pi*i*j*t):")
    print("      integral F_S(t) |q(t)|^2 dt >= 0")
    print("  Equivalently T_K=(hat F_S(i-j)) is PSD.")
    print("  A negative eigenvalue of T_K is a rigorous dual obstruction to covering,")
    print("  up to numerical robustness of the displayed eigenvalue.")


LENS_VECTORS: dict[str, tuple[int, ...]] = {
    "Toeplitz_PSD_dual": (6, 6, 6, 5, 4, 6),
    "low_AP_harmonics": (5, 5, 4, 4, 3, 5),
    "apex_14_modes": (4, 5, 5, 5, 4, 5),
    "C27_K33_modes": (4, 4, 6, 6, 5, 5),
    "few_apex_lift_modes": (5, 4, 5, 5, 6, 5),
    "boundary_moment_packets": (5, 3, 5, 6, 6, 6),
    "exact_interval_front": (4, 6, 5, 4, 4, 6),
    "raw_runner_tournament": (1, 1, 1, 1, 1, 1),
}


def beats(a: str, b: str) -> bool:
    va = LENS_VECTORS[a]
    vb = LENS_VECTORS[b]
    score = sum(1 if x > y else -1 if x < y else 0 for x, y in zip(va, vb))
    if score:
        return score > 0
    return list(LENS_VECTORS).index(a) < list(LENS_VECTORS).index(b)


def tournament_fingerprint() -> dict[str, object]:
    names = list(LENS_VECTORS)
    outdeg = Counter()
    c3 = 0
    for a in names:
        outdeg[a] = sum(1 for b in names if a != b and beats(a, b))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            for k, c in enumerate(names):
                if i < j < k:
                    if beats(a, b) and beats(b, c) and beats(c, a):
                        c3 += 1
                    if beats(a, c) and beats(c, b) and beats(b, a):
                        c3 += 1
    return {
        "score_hist": dict(sorted(Counter(outdeg.values()).items())),
        "c3": c3,
        "order": " -> ".join(sorted(names, key=lambda x: (-outdeg[x], names.index(x)))),
    }


def print_results(results: list[DualResult], degrees: tuple[int, ...], tol: float) -> None:
    print("[2] Scan summary")
    print(f"  rows audited={len(results)}")
    print(f"  degrees={','.join(map(str, degrees))}")
    print(f"  negativity tolerance={tol:g}")
    print(f"  dual-certified rows={sum(r.first_negative_degree is not None for r in results)}")
    print(f"  PSD-through-max rows={sum(r.first_negative_degree is None for r in results)}")
    print(f"  status counts={dict(sorted(Counter(r.exact.status for r in results).items()))}")
    print(f"  family counts={dict(sorted(Counter(r.row.family for r in results).items()))}")
    print(f"  first-negative degree counts={dict(sorted(Counter(r.first_negative_degree or 'none' for r in results).items(), key=lambda kv: str(kv[0])))}")
    print(f"  dominant-band counts={dict(sorted(Counter(r.dominant_band for r in results).items()))}")


def print_representatives(results: list[DualResult], limit: int) -> None:
    print("[3] Representative rows")
    ordered = sorted(results, key=lambda r: (r.first_negative_degree is None, r.first_negative_degree or 10**9, r.min_eig, r.row.name))
    print("  strongest / earliest dual certificates:")
    for r in ordered[:limit]:
        first = r.first_negative_degree if r.first_negative_degree is not None else "-"
        energies = ", ".join(f"{k}:{v:.3f}" for k, v in r.band_energy[:3])
        print(
            f"    {r.row.name:34s} qdiv={r.row.qdiv:<3d} status={r.exact.status:14s} "
            f"mu={fmt_frac(r.exact.mass):>10s} firstK={str(first):>3s} "
            f"minEig={r.min_eig:+.6e}@{r.min_eig_degree:<3d} band={r.dominant_band:18s} {energies}"
        )
    print("  PSD-through-max boundary/invisible rows:")
    invisible = [r for r in results if r.first_negative_degree is None]
    for r in invisible[:limit]:
        energies = ", ".join(f"{k}:{v:.3f}" for k, v in r.band_energy[:3])
        print(
            f"    {r.row.name:34s} qdiv={r.row.qdiv:<3d} status={r.exact.status:14s} "
            f"mu={fmt_frac(r.exact.mass):>10s} minEig={r.min_eig:+.6e}@{r.min_eig_degree:<3d} "
            f"band={r.dominant_band:18s} {energies}"
        )


def print_theorem_readout(results: list[DualResult]) -> None:
    tight = [r for r in results if r.exact.status == "boundary_only"]
    positive_cert = [r for r in results if r.exact.status == "positive_open" and r.first_negative_degree is not None]
    positive_invisible = [r for r in results if r.exact.status == "positive_open" and r.first_negative_degree is None]
    print("[4] Theorem-facing readout")
    print(f"  boundary-only rows certified negative={sum(r.first_negative_degree is not None for r in tight)}/{len(tight)}")
    print(f"  positive-open rows certified negative={len(positive_cert)}")
    print(f"  positive-open rows still PSD-through-max={len(positive_invisible)}")
    if positive_invisible:
        print("  invisible positive-open rows show this dual lens is incomplete:")
        for r in positive_invisible[:8]:
            print(f"    {r.row.name}: mu={fmt_frac(r.exact.mass)}, qdiv={r.row.qdiv}, minEig={r.min_eig:+.3e}")
    print("  Proposed new proof route:")
    print("    prove a structured Toeplitz negativity theorem for every qdiv>14")
    print("    covering residual after AP/GW boundary atoms and K33/state-lift")
    print("    exceptions are removed; use packets only after the harmonic dual lens")
    print("    identifies the missing frequency band.")


def print_tournament_analysis() -> None:
    fp = tournament_fingerprint()
    print("[5] Tournament Analysis")
    print("  vertices are harmonic proof lenses, not runners or packets.")
    print("  pairwise observable:")
    print("    which lens preserves the implication cover => F_S>=0 => Toeplitz PSD")
    print("    while retaining enough labelled information for known exits.")
    print("  tie Hamiltonian path:")
    print("    " + " -> ".join(LENS_VECTORS))
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']}")
    print(f"  score order: {fp['order']}")


def parse_degrees(text: str) -> tuple[int, ...]:
    return tuple(sorted({int(x) for x in text.split(",") if x.strip()}))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=120)
    parser.add_argument("--two-swap-limit", type=int, default=32)
    parser.add_argument("--max-generated", type=int, default=40)
    parser.add_argument("--degrees", default="16,24,32,40,48,56,64,72,80,90")
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--representatives", type=int, default=12)
    args = parser.parse_args()

    degrees = parse_degrees(args.degrees)
    rows = unique_rows(curated_rows() + generated_rows(args.single_limit, args.two_swap_limit, args.max_generated))
    results = [scan_row(row, degrees, args.tol) for row in rows]

    print("LRC14 Fourier-Toeplitz dual scout")
    print("=" * 78)
    print(
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"max_generated={args.max_generated}"
    )
    print_assumption_challenge()
    print_dual_lemma()
    print_results(results, degrees, args.tol)
    print_representatives(results, args.representatives)
    print_theorem_readout(results)
    print_tournament_analysis()


if __name__ == "__main__":
    main()
