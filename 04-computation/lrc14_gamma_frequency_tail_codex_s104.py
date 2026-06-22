#!/usr/bin/env python3
"""
codex-S104: same-frequency additive-energy tail for the LRC(14) cover atom.

KPS S31k proved the first Fourier-frequency additive-energy coefficient

    Gamma_k(1) = sum_A (-1)^|A| |hat f_A(1)|^4 (1-|A|/7)^(k-4)

is positive for k=8..12.  This script asks whether that is the first term of a
clean convergent packet:

    Gamma_k^sf = sum_{m>=1, 7 not divide m} Gamma_k(m),
    Gamma_k(m) = sum_A (-1)^|A| |hat f_A(m)|^4 (1-|A|/7)^(k-4).

For fixed residue r mod 7, m^4 Gamma_k(m) is constant.  Hence the same-frequency
packet has an explicit 1/m^4 tail.  This isolates the additive-energy part of
the higher resonance error from the support >= 5 hidden-fold / shape terms.

Tournament Analysis note:
  vertices are proof lenses, not runners:
    m=1_lead, same_frequency_tail, hidden_fold_shape, support_cycle_tail.
  The observable is (explained_AP_deviation, monotonicity_payoff, remaining
  residual, formal_tail_cost).  This script tests the first two lenses and
  leaves the residual for the octahedral / hidden-sum machinery.
"""

from __future__ import annotations

import cmath
import math
from collections import Counter
from fractions import Fraction
from itertools import chain, combinations


SECTORS = tuple(range(1, 7))
TWO_PI = 2.0 * math.pi


def subsets(items):
    return chain.from_iterable(combinations(items, r) for r in range(len(items) + 1))


def gamma_hat_prefactor(m: int) -> complex:
    if m % 7 == 0:
        return 0j
    return (1.0 - cmath.exp(-2j * math.pi * m / 7.0)) / (2j * math.pi * m)


def sector_sum(A: tuple[int, ...], m: int) -> complex:
    return sum(cmath.exp(-2j * math.pi * m * j / 7.0) for j in A)


def gamma_k_m(k: int, m: int) -> float:
    """Same-frequency additive-energy coefficient at Fourier frequency m."""
    if m % 7 == 0:
        return 0.0
    pref = gamma_hat_prefactor(m)
    total = 0.0
    for A in subsets(SECTORS):
        a = len(A)
        hat = -pref * sector_sum(A, m)
        total += ((-1) ** a) * (abs(hat) ** 4) * ((1.0 - a / 7.0) ** (k - 4))
    return float(total)


def residue_constant(k: int, r: int) -> float:
    """C_{k,r} with Gamma_k(m)=C_{k,r}/m^4 for m == r mod 7."""
    assert 1 <= r <= 6
    return gamma_k_m(k, r) * (r**4)


def same_frequency_partial(k: int, H: int) -> float:
    return sum(gamma_k_m(k, m) for m in range(1, H + 1) if m % 7)


def residue_tail_sum(r: int, H: int, extra_terms: int = 2000) -> tuple[float, float]:
    """Return (sum, integral_upper_remainder) of sum_{m>H, m=r mod 7} 1/m^4."""
    q0 = 0
    while r + 7 * q0 <= H:
        q0 += 1
    q1 = q0 + extra_terms
    s = sum(1.0 / ((r + 7 * q) ** 4) for q in range(q0, q1 + 1))
    # decreasing positive f(q)=(r+7q)^-4, sum_{q>q1} f(q) <= int_{q1}^inf f(x) dx
    rem = 1.0 / (21.0 * ((r + 7 * q1) ** 3))
    return s, rem


def same_frequency_tail_bound(k: int, H: int) -> float:
    """Absolute tail bound for sum_{m>H} Gamma_k(m)."""
    total = 0.0
    for r in range(1, 7):
        s, rem = residue_tail_sum(r, H)
        total += abs(residue_constant(k, r)) * (s + rem)
    return total


def same_frequency_series(k: int, H: int = 5000) -> tuple[float, float]:
    partial = same_frequency_partial(k, H)
    return partial, same_frequency_tail_bound(k, H)


def p0_decorr(k: int) -> Fraction:
    total = Fraction(0)
    for A in subsets(SECTORS):
        a = len(A)
        total += ((-1) ** a) * Fraction(7 - a, 7) ** k
    return total


def sector_of(e: int, x: Fraction) -> int:
    y = e * x
    y -= y.numerator // y.denominator
    return (7 * y.numerator) // y.denominator


def p0_exact(E: list[int]) -> Fraction:
    speeds = sorted(set(int(e) for e in E if e != 0))
    bps = {Fraction(0), Fraction(1)}
    for e in speeds:
        ae = abs(e)
        for m in range(0, 7 * ae + 1):
            bps.add(Fraction(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = Fraction(0)
    target = set(SECTORS)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = {sector_of(e, mid) for e in speeds}
        if target <= hit:
            total += hi - lo
    return total


def additive_energy(E: list[int]) -> int:
    c = Counter(a + b for a in E for b in E)
    return sum(v * v for v in c.values())


def astar(E: list[int]) -> int:
    k = len(E)
    return additive_energy(E) - (2 * k * k - k)


def row_data(label: str, E: list[int], series_by_k: dict[int, float]) -> dict[str, float | int | str]:
    k = len(E)
    p0 = p0_exact(E)
    pd = p0_decorr(k)
    dev = float(p0 - pd)
    A = astar(E)
    g1 = gamma_k_m(k, 1)
    gsf = series_by_k[k]
    return {
        "label": label,
        "k": k,
        "E": tuple(E),
        "p0": float(p0),
        "p0_num": p0.numerator,
        "p0_den": p0.denominator,
        "decorr": float(pd),
        "dev": dev,
        "Astar": A,
        "lead": A * g1,
        "samefreq": A * gsf,
        "residual_after_samefreq": dev - A * gsf,
        "lead_fraction_of_dev": (A * g1 / dev) if dev else float("nan"),
        "samefreq_fraction_of_dev": (A * gsf / dev) if dev else float("nan"),
    }


def print_gamma_tables() -> dict[int, float]:
    print("=== Same-frequency additive-energy packet Gamma_k(m) ===")
    print("Gamma_k(m) = C_{k,r}/m^4 for m == r mod 7; all constants below are signed.")
    series_by_k: dict[int, float] = {}
    for k in range(8, 14):
        constants = [residue_constant(k, r) for r in range(1, 7)]
        partial_12 = same_frequency_partial(k, 12)
        tail_12 = same_frequency_tail_bound(k, 12)
        series, tail_5000 = same_frequency_series(k, 5000)
        series_by_k[k] = series
        pos_cert = partial_12 - tail_12
        print(f"\nk={k}")
        print("  C_{k,r}, r=1..6:")
        print("   ", " ".join(f"{c:+.9e}" for c in constants))
        print(f"  Gamma(1)={gamma_k_m(k,1):+.9e}")
        print(f"  partial H=12={partial_12:+.9e}; abs tail <= {tail_12:.3e}; positivity cert >= {pos_cert:+.9e}")
        print(f"  series H=5000={series:+.9e}; abs tail <= {tail_5000:.3e}; ratio series/Gamma1={series/gamma_k_m(k,1):.6f}")
    return series_by_k


def print_row_validation(series_by_k: dict[int, float]) -> None:
    rows = [
        ("AP-k8", list(range(8))),
        ("nearAP-k8-gap", [0, 1, 2, 3, 4, 5, 6, 8]),
        ("sidon-k8", [0, 1, 3, 7, 12, 20, 30, 44]),
        ("AP-k9", list(range(9))),
        ("nearAP-k9-gap", [0, 1, 2, 3, 4, 5, 6, 7, 9]),
        ("dilate2-AP-k9", [2 * x for x in range(9)]),
        ("sidon-k9", [0, 1, 3, 7, 12, 20, 30, 44, 60]),
        ("AP-k10", list(range(10))),
        ("nearAP-k10-gap", [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]),
        ("AP-k11", list(range(11))),
        ("AP-k12", list(range(12))),
        ("AP-k13", list(range(13))),
    ]
    print("\n=== Exact p0 validation against additive-energy packet ===")
    print(
        "label                 k  A*     p0             dev        lead       samefreq   resid_sf   sf/dev"
    )
    for label, E in rows:
        d = row_data(label, E, series_by_k)
        print(
            f"{label:20s} {d['k']:2d} {d['Astar']:5d} "
            f"{d['p0_num']}/{d['p0_den']:<8d} "
            f"{d['dev']:+.6f} {d['lead']:+.6f} {d['samefreq']:+.6f} "
            f"{d['residual_after_samefreq']:+.6f} {d['samefreq_fraction_of_dev']:+.3f}"
        )


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n == 0:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / math.sqrt(vx * vy)


def scan_bounded_bank(k: int, max_value: int, series_by_k: dict[int, float]) -> None:
    """Exact anchored bank scan for the residual-leak proof target."""
    pd = float(p0_decorr(k))
    gsf = series_by_k[k]
    AP = list(range(k))
    AP_A = astar(AP)
    AP_p0 = float(p0_exact(AP))
    AP_R = (AP_p0 - pd) - gsf * AP_A
    rows = []
    for tail in combinations(range(1, max_value + 1), k - 1):
        E = [0, *tail]
        A = astar(E)
        p0 = float(p0_exact(E))
        R = (p0 - pd) - gsf * A
        energy_gain = gsf * (AP_A - A)
        residual_leak = R - AP_R
        if energy_gain > 1e-15:
            ratio = residual_leak / energy_gain
        else:
            ratio = float("-inf")
        rows.append((p0, A, R, energy_gain, residual_leak, ratio, tuple(E)))

    rows_by_p0 = sorted(rows, reverse=True, key=lambda row: row[0])
    risky = sorted(
        [row for row in rows if row[3] > 1e-15 and row[4] > 0],
        reverse=True,
        key=lambda row: row[5],
    )
    As = [float(row[1]) for row in rows]
    p0s = [row[0] for row in rows]
    Rs = [row[2] for row in rows]
    print(f"\n=== Exact bounded residual-leak scan: k={k}, anchored max<= {max_value} ===")
    print(f"rows={len(rows)}  AP_A*={AP_A}  AP_p0={AP_p0:.9f}  AP_Rsf={AP_R:+.9f}")
    print(f"corr(A*, p0)={pearson(As,p0s):+.4f}  corr(A*, residual_after_sf)={pearson(As,Rs):+.4f}")
    print("top p0 rows:")
    for p0, A, R, energy_gain, residual_leak, ratio, E in rows_by_p0[:6]:
        print(f"  p0={p0:.9f} A*={A:4d} R={R:+.9f} E={E}")
    print("largest residual-leak / same-frequency-energy-gap ratios:")
    for p0, A, R, energy_gain, residual_leak, ratio, E in risky[:8]:
        ap_minus = AP_p0 - p0
        print(
            f"  ratio={ratio:.6f} AP-p0={ap_minus:+.9f} "
            f"Egap={energy_gain:.9f} Rleak={residual_leak:.9f} "
            f"A*={A:4d} E={E}"
        )
    violations = [row for row in risky if row[5] > 1.0 + 1e-12]
    print(f"residual-leak violations of AP dominance inequality: {len(violations)}")


def print_tournament_lens() -> None:
    print("\n=== Tournament-analysis lens ranking ===")
    lenses = [
        ("m=1_lead", 3, 3, 1, 1),
        ("same_frequency_tail", 3, 3, 1, 2),
        ("hidden_fold_shape", 2, 2, 3, 3),
        ("support_cycle_tail", 1, 1, 3, 4),
    ]
    print("vertices: lens, explained_AP, monotonicity, remaining_residual, formal_tail_cost")
    for lens in lenses:
        print(f"  {lens[0]:20s} {lens[1:]}")
    print("Hamiltonian path: m=1_lead -> same_frequency_tail -> hidden_fold_shape -> support_cycle_tail")
    print("challenged assumption: additive energy alone is not enough; shifted/hidden-fold rows keep shape data.")


def main() -> None:
    series_by_k = print_gamma_tables()
    print_row_validation(series_by_k)
    scan_bounded_bank(8, 14, series_by_k)
    scan_bounded_bank(9, 14, series_by_k)
    print_tournament_lens()
    print("\nInterpretation:")
    print("  The same-frequency additive-energy packet is positive for k=8..13 with a tiny explicit")
    print("  1/m^4 tail after H=12.  It strengthens the KPS m=1 monotonicity signal.")
    print("  However, AP residuals after this packet are not zero; the remaining proof obligation")
    print("  lives in hidden-fold / support-cycle shape terms, not in scalar additive energy alone.")


if __name__ == "__main__":
    main()
