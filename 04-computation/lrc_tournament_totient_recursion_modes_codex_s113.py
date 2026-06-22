#!/usr/bin/env python3
"""Totient/coprime-density audit for the three tournament recursion modes.

The three tournament recurrences are cell-address recurrences:

  full:      A+B+C-D-E-F+G
  even half: A+B-C
  odd half:  A+B-C+D-E-F+G

They are exact for the full and half tiling cell counts, not for multiplicative
objects such as phi(n), phi(n)/n, or exact-period LRC witness packets.  This
script treats that failure as useful signal: the signed residual is an
"Euler-factor curvature" showing which coprime-density labels must survive
before LRC exact-period packets are scalarized.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, gcd


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    m = n
    while d * d <= m:
        while m % d == 0:
            out[d] = out.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


def fmt_factor(n: int) -> str:
    f = factor(n)
    if not f:
        return "1"
    parts = []
    for p in sorted(f):
        e = f[p]
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(parts)


def phi(n: int) -> int:
    ans = n
    for p in factor(n):
        ans = ans // p * (p - 1)
    return ans


def unit_density(n: int) -> Fraction:
    return Fraction(phi(n), n)


def ap_endpoint_density(n: int) -> Fraction:
    if n <= 1:
        return Fraction(0)
    return Fraction(phi(n), n - 1)


def full_cells(n: int) -> int:
    return comb(n - 1, 2) if n >= 2 else 0


def half_cells(n: int) -> int:
    return ((n - 1) * (n - 1)) // 4 if n >= 2 else 0


def slots(mode: str, n: int) -> list[tuple[str, int, int]]:
    """Return (label, sign, subtournament_size) for the mode at size n."""
    if mode == "full":
        return [
            ("A", +1, n - 1),
            ("B", +1, n - 1),
            ("C", +1, n - 1),
            ("D", -1, n - 2),
            ("E", -1, n - 2),
            ("F", -1, n - 2),
            ("G", +1, n - 3),
        ]
    if mode == "even_half":
        return [("A", +1, n - 1), ("B", +1, n - 1), ("C", -1, n - 2)]
    if mode == "odd_half":
        return [
            ("A", +1, n - 1),
            ("B", +1, n - 1),
            ("C", -1, n - 2),
            ("D", +1, n - 2),
            ("E", -1, n - 3),
            ("F", -1, n - 3),
            ("G", +1, n - 4),
        ]
    raise ValueError(mode)


def slot_expr(mode: str, n: int) -> str:
    return " ".join(
        f"{'+' if s > 0 else '-'}{lab}:{m}({fmt_factor(m)})" for lab, s, m in slots(mode, n)
    ).lstrip("+")


def op_value(mode: str, n: int, func) -> Fraction:
    total = Fraction(0)
    for _lab, sign, size in slots(mode, n):
        total += sign * Fraction(func(size))
    return total


def residual(mode: str, n: int, func) -> Fraction:
    return Fraction(func(n)) - op_value(mode, n, func)


def exponent_curvature(mode: str, n: int) -> dict[int, int]:
    """Prime-exponent residual: v_p(n) - signed sum v_p(slot sizes)."""
    cur = Counter(factor(n))
    for _lab, sign, size in slots(mode, n):
        for p, e in factor(size).items():
            cur[p] -= sign * e
    return {p: cur[p] for p in sorted(cur) if cur[p] != 0}


def divisor_sum_check(limit: int) -> None:
    print("Totient copy law: n = sum_{d|n} phi(d)")
    failures = []
    for n in range(1, limit + 1):
        s = sum(phi(d) for d in range(1, n + 1) if n % d == 0)
        if s != n:
            failures.append((n, s))
    print(f"  checked n<= {limit}: failures={len(failures)}")
    print()


def exactness_checks(limit: int) -> None:
    print("Cell-count exactness versus multiplicative curvature")
    full_fail = []
    even_fail = []
    odd_fail = []
    for n in range(5, limit + 1):
        if residual("full", n, full_cells) != 0:
            full_fail.append(n)
    for n in range(4, limit + 1, 2):
        if residual("even_half", n, half_cells) != 0:
            even_fail.append(n)
    for n in range(5, limit + 1, 2):
        if residual("odd_half", n, half_cells) != 0:
            odd_fail.append(n)
    print(f"  full recurrence exact for full cells through {limit}: {not full_fail}")
    print(f"  even half exact for half cells through {limit}: {not even_fail}")
    print(f"  odd half exact for half cells through {limit}: {not odd_fail}")
    for mode, ns in (("full", range(5, 13)), ("even_half", range(4, 14, 2)), ("odd_half", range(5, 15, 2))):
        nonzero = sum(1 for n in ns if residual(mode, n, phi) != 0)
        print(f"  phi residuals nonzero for {mode} sample: {nonzero}/{len(list(ns))}")
    print()


def spotlight() -> None:
    print("LRC14-adjacent slot tables")
    rows = [
        ("full", 14),
        ("even_half", 14),
        ("full", 15),
        ("odd_half", 15),
    ]
    for mode, n in rows:
        print(f"{mode:9s} n={n}: {slot_expr(mode, n)}")
        for name, func in (
            ("phi", phi),
            ("rho=phi/n", unit_density),
            ("APdens=phi/(n-1)", ap_endpoint_density),
        ):
            r = residual(mode, n, func)
            print(f"  residual {name:17s}: {r} ({float(r):+.6f})")
        print(f"  prime-exponent curvature: {exponent_curvature(mode, n)}")
    print()


def largest_residuals(limit: int) -> None:
    print("Largest coprime-density residuals, normalized by exact value")
    records: list[tuple[Fraction, str, int, Fraction, Fraction, dict[int, int]]] = []
    by_mode: dict[str, list[tuple[Fraction, str, int, Fraction, Fraction, dict[int, int]]]] = {
        "full": [],
        "even_half": [],
        "odd_half": [],
    }
    for mode in ("full", "even_half", "odd_half"):
        ns = range(5, limit + 1)
        if mode == "even_half":
            ns = range(4, limit + 1, 2)
        elif mode == "odd_half":
            ns = range(5, limit + 1, 2)
        for n in ns:
            val = unit_density(n)
            r = residual(mode, n, unit_density)
            if val:
                rec = (abs(r) / val, mode, n, r, val, exponent_curvature(mode, n))
                records.append(rec)
                by_mode[mode].append(rec)
    records.sort(reverse=True, key=lambda x: (x[0], x[2]))
    print("  global top residuals:")
    for score, mode, n, r, val, curv in records[:8]:
        print(
            f"  {mode:9s} n={n:2d} factor={fmt_factor(n):12s} "
            f"rho={val} residual={r} rel={score} curvature={curv}"
        )
    print("  per-mode top residuals:")
    for mode in ("full", "even_half", "odd_half"):
        by_mode[mode].sort(reverse=True, key=lambda x: (x[0], x[2]))
        print(f"    {mode}:")
        for score, _mode, n, r, val, curv in by_mode[mode][:4]:
            print(
                f"      n={n:2d} factor={fmt_factor(n):12s} "
                f"rho={val} residual={r} rel={score} curvature={curv}"
            )
    print()


def multiplicativity_defects(limit: int) -> None:
    print("Multiplicativity defects for unit-density rho=phi/n")
    coprime_fail = 0
    noncoprime_examples: list[tuple[Fraction, int, int, Fraction, Fraction, Fraction]] = []
    for a in range(2, limit + 1):
        for b in range(2, limit + 1):
            defect = unit_density(a * b) - unit_density(a) * unit_density(b)
            if gcd(a, b) == 1:
                if defect != 0:
                    coprime_fail += 1
            elif defect != 0:
                noncoprime_examples.append((abs(defect), a, b, defect, unit_density(a * b), unit_density(a) * unit_density(b)))
    noncoprime_examples.sort(reverse=True)
    print(f"  coprime multiplicativity failures for a,b<= {limit}: {coprime_fail}")
    print("  largest non-coprime defects:")
    for _absd, a, b, defect, rho_ab, prod in noncoprime_examples[:8]:
        print(
            f"    a={a:2d}({fmt_factor(a):6s}) b={b:2d}({fmt_factor(b):6s}) "
            f"rho(ab)={rho_ab} rho(a)rho(b)={prod} defect={defect}"
        )
    print()


def mode_tournament() -> None:
    print("Tournament Analysis over proof carriers")
    features = {
        "exact_period_phi_packets": {
            "coprime_density",
            "euler_factor_labels",
            "CRT_multiplicative",
            "LRC_witness_packets",
            "not_cell_affine",
        },
        "euler_curvature_residual": {
            "coprime_density",
            "euler_factor_labels",
            "signed_boundary",
            "slot_sizes",
            "not_cell_affine",
        },
        "three_recursion_slot_atlas": {
            "signed_boundary",
            "slot_sizes",
            "mirror_parity",
            "full_half_modes",
        },
        "half_tiling_fixed_line": {
            "mirror_parity",
            "fixed_line",
            "slot_sizes",
            "full_half_modes",
        },
        "scalar_cell_count": {
            "cell_affine",
            "signed_boundary",
            "full_half_modes",
        },
        "raw_subtournament_size": {"slot_sizes"},
        "raw_runner_vertices": set(),
    }
    names = list(features)
    adj = {x: set() for x in names}
    scores = Counter({name: 0 for name in names})
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            key_a = (len(features[a]), len(features[a] & features[b]), -i)
            key_b = (len(features[b]), len(features[a] & features[b]), -names.index(b))
            if key_a >= key_b:
                adj[a].add(b)
                scores[a] += 1
            else:
                adj[b].add(a)
                scores[b] += 1
    cycles3 = 0
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            for k, c in enumerate(names):
                if k <= j:
                    continue
                if b in adj[a] and c in adj[b] and a in adj[c]:
                    cycles3 += 1
                if c in adj[a] and b in adj[c] and a in adj[b]:
                    cycles3 += 1
    path = sorted(names, key=lambda x: (scores[x], len(features[x]), x), reverse=True)
    print("  vertices are proof carriers, not runners/arcs")
    print("  pairwise observable=(#labels preserved, overlap, declaration order)")
    print(f"  scores={dict(scores)}")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={cycles3}")
    print("  Hamiltonian path=" + " > ".join(path))
    print()


def main() -> None:
    print("Totient/coprime-density audit for tournament recursion modes")
    print("=" * 78)
    divisor_sum_check(120)
    exactness_checks(80)
    spotlight()
    largest_residuals(80)
    multiplicativity_defects(40)
    mode_tournament()
    print("Synthesis:")
    print("  The three A..G recurrences are exact address recurrences for cell carriers.")
    print("  Totient and coprime density are multiplicative exact-period packet counts.")
    print("  Their recurrence residual is not error; it is Euler-factor curvature.")
    print("  LRC14 should keep this curvature as labelled CRT/chi7/coimage data before")
    print("  applying scalar cap, floor, or Fejer estimates.")


if __name__ == "__main__":
    main()
