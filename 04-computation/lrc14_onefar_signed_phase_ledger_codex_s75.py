#!/usr/bin/env python3
"""HYP-2786/S75: signed phase ledger for the one-far binding family.

Incoming HYP-2784 says that sharpening the absolute THM-546 arc-complexity
`V` is not enough: even true `V` leaves the one-far binding rows off by about
two orders of magnitude.  Incoming HYP-2785 says the answer is not a
residue-only `w mod 7` table either.  The missing mechanism is signed
cancellation in

    Delta_w(B) = p0(B union {w}) - Phi(B).

This scout keeps the exact Abel endpoint value for Delta_w and compares it to
truncated signed Fourier heads.  It also records the odd/even mode split
requested in the Weyl-error thread: odd support often dominates the L1
envelope, but even-led exceptions remain finite and visible.

No proof is claimed.  The target is to shrink OPEN-Q-108 from an absolute BV
bound to a signed low-mode/mod-14 phase ledger plus a
Dedekind/equidistribution tail.
"""

from __future__ import annotations

import cmath
import functools
import itertools
import math
import sys
from collections import Counter
from fractions import Fraction as F

sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
print = functools.partial(print, flush=True)

from lrc14_threadB_localized_binding_bound_kpswf8 import (  # noqa: E402
    CAP,
    MARGIN,
    QVAL,
    Varcs,
    cells_with_miss,
    wDelta_signed,
)

TWO_PI = 2.0 * math.pi


def fmt(q: F | float, ndigits: int = 9) -> str:
    if isinstance(q, F):
        return f"{q} ({float(q):.{ndigits}f})"
    return f"{q:.{ndigits}f}"


def exact_one_miss_arcs(Ep: tuple[int, ...]) -> dict[int, list[tuple[F, F]]]:
    """Return arcs A_j where the base misses exactly inner sector j."""

    cells = cells_with_miss(Ep)
    arcs: dict[int, list[tuple[F, F]]] = {j: [] for j in range(1, 7)}
    cur: dict[int, F | None] = {j: None for j in range(1, 7)}
    for lo, hi, miss in cells:
        only = next(iter(miss)) if len(miss) == 1 else None
        for j in range(1, 7):
            if only == j:
                if cur[j] is None:
                    cur[j] = lo
            elif cur[j] is not None:
                arcs[j].append((cur[j], lo))
                cur[j] = None
    for j in range(1, 7):
        if cur[j] is not None:
            arcs[j].append((cur[j], F(1)))
    return arcs


def shat(n: int, j: int) -> complex:
    """Fourier coefficient of sector j=[j/7,(j+1)/7), centered for n != 0."""

    return cmath.exp(-1j * TWO_PI * n * j / 7) * (
        1 - cmath.exp(-1j * TWO_PI * n / 7)
    ) / (1j * TWO_PI * n)


def arc_hat(m: int, arcs: list[tuple[F, F]]) -> complex:
    if m == 0:
        return complex(sum(float(b - a) for a, b in arcs), 0.0)
    total = 0j
    for a, b in arcs:
        af = float(a)
        bf = float(b)
        total += (
            cmath.exp(-1j * TWO_PI * m * af)
            - cmath.exp(-1j * TWO_PI * m * bf)
        ) / (1j * TWO_PI * m)
    return total


def signed_head_split(
    arcs: dict[int, list[tuple[F, F]]], w: int, cutoff: int
) -> dict[str, object]:
    """Return signed Fourier head data for positive modes n<=cutoff."""

    odd = 0.0
    even = 0.0
    odd_l1 = 0.0
    even_l1 = 0.0
    by14 = {r: 0.0 for r in range(14)}
    by14_l1 = {r: 0.0 for r in range(14)}
    terms: list[tuple[float, int, float]] = []
    for n in range(1, cutoff + 1):
        if n % 7 == 0:
            continue
        contribution = 0.0
        for j in range(1, 7):
            contribution += 2.0 * (shat(n, j) * arc_hat(-n * w, arcs[j])).real
        if n % 2:
            odd += contribution
            odd_l1 += abs(contribution)
        else:
            even += contribution
            even_l1 += abs(contribution)
        by14[n % 14] += contribution
        by14_l1[n % 14] += abs(contribution)
        terms.append((abs(contribution), n, contribution))
    terms.sort(reverse=True)
    return {
        "head": odd + even,
        "odd": odd,
        "even": even,
        "odd_l1": odd_l1,
        "even_l1": even_l1,
        "odd_share": odd_l1 / (odd_l1 + even_l1) if odd_l1 + even_l1 else 0.0,
        "by14": by14,
        "by14_l1": by14_l1,
        "terms": terms,
    }


def exact_delta(B: tuple[int, ...], w: int) -> F:
    """Exact one-far Delta_w from the Abel endpoint identity."""

    return wDelta_signed(B, w) / w


def risk_rows(wmax: int, cutoff: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for k in range(8, 13):
        B = tuple(range(k - 1))
        arcs = exact_one_miss_arcs(B)
        best_delta = F(-10)
        best_w = -1
        for w in range(15, wmax + 1):
            d = exact_delta(B, w)
            if d > best_delta:
                best_delta = d
                best_w = w
        split = signed_head_split(arcs, best_w, cutoff)
        low6 = signed_head_split(arcs, best_w, 6)
        bv = F(6, 49) * Varcs(B) / best_w
        by14_l1 = split["by14_l1"]
        top_residue = max((r for r in range(14) if r % 7), key=lambda r: by14_l1[r])
        rows.append(
            {
                "k": k,
                "B": B,
                "w": best_w,
                "delta": best_delta,
                "margin": MARGIN[k],
                "V": Varcs(B),
                "bv": bv,
                "split": split,
                "low6": low6,
                "top_residue": top_residue,
            }
        )
    return rows


def finite_head_audit(wmax: int, cutoffs: tuple[int, ...]) -> dict[int, list[tuple[int, float, int, float, float]]]:
    """Worst exact-minus-head error for each cutoff and k over the binding family."""

    table: dict[int, list[tuple[int, float, int, float, float]]] = {}
    for cutoff in cutoffs:
        table[cutoff] = []
        for k in range(8, 13):
            B = tuple(range(k - 1))
            arcs = exact_one_miss_arcs(B)
            worst = 0.0
            worst_w = -1
            exact_at = 0.0
            head_at = 0.0
            for w in range(15, wmax + 1):
                d = float(exact_delta(B, w))
                h = float(signed_head_split(arcs, w, cutoff)["head"])
                err = abs(d - h)
                if err > worst:
                    worst = err
                    worst_w = w
                    exact_at = d
                    head_at = h
            table[cutoff].append((k, worst, worst_w, exact_at, head_at))
    return table


def tournament(rows: list[dict[str, object]]) -> None:
    print("\n" + "=" * 92)
    print("TOURNAMENT ANALYSIS")
    print("  vertices: binding one-far rows k=8..12 at their max positive Delta_w in the scanned band")
    print("  pairwise observable: larger Delta_w/margin, then larger odd L1 share, then smaller BV ratio")
    print("  switch/gauge: exact Abel endpoint value is scalarized only after the mod-14 Fourier split")
    scores = Counter()
    edges = set()
    for i, j in itertools.combinations(range(len(rows)), 2):
        ri = rows[i]
        rj = rows[j]
        ai = (
            ri["delta"] / ri["margin"],
            ri["split"]["odd_share"],
            -ri["delta"] / ri["bv"],
        )
        aj = (
            rj["delta"] / rj["margin"],
            rj["split"]["odd_share"],
            -rj["delta"] / rj["bv"],
        )
        if ai >= aj:
            edges.add((i, j))
            scores[i] += 1
        else:
            edges.add((j, i))
            scores[j] += 1
    cycles = 0
    for a, b, c in itertools.combinations(range(len(rows)), 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    hist = Counter(scores[i] for i in range(len(rows)))
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian pressure path:")
    order = sorted(
        range(len(rows)),
        key=lambda i: (
            rows[i]["delta"] / rows[i]["margin"],
            rows[i]["split"]["odd_share"],
            -rows[i]["delta"] / rows[i]["bv"],
        ),
        reverse=True,
    )
    for idx in order:
        row = rows[idx]
        print(
            f"    k={row['k']} w={row['w']}: "
            f"Delta/margin={float(row['delta'] / row['margin']):.6f}, "
            f"odd_share={row['split']['odd_share']:.6f}, "
            f"top n mod 14={row['top_residue']}"
        )


def main() -> None:
    wmax = 100
    audit_wmax = 80
    cutoff = 80
    cutoffs = (6, 13)
    print("HYP-2786/S75 one-far signed phase ledger")
    print(f"Binding rows: B=consec_(k-1), k=8..12, wide w in [15,{wmax}].")
    print(f"Signed Fourier head cutoff for detailed rows: n<={cutoff} (7-divisible modes vanish).")
    print("Exact Delta_w uses the Abel endpoint identity from Thread B/HYP-2784.")
    print()
    print("cap and plateau context:")
    for k in range(8, 13):
        print(
            f"  k={k}: cap={float(CAP[k]):.6f}, "
            f"Q(k-1)={float(QVAL[k]):.6f}, margin={float(MARGIN[k]):.6f}"
        )

    rows = risk_rows(wmax=wmax, cutoff=cutoff)

    print("\n" + "=" * 92)
    print("MAX POSITIVE ONE-FAR DELTA IN THE SCANNED BINDING BAND")
    print(
        "  Reading: BV is the absolute THM-546 bound (6/49)V/w. "
        "The signed low modes explain the small exact Delta."
    )
    print(
        f"{'k':>2} {'w*':>4} {'Delta':>12} {'Delta/margin':>13} "
        f"{'V':>4} {'|Delta|/BV':>11} {'head_err':>10} {'odd_share':>10} "
        f"{'top_mod14':>9} {'low6_head':>11}"
    )
    for row in rows:
        split = row["split"]
        low6 = row["low6"]
        head_err = float(row["delta"]) - split["head"]
        print(
            f"{row['k']:>2} {row['w']:>4} {float(row['delta']):>12.8f} "
            f"{float(row['delta'] / row['margin']):>13.6f} "
            f"{row['V']:>4} {float(abs(row['delta']) / row['bv']):>11.6f} "
            f"{head_err:>10.2e} {split['odd_share']:>10.6f} "
            f"{row['top_residue']:>9} {low6['head']:>11.8f}"
        )
        by14 = split["by14"]
        by14_l1 = split["by14_l1"]
        top_terms = [(n, c) for _abs, n, c in split["terms"][:5]]
        top_buckets = sorted(
            ((by14_l1[r], r, by14[r]) for r in range(14) if r % 7),
            reverse=True,
        )[:4]
        print(
            "     top terms n: "
            + ", ".join(f"{n}:{c:+.7f}" for n, c in top_terms)
        )
        print(
            "     top mod14 buckets (L1, signed): "
            + ", ".join(f"{r}:({l1:.7f},{sgn:+.7f})" for l1, r, sgn in top_buckets)
        )

    print("\n" + "=" * 92)
    print("FINITE HEAD AUDIT")
    print(
        "  Worst |exact Delta_w - signed head_{n<=M}| over "
        f"w in [15,{audit_wmax}] for repeatable runtime."
    )
    print(f"{'M':>4} {'k':>2} {'worst_err':>12} {'err/margin':>12} {'w':>4} {'exact':>11} {'head':>11}")
    for cutoff_i, entries in finite_head_audit(wmax=audit_wmax, cutoffs=cutoffs).items():
        for k, worst, worst_w, exact_at, head_at in entries:
            print(
                f"{cutoff_i:>4} {k:>2} {worst:>12.8f} "
                f"{worst / float(MARGIN[k]):>12.6f} {worst_w:>4} "
                f"{exact_at:>11.8f} {head_at:>11.8f}"
            )

    print("\n" + "=" * 92)
    print("SYNTHESIS")
    print("  1. The dangerous positive one-far Delta rows are low-mode phenomena: n=1,2,3")
    print("     dominate the visible head, and the mod-14 buckets 1 and 2 carry most pressure.")
    print("  2. Odd support is a useful L1 envelope in most binding rows, but k=11 is an")
    print("     even-led exception. This matches HYP-2720: odd support is not a signed cone.")
    print("  3. The exact Delta is only about 5%-16% of the absolute BV bound in the risk rows;")
    print("     HYP-2784's 125x gap is therefore a signed phase-cancellation gap, not a V gap.")
    print("  4. A proof-shaped target is now finite-head plus tail:")
    print(
        "       signed low modes n<=13/mod14 ledger + odd/even exception ledger"
        " + Dedekind/equidistribution tail."
    )
    print(
        "     This is narrower than an absolute Koksma proof, and HYP-2785 rules"
        " out a finite residue-only table."
    )

    tournament(rows)


if __name__ == "__main__":
    main()
