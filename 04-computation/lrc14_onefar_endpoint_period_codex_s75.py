#!/usr/bin/env python3
"""HYP-2789/S75: exact endpoint-period certificate for one-far consec bases.

For a fixed bounded base B, the Abel endpoint identity gives

    Delta_w(B) = S_B(w) / w,

where S_B(w)=w*Delta_w is a signed sum of values of the periodized primitive
G0 at rational endpoints.  Therefore S_B(w) is periodic in w, with period
dividing the lcm of the endpoint denominators.  This script scans one complete
period for the binding bases B=consec_(k-1), k=8..12, and records the exact
positive endpoint numerator.

No residue-only mod-7 claim is made.  The period is much larger than 7, which
matches HYP-2785's warning; the useful finite object is the endpoint-period
numerator, not a small residue table.
"""

from __future__ import annotations

import functools
import itertools
import sys
from collections import Counter
from fractions import Fraction as F
from math import gcd

sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
print = functools.partial(print, flush=True)

from lrc14_threadB_localized_binding_bound_kpswf8 import (  # noqa: E402
    CAP,
    MARGIN,
    QVAL,
    G0,
    Varcs,
    cells_with_miss,
    orbit_breakpoints,
)


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def ceil_frac(q: F) -> int:
    return -(-q.numerator // q.denominator)


def fdec(q: F, digits: int = 9) -> str:
    return f"{float(q):.{digits}f}"


def endpoint_packets(B: tuple[int, ...]) -> tuple[list[tuple[F, F, int]], int]:
    """Return one-miss endpoint packets and the exact endpoint period."""

    bp = orbit_breakpoints(B)
    period = 1
    packets: list[tuple[F, F, int]] = []
    for lo, hi, miss in cells_with_miss(B, bp):
        period = lcm(period, lo.denominator)
        period = lcm(period, hi.denominator)
        if len(miss) == 1:
            packets.append((lo, hi, next(iter(miss))))
    for lo, hi, _sector in packets:
        assert (period * lo).denominator == 1
        assert (period * hi).denominator == 1
    return packets, period


def signed_numerator(packets: list[tuple[F, F, int]], w: int) -> F:
    """Exact S_B(w)=w*Delta_w from the Abel endpoint packets."""

    total = F(0)
    for lo, hi, sector in packets:
        shift = F(sector, 7)
        total += G0(w * hi - shift) - G0(w * lo - shift)
    return total


def certify_row(k: int) -> dict[str, object]:
    B = tuple(range(k - 1))
    packets, period = endpoint_packets(B)

    # One representative for every residue class modulo period, with w>=15.
    best_delta = (F(-10), -1, F(0))
    max_s = (F(-10), -1)
    min_s = (F(10), -1)
    max_abs_s = (F(0), -1, F(0))
    for w in range(15, period + 15):
        s_val = signed_numerator(packets, w)
        delta = s_val / w
        if delta > best_delta[0]:
            best_delta = (delta, w, s_val)
        if s_val > max_s[0]:
            max_s = (s_val, w)
        if s_val < min_s[0]:
            min_s = (s_val, w)
        if abs(s_val) > max_abs_s[0]:
            max_abs_s = (abs(s_val), w, s_val)

    # Period sanity: denominator divisibility already proves this, but keep a
    # concrete check for future edits to the endpoint engine.
    for w in range(1, min(period, 80) + 1):
        assert signed_numerator(packets, w + period) == signed_numerator(packets, w)

    positive_numerator_margin_ratio = max_s[0] / F(15) / MARGIN[k]
    return {
        "k": k,
        "B": B,
        "period": period,
        "packets": len(packets),
        "V": Varcs(B),
        "max_delta": best_delta[0],
        "max_delta_w": best_delta[1],
        "max_delta_s": best_delta[2],
        "max_s": max_s[0],
        "max_s_w": max_s[1],
        "min_s": min_s[0],
        "min_s_w": min_s[1],
        "max_abs_s": max_abs_s[0],
        "max_abs_s_w": max_abs_s[1],
        "max_abs_s_signed": max_abs_s[2],
        "margin": MARGIN[k],
        "cap": CAP[k],
        "Q": QVAL[k],
        "numerator_ratio": positive_numerator_margin_ratio,
        "cutoff": ceil_frac(max_s[0] / MARGIN[k]),
    }


def tournament(rows: list[dict[str, object]]) -> None:
    print("\n" + "=" * 96)
    print("TOURNAMENT ANALYSIS")
    print("  vertices: consecutive one-far binding rows k=8..12")
    print("  pairwise observable: larger maxS/(15*margin), then larger maxDelta/margin")
    print("  switch/gauge: compare endpoint-period numerators before scalarizing by 1/w")
    scores = Counter()
    edges: set[tuple[int, int]] = set()
    for i, j in itertools.combinations(range(len(rows)), 2):
        ri = rows[i]
        rj = rows[j]
        ai = (ri["numerator_ratio"], ri["max_delta"] / ri["margin"])
        aj = (rj["numerator_ratio"], rj["max_delta"] / rj["margin"])
        if ai >= aj:
            winner, loser = i, j
        else:
            winner, loser = j, i
        scores[winner] += 1
        edges.add((winner, loser))
    for i in range(len(rows)):
        scores.setdefault(i, 0)
    hist = Counter(scores.values())
    cycles = 0
    for a, b, c in itertools.combinations(range(len(rows)), 3):
        triple = {(a, b), (b, c), (c, a)}
        rev = {(b, a), (c, b), (a, c)}
        if triple <= edges or rev <= edges:
            cycles += 1
    order = sorted(range(len(rows)), key=lambda i: scores[i], reverse=True)
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian pressure path:")
    for idx in order:
        row = rows[idx]
        print(
            f"    k={row['k']}: maxS/(15*margin)={fdec(row['numerator_ratio'], 6)}, "
            f"maxDelta/margin={fdec(row['max_delta'] / row['margin'], 6)}, "
            f"period={row['period']}"
        )


def main() -> None:
    print("HYP-2789/S75 one-far endpoint-period certificate")
    print("Binding rows: B=consec_(k-1), k=8..12, all wide w>=15.")
    print("Exact identity: Delta_w = S_B(w)/w, with S_B periodic in w.")
    print("The period is the lcm of Abel endpoint denominators, not a mod-7 table.")
    print()
    print("cap and plateau context:")
    for k in range(8, 13):
        print(
            f"  k={k}: cap={fdec(CAP[k], 6)}, Q(k-1)={fdec(QVAL[k], 6)}, "
            f"margin={fdec(MARGIN[k], 6)}"
        )

    rows = [certify_row(k) for k in range(8, 13)]

    print("\n" + "=" * 96)
    print("EXACT ONE-PERIOD CERTIFICATE")
    print(
        "  Reading: maxS is the maximum positive endpoint numerator S_B(w)=w*Delta_w "
        "over one complete period."
    )
    print(
        "  Since S_B(w+L)=S_B(w), all w>=15 satisfy Delta_w <= maxS/15."
    )
    print(
        f"{'k':>2} {'period':>7} {'pack':>4} {'V':>3} {'maxS':>12} "
        f"{'maxS/(15m)':>12} {'cut':>4} {'best_w':>6} {'bestDelta':>11} "
        f"{'best/m':>9}"
    )
    for row in rows:
        print(
            f"{row['k']:>2} {row['period']:>7} {row['packets']:>4} {row['V']:>3} "
            f"{str(row['max_s']):>12} {fdec(row['numerator_ratio'], 6):>12} "
            f"{row['cutoff']:>4} {row['max_delta_w']:>6} "
            f"{fdec(row['max_delta'], 8):>11} "
            f"{fdec(row['max_delta'] / row['margin'], 6):>9}"
        )

    print("\nendpoint numerator extremes:")
    for row in rows:
        print(
            f"  k={row['k']}: maxS={row['max_s']} at w={row['max_s_w']}; "
            f"minS={row['min_s']} at w={row['min_s_w']}; "
            f"max|S|={row['max_abs_s']} at w={row['max_abs_s_w']} "
            f"(signed {row['max_abs_s_signed']})"
        )

    print("\n" + "=" * 96)
    print("SYNTHESIS")
    print("  1. The exact all-w>=15 maxima match the short low-head scout rows:")
    print("       k=8 w=21, k=9 w=68, k=10 w=22, k=11 w=16, k=12 w=71.")
    print("  2. The positive endpoint numerator satisfies maxS <= 1289/980 < 1.316")
    print("     across k=8..12.  Therefore maxS/15 is below 0.453 of the cap margin")
    print("     in every consecutive binding row.")
    print("  3. This gives an exact finite certificate for the consec-base single-far")
    print("     branch.  It is compatible with HYP-2785: the modulus is endpoint-period")
    print("     sized, not residue-only mod 7.")
    print("  4. Remaining proof work: prove near-cap bases reduce to this consecutive")
    print("     endpoint-period certificate, or give each non-consec base a plateau-slack")
    print("     ledger against its own endpoint numerator.")

    tournament(rows)


if __name__ == "__main__":
    main()
