#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Index-alignment audit for the genuine-wide doublet family.

Incoming HYP-2794 correctly identifies the actual-size doublet family

    E_N(M) = {0,1,...,N-3} union {M,M+1},

but the companion signed-bound script prints labels ``k`` for rows built as
``range(k-1) union {M,M+1}``, which have size ``k+1``.  This script records the
exact cap/Q alignment so follow-up finite checks use the correct row size.
"""
from __future__ import annotations

import functools
import sys

sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
print = functools.partial(print, flush=True)

from lrc14_threadA_regime_dichotomy_kpswf8 import CAP, QVAL, p0_fast


def actual_doublet(size: int, M: int) -> tuple[int, ...]:
    return tuple(list(range(size - 2)) + [M, M + 1])


def printed_label_doublet(k_label: int, M: int) -> tuple[int, ...]:
    return tuple(list(range(k_label - 1)) + [M, M + 1])


def scan_family(size: int, lo: int = 15, hi: int = 120):
    vals = {M: p0_fast(actual_doublet(size, M)) for M in range(lo, hi + 1)}
    M_star = max(vals, key=vals.get)
    return M_star, vals[M_star]


def main() -> None:
    print("=" * 88)
    print("LRC14 doublet index alignment audit (codex-S77)")
    print("=" * 88)
    print("Actual row-size convention: E_N(M)={0,...,N-3} union {M,M+1}.")
    print()
    print("Actual-size doublet maxima over 15<=M<=120:")
    print(f"{'N':>3} {'|E|':>3} {'M*':>4} {'p0':>18} {'cap_N-p0':>14} {'Q(N-1)-p0':>14}")
    for N in range(8, 13):
        M_star, val = scan_family(N)
        print(
            f"{N:>3} {len(actual_doublet(N, M_star)):>3} {M_star:>4} "
            f"{str(val):>18} {str(CAP[N]-val):>14} {str(QVAL[N]-val):>14}"
        )
    print()
    print("Incoming signed-bound script label audit:")
    print("  It builds range(k-1) union {M,M+1}, which has actual size k+1.")
    print(f"{'label k':>7} {'actual N':>8} {'M*':>4} {'p0':>18} {'cap_k-p0':>14} {'cap_N-p0':>14}")
    for k_label in range(8, 13):
        vals = {
            M: p0_fast(printed_label_doublet(k_label, M))
            for M in range(15, 121)
        }
        M_star = max(vals, key=vals.get)
        val = vals[M_star]
        N = len(printed_label_doublet(k_label, M_star))
        cap_actual = CAP.get(N)
        actual_gap = cap_actual - val if cap_actual is not None else "n/a"
        print(
            f"{k_label:>7} {N:>8} {M_star:>4} {str(val):>18} "
            f"{str(CAP[k_label]-val):>14} {str(actual_gap):>14}"
        )
    print()
    print("Verdict: use actual size N for cap_N and Q(N-1).")
    print("The cap margins only get larger after correcting the shifted labels for N=9..12.")


if __name__ == "__main__":
    main()
