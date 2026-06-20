#!/usr/bin/env python3
"""HYP-2674 scout: uniform Delta tail as same-sign packet alignment.

After HYP-2653d, the far-element target is no longer a bounded w*Delta_w
constant.  The target is a uniform cap on Delta_w itself after a finite
max(E') cutoff.  This script asks what the corrected Delta_w object looks like
inside the dyadic family that all recent routes keep rediscovering.

Main observation: the dangerous finite rows are not hidden by cancellation.
Their one-missed-sector packets are all same-sign.  The proof target becomes:
classify the finite same-sign packet alignments, then prove the remaining tail
has enough sign changes or small packet mass to stay below cap_k-Q(k-1).
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from math import gcd


F = Fraction
CAP9_MINUS_Q8 = F(129643, 980980)


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def missed_dist(E: tuple[int, ...]) -> list[F]:
    E = tuple(sorted(set(E)))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(bps)
    p = [F(0)] * 7
    for lo, hi in zip(bps, bps[1:]):
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in E:
            v = frac_part(e * mid)
            hit.add((v.numerator * 7) // v.denominator)
        missed = sum(1 for j in range(1, 7) if j not in hit)
        p[missed] += hi - lo
    return p


def phi(E: tuple[int, ...]) -> F:
    p = missed_dist(E)
    return p[0] + F(1, 7) * p[1]


def p0(E: tuple[int, ...]) -> F:
    return missed_dist(E)[0]


def g0(y: F) -> F:
    y = frac_part(y)
    return y * F(6, 7) if y < F(1, 7) else F(6, 49) - (y - F(1, 7)) * F(1, 7)


def one_missed_cells(Ep: tuple[int, ...]) -> list[tuple[F, F, int]]:
    Ep = tuple(sorted(set(Ep)))
    bps = {F(0), F(1)}
    for e in Ep:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(bps)
    cells = []
    for lo, hi in zip(bps, bps[1:]):
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in Ep:
            v = frac_part(e * mid)
            hit.add((v.numerator * 7) // v.denominator)
        missed = [j for j in range(1, 7) if j not in hit]
        if len(missed) == 1:
            cells.append((lo, hi, missed[0]))
    return cells


def packet_delta(Ep: tuple[int, ...], w: int) -> tuple[F, dict[int, F], int]:
    packets = {s: F(0) for s in range(1, 7)}
    cells = one_missed_cells(Ep)
    for a, b, s in cells:
        packets[s] += g0(w * b - F(s, 7)) - g0(w * a - F(s, 7))
    return sum(packets.values()) / w, packets, len(cells)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.6f})"


def sign_word(packets: dict[int, F]) -> str:
    out = []
    for s in range(1, 7):
        v = packets[s]
        out.append("+" if v > 0 else "-" if v < 0 else "0")
    return "".join(out)


def summarize_row(name: str, Ep: tuple[int, ...], w: int) -> tuple[F, str]:
    delta_direct = p0(tuple(sorted(Ep + (w,)))) - phi(Ep)
    delta_packet, packets, n_cells = packet_delta(Ep, w)
    if delta_direct != delta_packet:
        raise AssertionError((name, delta_direct, delta_packet))
    p1 = missed_dist(Ep)[1]
    pos = sum(v for v in packets.values() if v > 0)
    neg = -sum(v for v in packets.values() if v < 0)
    print(name)
    print(f"  E'={Ep}, w={w}, max(E')={max(Ep)}")
    print(f"  Delta_w                    = {fmt(delta_direct)}")
    print(f"  margin_9 - Delta_w         = {fmt(CAP9_MINUS_Q8 - delta_direct)}")
    print(f"  p1(E')                     = {fmt(p1)}")
    print(f"  Delta_w/p1                 = {fmt(delta_direct / p1) if p1 else 'undefined'}")
    print(f"  w*Delta_w diagnostic       = {fmt(delta_direct * w)}")
    print(f"  packet sign word s=1..6    = {sign_word(packets)}")
    print(f"  packet positive/negative   = +{fmt(pos)} / -{fmt(neg)}")
    print(f"  one-missed cells           = {n_cells}")
    print("  packets:")
    for s in range(1, 7):
        print(f"    missed sector {s}: {fmt(packets[s])}")
    print()
    return delta_direct, sign_word(packets)


def dyadic_family(mult: int, s_max: int = 120) -> dict[str, tuple[F, int]]:
    best = (F(-10), -1)
    best_after_20 = (F(-10), -1)
    best_after_40 = (F(-10), -1)
    all_positive_best = None
    for s in range(3, s_max + 1):
        Ep = (0, 1, 2, 4, 8, 3 * s, 4 * s, 5 * s)
        w = mult * s
        delta, packets, _ = packet_delta(Ep, w)
        if delta > best[0]:
            best = (delta, s)
        if s > 20 and delta > best_after_20[0]:
            best_after_20 = (delta, s)
        if s > 40 and delta > best_after_40[0]:
            best_after_40 = (delta, s)
        if sign_word(packets) == "++++++":
            if all_positive_best is None or delta > all_positive_best[0]:
                all_positive_best = (delta, s)
    assert all_positive_best is not None
    return {
        "best": best,
        "best_after_20": best_after_20,
        "best_after_40": best_after_40,
        "all_positive_best": all_positive_best,
    }


def print_family_scan() -> None:
    print("dyadic family scan")
    print("  E_s={0,1,2,4,8,3s,4s,5s}; sampled exactly for s=3..120")
    for mult in (6, 10):
        data = dyadic_family(mult)
        print(f"  w={mult}s")
        for key, (delta, s) in data.items():
            print(
                f"    {key:17s}: s={s:3d}, Delta={fmt(delta)}, "
                f"margin_9-Delta={fmt(CAP9_MINUS_Q8-delta)}"
            )
    print()


def print_tournament_summary(rows: list[tuple[str, F]]) -> None:
    print("Tournament Analysis")
    print("  vertices: finite_B13_packet, dyadic_s4_alignment, nonshell_warning,")
    print("            dyadic_w6_tail_after20, dyadic_w10_tail_after20")
    print("  observable: larger positive Delta_w points to the riskier proof obligation.")
    ordered = sorted(rows, key=lambda item: item[1], reverse=True)
    scores = {name: len(rows) - 1 - rank for rank, (name, _delta) in enumerate(ordered)}
    print("  Hamiltonian path by positive Delta_w:")
    print("    " + " > ".join(name for name, _delta in ordered))
    print("  score histogram:", sorted(scores.values()))
    print("  directed 3-cycles: 0 (risk order is transitive in this scout)")
    print("  challenged assumption: the uniform tail is a generic discrepancy problem.")
    print("    The risky rows are same-sign packet alignments; cancellation only needs")
    print("    to be proved after the finite alignment pocket is isolated.")
    print()


def main() -> None:
    print("HYP-2674 LRC14 uniform Delta packet-alignment scout")
    print(f"k=9 far-tail margin cap_9-Q(8) = {fmt(CAP9_MINUS_Q8)}")
    print()

    rows = []
    for name, Ep, w in [
        ("finite_B13_packet", (0, 1, 2, 4, 6, 7, 8, 10), 12),
        ("dyadic_s4_alignment", (0, 1, 2, 4, 8, 12, 16, 20), 24),
        ("nonshell_warning", (0, 2, 3, 5, 6, 15), 18),
    ]:
        delta, _word = summarize_row(name, Ep, w)
        rows.append((name, delta))

    print_family_scan()
    fam6 = dyadic_family(6)
    fam10 = dyadic_family(10)
    rows.extend(
        [
            ("dyadic_w6_tail_after20", fam6["best_after_20"][0]),
            ("dyadic_w10_tail_after20", fam10["best_after_20"][0]),
        ]
    )
    print_tournament_summary(rows)
    print("PASS: uniform Delta packet-alignment scout complete.")


if __name__ == "__main__":
    main()
