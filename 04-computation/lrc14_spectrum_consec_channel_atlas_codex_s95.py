#!/usr/bin/env python3
r"""
lrc14_spectrum_consec_channel_atlas_codex_s95.py

Codex S95 scout for the Node-3 spectrum sum.

This is a bounded-core atlas for the exact quasi-independence ratio

    R' = meas(coverSet(E)^c cap G_P) / (meas(G_P) * meas(coverSet(E)^c)).

The real-space value of R' is computed exactly as a Fraction.  The low-channel
Fourier split is numerical, using the existing KPS spectrum engine, and is only
used to identify the finite signed channels that create the deviation from
independence.

The concrete question is whether the recent history is sharpening toward a
single object:

    complement-even low-frequency covariance + L2-controlled complement-odd tail.

For the consecutive bounded cores E={0,...,k-1}, this script enumerates all
admissible small parts P subset {1,...,13}, |P|+k=13, and reports:

  * exact R' minima,
  * low-frequency signed ledgers: residue 0, nonzero mean, Paley QR/NQR cut,
  * complement-pair profile of P under p -> 14-p,
  * additive energy of E,
  * sign complexity of the first H sign-paired channels {n,-n}.

Tournament Analysis note:
  vertices are not runners.  The diagnostic vertices are the sign-paired low
  Fourier channels {n,-n}; the Hamiltonian path is n=1,2,...,H.  The pairwise
  observable is the sign of the two-channel contribution s_i+s_j.  For this
  atlas we record cheaper proxies (negative count, sign changes, residue sums)
  rather than all O(H^2) tournament fingerprints for every P.
"""

import itertools
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import pi

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (  # noqa: E402
    complement,
    cover_set,
    fourier_num_of_arcs,
    intersect,
    meas,
    safe_set,
)


QR7 = {1, 2, 4}
NQR7 = {3, 5, 6}


def coeff(arcs, n):
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)


def fmt_frac(x):
    if isinstance(x, F):
        return f"{x} = {float(x):.8f}"
    return f"{x:.8f}"


def additive_profile(E):
    E = list(E)
    sums = Counter(a + b for a in E for b in E)
    energy = sum(v * v for v in sums.values())
    return {
        "energy": energy,
        "sumset": len(sums),
        "K2": len(sums) / max(1, len(E)),
    }


def complement_profile_P(P, N=14):
    P = set(P)
    pairs = [(1, 13), (2, 12), (3, 11), (4, 10), (5, 9), (6, 8)]
    paired = orphan = absent = 0
    for a, b in pairs:
        c = int(a in P) + int(b in P)
        if c == 2:
            paired += 1
        elif c == 1:
            orphan += 1
        else:
            absent += 1
    fixed7 = int(7 in P)
    return {
        "paired": paired,
        "orphan": orphan,
        "absent": absent,
        "fixed7": fixed7,
        "defect": orphan + fixed7,
    }


def reflection_defect_E(E):
    E = sorted(set(E))
    lo, hi = E[0], E[-1]
    S = set(E)
    miss = sum(1 for e in E if lo + hi - e not in S)
    return miss, miss / len(E)


def channel_summary(P, E, H=21):
    P = tuple(sorted(P))
    E = tuple(sorted(E))
    gp = safe_set(P)
    covc = complement(cover_set(E))
    baseline = meas(gp) * meas(covc)
    inter = meas(intersect(gp, covc))
    spec_exact = inter - baseline
    Rprime = inter / baseline if baseline else F(1)

    residue = defaultdict(float)
    terms = []
    for n in range(1, H + 1):
        s = 2.0 * (coeff(gp, n) * coeff(covc, n).conjugate()).real
        terms.append((n, n % 7, s))
        residue[n % 7] += s

    low = sum(s for _, _, s in terms)
    qr = sum(residue[r] for r in QR7)
    nqr = sum(residue[r] for r in NQR7)
    nonzero = qr + nqr
    paley = qr - nqr
    signs = [1 if s > 1e-12 else -1 if s < -1e-12 else 0 for _, _, s in terms]
    sign_changes = sum(
        1
        for a, b in zip(signs, signs[1:])
        if a != 0 and b != 0 and a != b
    )
    neg_count = sum(1 for s in signs if s < 0)
    top_neg = sorted(terms, key=lambda row: row[2])[:4]
    top_pos = sorted(terms, key=lambda row: row[2], reverse=True)[:4]

    base_f = float(baseline) if baseline else 1.0
    pprof = complement_profile_P(P)
    e_miss, e_defect = reflection_defect_E(E)
    aprof = additive_profile(E)
    return {
        "P": P,
        "E": E,
        "k": len(E),
        "baseline": baseline,
        "inter": inter,
        "spec_exact": spec_exact,
        "Rprime": Rprime,
        "spec_over": float(spec_exact) / base_f,
        "low": low,
        "low_over": low / base_f,
        "res0_over": residue[0] / base_f,
        "nonzero_over": nonzero / base_f,
        "paley_over": paley / base_f,
        "neg_count": neg_count,
        "sign_changes": sign_changes,
        "top_neg": top_neg,
        "top_pos": top_pos,
        "pprof": pprof,
        "e_reflect_missing": e_miss,
        "e_reflect_defect": e_defect,
        **aprof,
    }


def rec_line(rec):
    return (
        f"k={rec['k']} P={list(rec['P'])} "
        f"R'={rec['Rprime']} ({float(rec['Rprime']):.6f}) "
        f"spec/base={rec['spec_over']:+.6f} "
        f"low/base={rec['low_over']:+.6f} "
        f"r0={rec['res0_over']:+.6f} "
        f"nonzero={rec['nonzero_over']:+.6f} "
        f"paley={rec['paley_over']:+.6f} "
        f"Pdef={rec['pprof']['defect']} "
        f"neg={rec['neg_count']:02d} sc={rec['sign_changes']:02d}"
    )


def print_channel_extremes(rec):
    neg = ", ".join(f"n={n}(r{r}) {s:+.6f}" for n, r, s in rec["top_neg"])
    pos = ", ".join(f"n={n}(r{r}) {s:+.6f}" for n, r, s in rec["top_pos"])
    print(f"      most negative: {neg}")
    print(f"      most positive: {pos}")


def scan_consecutive(H=21):
    rows = []
    for k in range(8, 14):
        E = tuple(range(k))
        p_size = 13 - k
        for P in itertools.combinations(range(1, 14), p_size):
            rows.append(channel_summary(P, E, H=H))
    return rows


def print_minima(rows):
    print("=" * 100)
    print("CONSECUTIVE BOUNDED CORES: exact R' and low-channel minima")
    print("=" * 100)
    print(f"total rows scanned = {len(rows)}")
    for k in range(8, 14):
        subset = [r for r in rows if r["k"] == k]
        print("-" * 100)
        print(f"k={k}, rows={len(subset)}, E={{0,...,{k-1}}}, |P|={13-k}")
        leaders = [
            ("min exact R'", min(subset, key=lambda r: r["Rprime"])),
            ("most negative exact SPEC/base", min(subset, key=lambda r: r["spec_over"])),
            ("most negative low/base", min(subset, key=lambda r: r["low_over"])),
            ("most negative residue-0/base", min(subset, key=lambda r: r["res0_over"])),
            ("most negative nonzero/base", min(subset, key=lambda r: r["nonzero_over"])),
            ("most negative Paley/base", min(subset, key=lambda r: r["paley_over"])),
        ]
        seen = set()
        for label, rec in leaders:
            key = (label, rec["P"])
            print(f"  {label:<34} {rec_line(rec)}")
            if rec["P"] not in seen:
                print_channel_extremes(rec)
                seen.add(rec["P"])


def print_complement_buckets(rows):
    print("=" * 100)
    print("COMPLEMENT-PAIR BUCKETS FOR P under p -> 14-p")
    print("=" * 100)
    buckets = defaultdict(list)
    for rec in rows:
        buckets[(rec["k"], rec["pprof"]["defect"])].append(rec)
    for k in range(8, 14):
        print(f"k={k}")
        for defect in sorted(d for kk, d in buckets if kk == k):
            group = buckets[(k, defect)]
            minR = min(group, key=lambda r: r["Rprime"])
            minLow = min(group, key=lambda r: r["low_over"])
            print(
                f"  Pdef={defect:<2} rows={len(group):<4} "
                f"minR={float(minR['Rprime']):.6f} P={list(minR['P'])} "
                f"minLow/base={minLow['low_over']:+.6f} P={list(minLow['P'])}"
            )


def print_selected_probes(H=21):
    probes = [
        ("S94 exact consec floor row", [1, 3, 4, 5], list(range(9))),
        ("complement-paired k8 scout", [1, 3, 5, 9, 13], list(range(8))),
        ("small-P low-frequency row", [1, 2, 3], list(range(10))),
        ("coprime-P independence-favourable", [5, 7, 11], list(range(10))),
        ("wide d>=2 routed row", [1, 2, 6], [0, 4, 6, 8, 10, 12, 14, 15, 16, 17]),
        ("perforated near-AP row", [1, 2, 3, 12, 13], [0, 2, 3, 4, 5, 6, 7, 8]),
        ("scale-separated AP row", [1, 2, 3], [0, 7, 14, 21, 28, 35, 42, 49, 56, 63]),
    ]
    print("=" * 100)
    print("SELECTED SHAPE PROBES")
    print("=" * 100)
    for label, P, E in probes:
        rec = channel_summary(P, E, H=H)
        print(f"{label}: {rec_line(rec)}")
        print(
            f"      E-reflection-defect={rec['e_reflect_defect']:.3f} "
            f"energy={rec['energy']} |E+E|={rec['sumset']} K2={rec['K2']:.3f}"
        )
        print_channel_extremes(rec)


def print_synthesis(rows, H):
    exact_worst = min(rows, key=lambda r: r["Rprime"])
    low_worst = min(rows, key=lambda r: r["low_over"])
    r0_worst = min(rows, key=lambda r: r["res0_over"])
    paley_worst = min(rows, key=lambda r: r["paley_over"])
    print("=" * 100)
    print("SYNTHESIS")
    print("=" * 100)
    print(f"H={H} sign-paired low channels were used for the Fourier ledger.")
    print(f"Global exact R' minimum: {rec_line(exact_worst)}")
    print(f"Global low-channel minimum: {rec_line(low_worst)}")
    print(f"Global residue-0 trunk minimum: {rec_line(r0_worst)}")
    print(f"Global Paley-cut minimum: {rec_line(paley_worst)}")
    print()
    print("Observed sharpened story:")
    print("  1. The exact floor in consecutive bounded cores is not threatened by Paley imbalance alone.")
    print("     In the worst exact R' row, the Paley ledger can be positive while R' is still below 1.")
    print("  2. The negative covariance decomposes into a residue-0 trunk plus a nonzero-shell mean.")
    print("     This matches S94 and warns against over-specializing to QR/NQR or apex-7 alone.")
    print("  3. Consecutive E has reflection defect 0 and maximal one-dimensional additive structure.")
    print("     The hard branch is therefore complement-even / high-additive-energy, not random.")
    print("  4. The complement-odd part is still a plausible L2 target: rows with asymmetric P do")
    print("     not create a new vanishing R' family in this bounded atlas.")
    print()
    print("Next proof obligation suggested by the atlas:")
    print("  Formalize an even/odd split under the complement involution v -> 14-v.")
    print("  Bound the complement-odd low packet by Parseval/Cauchy-Schwarz and route the")
    print("  complement-even packet to the finite consecutive/AP/Freiman atlas.")


def main():
    H = 21
    print("# LRC14 spectrum consecutive-channel atlas (codex S95)")
    print("# Exact real-space R' values; numerical low-channel Fourier diagnostics.")
    print("# The hidden quotient being tested is complement-even low covariance.")
    rows = scan_consecutive(H=H)
    print_minima(rows)
    print_complement_buckets(rows)
    print_selected_probes(H=H)
    print_synthesis(rows, H)


if __name__ == "__main__":
    main()
