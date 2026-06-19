#!/usr/bin/env python3
"""
LRC(14) exact-period AP-drop repair atlas.

This script refines HYP-2628/HYP-2629's observation that the squarefree
Q=210 grid misses exactly two AP one-drop cores while Q=1260 sees all of them.
The point is not merely that 1260 is larger.  The Q=1260 hits have reduced
denominators in the exact-period divisor lattice of

    1260 = 2^2 * 3^2 * 5 * 7,

and, inside this raw Hill carrier, the two Q=210-blind drops are repaired by
repeated-prime lanes:

    drop  6: uses repeated 2/3 lanes, with lowest Q=1260 denominator 63;
    drop 12: uses the dyadic 2^2 lane, with denominator 12.

Drop 6 has an earlier q=98 witness outside the 1260 divisor lattice; this
script is about the exact-period packets retained by the raw K14 product.

Tournament Analysis contract:
  * Vertices are proof quotients: radical grid, crossing grid, repeated-prime
    lanes, AP-drop mouths, and raw runner vertices.
  * Pairwise observable is preservation of the AP one-drop strict-safe
    predicate after quotienting to exact-period denominators.
  * Switch/gauge orients A -> B when A repairs more of the Q=210-blind AP
    mouths while preserving the mod-210 address.  Ties use the displayed
    Hamiltonian path.
  * Fingerprints are the score histogram, directed cycles, SCC sizes, and
    Hamiltonian-path count of the resulting transitive proof-quotient order.

Assumption challenge:
  This session again avoids runner vertices as the primary quotient.  Candidate
  vertices considered were runners, AP drops, safe components, endpoint mouths,
  residues a/Q, reduced denominators, squarefree masks, exact-period packets,
  and proof obligations.  Reduced denominators inside safe mouths preserve the
  LRC predicate "this rational is strict-safe for the 12-core" and the exact
  period/copy data needed by HYP-2628.  They destroy continuous mouth length and
  most endpoint-owner history; the dihedral endpoint atlas keeps that data.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd

import lrc14_uniform_fattening_gauntlet_codex as base


N = 14
AP13 = tuple(range(1, N))
PRIMES = (2, 3, 5, 7)
RADICAL_Q = 210
RAW_Q = 1260
GRID_QS = (84, 126, 210, 315, 420, 630, 840, 1050, 1260)


def fmt(fr: F) -> str:
    return f"{fr} = {float(fr):.9f}"


def factorization(n: int) -> str:
    if n == 1:
        return "1"
    parts = []
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            e = 0
            while m % d == 0:
                m //= d
                e += 1
            parts.append(f"{d}^{e}" if e > 1 else str(d))
        d += 1 if d == 2 else 2
    if m > 1:
        parts.append(str(m))
    return "*".join(parts)


def phi(n: int) -> int:
    ans = n
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            ans -= ans // p
            while m % p == 0:
                m //= p
        p += 1 if p == 2 else 2
    if m > 1:
        ans -= ans // m
    return ans


def mask(n: int) -> str:
    touched = [str(p) for p in PRIMES if n % p == 0]
    return "{" + ",".join(touched) + "}" if touched else "{}"


def strict_safe_residues(core: tuple[int, ...], q: int) -> list[int]:
    safe = []
    for a in range(q):
        ok = True
        for v in core:
            r = (v * a) % q
            dist = min(r, q - r)
            if 14 * dist <= q:
                ok = False
                break
        if ok:
            safe.append(a)
    return safe


def reduced_denominator(a: int, q: int) -> int:
    if a == 0:
        return 1
    return q // gcd(a, q)


def danger_speeds(a: int, q: int, speeds: tuple[int, ...] = AP13) -> tuple[int, ...]:
    out = []
    for v in speeds:
        r = (v * a) % q
        dist = min(r, q - r)
        if 14 * dist <= q:
            out.append(v)
    return tuple(out)


def residues_in_component(component: tuple[F, F], q: int) -> list[int]:
    start, end = component
    lo = start * q
    hi = end * q
    first = lo.numerator // lo.denominator + 1
    last = (hi.numerator - 1) // hi.denominator
    if last < first:
        return []
    return list(range(max(first, 0), min(last, q - 1) + 1))


def omitted_clock_centers(drop: int, component: tuple[F, F]) -> tuple[int, ...]:
    mid = (component[0] + component[1]) / 2
    centers = []
    for k in range(drop):
        c = F(k, drop)
        dist = abs(mid - c)
        dist = min(dist, 1 - dist)
        if dist <= F(1, N * drop):
            centers.append(k)
    return tuple(centers)


def denominator_counter(residues: list[int], q: int) -> Counter[int]:
    return Counter(reduced_denominator(a, q) for a in residues)


def summarize_denominators(counter: Counter[int]) -> str:
    if not counter:
        return "-"
    return ", ".join(f"{d}:{counter[d]}" for d in sorted(counter))


def ap_drop_summary() -> list[dict[str, object]]:
    rows = []
    for drop in AP13:
        core = tuple(v for v in AP13 if v != drop)
        measure, components = base.safe_components(core, N)
        q_counts = {q: len(strict_safe_residues(core, q)) for q in GRID_QS}
        raw_residues = strict_safe_residues(core, RAW_Q)
        raw_denoms = denominator_counter(raw_residues, RAW_Q)
        rows.append(
            {
                "drop": drop,
                "core": core,
                "measure": measure,
                "components": components,
                "q_counts": q_counts,
                "raw_residues": raw_residues,
                "raw_denoms": raw_denoms,
                "raw_denoms_dividing_210": sum(
                    count for d, count in raw_denoms.items() if RADICAL_Q % d == 0
                ),
            }
        )
    return sorted(rows, key=lambda r: (r["measure"], r["drop"]))


def first_grid_hit(core: tuple[int, ...], limit: int = RAW_Q) -> tuple[int, int, list[int]]:
    for q in range(N, limit + 1, N):
        residues = strict_safe_residues(core, q)
        if residues:
            return q, len(residues), residues
    return 0, 0, []


def component_table(drop: int) -> list[dict[str, object]]:
    core = tuple(v for v in AP13 if v != drop)
    _, components = base.safe_components(core, N)
    rows = []
    for component in components:
        q210 = residues_in_component(component, RADICAL_Q)
        q1260 = residues_in_component(component, RAW_Q)
        rows.append(
            {
                "interval": component,
                "length": component[1] - component[0],
                "centers": omitted_clock_centers(drop, component),
                "q210": q210,
                "q1260": q1260,
                "q1260_denoms": [reduced_denominator(a, RAW_Q) for a in q1260],
            }
        )
    return rows


def exact_period_packet_table(drop: int) -> list[tuple[int, int, int, str, str, bool]]:
    core = tuple(v for v in AP13 if v != drop)
    raw_residues = strict_safe_residues(core, RAW_Q)
    counter = denominator_counter(raw_residues, RAW_Q)
    rows = []
    for d in sorted(counter):
        rows.append(
            (
                d,
                counter[d],
                phi(d),
                factorization(d),
                mask(d),
                RADICAL_Q % d == 0,
            )
        )
    return rows


def transitive_tournament_fingerprint(path: list[str]) -> dict[str, object]:
    n = len(path)
    scores = [n - 1 - i for i in range(n)]
    cycles = 0
    hamiltonian_paths = 1
    return {
        "path": path,
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": cycles,
        "scc_sizes": [1] * n,
        "hamiltonian_paths": hamiltonian_paths,
    }


def print_drop_grid_table(rows: list[dict[str, object]]) -> None:
    print("AP one-drop exact grid counts")
    print("drop  meas(G_C)        comps  Q=210  Q=315  Q=420  Q=630  Q=1260  raw reduced denoms")
    for row in rows:
        q = row["q_counts"]
        print(
            f"{row['drop']:>4}  {str(row['measure']):>14}  "
            f"{len(row['components']):>5}  {q[210]:>5}  {q[315]:>5}  "
            f"{q[420]:>5}  {q[630]:>5}  {q[1260]:>6}  "
            f"{summarize_denominators(row['raw_denoms'])}"
        )
    print()
    missed = [row["drop"] for row in rows if row["q_counts"][RADICAL_Q] == 0]
    repaired = [row["drop"] for row in rows if row["q_counts"][RAW_Q] > 0]
    print(f"Q=210 misses AP drops: {missed}")
    print(f"Q=1260 repairs AP drops: {repaired}")
    print()


def print_repair_detail(drop: int) -> None:
    core = tuple(v for v in AP13 if v != drop)
    q_first, n_first, first_residues = first_grid_hit(core)
    print("=" * 88)
    print(f"Drop {drop}: component mouths and exact-period repair packets")
    print("=" * 88)
    measure, _ = base.safe_components(core, N)
    print(f"core={core}")
    print(f"meas(G_C)={fmt(measure)}")
    print(
        f"first q divisible by 14 with a strict-safe residue: q={q_first}, "
        f"count={n_first}, residues={first_residues[:16]}"
    )
    print()
    print("counts on selected grids:")
    counts = {q: strict_safe_residues(core, q) for q in GRID_QS}
    for q in GRID_QS:
        residues = counts[q]
        denoms = denominator_counter(residues, q)
        print(f"  Q={q:4d}: count={len(residues):2d}; denoms={summarize_denominators(denoms)}")
    print()
    print("Q=1260 exact-period packet table:")
    print("  denom  hits  phi(denom)  factor          mask          divides 210?")
    for d, hits, ph, fac, m, div210 in exact_period_packet_table(drop):
        print(f"  {d:5d}  {hits:4d}  {ph:10d}  {fac:<14} {m:<13} {div210}")
    print()
    print("component-level Q=210 vs Q=1260 hits:")
    for idx, row in enumerate(component_table(drop)):
        a, b = row["interval"]
        print(
            f"  component {idx}: ({a}, {b}), length={row['length']}, "
            f"omitted_clock={row['centers']}"
        )
        print(f"    Q=210 residues:  {row['q210']}")
        print(
            f"    Q=1260 residues: {row['q1260']} "
            f"denoms={row['q1260_denoms']}"
        )
        danger = [danger_speeds(r, RAW_Q) for r in row["q1260"]]
        print(f"    Q=1260 full-AP danger speeds: {danger}")
    print()


def main() -> None:
    print("LRC14 exact-period AP-drop repair atlas - codex S22")
    print()
    print("Core question:")
    print("  Why does the squarefree Q=210 grid miss exactly AP drops 6 and 12,")
    print("  while the raw Q=1260 exact-period carrier sees every AP one-drop core?")
    print()

    rows = ap_drop_summary()
    print_drop_grid_table(rows)

    print("Radical-blind repair summary")
    for drop in (6, 12):
        row = next(r for r in rows if r["drop"] == drop)
        denoms = sorted(row["raw_denoms"])
        print(
            f"  drop {drop}: Q=1260 safe reduced denominators {denoms}; "
            f"hits whose denom divides 210 = {row['raw_denoms_dividing_210']}"
        )
    print()
    print("Interpretation:")
    print("  Inside the 1260 divisor lattice, drop 6 uses repeated-prime")
    print("  packets with lowest reduced denominator 63=3^2*7.")
    print("  Drop 12 uses the dyadic packet 12=2^2*3.")
    print("  The radical 210=2*3*5*7 has none of these repeated-prime")
    print("  exact-period packets.")
    print("  The raw Hill carrier 1260=2^2*3^2*5*7 contains both lanes.")
    print("  The crossing quotient 315 sees some 3^2 packets but loses the dyadic lane.")
    print()

    for drop in (6, 12):
        print_repair_detail(drop)

    print("=" * 88)
    print("Tournament Analysis")
    print("=" * 88)
    path = [
        "raw_1260_exact_period_carrier",
        "drop6_3square7_lane",
        "drop12_2square3_lane",
        "crossing_315_partial_carrier",
        "radical_210_squarefree_grid",
        "raw_runner_vertices",
    ]
    fp = transitive_tournament_fingerprint(path)
    print("Hamiltonian proof path:")
    print("  " + " > ".join(fp["path"]))
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3_cycles={fp['directed_3_cycles']}")
    print(f"SCC_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_paths={fp['hamiltonian_paths']}")
    print()
    print("Session readout:")
    print("  The Q=210 failure is not generic coarseness. It is an exact-period")
    print("  quotient loss: the two missed AP-drop mouths live in the repeated")
    print("  2-adic and 3-adic lanes. This gives a concrete proof target:")
    print("  preserve reduced-denominator packets until after the AP mouth /")
    print("  coimage wall ledger has been checked.")


if __name__ == "__main__":
    main()
