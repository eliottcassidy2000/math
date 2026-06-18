#!/usr/bin/env python3
"""
lrc14_colored_discrepancy_bound_codex.py

Codex 2026-06-18: turn the LRC(14) S3 phase-color reservoir into a
deterministic finite-grid lower bound.

The phase-color script found the exact CRT object:

    witnesses at q=14V = colored grid hits in G_P cap C_b(E),
    x = (b+14t)/(14V), b=0,...,13.

This script tests the simplest possible discrepancy lemma.  If a color layer
is a union of m disjoint intervals, then the shifted mesh-1/V grid in that
color has at least

    V * measure(layer) - m

hits in the open layer.  Summing over colors gives

    actual_count(P,E,V) >= V * Sigma(P,E) - K(P,E),

where K is the total number of colored interval components.  This is crude but
rigorous once the interval representation is correct; it converts the
remaining CRT placement problem into:

    V > K(P,E) / Sigma(P,E).

Tournament Analysis declaration.
  Vertex set: phase colors b=0,...,13.
  Pairwise observable: boundary efficiency over the structured bank,
    aggregate_mass_b / max(1, aggregate_component_count_b).
  Gauge/tie path: larger efficiency wins; ties follow increasing color.
  Fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian paths.

Assumption challenge.
  I considered runners, residues, phase colors, interval components, boundary
  endpoints, and proof-obligation cutoffs.  Component counting preserves the
  finite grid-discrepancy predicate, but destroys which runner or relation
  created each endpoint.  The challenged assumption is that the final
  integer-placement step needs number theory first; at this level it may be a
  plain interval-component inequality plus a bounded tail.
"""

from __future__ import annotations

import itertools
import random
import sys
from collections import Counter
from fractions import Fraction as F
from math import ceil, gcd

import lrc14_global_threshold_ladder_codex as ladder
import lrc14_phase_color_reservoir_codex as pc

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

RNG = random.Random(20260618 + 2594)


def floor_frac(x: F) -> int:
    return x.numerator // x.denominator


def ceil_frac(x: F) -> int:
    return -((-x.numerator) // x.denominator)


def open_grid_count(intervals: list[tuple[F, F]], b: int, V: int) -> int:
    """Count t with (b+14t)/(14V) strictly inside the interval union."""
    total = 0
    for lo, hi in intervals:
        lower = (14 * V * lo - b) / 14
        upper = (14 * V * hi - b) / 14
        t0 = floor_frac(lower) + 1
        t1 = ceil_frac(upper) - 1
        t0 = max(t0, 0)
        t1 = min(t1, V - 1)
        if t0 <= t1:
            total += t1 - t0 + 1
    return total


def layer_data(P: tuple[int, ...], E: tuple[int, ...]) -> list[dict]:
    rows = []
    for b in range(14):
        comps = pc.color_components(P, E, b)
        rows.append(
            {
                "b": b,
                "measure": ladder.measure(comps),
                "components": len(comps),
                "max_width": max((hi - lo for lo, hi in comps), default=F(0)),
                "intervals": comps,
            }
        )
    return rows


def sigma_K(P: tuple[int, ...], E: tuple[int, ...]) -> tuple[F, int, list[dict]]:
    rows = layer_data(P, E)
    sig = sum((row["measure"] for row in rows), F(0))
    K = sum((row["components"] for row in rows), 0)
    return sig, K, rows


def cutoff_from(sig: F, K: int) -> int | None:
    if sig <= 0:
        return None
    return floor_frac(F(K, 1) / sig) + 1


def exact_open_count(P: tuple[int, ...], E: tuple[int, ...], V: int) -> int:
    return sum(open_grid_count(row["intervals"], row["b"], V) for row in layer_data(P, E))


def gcd_all(vals) -> int:
    g = 0
    for v in vals:
        g = gcd(g, v)
    return g


def random_shape(k: int, spread: int) -> tuple[int, ...]:
    body = sorted(RNG.sample(range(1, spread + 1), k - 1))
    return (0, *body)


def tail_bank(count: int = 220) -> list[tuple[tuple[int, ...], tuple[int, ...], str]]:
    out = []
    seen = set()
    attempts = 0
    while len(out) < count and attempts < 20 * count:
        attempts += 1
        k = RNG.randint(3, 13)
        psz = 13 - k
        P = tuple(sorted(RNG.sample(range(1, 14), psz)))
        spread = RNG.choice([3 * k, 5 * k, 8 * k, 13 * k, 21 * k, 34 * k])
        E = random_shape(k, spread)
        if gcd_all(P + E) != 1:
            continue
        key = (P, E)
        if key in seen:
            continue
        seen.add(key)
        out.append((P, E, f"k={k},spread={spread}"))
    return out


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = set()
    order: list[int] = []

    def dfs1(v: int) -> None:
        seen.add(v)
        for w, ok in enumerate(adj[v]):
            if ok and w not in seen:
                dfs1(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs1(v)
    seen.clear()
    sizes = []
    for start in reversed(order):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for w, ok in enumerate(radj[v]):
                if ok and w not in seen:
                    seen.add(w)
                    stack.append(w)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def hamiltonian_paths_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [Counter() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, cnt in list(dp[mask].items()):
            if not cnt:
                continue
            for w in range(n):
                if not (mask & (1 << w)) and adj[v][w]:
                    dp[mask | (1 << w)][w] += cnt
    return sum(dp[(1 << n) - 1].values())


def efficiency_tournament(mass: list[F], comps: list[int]) -> dict:
    eff = [mass[b] / max(1, comps[b]) for b in range(14)]
    adj = [[False] * 14 for _ in range(14)]
    for i in range(14):
        for j in range(14):
            if i == j:
                continue
            if eff[i] > eff[j] or (eff[i] == eff[j] and i < j):
                adj[i][j] = True
    scores = [sum(row) for row in adj]
    return {
        "efficiency": eff,
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles3": count_directed_3cycles(adj),
        "scc": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_paths_count(adj),
        "leaders": sorted(zip(scores, range(14), eff), reverse=True),
    }


def print_named() -> None:
    named = [
        ("quarter_min", (1, 2, 3), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        ("near_via_min", (1, 2, 3, 11), (0, 2, 3, 4, 5, 6, 7, 8, 10)),
        ("via_zero_k7", (1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8)),
        ("via_zero_k9", (1, 2, 3, 13), (0, 2, 3, 4, 5, 6, 7, 8, 9)),
        ("broad_1_90", (1, 2, 9), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
    ]
    print("\n" + "=" * 88)
    print("STEP 1. Named hard shapes: component penalty and finite witnesses")
    print("=" * 88)
    for label, P, E in named:
        sig, K, rows = sigma_K(P, E)
        cutoff = cutoff_from(sig, K)
        max_color_K = max(row["components"] for row in rows)
        max_width = max(row["max_width"] for row in rows)
        print(f"  {label}: P={P} E={E}")
        print(
            f"    Sigma={sig} ({float(sig):.6f})  K={K}  "
            f"cutoff V>K/Sigma => {cutoff}  max_color_K={max_color_K} "
            f"max_width={max_width}"
        )
        lifts = pc.candidate_lifts(P, E, max(700, (cutoff or 0) + 80))
        sample_lifts = lifts[:3] + lifts[-3:]
        seen = set()
        for V in sample_lifts:
            if V in seen:
                continue
            seen.add(V)
            S = pc.build_S(P, E, V)
            actual = pc.actual_crt_count(S, V)
            open_count = exact_open_count(P, E, V)
            crude = sig * V - K
            print(
                f"      V={V:4d}: actual={actual:5d} open={open_count:5d} "
                f"V*Sigma-K={crude} ({float(crude):8.2f}) "
                f"certified={crude > 0}"
            )


def scan_structured() -> tuple[list[F], list[int]]:
    print("\n" + "=" * 88)
    print("STEP 2. Structured bank: worst K/Sigma cutoff")
    print("=" * 88)
    agg_mass = [F(0) for _ in range(14)]
    agg_comps = [0 for _ in range(14)]
    global_worst = (F(0), None)
    global_max_K = (0, None)
    global_min_sig = (F(10), None)
    total = 0
    for k in range(3, 14):
        psz = 13 - k
        local_worst = (F(0), None)
        local_min_sig = (F(10), None)
        local_max_K = (0, None)
        for E in ladder.shapes_for_k(k):
            phase_sets = [pc.phase_safe_set(E, b) for b in range(14)]
            for P in ladder.powerset_P(psz):
                gp = ladder.safe_set(P)
                rows = []
                sig = F(0)
                K = 0
                for b in range(14):
                    comps = ladder.intersect(gp, phase_sets[b])
                    m = ladder.measure(comps)
                    c = len(comps)
                    sig += m
                    K += c
                    rows.append((m, c))
                    agg_mass[b] += m
                    agg_comps[b] += c
                total += 1
                ratio = F(K, 1) / sig if sig else F(10**9)
                if ratio > local_worst[0]:
                    local_worst = (ratio, (P, E, sig, K))
                if ratio > global_worst[0]:
                    global_worst = (ratio, (P, E, sig, K))
                if sig < local_min_sig[0]:
                    local_min_sig = (sig, (P, E, K))
                if sig < global_min_sig[0]:
                    global_min_sig = (sig, (P, E, K))
                if K > local_max_K[0]:
                    local_max_K = (K, (P, E, sig))
                if K > global_max_K[0]:
                    global_max_K = (K, (P, E, sig))
        cutoff = cutoff_from(local_worst[1][2], local_worst[1][3]) if local_worst[1] else None
        print(
            f"  k={k:2d}: worst K/Sigma={local_worst[0]} ({float(local_worst[0]):8.3f}) "
            f"cutoff={cutoff:5d} at P,E={local_worst[1][:2]} "
            f"minSigma={local_min_sig[0]} K_at_min={local_min_sig[1][2]} maxK={local_max_K[0]}"
        )
    print(f"\n  total structured cases={total}")
    print(
        f"  GLOBAL worst K/Sigma={global_worst[0]} ({float(global_worst[0]):.3f}) "
        f"at P,E={global_worst[1][:2]} Sigma={global_worst[1][2]} K={global_worst[1][3]}"
    )
    print(
        f"  GLOBAL min Sigma={global_min_sig[0]} ({float(global_min_sig[0]):.6f}) "
        f"at P,E={global_min_sig[1][:2]} K={global_min_sig[1][2]}"
    )
    print(
        f"  GLOBAL max K={global_max_K[0]} at P,E={global_max_K[1][:2]} "
        f"Sigma={global_max_K[1][2]}"
    )
    tour = efficiency_tournament(agg_mass, agg_comps)
    print("\n  Boundary-efficiency phase-color tournament:")
    print(f"    aggregate_mass={[str(x) for x in agg_mass]}")
    print(f"    aggregate_components={agg_comps}")
    print(f"    efficiency={[str(x) for x in tour['efficiency']]}")
    print(
        f"    score_hist={tour['score_hist']} cycles3={tour['cycles3']} "
        f"scc={tour['scc']} hamiltonian_paths={tour['hamiltonian_paths']}"
    )
    print(f"    leaders={[(s,b,str(e)) for s,b,e in tour['leaders'][:8]]}")
    return agg_mass, agg_comps


def scan_tail() -> None:
    print("\n" + "=" * 88)
    print("STEP 3. Random large-spread tail: cutoff compared with natural V scale")
    print("=" * 88)
    worst_cutoff_scale = (F(0), None)
    worst_ratio = (F(0), None)
    min_sig = (F(10), None)
    max_K = (0, None)
    rows = []
    for P, E, tag in tail_bank():
        sig, K, _ = sigma_K(P, E)
        if sig <= 0:
            ratio = F(10**9)
            cutoff = None
        else:
            ratio = F(K, 1) / sig
            cutoff = cutoff_from(sig, K)
        natural = max(E) + 14
        scale = F(cutoff or 10**9, natural)
        rows.append((scale, ratio, cutoff, natural, sig, K, P, E, tag))
        if scale > worst_cutoff_scale[0]:
            worst_cutoff_scale = (scale, rows[-1])
        if ratio > worst_ratio[0]:
            worst_ratio = (ratio, rows[-1])
        if sig < min_sig[0]:
            min_sig = (sig, rows[-1])
        if K > max_K[0]:
            max_K = (K, rows[-1])

    print(f"  tail cases={len(rows)}")
    print(
        f"  worst cutoff/(maxE+14)={worst_cutoff_scale[0]} "
        f"({float(worst_cutoff_scale[0]):.3f})"
    )
    print(f"    row={worst_cutoff_scale[1][6:9]} cutoff={worst_cutoff_scale[1][2]} natural={worst_cutoff_scale[1][3]} Sigma={worst_cutoff_scale[1][4]} K={worst_cutoff_scale[1][5]}")
    print(
        f"  worst K/Sigma={worst_ratio[0]} ({float(worst_ratio[0]):.3f})"
    )
    print(f"    row={worst_ratio[1][6:9]} cutoff={worst_ratio[1][2]} natural={worst_ratio[1][3]} Sigma={worst_ratio[1][4]} K={worst_ratio[1][5]}")
    print(
        f"  min Sigma={min_sig[0]} ({float(min_sig[0]):.6f})"
    )
    print(f"    row={min_sig[1][6:9]} cutoff={min_sig[1][2]} natural={min_sig[1][3]} K={min_sig[1][5]}")
    print(f"  max K={max_K[0]}")
    print(f"    row={max_K[1][6:9]} cutoff={max_K[1][2]} natural={max_K[1][3]} Sigma={max_K[1][4]}")

    print("\n  Five smallest cutoff/natural rows:")
    for scale, ratio, cutoff, natural, sig, K, P, E, tag in sorted(rows)[:5]:
        print(
            f"    scale={float(scale):.3f} cutoff={cutoff} natural={natural} "
            f"K/Sigma={float(ratio):.2f} Sigma={float(sig):.4f} K={K} {tag}"
        )
    print("\n  Five largest cutoff/natural rows:")
    for scale, ratio, cutoff, natural, sig, K, P, E, tag in sorted(rows)[-5:]:
        print(
            f"    scale={float(scale):.3f} cutoff={cutoff} natural={natural} "
            f"K/Sigma={float(ratio):.2f} Sigma={float(sig):.4f} K={K} {tag}"
        )


def scan_tail_exact_discrepancy() -> None:
    print("\n" + "=" * 88)
    print("STEP 4. Exact q=14V discrepancy on random covering tails")
    print("=" * 88)
    # Use an independent deterministic stream so this section is stable even
    # if the preceding component scan changes size.
    local = random.Random(20260618 + 2594 + 77)
    rows = []
    seen = set()
    attempts = 0
    while len(rows) < 160 and attempts < 8000:
        attempts += 1
        k = local.randint(3, 13)
        psz = 13 - k
        P = tuple(sorted(local.sample(range(1, 14), psz)))
        spread = local.choice([3 * k, 5 * k, 8 * k, 13 * k, 21 * k, 34 * k])
        body = tuple(sorted(local.sample(range(1, spread + 1), k - 1)))
        E = (0,) + body
        if gcd_all(P + E) != 1:
            continue
        key = (P, E)
        if key in seen:
            continue
        seen.add(key)
        sig, K, _ = sigma_K(P, E)
        found = 0
        for V in range(max(E) + 14, max(E) + 420):
            S = pc.build_S(P, E, V)
            if len(S) != 13 or min(V - e for e in E) <= 13:
                continue
            if gcd_all(S) != 1 or not pc.is_covering(S):
                continue
            actual = pc.actual_crt_count(S, V)
            expected = sig * V
            ratio = F(actual, 1) / expected if expected else F(0)
            deficit = expected - actual
            rows.append((ratio, deficit, V, sig, K, P, E, f"k={k},spread={spread}", actual))
            found += 1
            if found >= 2 or len(rows) >= 160:
                break

    if not rows:
        print("  no covering tail rows found")
        return
    min_ratio = min(rows, key=lambda r: r[0])
    max_deficit = max(rows, key=lambda r: r[1])
    min_actual = min(rows, key=lambda r: r[8])
    positive_deficits = [r[1] for r in rows if r[1] > 0]
    negative_deficits = [r[1] for r in rows if r[1] < 0]
    print(f"  exact covering rows={len(rows)}")
    print(f"  zero actual witnesses: {sum(1 for r in rows if r[8] == 0)}")
    print(f"  negative deficits (actual > V*Sigma): {len(negative_deficits)}")
    print(
        f"  min actual/(V*Sigma)={min_ratio[0]} ({float(min_ratio[0]):.6f}) "
        f"at V={min_ratio[2]} actual={min_ratio[8]} deficit={min_ratio[1]} "
        f"({float(min_ratio[1]):.3f})"
    )
    print(f"    row={min_ratio[5:8]} Sigma={min_ratio[3]} K={min_ratio[4]}")
    print(
        f"  max positive deficit V*Sigma-actual={max_deficit[1]} "
        f"({float(max_deficit[1]):.3f}) at V={max_deficit[2]} "
        f"ratio={float(max_deficit[0]):.6f} actual={max_deficit[8]}"
    )
    print(f"    row={max_deficit[5:8]} Sigma={max_deficit[3]} K={max_deficit[4]}")
    print(
        f"  min actual count={min_actual[8]} at V={min_actual[2]} "
        f"ratio={float(min_actual[0]):.6f} deficit={float(min_actual[1]):.3f}"
    )
    print(f"    row={min_actual[5:8]} Sigma={min_actual[3]} K={min_actual[4]}")
    if positive_deficits:
        print(
            f"  positive deficit range: min={min(positive_deficits)} "
            f"max={max(positive_deficits)}"
        )
    print("\n  Eight lowest-ratio exact rows:")
    for ratio, deficit, V, sig, K, P, E, tag, actual in sorted(rows)[:8]:
        print(
            f"    ratio={float(ratio):.6f} deficit={float(deficit):7.3f} "
            f"V={V:4d} actual={actual:4d} Sigma={float(sig):.4f} K={K:4d} {tag}"
        )


def verify_interval_lemma() -> None:
    print("\n" + "=" * 88)
    print("STEP 5. Random exact checks of actual >= open >= V*Sigma-K")
    print("=" * 88)
    named = [
        ((1, 2, 3), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        ((1, 2, 3, 11), (0, 2, 3, 4, 5, 6, 7, 8, 10)),
        ((1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8)),
        ((1, 2, 3, 13), (0, 2, 3, 4, 5, 6, 7, 8, 9)),
    ]
    failures = []
    tested = 0
    worst_slack = (F(10**9), None)
    for _ in range(120):
        P, E = RNG.choice(named)
        lifts = pc.candidate_lifts(P, E, 900)
        if not lifts:
            continue
        V = RNG.choice(lifts)
        sig, K, _ = sigma_K(P, E)
        S = pc.build_S(P, E, V)
        actual = pc.actual_crt_count(S, V)
        open_count = exact_open_count(P, E, V)
        lower = sig * V - K
        tested += 1
        slack = F(open_count, 1) - lower
        if slack < worst_slack[0]:
            worst_slack = (slack, (P, E, V, actual, open_count, sig, K, lower))
        if actual < open_count or F(open_count, 1) < lower:
            failures.append((P, E, V, actual, open_count, lower))
    print(f"  tested={tested} failures={len(failures)}")
    print(f"  worst open-(V*Sigma-K) slack={worst_slack[0]} ({float(worst_slack[0]):.3f})")
    print(f"    row={worst_slack[1]}")
    if failures:
        for row in failures[:5]:
            print("    FAILURE", row)


def main() -> None:
    print("=" * 88)
    print("LRC(14) colored discrepancy bound scout")
    print("=" * 88)
    print("Testing actual_count >= open_count >= V*Sigma(P,E)-K(P,E).")
    print("K is the total number of components over the 14 phase-color layers.")
    print_named()
    scan_structured()
    scan_tail()
    scan_tail_exact_discrepancy()
    verify_interval_lemma()
    print("\n" + "=" * 88)
    print("TAKEAWAY")
    print("=" * 88)
    print("  1. The exact phase-color reservoir admits a deterministic lower bound:")
    print("        actual_count(q=14V) >= V*Sigma(P,E) - K(P,E).")
    print("  2. This proves every fixed (P,E) shape for all V > K/Sigma, reducing")
    print("     the finite placement problem to a computable small-V check.")
    print("  3. The remaining proof burden is now sharper: prove a uniform Sigma")
    print("     floor and a relation-lattice/tail bound showing K/Sigma is below")
    print("     the available V scale, or enumerate the finite residue tail.")
    print("  4. LRC(14) is still not proved; the obstruction has been compressed")
    print("     from colored CRT placement to a boundary-complexity inequality.")


if __name__ == "__main__":
    main()
