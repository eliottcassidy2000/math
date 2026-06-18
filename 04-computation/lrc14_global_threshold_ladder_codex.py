#!/usr/bin/env python3
"""
lrc14_global_threshold_ladder_codex.py

Codex 2026-06-18: exact threshold-ladder scout for the LRC(14) S3
density-floor route.

The recent B(k)/G_P-intersection work exposed a naming trap:

    rho_alpha(P,E) = meas(G_P cap {x : maxgap({e*x:e in E}) > alpha})

The old via-max sufficient condition used alpha=2/7.  That uniform floor is
false.  The global witness criterion only needs alpha=1/7.  This script asks a
sharper question: how much slack alpha>1/7 survives in adversarial families,
and where exactly does the via-max anti-correlation first bind?

Tournament Analysis declaration.
  Vertex set chosen here: proof obligations / threshold levels, not runners.
  Pairwise observable: over the structured and random case bank, compare
    (zero count, exact minimum rho_alpha, exact minimum critical gap).
  Gauge: orient A -> B if A is more robust: fewer zeros first, then larger
    minimum rho, then larger critical slack.  Ties follow increasing alpha as
    the fixed Hamiltonian path.
  Fingerprints: score histogram, directed 3-cycles, SCC sizes, Hamiltonian path
    count.

Assumption challenge.
  I considered runner, gap, fixed section, safe-component, denominator-event,
  residue, Fourier-mode, relation-lattice, and proof-obligation vertices.  The
  proof-obligation quotient preserves the implication chain needed by the LRC
  predicate ("if this threshold reservoir survives, denominator forcing has
  slack"), but it destroys interval topology, individual blockers, and short
  relation-lattice structure.  The challenged assumption is that the symbol
  rho* can be used without specifying alpha: alpha=2/7 and alpha=1/7 are
  different mathematical objects.
"""

from __future__ import annotations

import itertools
import random
import sys
from collections import Counter, defaultdict, deque
from fractions import Fraction as F

H = F(1, 14)
RNG = random.Random(20260618)

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

THRESHOLDS: list[tuple[str, F]] = [
    ("global_1_7", F(1, 7)),
    ("slack_1_6", F(1, 6)),
    ("slack_3_14", F(3, 14)),
    ("quarter_1_4", F(1, 4)),
    ("near_via_4_15", F(4, 15)),
    ("via_2_7", F(2, 7)),
]


def frac_s(x: F) -> str:
    return str(x)


def merge(arcs: list[tuple[F, F]]) -> list[tuple[F, F]]:
    arcs = sorted(arcs)
    out: list[tuple[F, F]] = []
    for a, b in arcs:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def measure(arcs: list[tuple[F, F]]) -> F:
    return sum((b - a for a, b in arcs), F(0))


def intersect(A: list[tuple[F, F]], B: list[tuple[F, F]]) -> list[tuple[F, F]]:
    A = sorted(A)
    B = sorted(B)
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def complement(arcs: list[tuple[F, F]]) -> list[tuple[F, F]]:
    arcs = merge(arcs)
    out: list[tuple[F, F]] = []
    prev = F(0)
    for a, b in arcs:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        out.append((prev, F(1)))
    return out


def danger_arcs(u: int, h: F = H) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    for j in range(u):
        c = F(j, u)
        a = (c - h / u) % 1
        b = (c + h / u) % 1
        if a < b:
            out.append((a, b))
        else:
            out.append((a, F(1)))
            out.append((F(0), b))
    return out


_SAFE_CACHE: dict[tuple[int, ...], list[tuple[F, F]]] = {}


def safe_set(P: tuple[int, ...]) -> list[tuple[F, F]]:
    P = tuple(sorted(P))
    if P not in _SAFE_CACHE:
        if not P:
            _SAFE_CACHE[P] = [(F(0), F(1))]
        else:
            _SAFE_CACHE[P] = complement(merge([iv for u in P for iv in danger_arcs(u)]))
    return _SAFE_CACHE[P]


def order_breakpoints(E: tuple[int, ...]) -> list[F]:
    E = tuple(sorted(set(E)))
    bps = {F(0), F(1)}
    for a, b in itertools.combinations(E, 2):
        d = b - a
        for m in range(d + 1):
            bps.add(F(m, d))
    return sorted(x for x in bps if 0 <= x <= 1)


_GOOD_CACHE: dict[tuple[tuple[int, ...], F], list[tuple[F, F]]] = {}


def good_set_thr(E: tuple[int, ...], thr: F) -> list[tuple[F, F]]:
    E = tuple(sorted(set(E)))
    key = (E, thr)
    if key in _GOOD_CACHE:
        return _GOOD_CACHE[key]
    good: list[tuple[F, F]] = []
    bps = order_breakpoints(E)
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        pts = sorted(((e * xm) % 1, e) for e in E)
        order = [e for _, e in pts]
        floors = [int((e * xm) // 1) for e in order]
        for idx, e_cur in enumerate(order):
            f_cur = floors[idx]
            if idx < len(order) - 1:
                e_next = order[idx + 1]
                f_next = floors[idx + 1]
                wrap = F(0)
            else:
                e_next = order[0]
                f_next = floors[0]
                wrap = F(1)
            # gap = (e_next*x - f_next) - (e_cur*x - f_cur) + wrap
            A = F(e_next - e_cur)
            C = F(f_cur - f_next) + wrap
            if A == 0:
                if C > thr:
                    good.append((x0, x1))
                continue
            xb = (thr - C) / A
            if A > 0:
                lo, hi = max(x0, xb), x1
            else:
                lo, hi = x0, min(x1, xb)
            if lo < hi:
                good.append((lo, hi))
    _GOOD_CACHE[key] = merge(good)
    return _GOOD_CACHE[key]


def rho_alpha(P: tuple[int, ...], E: tuple[int, ...], alpha: F) -> F:
    return measure(intersect(safe_set(P), good_set_thr(E, alpha)))


def max_gap_on_gp(P: tuple[int, ...], E: tuple[int, ...]) -> F:
    """Supremum of maxgap({e*x}) over x in G_P, exact Fraction."""
    E = tuple(sorted(set(E)))
    gp = safe_set(P)
    best = F(0)
    bps = order_breakpoints(E)
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        cell = [(x0, x1)]
        sub = intersect(cell, gp)
        if not sub:
            continue
        xm = (x0 + x1) / 2
        pts = sorted(((e * xm) % 1, e) for e in E)
        order = [e for _, e in pts]
        floors = [int((e * xm) // 1) for e in order]
        gap_forms: list[tuple[F, F]] = []
        for idx, e_cur in enumerate(order):
            f_cur = floors[idx]
            if idx < len(order) - 1:
                e_next = order[idx + 1]
                f_next = floors[idx + 1]
                wrap = F(0)
            else:
                e_next = order[0]
                f_next = floors[0]
                wrap = F(1)
            gap_forms.append((F(e_next - e_cur), F(f_cur - f_next) + wrap))
        for lo, hi in sub:
            for A, C in gap_forms:
                if A > 0:
                    val = A * hi + C
                elif A < 0:
                    val = A * lo + C
                else:
                    val = C
                if val > best:
                    best = val
    return best


def powerset_P(psz: int) -> list[tuple[int, ...]]:
    return list(itertools.combinations(range(1, 14), psz))


def shapes_for_k(k: int) -> list[tuple[int, ...]]:
    """Adversarial bounded-spread shapes: consecutive, perforated, tails."""
    out: list[tuple[int, ...]] = [tuple(range(k))]

    base1 = list(range(k + 1))
    for drop in range(1, k + 1):
        e = tuple(x for x in base1 if x != drop)
        if len(e) == k and e[0] == 0:
            out.append(e)

    base2 = list(range(k + 2))
    two_drops = {
        (1, 2),
        (1, k),
        (1, k + 1),
        (2, k),
        (2, k + 1),
        (max(1, k - 1), k),
    }
    for d1, d2 in sorted(two_drops):
        if d1 == d2:
            continue
        e = tuple(x for x in base2 if x not in (d1, d2))
        if len(e) == k and e[0] == 0:
            out.append(e)

    # Relation-lattice stress: one or two far tail points glued to a small AP.
    for tail in (k + 1, k + 3, 2 * k, 3 * k):
        if tail > k - 2:
            e = tuple(list(range(k - 1)) + [tail])
            if len(set(e)) == k:
                out.append(e)
    for tail1, tail2 in ((2 * k, 2 * k + 1), (3 * k, 3 * k + 2)):
        if k >= 5:
            e = tuple(list(range(k - 2)) + [tail1, tail2])
            if len(set(e)) == k:
                out.append(e)

    return list(dict.fromkeys(out))


def sample_random_cases(n: int = 500) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    cases: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    seen = set()
    attempts = 0
    while len(cases) < n and attempts < 20 * n:
        attempts += 1
        k = RNG.randint(3, 13)
        psz = 13 - k
        P = tuple(sorted(RNG.sample(range(1, 14), psz)))
        spread = RNG.choice([k - 1, k, k + 1, k + 2, 2 * k, 3 * k])
        if spread < k - 1:
            continue
        body = tuple(sorted(RNG.sample(range(1, spread + 1), k - 1)))
        E = (0,) + body
        key = (P, E)
        if key not in seen:
            seen.add(key)
            cases.append(key)
    return cases


def update_metric(metrics, name: str, alpha: F, rho: F, case):
    m = metrics[name]
    m["tested"] += 1
    if rho == 0:
        m["zeros"] += 1
        if len(m["zero_examples"]) < 6:
            m["zero_examples"].append(case)
    if rho < m["min_rho"]:
        m["min_rho"] = rho
        m["min_case"] = case


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    n = len(adj)
    c = 0
    for a, b, cidx in itertools.combinations(range(n), 3):
        cyc1 = adj[a][b] and adj[b][cidx] and adj[cidx][a]
        cyc2 = adj[a][cidx] and adj[cidx][b] and adj[b][a]
        if cyc1 or cyc2:
            c += 1
    return c


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

    seen = set()
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
                if mask & (1 << w) == 0 and adj[v][w]:
                    dp[mask | (1 << w)][w] += cnt
    return sum(dp[(1 << n) - 1].values())


def obligation_tournament(metrics):
    names = [name for name, _ in THRESHOLDS]
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for i, ni in enumerate(names):
        for j, nj in enumerate(names):
            if i == j:
                continue
            ai = metrics[ni]
            aj = metrics[nj]
            left = (ai["zeros"], -ai["min_rho"])
            right = (aj["zeros"], -aj["min_rho"])
            if left < right:
                adj[i][j] = True
            elif left == right and i < j:
                adj[i][j] = True
    scores = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]
    return {
        "names": names,
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles3": count_directed_3cycles(adj),
        "scc": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_paths_count(adj),
        "leaders": sorted(zip(scores, names), reverse=True),
    }


def print_case_ladder(label: str, P: tuple[int, ...], E: tuple[int, ...]) -> None:
    print(f"\n{label}")
    print(f"  P={P}  E={E}  k={len(E)}  spread={max(E)-min(E)}")
    crit = max_gap_on_gp(P, E)
    print(f"  critical sup alpha on G_P = {crit} = {float(crit):.6f}")
    for name, alpha in THRESHOLDS:
        r = rho_alpha(P, E, alpha)
        print(
            f"    {name:14s} alpha={str(alpha):>5s}: "
            f"rho={str(r):>16s} = {float(r):.6f}"
        )


def main() -> None:
    print("=" * 88)
    print("LRC(14) global-threshold ladder: rho_alpha(P,E)")
    print("=" * 88)
    print("Definitions: G_P={x: ||p x||>=1/14 for all p in P};")
    print("rho_alpha(P,E)=meas(G_P cap {maxgap({e*x})>alpha}).")
    print("alpha=2/7 is the old via-max sufficient condition; alpha=1/7 is the")
    print("global witness threshold.  All decisions below use Fraction arithmetic.")

    print("\n" + "=" * 88)
    print("STEP 1. Explicit anti-correlation ladders")
    print("=" * 88)
    witness_cases = [
        ((1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8)),
        ((1, 2, 3, 12, 13), (0, 2, 3, 4, 5, 6, 7, 8)),
        ((1, 2, 3, 13), (0, 2, 3, 4, 5, 6, 7, 8, 9)),
        ((1, 2, 3), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        ((1, 2, 9), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
    ]
    for idx, (P, E) in enumerate(witness_cases, 1):
        print_case_ladder(f"  Case {idx}", P, E)

    print("\n" + "=" * 88)
    print("STEP 2. Consecutive clusters, all admissible P")
    print("=" * 88)
    all_metrics = {
        name: {
            "tested": 0,
            "zeros": 0,
            "zero_examples": [],
            "min_rho": F(10),
            "min_case": None,
        }
        for name, _ in THRESHOLDS
    }
    consecutive_summary = {}
    for k in range(3, 14):
        E = tuple(range(k))
        psz = 13 - k
        row = {}
        for name, alpha in THRESHOLDS:
            mn = (F(10), None)
            z = 0
            for P in powerset_P(psz):
                r = rho_alpha(P, E, alpha)
                update_metric(all_metrics, name, alpha, r, ("consecutive", P, E))
                if r == 0:
                    z += 1
                if r < mn[0]:
                    mn = (r, P)
            row[name] = (mn, z)
        consecutive_summary[k] = row
        print(f"  k={k:2d}, E=0..{k-1}:")
        for name, _alpha in THRESHOLDS:
            (mn, Pmin), z = row[name]
            print(f"    {name:14s} min={str(mn):>14s} = {float(mn):.6f} zeros={z:4d} at P={Pmin}")

    print("\n" + "=" * 88)
    print("STEP 3. Structured bounded-spread/perforated shapes, all admissible P")
    print("=" * 88)
    structured_count = 0
    for k in range(3, 14):
        psz = 13 - k
        Es = shapes_for_k(k)
        local_metrics = {
            name: {"min_rho": F(10), "min_case": None, "zeros": 0}
            for name, _ in THRESHOLDS
        }
        for E in Es:
            for P in powerset_P(psz):
                structured_count += 1
                for name, alpha in THRESHOLDS:
                    r = rho_alpha(P, E, alpha)
                    update_metric(all_metrics, name, alpha, r, ("structured", P, E))
                    if r == 0:
                        local_metrics[name]["zeros"] += 1
                    if r < local_metrics[name]["min_rho"]:
                        local_metrics[name]["min_rho"] = r
                        local_metrics[name]["min_case"] = (P, E)
        print(f"  k={k:2d}: shapes={len(Es):3d}, cases={len(Es) * len(powerset_P(psz)):7d}")
        for name, _alpha in THRESHOLDS:
            lm = local_metrics[name]
            print(
                f"    {name:14s} min={str(lm['min_rho']):>14s} = "
                f"{float(lm['min_rho']):.6f} zeros={lm['zeros']:5d}"
            )

    print("\n" + "=" * 88)
    print("STEP 4. Random bounded-spread stress")
    print("=" * 88)
    random_cases = sample_random_cases(500)
    random_metrics = {
        name: {"min_rho": F(10), "min_case": None, "zeros": 0}
        for name, _ in THRESHOLDS
    }
    for P, E in random_cases:
        for name, alpha in THRESHOLDS:
            r = rho_alpha(P, E, alpha)
            update_metric(all_metrics, name, alpha, r, ("random", P, E))
            if r == 0:
                random_metrics[name]["zeros"] += 1
            if r < random_metrics[name]["min_rho"]:
                random_metrics[name]["min_rho"] = r
                random_metrics[name]["min_case"] = (P, E)
    print(f"  tested random cases: {len(random_cases)}")
    for name, _alpha in THRESHOLDS:
        rm = random_metrics[name]
        print(
            f"    {name:14s} min={str(rm['min_rho']):>14s} = "
            f"{float(rm['min_rho']):.6f} zeros={rm['zeros']:4d} at {rm['min_case']}"
        )

    print("\n" + "=" * 88)
    print("STEP 5. Aggregate threshold floors and zero examples")
    print("=" * 88)
    print(f"  total structured cases scanned: {structured_count}")
    print("  Aggregate rho_alpha metrics:")
    for name, alpha in THRESHOLDS:
        m = all_metrics[name]
        print(
            f"    {name:14s} alpha={str(alpha):>5s} tested={m['tested']:7d} "
            f"zeros={m['zeros']:6d} min={str(m['min_rho']):>16s} = {float(m['min_rho']):.6f}"
        )
        print(f"      min_case={m['min_case']}")
        if m["zero_examples"]:
            print(f"      first_zero_examples={m['zero_examples']}")

    print("\n  Critical-threshold probe on minima and zero examples:")
    probe_cases = []
    for P, E in witness_cases:
        probe_cases.append(("explicit", P, E))
    for name, _alpha in THRESHOLDS:
        m = all_metrics[name]
        if m["min_case"] is not None:
            _tag, P, E = m["min_case"]
            probe_cases.append((f"min_{name}", P, E))
        for zcase in m["zero_examples"]:
            _tag, P, E = zcase
            probe_cases.append((f"zero_{name}", P, E))
    dedup = []
    seen = set()
    for tag, P, E in probe_cases:
        key = (P, E)
        if key not in seen:
            seen.add(key)
            dedup.append((tag, P, E))
    crit_probe = (F(10), None)
    for tag, P, E in dedup:
        crit = max_gap_on_gp(P, E)
        if crit < crit_probe[0]:
            crit_probe = (crit, (tag, P, E))
        print(f"    {tag:18s} critical={crit} = {float(crit):.6f} P={P} E={E}")
    print(f"  probe critical floor={crit_probe[0]} = {float(crit_probe[0]):.6f} at {crit_probe[1]}")

    print("\n" + "=" * 88)
    print("STEP 6. Proof-obligation tournament")
    print("=" * 88)
    tour = obligation_tournament(all_metrics)
    print("  vertices:", ", ".join(tour["names"]))
    print(f"  score_hist={tour['score_hist']} cycles3={tour['cycles3']} "
          f"scc={tour['scc']} hamiltonian_paths={tour['hamiltonian_paths']}")
    print(f"  leaders={tour['leaders']}")

    print("\n" + "=" * 88)
    print("STEP 7. Component-width diagnostic for finite denominator placement")
    print("=" * 88)
    component_cases = [
        ("quarter_min", (1, 2, 3), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10), F(1, 4)),
        ("near_via_min", (1, 2, 3, 11), (0, 2, 3, 4, 5, 6, 7, 8, 10), F(4, 15)),
        ("via_zero_k7_at_4_15", (1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8), F(4, 15)),
        ("via_zero_k9_at_4_15", (1, 2, 3, 13), (0, 2, 3, 4, 5, 6, 7, 8, 9), F(4, 15)),
    ]
    for label, P, E, alpha in component_cases:
        both = intersect(safe_set(P), good_set_thr(E, alpha))
        widths = [b - a for a, b in both]
        total = measure(both)
        max_width = max(widths) if widths else F(0)
        min_width = min(widths) if widths else F(0)
        print(f"  {label}: alpha={alpha} P={P} E={E}")
        print(
            f"    measure={total} = {float(total):.6f} components={len(both)} "
            f"max_width={max_width} min_width={min_width}"
        )
        print(f"    first_components={both[:6]}")

    print("\n" + "=" * 88)
    print("TAKEAWAY")
    print("=" * 88)
    print("  1. The via-max alpha=2/7 density floor is exactly the wrong target:")
    print("     explicit admissible anti-correlation cases have rho_{2/7}=0.")
    print("  2. The global alpha=1/7 object is much fatter in every case scanned;")
    print("     this is the correct reservoir for a direct LRC witness.")
    print("  3. The critical-gap scan suggests a stronger open subtarget:")
    print("     prove a uniform critical slack sup_{x in G_P} maxgap(E*x) >= beta")
    print("     with beta>1/7, and then separately prove a density floor below beta.")
    print("     The tested anti-correlations bind only at the via boundary, not at")
    print("     the global threshold.")
    print("  4. This does not prove LRC(14).  It moves the proof target from")
    print("     'rho*>=c0' to the precise two-part statement: global/slack reservoir")
    print("     plus Weyl/CRT placement cannot be fully covered by conditional blockers.")


if __name__ == "__main__":
    main()
